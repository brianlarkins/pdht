/* CHT-MPI, version 1.0. An MPI-based Chunks and Tasks library
 *                       implementation.
 * For copyright and license information, see below under "Copyright and
 * license".
 * 
 * Primary academic reference: 
 * Chunks and Tasks: A programming model for parallelization of dynamic
 * algorithms,
 * Emanuel H. Rubensson and Elias Rudberg,
 * Parallel Computing 00, 00 (2013),
 * <http://dx.doi.org/10.1016/j.parco.2013.09.006>
 * 
 * For further information about Chunks and Tasks, see
 * <http://www.chunks-and-tasks.org>.
 * 
 * === Copyright and license ===
 * 
 * Copyright (c) 2009-2014 Emanuel H. Rubensson and Elias Rudberg. All
 *                         rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * 
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer listed
 *   in this license in the documentation and/or other materials provided
 *   with the distribution.
 * 
 * - Neither the name of the copyright holders nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * The copyright holders provide no reassurances that the source code
 * provided does not infringe any patent, copyright, or any other
 * intellectual property rights of third parties. The copyright holders
 * disclaim any liability to any recipient for claims brought against
 * recipient by any third party for infringement of that parties
 * intellectual property rights.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "TaskSchedulerService_parent.h"
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include "services/services_utils.h"

namespace cht {
  namespace TaskSchedulerService {

    void Parent::start_derived() {
      AccessKey key(this);
      OS.outputInfo("This is TaskSchedulerService::start() (parent)");
      // Set task type string <-> int maps and send to workers
      std::list<std::string> strList;
      cht::obj_factory<Task>::instance().getObjectTypeIDStrings(strList);
      cht::Service::populateMapsForGivenStrList(strList,
						taskTypeIDToIntIDMap,
						intIDToTaskTypeIDMap);
      cht::Service::sendStrBufToWorkers(strList, 
					TAG_TaskTypeID_map,
					key.comm_to_workers(),
					MW);      

      // Send parameters to workers.
      int intVector[8];
      intVector[0] = n_worker_threads;
      intVector[1] = n_prefetch_tasks;
      intVector[2] = (int)stealing_disabled_flag;
      intVector[3] = (int)steal_only_from_0_flag;
      intVector[4] = status_report_file_interval_seconds;
      intVector[5] = (int)do_mutex_lock_checking;
      intVector[6] = child_chunk_fetch_level;
      intVector[7] = status_report_interval_milliseconds;

      for(int i = 0; i < key.n_workers(); i++)
	MW._MPI_Send(intVector, 8*sizeof(int), 
		    MPI_UNSIGNED_CHAR, i, 
		    TAG_Thread_params, *key.comm_to_workers());  
    }

    void Parent::stop_derived() {
      AccessKey key(this);
      OS.outputInfo("This is TaskSchedulerService::stop() (parent)");
      for(int i = 0; i < key.n_workers(); i++)
	MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, i, 
		    TAG_Workers_finished, *key.comm_to_workers(), MPI_STATUS_IGNORE);
      for(int i = 0; i < key.n_workers(); i++)
	MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, i, 
		    TAG_Workers_finished, *key.comm_to_workers());
    }

    // FIXME: executeMotherTask should be template function 
    ChunkID Parent::
    executeMotherTask(std::string taskTypeID,
		      std::vector<ChunkID> const & taskInputChunks) {
      double starttime = cht::get_wall_seconds();
      // Create new TaskInfoStruct object.
      int taskTypeInt = taskTypeIDToIntIDMap[taskTypeID];
      int ownerRank = -1; // Use ownerRank=-1 for parent.
      TaskID taskID(ownerRank, ++idCounter, taskTypeInt); 
      int ancestorDepth = 0;
      TaskInfoStruct* taskInfo = new TaskInfoStruct(taskID, TASK_ID_NULL, NULL, taskTypeID, taskInputChunks, ownerRank, ancestorDepth);
      AccessKey key(this);
      // Send task info to a worker
      std::stringstream ss1;
      ss1 << taskInfo->str() 
	  << " execute_mother_task  taskType = '" << taskTypeID << "'.";
      OS.outputInfo(ss1.str());
      // Decide which worker to send the task to.
      int rankToSendTo = 0;
      if(send_to_separate_workers_flag) {
	AccessKey key(this);
	rankToSendTo = workerRankCounter % key.n_workers();
	workerRankCounter++;
      }
      std::stringstream ss1b;
      ss1b << "Sending task info to worker with rank " << rankToSendTo << "...";
      OS.outputInfo(ss1b.str());
      int buffer_length_needed = taskInfo->pack_size();
      std::vector<char> msgBuf(buffer_length_needed);
      taskInfo->pack(&msgBuf[0]);
      MW._MPI_Send(&msgBuf[0], buffer_length_needed, MPI_UNSIGNED_CHAR, rankToSendTo, TAG_Mother_task_info, *key.comm_to_workers());
      OS.outputInfo("MPI_Send for sending task info to a worker finished.");
      OS.outputInfo("Now waiting for TAG_Mother_task_finished message from that worker.");
      // Note: we could just use MPI_Recv here, but that turned out to be
      // consuming a lot of cpu time for some MPI implementations,
      // e.g. openmpi/1.3.4 on isis. To be sure that the parent process
      // really do not use any cpu time, we use a loop with sleep() and
      // Iprobe() instead.
      while(1) {
	int flag;
	MW._MPI_Iprobe(rankToSendTo, TAG_Mother_task_finished, *key.comm_to_workers(), &flag, MPI_STATUS_IGNORE);
	if(flag)
	  break;
	int one_milli_second = 1000;
	usleep(100*one_milli_second);
      }
      ChunkID result_chunkID;
      MW._MPI_Recv(&result_chunkID, sizeof(ChunkID), MPI_UNSIGNED_CHAR, rankToSendTo, TAG_Mother_task_finished, *key.comm_to_workers(), MPI_STATUS_IGNORE);
      OS.outputInfo("cht::TaskSchedulerService::Parent::executeMotherTask done!");
      double finishtime = cht::get_wall_seconds();
      std::stringstream ss;
      ss << "executeMotherTask for '" << taskTypeID << "' took " << (finishtime - starttime) << " wall seconds.";
      OS.outputInfo(ss.str());
      return result_chunkID;
    }

    void Parent::setNoOfWorkerThreads(int n_worker_threads_) {
      // This routine may only be called when the service is stopped,
      // since the parameters are transferred to workers in start(). We
      // check that service is stopped by demanding that comm_to_workers
      // pointer is NULL.
      if ( serviceIsRunning() )
	throw std::runtime_error("Error! cht::TaskSchedulerService::Parent::setNoOfWorkerThreads() called while service running.");
      // Be aware that service in principle could start to run while
      // setting the parameters below...
      n_worker_threads = n_worker_threads_;
    }

    void Parent::setNoOfPrefetchTasks(int n_prefetch_tasks_) {
      // This routine may only be called when the service is stopped,
      // since the parameters are transferred to workers in start(). We
      // check that service is stopped by demanding that comm_to_workers
      // pointer is NULL.
      if ( serviceIsRunning() )
	throw std::runtime_error("Error! cht::TaskSchedulerService::Parent::setNoOfPrefetchTasks() called while service running.");
      // Be aware that service in principle could start to run while
      // setting the parameters below...
      n_prefetch_tasks = n_prefetch_tasks_;
    }

    void Parent::setChildChunkFetchLevel(int child_chunk_fetch_level_) {
      // This routine may only be called when the service is stopped,
      // since the parameters are transferred to workers in start(). We
      // check that service is stopped by demanding that comm_to_workers
      // pointer is NULL.
      if ( serviceIsRunning() )
	throw std::runtime_error("Error! cht::TaskSchedulerService::Parent::setChildChunkFetchLevel() called while service running.");
      // Be aware that service in principle could start to run while
      // setting the parameters below...
      child_chunk_fetch_level = child_chunk_fetch_level_;
    }

    void Parent::setDebugParams(bool stealing_disabled_flag_,
				bool steal_only_from_0_flag_,
				bool send_to_separate_workers_flag_,
				int status_report_file_interval_seconds_,
				int status_report_interval_milliseconds_,
				bool do_mutex_lock_checking_) {
      // This routine may only be called when the service is stopped,
      // since the parameters are transferred to workers in start(). We
      // check that service is stopped by demanding that comm_to_workers
      // pointer is NULL.
      if ( serviceIsRunning() )
	throw std::runtime_error("Error! cht::TaskSchedulerService::Parent::setDebugParams() called while service running.");
      // Be aware that service in principle could start to run while
      // setting the parameters below...
      status_report_file_interval_seconds = status_report_file_interval_seconds_;
      status_report_interval_milliseconds = status_report_interval_milliseconds_;
      stealing_disabled_flag = stealing_disabled_flag_;
      steal_only_from_0_flag = steal_only_from_0_flag_;
      send_to_separate_workers_flag = send_to_separate_workers_flag_;
      do_mutex_lock_checking = do_mutex_lock_checking_;
    }

    // Default values for n threads 0 here, indicating that workers
    // should decide how many threads to use, can be changed in setParams
    // function.
    Parent::Parent() 
      : OS( OutputService::Parent::instance() ),
	n_worker_threads(0), 
	n_prefetch_tasks(3), // Default n_prefetch_tasks value set here
	child_chunk_fetch_level(-1), // Default child_chunk_fetch_level value set here
	stealing_disabled_flag(false),
	steal_only_from_0_flag(false),
	send_to_separate_workers_flag(false),
	status_report_file_interval_seconds(0),
	status_report_interval_milliseconds(0),
	do_mutex_lock_checking(false),
	workerRankCounter(0)
    {
  
    }


  }; // end namespace
}; // end namespace
