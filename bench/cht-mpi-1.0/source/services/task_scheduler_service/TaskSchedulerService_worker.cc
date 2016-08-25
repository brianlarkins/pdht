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
#include "TaskSchedulerService_worker.h"
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include "services/services_utils.h"
#include "services/chunk_obj_service/ChunkObjService_worker.h"
#include "services/chunk_obj_cache_service/ChunkObjCacheService_worker.h"

static std::string GetFileAndLineString(const char* file, int line) {
  std::stringstream ss;
  ss << "File: " << file << " Line: " << line;
  return ss.str();
}

#define LockMutex_macro LockMutex(GetFileAndLineString(__FILE__, __LINE__));

#define DebugOutput_macro  {   usleep(44); std::cout << "DebugOutput: pthread_self = " << pthread_self() << "   " << GetFileAndLineString(__FILE__, __LINE__) << std::endl; usleep(44); } 

namespace cht {
  namespace TaskSchedulerService {

    static unsigned int detect_number_of_cores() {
      const char *env = getenv("OMP_NUM_THREADS");
      int defThreads = 1;
      if ( !(env && (defThreads=atoi(env)) > 0) ) {
	FILE *f = fopen("/proc/cpuinfo", "rt");
	if(f) {
	  char line[256];
	  defThreads = 0;
	  while(fgets(line, sizeof(line), f))
	    if(strncmp(line, "processor", 9) == 0)
	      defThreads++;
	  fclose(f);
	  /* Protect against case when /proc/cpuinfo exits but
	     contains garbage. Unlikely but possible. */
	  if(defThreads == 0)
	    defThreads = 1;
	}
      }
      if (defThreads <= 0)
	return 1;
      return (unsigned int)defThreads;
    }  


    void* global_thread_func_main_worker(void* arg)
    {
      Worker* p = (Worker*) arg;
      try {
	p->WorkerFunc();
      } 
      catch (std::runtime_error e) {
	std::cerr << "Error! Exception std::runtime_error caught in TaskSchedulerService "
	  "global_thread_func_main_worker()." << std::endl << "what() = '" << e.what() << "'" << std::endl;
        p->call_mpi_abort();
      }
      catch (std::exception & e) {
	std::cerr << "Error! Exception std::exception caught in TaskSchedulerService "
	  "global_thread_func_main_worker()." << std::endl << "what() = '" << e.what() << "'" << std::endl;
        p->call_mpi_abort();
      }
      catch ( ... ) {
	std::cerr << "Error! Exception caught in TaskSchedulerService global_thread_func_main_worker()." << std::endl;
        p->call_mpi_abort();
      }
      return NULL;
    }

    void* global_thread_func_request_handler(void* arg)
    {
      Worker* p = (Worker*) arg;
      try {
	p->RequestHandlerFunc();
      }
      catch (std::runtime_error e) {
	std::cerr << "Error! Exception std::runtime_error caught in TaskSchedulerService "
	  "global_thread_func_request_handler()." << std::endl << "what() = '" << e.what() << "'" << std::endl;
        p->call_mpi_abort();
      }
      catch (std::exception & e) {
	std::cerr << "Error! Exception std::exception caught in TaskSchedulerService "
	  "global_thread_func_request_handler()." << std::endl << "what() = '" << e.what() << "'" << std::endl;
        p->call_mpi_abort();
      }
      catch ( ... ) {
	std::cerr << "Error! Exception caught in TaskSchedulerService global_thread_func_request_handler()." << std::endl;
        p->call_mpi_abort();
      }
      return NULL;
    }

    void* global_thread_func_data_fetcher(void* arg)
    {
      Worker* p = (Worker*) arg;
      try {
	p->DataFetcherFunc();
      }
      catch (std::runtime_error e) {
	std::cerr << "Error! Exception std::runtime_error caught in TaskSchedulerService "
	  "global_thread_func_data_fetcher()." << std::endl << "what() = '" << e.what() << "'" << std::endl;
        p->call_mpi_abort();
      }
      catch (std::exception & e) {
	std::cerr << "Error! Exception std::exception caught in TaskSchedulerService "
	  "global_thread_func_data_fetcher()." << std::endl << "what() = '" << e.what() << "'" << std::endl;
        p->call_mpi_abort();
      }
      catch ( ... ) {
	std::cerr << "Error! Exception caught in TaskSchedulerService global_thread_func_data_fetcher()." << std::endl;
        p->call_mpi_abort();
      }
      return NULL;
    }

    void* global_thread_func_monitoring(void* arg)
    {
      Worker* p = (Worker*) arg;
      try {
	p->MonitoringFunc();
      }
      catch (std::runtime_error e) {
	std::cerr << "Error! Exception std::runtime_error caught in TaskSchedulerService "
	  "global_thread_func_monitoring()." << std::endl << "what() = '" << e.what() << "'" << std::endl;
        p->call_mpi_abort();
      }
      catch (std::exception & e) {
	std::cerr << "Error! Exception std::exception caught in TaskSchedulerService "
	  "global_thread_func_monitoring()." << std::endl << "what() = '" << e.what() << "'" << std::endl;
        p->call_mpi_abort();
      }
      catch ( ... ) {
	std::cerr << "Error! Exception caught in TaskSchedulerService global_thread_func_monitoring()." << std::endl;
        p->call_mpi_abort();
      }
      return NULL;
    }

    void* global_thread_func_for_mutex_lock(void* arg) {
      Worker* p = (Worker*) arg;
      p->MutexLockThreadFunc();
      return NULL;
    }

    void Worker::LockMutex(std::string debugInfoStr) {
      mutex.lock();
      mutexIsLocked = true;
      mutexLockDebugStr = debugInfoStr;
    }

    void Worker::UnlockMutex() {
      assert(mutexIsLocked);
      mutexLockDebugStr = "";
      mutexIsLocked = false;
      mutex.unlock();
    }


    void Worker::call_mpi_abort() {
      AccessKey key(this);
      MW._MPI_Abort(*key.comm_to_parent(), -1);
    }

    void Worker::start_derived() {
      AccessKey key(this);
      OS.outputInfo("This is TaskSchedulerService::start() (worker)");
      // Get task type string <-> int maps from parent
      std::list<std::string> strList;
      cht::Service::receiveStrBufFromParent(strList, 
					    TAG_TaskTypeID_map,
					    key.comm_to_parent(),  
					    MW);
      cht::Service::populateMapsForGivenStrList(strList,
						taskTypeIDToIntIDMap,
						intIDToTaskTypeIDMap);
      // Check that tasktype ID strings match the ones in local
      // Task factory
      std::list<std::string> localStrList;
      cht::obj_factory<Task>::instance().getObjectTypeIDStrings(localStrList);
      cht::Service::checkStrListAgainstMap(localStrList, taskTypeIDToIntIDMap);
      
      // Create two threads here: main work thread and request handler
      // thread.
      finishFlagForThreadCommunicationWorker = false;
      steal_attempt_in_progress = false;

      // Initialize random numbers using own rank, to avoid all workers
      // drawing the same random numbers.
      srand(key.my_rank());

      // Receive parameters from parent.
      int intVector[8];
      MW._MPI_Recv(intVector, 8*sizeof(int), MPI_UNSIGNED_CHAR, 0, 
		  TAG_Thread_params, *key.comm_to_parent(), MPI_STATUS_IGNORE);
      n_worker_threads = intVector[0];
      n_prefetch_tasks = intVector[1];
      stealing_disabled_flag = (bool)intVector[2];
      steal_only_from_0_flag = (bool)intVector[3];
      status_report_file_interval_seconds = intVector[4];
      do_mutex_lock_checking = (bool)intVector[5];
      child_chunk_fetch_level = intVector[6];
      status_report_interval_milliseconds = intVector[7];
      if(n_worker_threads < 1) {
	// Set n_worker_threads to number of cores minus 1 (or to 1 if there is only one core)
	int n_cores = detect_number_of_cores();
	n_worker_threads = n_cores > 1 ? n_cores - 1 : 1 ;
      }

      if (threadGroupWorkers || threadRequestHandler)
	throw std::runtime_error("threadGroupWorkers || threadRequestHandler in "
				 "cht::TaskSchedulerService::Worker::start_derived()");
      // threadGroupWorkers and threadRequestHandler are not
      // protected, only adressed in calls to new and delete
      threadGroupWorkers = new Threads::ThreadGroup("TS-worker", global_thread_func_main_worker, this,
						    n_worker_threads, n_worker_threads);
      threadRequestHandler = new Threads::Thread("TS-request", global_thread_func_request_handler, this);
      threadDataFetcher = new Threads::Thread("TS-fetcher", global_thread_func_data_fetcher, this);
      threadMonitoring = new Threads::Thread("TS-monitoring", global_thread_func_monitoring, this);
    }

    void Worker::stop_derived() {
      AccessKey key(this);
      OS.outputInfo("This is TaskSchedulerService::stop() (worker)");
      // Set flag to notify other threads that they should exit now.
      LockMutex_macro;
      finishFlagForThreadCommunicationWorker = true;
      data_fetcher_cond.signal();
      UnlockMutex();
      OS.outputInfo("TaskSchedulerService::stop() (worker) now waiting for threads to finish.");
      delete threadGroupWorkers;
      threadGroupWorkers = NULL;
      // signal to parent: my workers are ready
      MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, TAG_Workers_finished, *key.comm_to_parent());  
      // Request handler will receive signal from parent that all other workers are ready
      delete threadRequestHandler;
      threadRequestHandler = NULL;
      OS.outputInfo("TaskSchedulerService::stop() (worker) done!");
    }


    int Worker::GetTotNoOfTasks()
    {
      assert(mutexIsLocked);
      int sum = 0;
      sum += list_waiting_before_transaction.size();
      sum += list_pending_for_transaction.size();
      sum += list_doing_transaction_nonleaf_1.size();
      sum += list_doing_transaction_nonleaf_2.size();
      sum += list_doing_transaction_leaf.size();
      sum += list_waiting_before_finalize.size();
      sum += list_pending_for_finalize.size();
      sum += list_running_finalize.size();
      sum += list_stolen.size();
      return sum;
    }

    bool Worker::ExistsInMap(const TaskListMap & map, TaskID id)
    {
      assert(mutexIsLocked);
      if(map.find(id) != map.end())
	return true;
      return false;
    }

    void Worker::ReinsertTaskToPendingForTransactionLists(TaskInfoStruct* taskInfo) {
      assert(mutexIsLocked);
      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
      // Verify that this taskInfo does not already exist in maps.
      assert(helper_list_pending_for_transaction.find(taskInfo->get_id()) == helper_list_pending_for_transaction.end());
      assert(list_pending_for_transaction.find(TaskInfoPtr(*taskInfo)) == list_pending_for_transaction.end());
      // insert
      helper_list_pending_for_transaction.insert(TaskListMap::value_type(taskInfo->get_id(), taskInfo));
      list_pending_for_transaction.insert(TaskListSortedMap::value_type(TaskInfoPtr(*taskInfo), taskInfo));
      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
      // Let data fetcher know that we have added something to the pending_for_transaction list
      data_fetcher_cond.signal();
    }

    void Worker::AddTaskToPendingForTransactionLists(TaskInfoStruct* taskInfo) {
      AccessKey key(this);
      // It is assumed here that key outlives the following parameter copies
      int n_workers = key.n_workers();
      int my_rank = key.my_rank();

      assert(mutexIsLocked);
      assert(taskInfo != NULL);
      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());

      // Verify that this taskInfo does not already exist in maps.
      assert(helper_list_pending_for_transaction.find(taskInfo->get_id()) == helper_list_pending_for_transaction.end());
      assert(list_pending_for_transaction.find(TaskInfoPtr(*taskInfo)) == list_pending_for_transaction.end());
      
      // Before inserting into pending_for_transaction lists, populate ownerList. Note that this must be done before inserting into list, sine the ownerList is used by the < operator.
      // To populate ownerList, check which workers are the largest owners for the input chunks for this task.
      taskInfo->ownerList.clear();
      taskInfo->totInputDataSize = 0;
      std::vector<ChunkID> & inputChunkIDs = taskInfo->inputChunks;
      for(int i = 0; i < inputChunkIDs.size(); i++) {
	if(inputChunkIDs[i] != CHUNK_ID_NULL) {
	  taskInfo->totInputDataSize += inputChunkIDs[i].totalDeepSize;
	  for(int k = 0; k < ChunkID::MAX_NO_OF_OWNERS_IN_LIST; k++) {
	    int currRank = inputChunkIDs[i].majorOwners[k].rank;
	    size_t currSize = inputChunkIDs[i].majorOwners[k].size;
	    if(currRank == -1)
	      continue;
	    // Check if currRank already present in ownerList.
	    std::list<ChunkID::MajorOwnerInfo>::iterator it = taskInfo->ownerList.begin();
	    while(it != taskInfo->ownerList.end()) {
	      if(it->rank == currRank) {
		it->size += currSize;
		break;
	      }
	      it++;
	    }
	    if(it == taskInfo->ownerList.end()) {
	      // Not found
	      taskInfo->ownerList.push_back(inputChunkIDs[i].majorOwners[k]);
	    }
	  }
	}
      } // end for i
      // OK, now we have populated the ownerList for this task.
      // Check if this worker is among the major owners.
      taskInfo->currWorkerOwnedSize = 0;
      std::list<ChunkID::MajorOwnerInfo>::iterator it = taskInfo->ownerList.begin();
      while(it != taskInfo->ownerList.end()) {
	assert(it->rank >= 0 && it->rank < n_workers);
	if(it->rank == my_rank)
	  taskInfo->currWorkerOwnedSize = it->size;
	it++;
      }

      helper_list_pending_for_transaction.insert(TaskListMap::value_type(taskInfo->get_id(), taskInfo));
      list_pending_for_transaction.insert(TaskListSortedMap::value_type(TaskInfoPtr(*taskInfo), taskInfo));
      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
      // Find out if fallback execute should be used 
      int nInpChunks = taskInfo->inputChunks.size();
      std::vector<bool> input_type_is_ChunkID(nInpChunks);
      for(int ind = 0;ind < nInpChunks; ind ++) {
	// If any of the ChunkIDs for input parameters that do
	// not have the type "ChunkID" are CHUNK_ID_NULL we
	// want to use the fallback execute() that the user
	// hopefully has implemented in this case. 
	input_type_is_ChunkID[ind] = 
	  cht::arg_manager<Task>::instance().check(taskInfo->get_taskTypeID(), "ChunkID", ind);
	if ( !input_type_is_ChunkID[ind] ) {
	  if ( taskInfo->inputChunks[ind] == CHUNK_ID_NULL )
	    taskInfo->fallback_execute_should_be_used = true; 
	  else
	    if (!cht::arg_manager<Task>::instance().check(taskInfo->get_taskTypeID(), 
							  ChunkObjService::Worker::instance().getChunkTypeIDStr(taskInfo->inputChunks[ind].chunkTypeID), 
							  ind)) {
	      std::cerr << "Error: Wrong chunk type in input to task. taskTypeID = " << taskInfo->get_taskTypeID() << 
		" ind = " << ind << " ChunkTypeIDStr = " << ChunkObjService::Worker::instance().getChunkTypeIDStr(taskInfo->inputChunks[ind].chunkTypeID) << std::endl;
	      throw std::runtime_error("Wrong chunk type in input to task"); 
	    }
	  // FIXME: Better error message! Print task type, index, and chunk type 
	  // ELIAS NOTE 2013-10-09: Really important to have a better error message here. Very annoying for the user otherwise.
	}
      }
      // Try to populate taskInfo->input_chunks with local chunks
      bool all_chunks_found_locally = true;
      if ( taskInfo->fallback_execute_should_be_used == false ) {
	// Ok, the fallback execute will not be used
	taskInfo->input_chunks.resize(nInpChunks);
	for(int ind = 0;ind < nInpChunks; ind ++) {
	  if ( !input_type_is_ChunkID[ind] ) {
	    if ( taskInfo->inputChunks[ind] == CHUNK_ID_NULL )
	      throw std::runtime_error("taskInfo->inputChunks[ind] == CHUNK_ID_NULL when fallback execute will not be used?!?");
	    // If the chunk does not reside locally the following call should not do anything
	    if (!ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().getChunkIfLocal(taskInfo->inputChunks[ind], taskInfo->input_chunks[ind]))
	      all_chunks_found_locally = false;
	  }
	} // end for
      }
      pendingTaskCount_tot++;
      if ( taskInfo->fallback_execute_should_be_used || all_chunks_found_locally ) {
	taskInfo->ready_to_run = true;
 	pendingTaskCount_readyToRun++;
      }
      else
	// Let data fetcher know that we have added tasks that require remote chunks to the pending_for_transaction list 
	data_fetcher_cond.signal();
    }

    void Worker::RemoveTaskFromPendingForTransactionLists(TaskInfoStruct* taskInfo) {
      assert(mutexIsLocked);
      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());

      int list_pending_for_transaction_size_before = list_pending_for_transaction.size();
      int helper_list_pending_for_transaction_size_before = helper_list_pending_for_transaction.size();

      if (!helper_list_pending_for_transaction.erase(taskInfo->get_id()))
	throw std::runtime_error("Failed to erase element in helper_list_pending_for_transaction");
      if (!list_pending_for_transaction.erase(TaskInfoPtr(*taskInfo)))
	throw std::runtime_error("Failed to erase element in list_pending_for_transaction");

      if(list_pending_for_transaction.size() != helper_list_pending_for_transaction.size()) {
	std::cerr << "list_pending_for_transaction.size()        = " << list_pending_for_transaction.size()        << std::endl;
	std::cerr << "helper_list_pending_for_transaction.size() = " << helper_list_pending_for_transaction.size() << std::endl;
	std::cerr << "list_pending_for_transaction_size_before        = " << list_pending_for_transaction_size_before        << std::endl;
	std::cerr << "helper_list_pending_for_transaction_size_before = " << helper_list_pending_for_transaction_size_before << std::endl;
	for(TaskListMap::iterator it = helper_list_pending_for_transaction.begin(); it != helper_list_pending_for_transaction.end(); it++)
	  std::cerr << "helper_list_pending_for_transaction element: " << (it->second)->get_id().str() << std::endl;
      }

      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
      // Let data fetcher know that we have removed a task from the pending_for_transaction list
      data_fetcher_cond.signal();
    }

    void Worker::RemoveTaskFromWaitLists(TaskID tid) {
      assert(mutexIsLocked);
      TaskIDList tasksToMove;
      TaskListMap::iterator it;
      for ( it=list_waiting_before_transaction.begin() ; it != list_waiting_before_transaction.end(); it++ ) {
	TaskInfoStruct* taskInfo = (*it).second;
	assert((*it).first == taskInfo->get_id());
	taskInfo->waitList.remove(tid);
	if(taskInfo->waitList.empty()) {
	  // OK, waitList is empty for that task. We can move it from list_waiting_before_transaction to list_pending_for_transaction.
	  tasksToMove.push_back(taskInfo->get_id());
	}
      }
      TaskIDList::const_iterator it2;
      for( it2 = tasksToMove.begin() ; it2 != tasksToMove.end(); it2++ ) {
	TaskListMap::const_iterator itTmp = list_waiting_before_transaction.find(*it2);
	TaskInfoStruct* taskInfo = (*itTmp).second;
	setTaskInputChunkIDs(*taskInfo);
	list_waiting_before_transaction.erase(taskInfo->get_id());
	assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
	AddTaskToPendingForTransactionLists(taskInfo);
      }
    }

    // NOTE: mutex should be locked when calling this function.
    void Worker::RemoveTaskFromChildListOfCreator(TaskInfoStruct* taskInfo) {
      assert(mutexIsLocked);
      TaskID creatorTaskID = taskInfo->get_creator();
      // Look for creator in various lists where it could be
      bool creatorFound = false;
      TaskListMap::iterator it;
      if(!creatorFound) {
	it = list_doing_transaction_nonleaf_2.find(creatorTaskID);
	if(it != list_doing_transaction_nonleaf_2.end()) {
	  creatorFound = true;
	}
      }
      if(!creatorFound) {
	it = list_waiting_before_finalize.find(creatorTaskID);
	if(it != list_waiting_before_finalize.end()) {
	  creatorFound = true;
	}
      }
      if(!creatorFound) {
	// This could happen if this task was stolen
	// FIXME: verify here that task was stolen
	return;
      }
      TaskInfoStruct* creatorTaskInfo = (*it).second;
      creatorTaskInfo->childList.remove(taskInfo->get_id());
    }

    bool Worker::UpdateLists()
    {
      assert(mutexIsLocked);
      TaskListMap::iterator it;
      // Check if any task can be moved from list_waiting_before_finalize
      // to list_pending_for_finalize.
      for ( it=list_waiting_before_finalize.begin() ; it != list_waiting_before_finalize.end(); it++ )
	{
	  TaskInfoStruct* taskInfo = (*it).second;
	  // Check if all children have been removed for this task.
	  int noOfChildTasks = taskInfo->childList.size();
	  // Check if all registerChunk and copyChunk operations are finished for this task.
	  assert(taskInfo->get_task() != 0);
	  bool allRegisterAndCopyChunkOpsHaveFinished = taskInfo->get_task()->internal_allRegisterAndCopyChunkOpsHaveFinished();
	  if(noOfChildTasks == 0 && allRegisterAndCopyChunkOpsHaveFinished)
	    {
	      std::stringstream s;
	      s << "Moving to list_pending_for_finalize, taskInfo->id = " << taskInfo->get_id().str();
	      OS.outputDebug(s.str());
	      list_waiting_before_finalize.erase(taskInfo->get_id());
	      list_pending_for_finalize.insert(TaskListMap::value_type(taskInfo->get_id(), taskInfo));
	      return true;
	    }
	}
      return false;
    }


    void Worker::resetStatistics() {
      LockMutex_macro;
      timerForEverything.reset();
      taskStatistics.clear();
      stealAttemptStatistics.clear();
      stealAttemptSendStatistics.clear();
      RequestHandlerFuncStatisticsParent.clear();
      RequestHandlerFuncStatisticsWorker.clear();
      accumulatedIdleTime = 0;
      pendingTaskCount_tot = 0;
      pendingTaskCount_readyToRun = 0;
      transactionCount = 0;
      noOfStealsFromThisWorker = 0;
      noOfStealsByThisWorker = 0;
      stealStatisticsList.clear();
      OS.outputInfo("Status report: resetStatistics ----------");
      UnlockMutex();
    }


    static void print_stats(std::stringstream & ss, 
			    cht::simple_statistics & stats,
			    const char* name) {
      const int w = 13;
      const int p = 3;
      ss << name << ":\n"
	 << name << " count    = " << stats.count << "\n"
	 << name << " time_min = " 
	 << std::fixed << std::setw(w) << std::setprecision(6) << stats.time_min << "\n"
	 << name << " time_max = " 
	 << std::fixed << std::setw(w) << std::setprecision(6) << stats.time_max << "\n";
      if(stats.count > 0)
	ss << name << " time avg = " 
	   << std::fixed << std::setw(w) << std::setprecision(6) 
	   << stats.time_tot / stats.count << "\n";
    }


    void Worker::reportStatistics(std::string messageHeader) {
      const int w = 13;
      const int p = 3;
      std::stringstream ss;
      ss << messageHeader << "\n";
      LockMutex_macro;
      ss << "maxNoOfSimultaneousTasks = " << maxNoOfSimultaneousTasks << "\n";
      double readyToRunPercentage = 100.0 * (double)pendingTaskCount_readyToRun / pendingTaskCount_tot;
      ss << "pendingTaskCount_tot = " << pendingTaskCount_tot << ", pendingTaskCount_readyToRun = " << pendingTaskCount_readyToRun << " ==> " << readyToRunPercentage << "% ready-to-run\n";

      double sum_alltasks_wall = 0;
      double sum_alltasks_wall_work = 0;
      for(TaskStatsMap::iterator it = taskStatistics.begin(); it != taskStatistics.end(); it++) {
	if( ((*it).second.count_execute_leaf + (*it).second.count_execute_nonleaf) != (*it).second.count_finalize ) {
	  // This is not necessarily an error, may happen if task was stolen while execute was in progress.
	  //	  throw std::runtime_error("Error in cht::TaskSchedulerService:Worker::reportStatistics(): count mismatch.");
	}
	std::string taskTypeID = (*it).first;
	ss << "Task type '" << taskTypeID 
	   << "', count_leaf = " << (*it).second.count_execute_leaf 
	   << ", count_nonleaf = " << (*it).second.count_execute_nonleaf 
	   << "\n";
	ss << "                                             leaf         non-leaf              \n";
	ss << "                                             execute()    execute()   finalize()\n";
	for (int i = 0; i < Threads::ThreadGroup::numberOfPossibleStates; i++) {
	  ss << Threads::ThreadGroup::get_thread_state_name(i) << " (wall time) :"
	     << std::fixed << std::setw(w) << std::setprecision(p) 
	     << (*it).second.stats_execute_leaf[i].time_acc 
	     << std::fixed << std::setw(w) << std::setprecision(p) 
	     << (*it).second.stats_execute_nonleaf[i].time_acc 
	     << std::fixed << std::setw(w) << std::setprecision(p) 
	     << (*it).second.stats_finalize[i].time_acc
	     << "\n";
	  ss << Threads::ThreadGroup::get_thread_state_name(i) << " (  count  ) :"
	     << std::fixed << std::setw(w) << std::setprecision(p) 
	     << (*it).second.stats_execute_leaf[i].counter 
	     << std::fixed << std::setw(w) << std::setprecision(p) 
	     << (*it).second.stats_execute_nonleaf[i].counter 
	     << std::fixed << std::setw(w) << std::setprecision(p) 
	     << (*it).second.stats_finalize[i].counter
	     << "\n";
	  sum_alltasks_wall     += (*it).second.stats_execute_leaf[i].time_acc;
	  sum_alltasks_wall     += (*it).second.stats_execute_nonleaf[i].time_acc;
	  sum_alltasks_wall     += (*it).second.stats_finalize[i].time_acc;
	}
	ss << "execute_leaf_working_time_min = " 
	   << std::fixed << std::setw(w) << std::setprecision(p) 
	   << (*it).second.execute_leaf_working_time_min << "\n"
	   << "execute_leaf_working_time_max = " 
	   << std::fixed << std::setw(w) << std::setprecision(p) 
	   << (*it).second.execute_leaf_working_time_max << "\n";
	if((*it).second.count_execute_leaf > 0) {
	  double execute_leaf_working_time_average = 
	    (*it).second.stats_execute_leaf[Threads::ThreadGroup::working].time_acc / 
	    (*it).second.count_execute_leaf;
	  ss << "                      average = " 
	     << std::fixed << std::setw(w) << std::setprecision(p) 
	     << execute_leaf_working_time_average << "\n";
	}
	else
	  ss << "                      average = \n";
	sum_alltasks_wall_work     += (*it).second.stats_execute_leaf    [Threads::ThreadGroup::working].time_acc;
	sum_alltasks_wall_work     += (*it).second.stats_execute_nonleaf [Threads::ThreadGroup::working].time_acc;
	sum_alltasks_wall_work     += (*it).second.stats_finalize        [Threads::ThreadGroup::working].time_acc;
      }
  
      double totalElapsedWallTime = timerForEverything.get_elapsed_wall_seconds();

      ss << "Threads: "
	 << "  n_worker_threads = " << std::setw(2) << n_worker_threads
	 << "\n";  
      ss << "Tot time in execute() and finalize() (including blocking comm operations):" 
	 << std::fixed << std::setw(w) << std::setprecision(p) << sum_alltasks_wall
	 << "\n";
      ss << "Tot time working in execute() and finalize():" << std::fixed << std::setw(w) << std::setprecision(p) << sum_alltasks_wall_work 
	 << " --> worker active percentage: " 
	 << std::fixed << std::setw(6) << std::setprecision(2)
	 << ((sum_alltasks_wall_work / n_worker_threads) / totalElapsedWallTime)*100 << "%\n";
      ss << "Accumulated idle time                       :" << std::fixed << std::setw(w) << std::setprecision(p) << accumulatedIdleTime 
	 << " --> worker idle percentage  : " 
	 << std::fixed << std::setw(6) << std::setprecision(2)
	 << ((accumulatedIdleTime / n_worker_threads) / totalElapsedWallTime)*100 << "%\n";

      ss << "Total time used:\n";
      ss << "   wall       :      "
	 << std::fixed << std::setw(w) << std::setprecision(p) 
	 << timerForEverything.get_elapsed_wall_seconds() << "\n";
      ss << "   cpu_user   :      "
	 << std::fixed << std::setw(w) << std::setprecision(p) 
	 << timerForEverything.get_elapsed_cpu_seconds_user() << "\n";
      ss << "   cpu_system :      "
	 << std::fixed << std::setw(w) << std::setprecision(p) 
	 << timerForEverything.get_elapsed_cpu_seconds_system() << "\n";

      print_stats(ss, stealAttemptStatistics, "stealAttemptStatistics");
      print_stats(ss, stealAttemptSendStatistics, "stealAttemptSendStatistics");
      print_stats(ss, RequestHandlerFuncStatisticsParent, "RequestHandlerFuncStatisticsParent");
      print_stats(ss, RequestHandlerFuncStatisticsWorker, "RequestHandlerFuncStatisticsWorker");
      ss << "n_prefetch_tasks : " << n_prefetch_tasks << "\n";
      ss << "child_chunk_fetch_level : " << child_chunk_fetch_level << "\n";

      ss << "noOfStealsFromThisWorker = " << noOfStealsFromThisWorker << "\n";
      ss << "noOfStealsByThisWorker   = " << noOfStealsByThisWorker << "\n";
      int stealNumber = 0;
      for(std::list<StealInfoStruct>::const_iterator it = stealStatisticsList.begin(); it != stealStatisticsList.end(); it++) {
	stealNumber++;
	ss << "Steal " << stealNumber << " : victim=" << it->victimRank 
	   << ", transactionCount=" << it->transactionCountWhenStealOccurred 
	   << ", ancestorDepth=" << it->ancestorDepthOfStolenTask 
	   << ", wallTime=" << it->wallTimeWhenStealOccurred 
	   << ", onlyOwn=" << it->stealOnlyOwnStuff
	   << ", taskType" << it->taskTypeString << "\n";
      }
      ss << "Final transactionCount = " << transactionCount << "\n";

      OS.outputInfo(ss.str());
      OS.outputInfo("Status report: reportStatistics ----------");
  
      UnlockMutex();

      // FIXME: OUTPUT THREAD STATISTICS HERE

    }      




    void Worker::updateStatistics(const TaskInfoStruct & taskInfo,
						       cht::vector<Threads::ThreadGroup::State_stats> const & threadStats_start,
						       cht::vector<Threads::ThreadGroup::State_stats> const & threadStats_stop,
						       std::string execute_or_finalize_type) {
      std::string taskTypeString(taskInfo.get_taskTypeID());
      TaskStatsMap::iterator it = taskStatistics.find(taskTypeString);
      if(it == taskStatistics.end()) {
	TaskTypeStatistics stats;
	taskStatistics.insert(TaskStatsMap::value_type(taskTypeString, stats));
	it = taskStatistics.find(taskTypeString);
      }
      TaskTypeStatistics & stats = (*it).second;
      int* counterPtr = NULL; // to be set below
      Threads::ThreadGroup::State_stats* statsListPtr = NULL; // to be set below
      if(execute_or_finalize_type == "execute_leaf") {
	counterPtr = &stats.count_execute_leaf;
	statsListPtr = stats.stats_execute_leaf;
	double timeTakenInWorkingState = 
	  threadStats_stop[Threads::ThreadGroup::working].time_acc - 
	  threadStats_start[Threads::ThreadGroup::working].time_acc;
	if(*counterPtr == 0) {
	  stats.execute_leaf_working_time_min = timeTakenInWorkingState;
	  stats.execute_leaf_working_time_max = timeTakenInWorkingState;
	}
	else {
	  stats.execute_leaf_working_time_min = (stats.execute_leaf_working_time_min < timeTakenInWorkingState) ? 
	    stats.execute_leaf_working_time_min : timeTakenInWorkingState;
	  stats.execute_leaf_working_time_max = (stats.execute_leaf_working_time_max > timeTakenInWorkingState) ? 
	    stats.execute_leaf_working_time_max : timeTakenInWorkingState;
	}
      }
      else if(execute_or_finalize_type == "execute_nonleaf") {
	counterPtr = &stats.count_execute_nonleaf;
	statsListPtr = stats.stats_execute_nonleaf;
      }
      else if(execute_or_finalize_type == "finalize") {
	counterPtr = &stats.count_finalize;
	statsListPtr = stats.stats_finalize;
      }
      else 
	throw std::runtime_error("Error in cht::TaskSchedulerService:Worker::updateStatistics().");
      (*counterPtr)++;
      for (int i = 0; i < Threads::ThreadGroup::numberOfPossibleStates; ++i) {
	statsListPtr[i].counter  += (threadStats_stop[i].counter - threadStats_start[i].counter);
	statsListPtr[i].time_acc += (threadStats_stop[i].time_acc - threadStats_start[i].time_acc);
      }
    }


    void Worker::outputListContents(std::ostream & stream, const char* s, Worker::TaskListMap & theList) {
      stream << "outputListContents for list " << s << ", size = " << theList.size() << " : ";
      TaskListMap::iterator it;
      for ( it=theList.begin() ; it != theList.end(); it++ ) {
	TaskInfoStruct* taskInfo = (*it).second;
	stream << taskInfo->get_id().str() << " ";
      }
      stream << std::endl;
    }

    /* This function is good to have for debugging, if the CHT calculation hangs. */
    void Worker::writeStatusReportFile(int my_rank) {
      assert(mutexIsLocked);
      // Create filename based on MPI rank.
      std::stringstream ss;
      ss << "/scratch/cht_worker_status_report_rank_" << my_rank << ".txt";
      std::ofstream f;
      f.open(ss.str().c_str());
      f << "Status report for CHT worker with rank " << my_rank << "." << std::endl;
      outputListContents(f, "list_waiting_before_transaction", list_waiting_before_transaction);
      outputListContents(f, "helper_list_pending_for_transaction", helper_list_pending_for_transaction);
      outputListContents(f, "list_doing_transaction_nonleaf_1", list_doing_transaction_nonleaf_1);
      outputListContents(f, "list_doing_transaction_nonleaf_2", list_doing_transaction_nonleaf_2);
      outputListContents(f, "list_doing_transaction_leaf", list_doing_transaction_leaf);
      outputListContents(f, "list_waiting_before_finalize", list_waiting_before_finalize);
      outputListContents(f, "list_pending_for_finalize", list_pending_for_finalize);
      outputListContents(f, "list_running_finalize", list_running_finalize);
      outputListContents(f, "list_finished", list_finished);
      outputListContents(f, "list_stolen", list_stolen);
      f.close();
    }

    void Worker::doStatusReport() {
      assert(mutexIsLocked);
      // Check how many tasks are currently executing.
      int nExecuting = 0;
      int nExecuted = 0;
      int nReadyToRun = 0;
      TaskInfoStruct* taskInfo = NULL;
      TaskListSortedMap::const_iterator it = list_pending_for_transaction.begin();
      while(it != list_pending_for_transaction.end()) {
	taskInfo = it->second;
	if(taskInfo->executing)
	  nExecuting++;
	if(taskInfo->executed)
	  nExecuted++;
	if(taskInfo->ready_to_run)
	  nReadyToRun++;
	it++;
      }
      int nInPendingList = list_pending_for_transaction.size();
      int nInWaitingList = list_waiting_before_transaction.size();
      double currTime = timerForEverything.get_elapsed_wall_seconds();
      std::stringstream s;
      s << "Status report:" 
	<< " nInWaitingList=" << nInWaitingList
	<< " nInPendingList=" << nInPendingList
	<< " nExecuting=" << nExecuting
	<< " nExecuted=" << nExecuted
	<< " nReadyToRun=" << nReadyToRun
	<< " currTime=" << currTime;
      OS.outputInfo(s.str());
    }

    void Worker::sleepAndAddToAccumulatedIdleTime(int microsecondsToSleep) {
      double idleTimeStart = cht::get_wall_seconds();
      usleep(microsecondsToSleep);
      LockMutex_macro;
      accumulatedIdleTime += cht::get_wall_seconds() - idleTimeStart;
      UnlockMutex();      
    }

    void Worker::GetHelperRanks(const TaskListMap & taskListMap, std::list<int> & resultList) {
      TaskListMap::const_iterator it = taskListMap.begin();
      while(it != taskListMap.end()) {
	TaskInfoStruct* taskInfo = it->second;
	int currHelperRank = taskInfo->helperWorkerRank;
	// Check if currHelperRank is already present in resultList.
	std::list<int>::const_iterator it2 = resultList.begin();
	while(it2 != resultList.end()) {
	  if(*it2 == currHelperRank)
	    break;
	  it2++;
	}
	if(it2 == resultList.end())
	  resultList.push_back(currHelperRank);
	it++;
      }
    }


    /*
      WorkerFunc() description: 

      This routine is called by each of the worker nodes. Its main purpose
      is to perform the actual work by calling execute() and finalize()
      for different tasks. It also makes sure the different tasks are
      moved between the different lists when appropriate.
  
      At the same time as WorkerFunc() is running, RequestHandlerFunc() is
      also running but by another thread. Both of them access the
      different task lists, and they also communicate via certain flag
      variables. Conflicting access to such shared data is avoided using a
      mutex.

      Note that this function does NOT use any MPI receive calls; that is
      handled by the RequestHandlerFunc() thread.
    */
    void Worker::WorkerFunc()
    {
      AccessKey key(this);

      // It is assumed here that key outlives the following parameter copies
      int n_workers = key.n_workers();
      int my_rank = key.my_rank();
      MPI_Comm* comm_worker_world = key.comm_worker_world();
      MPI_Comm* comm_to_parent = key.comm_to_parent();

      // If this is worker with rank 0, we should receive the mother task
      // from the parent. However, that is done by the other thread. Here,
      // we do not care about that.

      double updateListsSecondsSum = 0;
      double mutexLockSecondsSum = 0;

      const int worker_threads_sleep_microseconds = 2000;

      while(1)
	{
	  double secondsLockMutexStart = cht::get_wall_seconds();
	  LockMutex_macro;
	  assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
	  mutexLockSecondsSum += cht::get_wall_seconds() - secondsLockMutexStart;

	  // Check if we are done.
	  if(finishFlagForThreadCommunicationWorker)
	    {
	      // Done!
	      UnlockMutex();
	      std::stringstream s1;
	      s1 << "cht::TaskSchedulerService:Worker::WorkerFunc() ending, \n"
		 << "updateListsSecondsSum = " << updateListsSecondsSum << "\n"
		 << "mutexLockSecondsSum   = " << mutexLockSecondsSum   << "\n";
	      OS.outputInfo(s1.str());
	      break;
	    }

	  // Check if any task list needs updating.
	  double secondsUpdateListsStart = cht::get_wall_seconds();
	  while(UpdateLists()) {} 
	  updateListsSecondsSum += cht::get_wall_seconds() - secondsUpdateListsStart;

	  // Update maxNoOfSimultaneousTasks is needed.
	  int totNoOfTasks = GetTotNoOfTasks();
	  maxNoOfSimultaneousTasks = (totNoOfTasks > maxNoOfSimultaneousTasks) ? 
	    totNoOfTasks : maxNoOfSimultaneousTasks;

	  // Check if any task can be moved from list_doing_transaction_nonleaf_1 to list_doing_transaction_nonleaf_2.
	  {
	    TaskListMap::iterator it;	
	    bool didSomething = false;
	    for ( it=list_doing_transaction_nonleaf_1.begin() ; it != list_doing_transaction_nonleaf_1.end(); it++ ) {
	      TaskInfoStruct* taskInfo = (*it).second;
	      assert(taskInfo->get_task() != 0);
	      if(taskInfo->chunkOpsStarted) {
		if(taskInfo->get_task()->internal_allRegisterAndCopyChunkOpsHaveFinished()) {
		  list_doing_transaction_nonleaf_1.erase(taskInfo->get_id());
		  list_doing_transaction_nonleaf_2.insert(TaskListMap::value_type(taskInfo->get_id(), taskInfo));
		  UnlockMutex();
		  taskInfo->get_task()->internal_registerTasksInList();
		  LockMutex_macro;
		  list_doing_transaction_nonleaf_2.erase(taskInfo->get_id());
		  list_waiting_before_finalize.insert(TaskListMap::value_type(taskInfo->get_id(), taskInfo));
		  UnlockMutex();
		  didSomething = true;
		  break;
		}
	      }
	    } // end for
	    if(didSomething)
	      continue;
	  }

	  // Check if there is any job pending for finalize.
	  if(list_pending_for_finalize.size() > 0)
	    {
	      // Take one job in list_pending_for_finalize and move it to 
	      // list_running_finalize and run finalize() for that task.
	      TaskListMap::iterator it = list_pending_for_finalize.begin();
	      TaskInfoStruct* taskInfo = (*it).second;
	      if(taskInfo->get_id() != (*it).first)
		throw std::runtime_error("Error!");
	      list_pending_for_finalize.erase(taskInfo->get_id());
	      list_running_finalize.insert(TaskListMap::value_type(taskInfo->get_id(), taskInfo));
	      UnlockMutex();
	      // Run finalize().
	      std::string taskTypeString(taskInfo->get_taskTypeID());
	      std::stringstream s;
	      s << "NOT calling finalize() for taskInfo->id = " << taskInfo->get_id().str() << " of task type '" << taskTypeString << "'";
	      OS.outputDebug(s.str());
	      double timer = cht::get_wall_seconds();
	      cht::vector<Threads::ThreadGroup::State_stats> threadStats_start;
	      Threads::getThreadStats(threadStats_start);
	      if (taskInfo->outputChunkOrTask_is_null())
		taskInfo->set_outputChunk(CHUNK_ID_NULL);
	      else if (taskInfo->get_outputChunkOrTask().is_chunkID())
		taskInfo->set_outputChunk(taskInfo->get_outputChunkOrTask().get_chunkID());
	      else {
		assert(taskInfo->get_outputChunkOrTask().is_taskID());
		// ELIAS NOTE 2012-02-04: It is important to lock mutex before calling getTaskResult() here,
		// since other threads may be in the process of modifying the task lists used by getTaskResult().
		LockMutex_macro;
		ChunkID resultChunkID = getTaskResult(taskInfo->get_outputChunkOrTask().get_taskID());
		UnlockMutex();
		taskInfo->set_outputChunk(resultChunkID);
	      }
	      cht::vector<Threads::ThreadGroup::State_stats> threadStats_stop;
	      Threads::getThreadStats(threadStats_stop);
	      LockMutex_macro;
	      updateStatistics(*taskInfo, threadStats_start, threadStats_stop, "finalize");
	      UnlockMutex();
	      std::stringstream s2;
	      s2 << "finalize() NOT completed for taskInfo->id = " << taskInfo->get_id().str() << " of task type '" << taskTypeString 
		 << "' wall seconds used: " << cht::get_wall_seconds() - timer;
	      OS.outputDebug(s2.str());
	      LockMutex_macro;
	      // Delete all temporary chunks owned by the task. This
	      // can be done now since all child tasks are finished.
	      // Since the deleteChunk call can take some time, we do not want to do that while the mutex is locked.
	      // So we create a list of ChunkIDs that should be deleted, then unlock the mutex and then call deleteChunk.
	      {
		std::list<ChunkID> chunks_to_delete;
		std::list<ID> const & tmp_chunks_to_be_deleted = taskInfo->getTemporaryChunks();
		std::list<ID>::const_iterator ch_it;
		for(ch_it = tmp_chunks_to_be_deleted.begin(); ch_it != tmp_chunks_to_be_deleted.end(); ch_it++) {
		  ChunkID cid_tmp;
		  if (ch_it->is_chunkID())
		    cid_tmp = ch_it->get_chunkID();
		  else {
		    assert(ch_it->is_taskID());
		    cid_tmp = getTaskResult(ch_it->get_taskID());
		  }
		  chunks_to_delete.push_back(cid_tmp);
		}
		taskInfo->clearTemporaryChunkList();
		// Unlock mutex before doing deleteChunk calls, since that may take some time.
		UnlockMutex();
		std::list<ChunkID>::const_iterator cid_it;
		for(cid_it = chunks_to_delete.begin(); cid_it != chunks_to_delete.end(); cid_it++)
		  ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().deleteChunk(*cid_it);
		LockMutex_macro;
	      }
	      // Check if task has children in list_finished and
	      // delete them.  This can be done since now there should
	      // not be any tasks left depending on the output of
	      // these children tasks.
	      {		
		std::list<TaskID> task_children;
		TaskListMap::const_iterator fi_it = list_finished.begin();
		while ( fi_it != list_finished.end() ) {
		  if ( fi_it->second->get_creator() == taskInfo->get_id() )
		    task_children.push_back(fi_it->first);
		  fi_it++;
		}
		std::list<TaskID>::iterator ch_it;
		for(ch_it = task_children.begin(); ch_it != task_children.end(); ch_it++) {
		  delete list_finished[*ch_it];
		  list_finished.erase(*ch_it);
		}
	      }
	      list_running_finalize.erase(taskInfo->get_id());
	      bool isMotherTask = (taskInfo->get_creator() == TASK_ID_NULL);
	      // Check if task was stolen, and delete taskInfo if that is the case.
	      if(taskInfo->get_stolen())
		{
		  // Send message to owner to report that the task is completed. 
		  std::stringstream s;
		  s.str("");
		  s << "Status report: returning stolen task " << taskInfo->get_id().str() << ".";
		  OS.outputInfo(s.str());
		  size_t msgSize = sizeof(TaskID) + sizeof(ChunkID);
		  std::vector<char> buffer(msgSize);
		  char * p = &buffer[0];
		  memcpy(p, &taskInfo->get_id(), sizeof(TaskID));
		  p += sizeof(TaskID);
		  memcpy(p, &taskInfo->get_outputChunk(), sizeof(ChunkID));
		  int tag = TAG_Stolen_task_finished;
		  MW._MPI_Send(&buffer[0], msgSize, MPI_UNSIGNED_CHAR, taskInfo->get_ownerWorkerRank(), tag, *comm_worker_world);
		  OS.outputDebug("deleting stolen taskinfo object...");
		  // Before deleting taskInfo object, check if this was the
		  // mother task. 
		  if(isMotherTask)
		    throw std::runtime_error("Error: mother task was stolen, should not happen.");
		  delete taskInfo;
		  taskInfo = NULL;
		  OS.outputDebug("deleted stolen taskinfo object after finalize() completed.");
		}
	      else {
		list_finished.insert(TaskListMap::value_type(taskInfo->get_id(), taskInfo));
		RemoveTaskFromWaitLists(taskInfo->get_id());
		RemoveTaskFromChildListOfCreator(taskInfo);
	      }
	      UnlockMutex();
	      if(isMotherTask)
		{
		  OS.outputInfo("Mother task finished!");
		  ChunkID cid_tmp = taskInfo->get_outputChunk();
		  MW._MPI_Send(&cid_tmp, sizeof(ChunkID), MPI_UNSIGNED_CHAR, 0, TAG_Mother_task_finished, *comm_to_parent);
		  LockMutex_macro;
		  list_finished.erase(taskInfo->get_id());
		  delete taskInfo;
		  thisWorkerHasMotherTask = false;
		  UnlockMutex();
		}
	      continue;
	    }

	  // Check if we can start the transaction for any task in the pending for transaction list.
	  assert(mutexIsLocked);
	  assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
          if(list_pending_for_transaction.size() > 0)
	    {
	      // Start from end of sorted list and try to do something about that task. Look for a task that is executed.
	      // We stop if we find a task that needs to be executed, since then that should be done before considering doing transaction for other tasks.
	      int noOfOngoingNonLeafTransactions = list_doing_transaction_nonleaf_1.size() + list_doing_transaction_nonleaf_2.size();
	      const int maxNoOfNonLeafTransactions = 1;
	      if(noOfOngoingNonLeafTransactions < maxNoOfNonLeafTransactions) {
		TaskInfoStruct* taskInfo = NULL;
		TaskListSortedMap::reverse_iterator it = list_pending_for_transaction.rbegin();
		while(it != list_pending_for_transaction.rend()) {
		  taskInfo = it->second;
		  if(taskInfo->executed) // in this case we should consider starting the transaction for this task.
		    break;
		  if(taskInfo->executed == false && taskInfo->executing == false) // in this case we should break since this task needs to be executed.
		    break;
		  it++;
		}
		if(it != list_pending_for_transaction.rend()) {
		  if(taskInfo->executed) {
		    // OK, now taskInfo points to a task that is executed. Check if we can perform the transaction for this task.
		    // To avoid unrolling parts of the task tree too early, we first check that any tasks that are currently executing have been executing for some time,
		    // since otherwise it may happen that some currently executing task could produce new tasks that should be processed before this task.
		    // Check how many tasks are currently executing. FIXME: use variable to keep track of this instead of going through list each time.
		    int nExecuting = 0;
		    TaskListSortedMap::const_iterator it = list_pending_for_transaction.begin();
		    while(it != list_pending_for_transaction.end()) {
		      if(it->second->executing)
			nExecuting++;
		      it++;
		    }
		    // Check how long time has passed since the last time we started execute() for some task.
		    double secondsElapsed = cht::get_wall_seconds() - timeWhenTaskWasLastSetAsExecuting;
		    if(nExecuting == 0 || secondsElapsed > 0.002) { // FIXME: DO NOT USE HARD-CODED VALUE HERE
		      // OK, proceed with the transaction.
		      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
		      RemoveTaskFromPendingForTransactionLists(taskInfo);
		      list_doing_transaction_nonleaf_1.insert(TaskListMap::value_type(taskInfo->get_id(), taskInfo));
		      transactionCount++;
		      UnlockMutex();
		      assert(taskInfo->get_task() != 0);
		      taskInfo->get_task()->internal_startRegisterAndCopyChunkOps();
		      LockMutex_macro;
		      taskInfo->chunkOpsStarted = true;
		      UnlockMutex();
		      continue;
		    }
		  }
		}
	      }
	    }

	  // Check if we can execute any task in the pending for transaction list.
	  assert(mutexIsLocked);
	  assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
          if(list_pending_for_transaction.size() > 0)
	    {
	      // Start from end of sorted list and try to do something about that task.
	      // Look for a task that is not executing and that has not been executed earlier and that is ready to run.
	      TaskInfoStruct* taskInfo = NULL;
	      TaskListSortedMap::reverse_iterator it = list_pending_for_transaction.rbegin();
	      while(it != list_pending_for_transaction.rend()) {
		taskInfo = it->second;
		if(taskInfo->executing == false && taskInfo->executed == false && taskInfo->ready_to_run)
		  break;
		it++;
	      }
	      if(it == list_pending_for_transaction.rend()) {
		// All tasks in list are already executing.
		UnlockMutex();
		sleepAndAddToAccumulatedIdleTime(worker_threads_sleep_microseconds);
		continue;
	      }
	      // At this point we know that taskInfo is a task that is not already executing.
	      // Set the "executing" flag so that other threads know they are not allowed to touch this task.
	      taskInfo->executing = true;
	      timeWhenTaskWasLastSetAsExecuting = cht::get_wall_seconds();
	      UnlockMutex();
	      // Run execute().
	      // NOTE:
	      // cht::obj_factory<Task>::instance().createObject() is
	      // creating a new task object using the new operator, so
	      // we must remember to call delete for the task pointer.
	      taskInfo->set_task( cht::obj_factory<Task>::instance().
				  createObject(taskInfo->get_taskTypeID()) );
	      std::stringstream s;
	      s << taskInfo->str() << " calling_execute";
	      OS.outputDebug(s.str());
	      double timer = cht::get_wall_seconds();
	      cht::vector<Threads::ThreadGroup::State_stats> threadStats_start;
	      Threads::getThreadStats(threadStats_start);
	      size_t nInpChunks = taskInfo->inputChunks.size();
	      assert(taskInfo->get_task() != 0);
	      taskInfo->get_task()->internal_setInputChunkIDs( taskInfo->inputChunks );
	      if (!taskInfo->fallback_execute_should_be_used) {
		taskInfo->get_task()->internal_setInputChunks(taskInfo->input_chunks);
		std::vector<BaseObj const *> inpChunksAndChunkIDs(nInpChunks);
		for(int ind = 0;ind < nInpChunks; ind ++)  {
		  if (cht::arg_manager<Task>::instance().check(taskInfo->get_taskTypeID(), "ChunkID", ind))
		    inpChunksAndChunkIDs[ind] = &taskInfo->inputChunks[ind];
		  else
		    inpChunksAndChunkIDs[ind] = taskInfo->input_chunks[ind].get();
		}
		taskInfo->get_task()->internal_setinputChunksAndChunkIDs( inpChunksAndChunkIDs );
	      }
	      ID taskResult;
	      try {
		taskResult = taskInfo->get_task()->internal_execute(taskInfo->get_id(), 
								    taskInfo->fallback_execute_should_be_used);
	      }
	      catch (std::exception & e) {
		std::cerr << "Error: std::exception caught when attempting to execute task of type '" << taskInfo->get_taskTypeID() << "'." << std::endl << "what() = '" << e.what() << "'" << std::endl;
		throw e;
	      }
	      // OK, now we have taskResult. We want to check that it seems reasonable.
	      if(taskResult.is_chunkID()) {
		// Check that the chunk type is correct.
		// We can only check output type if the declared output type is not "ChunkID".
		if(cht::arg_manager<Task>::instance().checkOutputType(taskInfo->get_taskTypeID(), "ChunkID") == false) {
		  std::string returnedChunkType = ChunkObjService::Worker::instance().getChunkTypeIDStr(taskResult.get_chunkID().chunkTypeID);
		  if(cht::arg_manager<Task>::instance().checkOutputType(taskInfo->get_taskTypeID(), returnedChunkType) == false) {
		    std::cerr << "Wrong output chunk type returned by task type " << taskInfo->get_taskTypeID() << "." << std::endl;
		    throw std::runtime_error("Wrong output chunk type.");
		  }
		}
	      }
	      else if(taskResult.is_taskID()) {
		// FIXME: Check that the output chunk type of the taskResult task type is correct.
	      }
	      LockMutex_macro;
	      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
	      taskInfo->get_task()->internal_get_child_task_ids(taskInfo->childList);
	      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
	      UnlockMutex();
	      // inpChunksAndChunkIDs is only used temporarily, we
	      // clear it immediately after execute (together with the
	      // other input vectors) :
	      taskInfo->get_task()->internal_clearInputVectors();
	      // Check that taskResult is persistent
	      {
		std::list<ID> const & tmp_chunks = taskInfo->getTemporaryChunks();
		std::list<ID>::const_iterator tmp_ch_it = tmp_chunks.begin();
		while (tmp_ch_it != tmp_chunks.end()) {
		  if (*tmp_ch_it == taskResult)
		    throw std::runtime_error("Output from task must be persistent!");
		  ++tmp_ch_it;
		}
	      }
	      taskInfo->set_outputChunkOrTask(taskResult);
	      cht::vector<Threads::ThreadGroup::State_stats> threadStats_stop;
	      Threads::getThreadStats(threadStats_stop);
	      // Now we need to check if the task has been stolen while this thread was doing execute for the task.
	      LockMutex_macro;
	      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
	      if(ExistsInMap(list_stolen, taskInfo->get_id())) {
		// Task was stolen. In this case we do nothing with the results of execute, the thief worker will take care of it.
		taskInfo->executing = false;
		taskInfo->executed = false;
		UnlockMutex();
		continue;
	      }
	      bool was_leaf_task = false;
	      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
	      assert(taskInfo->get_task() != 0);
	      if(taskInfo->get_task()->internal_getNoOfTasksToRegister() == 0)
		was_leaf_task= true;
	      if(was_leaf_task) {
		/* The task we just executed was a leaf task. */
		assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
		updateStatistics(*taskInfo, threadStats_start, threadStats_stop, "execute_leaf");
		assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
		RemoveTaskFromPendingForTransactionLists(taskInfo);
		list_doing_transaction_leaf.insert(TaskListMap::value_type(taskInfo->get_id(), taskInfo));
		transactionCount++;
		UnlockMutex();
		// Start any createChunk and copyChunk operations and
		// then place the task in
		// list_waiting_before_finalize.
		assert(taskInfo->get_task() != 0);
		taskInfo->get_task()->internal_startRegisterAndCopyChunkOps();
		LockMutex_macro;
		list_doing_transaction_leaf.erase(taskInfo->get_id());
		list_waiting_before_finalize.insert(TaskListMap::value_type(taskInfo->get_id(), taskInfo));
		taskInfo->executing = false;
		taskInfo->executed = true;
		UnlockMutex();
		continue;
	      }
	      else {
		/* The task we just executed was not a leaf task. */
		assert(taskInfo->get_task() != 0);
		taskInfo->executing = false;
		taskInfo->executed = true;
		assert(taskInfo->get_task() != 0);
		updateStatistics(*taskInfo, threadStats_start, threadStats_stop, "execute_nonleaf");
		assert(taskInfo->get_task() != 0);
		UnlockMutex();
		continue;
	      }
	    }

	  // Check if it is time to make a steal attempt.
	  int noOfOngoingNonLeafTransactions = list_doing_transaction_nonleaf_1.size() + list_doing_transaction_nonleaf_2.size();
	  int noOfOngoingLeafTransactions = list_doing_transaction_leaf.size();
	  int nStolenTasks = list_stolen.size();
	  if(   !steal_attempt_in_progress 
		&& !stealing_disabled_flag
		&& n_workers > 1
		&& noOfOngoingNonLeafTransactions == 0
		&& noOfOngoingLeafTransactions == 0
		&& list_pending_for_transaction.size() == 0 )
	    {
	      // Check if there are any tasks in the wait list (if so, that means this worker is stuck because some other worker stole part of the work from this worker).
	      int rank_to_steal_from = -1;
	      bool stealOnlyFromWithinOwnTasks = false;
	      if(nStolenTasks > 0) {
		// OK, there is at least one stolen task. Now create a corresponding list of worker ranks.
		stealOnlyFromWithinOwnTasks = true;
		std::list<int> helperWorkerRanks;
		GetHelperRanks(list_stolen, helperWorkerRanks);
		int nHelpers = helperWorkerRanks.size();
		int randInt = (int)((nHelpers) * ( (double)rand() / RAND_MAX ));
		if(!(randInt >= 0 && randInt < nHelpers))
		  std::cerr << "nHelpers = " << nHelpers << "   randInt = " << randInt << std::endl;
		assert(randInt >= 0 && randInt < nHelpers);
		int count = 0;
		std::list<int>::const_iterator it = helperWorkerRanks.begin();
		rank_to_steal_from = -1;
		while(it != helperWorkerRanks.end()) {
		  if(count == randInt)
		    rank_to_steal_from = *it;
		  count++;
		  it++;
		}
	      } // end if there is some stolen task
	      else {
		// No stolen task. In this case we should steal from a random worker.
		int randInt = (int)((n_workers-1) * ( (double)rand() / RAND_MAX ));
		rank_to_steal_from = randInt;
		if(rank_to_steal_from >= my_rank)
		  rank_to_steal_from++;
	      }
	      if(rank_to_steal_from < 0 || rank_to_steal_from >= n_workers || rank_to_steal_from == my_rank)
		throw std::runtime_error("Error: (rank_to_steal_from < 0 || rank_to_steal_from >= n_workers || rank_to_steal_from == my_rank).");
	      std::stringstream s;
	      s << "Sending task steal attempt message, rank_to_steal_from = " << rank_to_steal_from;
	      OS.outputDebug(s.str());
	      steal_attempt_start_time = cht::get_wall_seconds();
	      double sendStartTime = cht::get_wall_seconds();
	      MW._MPI_Send(&stealOnlyFromWithinOwnTasks, sizeof(bool), MPI_UNSIGNED_CHAR, 
			  rank_to_steal_from, TAG_Steal_task_attempt, *comm_worker_world);
	      stealAttemptSendStatistics.add_time(cht::get_wall_seconds() - sendStartTime);
	      // Set flag to make sure no more steal attempts are
	      // initiates until this one is complete.
	      steal_attempt_in_progress = true;
	      currStealAttemptIsOnlyForOwnStuff = stealOnlyFromWithinOwnTasks;
	      UnlockMutex();
	      continue;
	    } // end if consider making a steal attempt.

	  // If we get to this point in the loop it means there was
	  // nothing at all to do. In this case we wait a little while.
	  UnlockMutex();
	  sleepAndAddToAccumulatedIdleTime(worker_threads_sleep_microseconds);
	}  // end while main worker loop
    } // end Worker::WorkerFunc()


    const int TASK_MESSAGE_MAX_SIZE = 1000;

    void Worker::HandleMessageFromParent(int tag, int msgSize, char* recvBufParent, MPI_Comm* comm_to_parent, MPI_Request & request) {
      if(tag == TAG_Mother_task_info) {
	// OK, we have a TAG_Mother_task_info message from parent.
	// Determine message length.
	if(msgSize > TASK_MESSAGE_MAX_SIZE)
	  throw std::runtime_error("Error: (msgSize > TASK_MESSAGE_MAX_SIZE).");
	OS.outputInfo("Received task info from parent!");
	// OK, now we have the task info in message buffer. We need to
	// create a corresponding taskInfo object and add it to the
	// proper list.
	TaskInfoStruct* taskInfo = new TaskInfoStruct;
	taskInfo->unpack(recvBufParent, msgSize);
	LockMutex_macro;
	assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
	AddTaskToPendingForTransactionLists(taskInfo);
	std::string taskTypeString(taskInfo->get_taskTypeID());
	thisWorkerHasMotherTask = true;
	std::stringstream s1;
	s1 << taskInfo->str() << " received_from_parent taskType = '" << taskTypeString << "'";
	std::stringstream s2;
	s2 << "Status report: received mother task " << taskInfo->get_id().str() << " now! totInputDataSize(MB) = " 
	   << taskInfo->totInputDataSize/1000000 << ", currWorkerOwnedSize(MB) = " << taskInfo->currWorkerOwnedSize/1000000 << ", taskType = " << taskInfo->get_taskTypeID();
	UnlockMutex();
	// ELIAS NOTE 2013-10-24: At this point, the mutex has been
	// unlocked, so we must not access the taskInfo pointer
	// anymore. If we access the taskInfo pointer here, it may
	// happen that another thread has already processed the task
	// and deleted the object, then we could get a crash. We did
	// get segmentation fault crashes due to this when running the
	// cht_scf_test test program on Tintin. The problem was fixed
	// by doing all accesses to the taskInfo pointer above, before
	// the call to UnlockMutex().
	OS.outputInfo(s1.str());
	OS.outputInfo(s2.str());
	MW._MPI_Irecv(recvBufParent, TASK_MESSAGE_MAX_SIZE, MPI_UNSIGNED_CHAR, 
		     MPI_ANY_SOURCE, MPI_ANY_TAG, *comm_to_parent, &request);
      }
      else {
	// This point should never be reached.
	throw std::runtime_error("Error: unknown tag received from parent.");
      }
    }


    void Worker::HandleStealAttempt(int rank_who_wants_to_steal, bool stealOnlyFromWithinOwnTasks, MPI_Comm* comm_worker_world) {
      // Check if there is any task to steal.
      TaskInfoStruct* taskInfoToSteal = NULL;
      LockMutex_macro;
      // We only allow stealing from list_pending_for_execute.
      assert(  list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
      if(list_pending_for_transaction.size() > 1)
	{
	  // OK, there is maybe something to steal, at least
	  // the pending list is not empty.  We want to steal
	  // from the top of the list but only if the task is
	  // not stolen, and is not the mother task.
	  TaskInfoStruct* taskInfoTmp = list_pending_for_transaction.begin()->second; 
	  assert(taskInfoTmp != NULL);
	  if(stealOnlyFromWithinOwnTasks) {
	    // Look for a task with a top ancestor that matches the requestor rank.
	    TaskListSortedMap::const_iterator it = list_pending_for_transaction.begin();
	    while(it != list_pending_for_transaction.end()) {
	      taskInfoTmp = it->second;
	      assert(taskInfoTmp != NULL);
	      TaskInfoPtr taskInfoPtr(*taskInfoTmp);
	      std::list<TaskInfoStruct const *> ancestorList;
	      taskInfoPtr.get_ancestor_list(*taskInfoTmp, ancestorList);
	      TaskInfoStruct const * topAncestor = ancestorList.front();
	      TaskID topAncestorID = topAncestor->get_id();
	      if(topAncestorID.ownerRank == rank_who_wants_to_steal)
		break;
	      it++;
	    }
	    if(it == list_pending_for_transaction.end()) {
	      // No task with suitable ancestor found.
	      taskInfoTmp = NULL;
	    }
	  }
	  if(taskInfoTmp != NULL) {
	    // Check that this task was not stolen from another worker and that it is not the mother task
	    if(taskInfoTmp->get_stolen() == false && taskInfoTmp->get_creator() != TASK_ID_NULL && taskInfoTmp->executing == false) 
	      taskInfoToSteal = taskInfoTmp;
	  }
	}
      if(taskInfoToSteal == NULL)
	{
	  // Nothing to steal.
	  std::stringstream s;
	  s << "Denying steal attempt, rank_who_wants_to_steal = " << rank_who_wants_to_steal;
	  OS.outputDebug(s.str());
	  // Send TAG_Steal_failed message back.
	  MW._MPI_Send(0, 0, MPI_UNSIGNED_CHAR, rank_who_wants_to_steal, 
		      TAG_Steal_failed, *comm_worker_world);
	  UnlockMutex();
	}
      else {
	// OK, there is something to steal!
	bool doSteal = true;
	// Now check if the thief worker requested to steal only from within its own tasks.
	if(stealOnlyFromWithinOwnTasks) {
	  TaskInfoPtr taskInfoPtr(*taskInfoToSteal);
	  std::list<TaskInfoStruct const *> ancestorList;
	  taskInfoPtr.get_ancestor_list(*taskInfoToSteal, ancestorList);
	  TaskInfoStruct const * topAncestor = ancestorList.front();
	  TaskID topAncestorID = topAncestor->get_id();
	  if(topAncestorID.ownerRank != rank_who_wants_to_steal) {
	    // Nothing to steal.
	    std::stringstream s;
	    s << "Denying steal attempt because of wrong topAncestorID.ownerRank, rank_who_wants_to_steal = " << rank_who_wants_to_steal;
	    OS.outputDebug(s.str());
	    // Send TAG_Steal_failed message back.
	    MW._MPI_Send(0, 0, MPI_UNSIGNED_CHAR, rank_who_wants_to_steal, 
			TAG_Steal_failed, *comm_worker_world);
	    doSteal = false;
	    UnlockMutex();
	  }
	}
	if(doSteal) {
	  std::stringstream s;
	  s << taskInfoToSteal->str() << " sending_task_to_thief "
	    << "rank_who_wants_to_steal = " << rank_who_wants_to_steal;
	  OS.outputInfo(s.str());
	  // Remove a task from the list, pack it and send it.
	  // Now we have the taskInfo object. Remove it from list.
	  assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
	  RemoveTaskFromPendingForTransactionLists(taskInfoToSteal);
	  // Mark it as stolen.
	  taskInfoToSteal->set_stolen(true);
	  taskInfoToSteal->helperWorkerRank = rank_who_wants_to_steal;
	  // Add it to list_stolen where it will stay until we get a
	  // message saying the task has been completed.
	  list_stolen.insert(TaskListMap::value_type(taskInfoToSteal->get_id(), taskInfoToSteal));
	  // Pack and send task.
	  OS.outputInfo("Packing and sending stolen task.");
	  int buffer_length_needed = taskInfoToSteal->pack_size();
	  std::vector<char> buf(buffer_length_needed);
	  taskInfoToSteal->pack(&buf[0]);
	  MW._MPI_Send(&buf[0], buffer_length_needed, MPI_UNSIGNED_CHAR, 
		      rank_who_wants_to_steal, TAG_Stolen_task_info, *comm_worker_world);
	  OS.outputInfo("Stolen task sent!");
	  // Let data fetcher know that we have removed a task from the pending_for_execute list
	  data_fetcher_cond.signal();
	  noOfStealsFromThisWorker++;
	  UnlockMutex();
	} // end if doSteal
      } // end else there is something to steal
    }


    void Worker::HandleStolenTaskMessage(const char* recvBufWorker, int msgSize, int victimRank) {
      // Stolen task received!
      TaskInfoStruct* taskInfo = new TaskInfoStruct;
      taskInfo->unpack(recvBufWorker, msgSize);
      std::stringstream s;
      s << "Successfully stealing from rank " << victimRank;
      OS.outputInfo(s.str());
      s.str("");
      s << taskInfo->str() << " stolen_task_received "
	<< "ownerRank = " << taskInfo->get_ownerWorkerRank()
	<< " taskType = '" << taskInfo->get_taskTypeID() << "'";
      // Check that taskInfo->id.str() matches string from map
      if ( taskInfo->get_taskTypeID() != intIDToTaskTypeIDMap[taskInfo->get_id().taskTypeID])
	throw std::runtime_error("taskInfo->taskTypeID != intIDToTaskTypeIDMap[taskInfo->get_id().taskTypeID] in Worker::RequestHandlerFunc().");
      OS.outputInfo(s.str());
      LockMutex_macro;
      stealAttemptStatistics.add_time(cht::get_wall_seconds() - steal_attempt_start_time);
      steal_attempt_in_progress = false;
      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
      AddTaskToPendingForTransactionLists(taskInfo);
      noOfStealsByThisWorker++;
      StealInfoStruct stealInfo;
      stealInfo.victimRank = victimRank;
      stealInfo.transactionCountWhenStealOccurred = transactionCount;
      stealInfo.ancestorDepthOfStolenTask = taskInfo->get_ancestorDepth();
      stealInfo.wallTimeWhenStealOccurred = timerForEverything.get_elapsed_wall_seconds();
      stealInfo.stealOnlyOwnStuff = currStealAttemptIsOnlyForOwnStuff;
      stealInfo.taskTypeString = taskInfo->get_taskTypeID();
      stealStatisticsList.push_back(stealInfo);
      s.str("");
      s << "Status report: stealing task " << taskInfo->get_id().str() << " now!";
      OS.outputInfo(s.str());
      UnlockMutex();
    }
    

    void Worker::HandleStolenTaskFinishedMessage(const char* recvBufWorker, int msgSize) {
      int msgSizeExpected = sizeof(TaskID) + sizeof(ChunkID);
      if(msgSize != msgSizeExpected)
	throw std::runtime_error("Error in HandleStolenTaskFinishedMessage: (msgSize != msgSizeExpected).");
      TaskID stolenTaskID;
      const char * p = recvBufWorker;
      memcpy(&stolenTaskID, p, sizeof(TaskID));
      p += sizeof(TaskID);
      ChunkID cid_tmp;
      memcpy(&cid_tmp, p, sizeof(ChunkID));
      LockMutex_macro;
      list_stolen[stolenTaskID]->set_outputChunk(cid_tmp);
      // Great, a stolen task has been finished!
      std::stringstream s;
      s << list_stolen[stolenTaskID]->str() << " stolen_task_finished ";
      UnlockMutex();
      OS.outputInfo(s.str());
      // Now we would like to move this task from list_stolen to list_finished, but if task is executing we must wait for that to finish first.
      while(1) {
	bool breakNow = false;
	LockMutex_macro;
	if(list_stolen[stolenTaskID]->executing == false)
	  breakNow = true;
	UnlockMutex();
	if(breakNow)
	  break;
	usleep(2000);
      }
      LockMutex_macro;
      list_finished.insert(TaskListMap::value_type(stolenTaskID, 
						   list_stolen[stolenTaskID]));
      list_stolen.erase(stolenTaskID);
      RemoveTaskFromWaitLists(stolenTaskID);
      RemoveTaskFromChildListOfCreator(list_finished[stolenTaskID]);
      UnlockMutex();
      s.str("");
      s << "Stolen task with stolenTaskID = " << stolenTaskID.str() << " now removed from list_stolen.";
      OS.outputInfo(s.str());
    }


    /*
      RequestHandlerFunc() description: 

      This routine is called by each of the worker nodes. Its purpose is
      to receive any incoming messages, and act upon them by sometimes
      sending messages in return, and sometimes modifying the different
      task lists.
  
      The two functions WorkerFunc() and RequestHandlerFunc() are both
      running at the same time, but in different threads. Both of them
      access the different task lists, and they also communicate via
      certain flag variables. Conflicting access to such shared data is
      avoided using a mutex.

      Note that RequestHandlerFunc() handles all MPI receive calls; some
      of the received messages are the result of send operations issued by
      the WorkerFunc() thread.
    */
    void Worker::RequestHandlerFunc()
    {
      AccessKey key(this);
      // It is assumed here that key outlives the following parameter copies
      int n_workers = key.n_workers();
      int my_rank = key.my_rank();
      MPI_Comm* comm_worker_world = key.comm_worker_world();
      MPI_Comm* comm_to_parent = key.comm_to_parent();
  
      std::stringstream s;
      s << "Beginning of cht::TaskSchedulerService::Worker::RequestHandlerFunc, my_rank = " << my_rank;
      OS.outputInfo(s.str());

      const int nRequestsTot = 2;
      MPI_Request requestList[nRequestsTot]; // the order is: parent is [0], workers is [1] 
      char recvBufParent[TASK_MESSAGE_MAX_SIZE];
      char recvBufWorker[TASK_MESSAGE_MAX_SIZE];
      MW._MPI_Irecv(recvBufParent, TASK_MESSAGE_MAX_SIZE, MPI_UNSIGNED_CHAR, 
		   MPI_ANY_SOURCE, MPI_ANY_TAG, *comm_to_parent, &requestList[0]);
      MW._MPI_Irecv(recvBufWorker, TASK_MESSAGE_MAX_SIZE, MPI_UNSIGNED_CHAR, 
		   MPI_ANY_SOURCE, MPI_ANY_TAG, *comm_worker_world, &requestList[1]);
  
      while(1)
	{
	  // Wait for any message from parent or workers.
	  MPI_Status status;
	  int index = -1;
	  MW._MPI_Waitany(nRequestsTot, requestList, &index, &status);
	  // Now something has arrived!
	  // Check index to see if the arrived message is from parent or from another worker.
	  if(index == 0) {
	    // Message from parent.
	    double startTime = cht::get_wall_seconds();
	    if(status.MPI_SOURCE != 0)
	      throw std::runtime_error("Error: (status.MPI_SOURCE != 0).");
	    if(status.MPI_TAG == TAG_Workers_finished) {
	      // signal from parent that all other workers are ready
	      OS.outputInfo("cht::TaskSchedulerService:Worker::RequestHandlerFunc() breaking.");
	      MW._MPI_Irecv(recvBufParent, TASK_MESSAGE_MAX_SIZE, MPI_UNSIGNED_CHAR, 
			   MPI_ANY_SOURCE, MPI_ANY_TAG, *comm_to_parent, &requestList[0]);
	      break;
	    }
	    int msgSize;
	    MW._MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &msgSize);
	    HandleMessageFromParent(status.MPI_TAG, msgSize, recvBufParent, comm_to_parent, requestList[0]);
	    LockMutex_macro;
	    RequestHandlerFuncStatisticsParent.add_time(cht::get_wall_seconds() - startTime);
	    UnlockMutex();
	  } // end if message from parent
	  else {
	    // Message from another worker.
	    double startTime = cht::get_wall_seconds();
	    if(index != 1)
	      throw std::runtime_error("Error: (index != 1) for worker case.");
	    if(status.MPI_TAG == TAG_Steal_task_attempt) {
	      // Determine message length.
	      int msgSize;
	      MW._MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &msgSize);
	      if(msgSize != sizeof(bool))
		throw std::runtime_error("Error: (msgSize != sizeof(bool)) for TAG_Steal_task_attempt.");
	      bool stealOnlyFromWithinOwnTasks;
	      memcpy(&stealOnlyFromWithinOwnTasks, recvBufWorker, sizeof(bool));
	      // OK, we have a steal attempt.
	      int rank_who_wants_to_steal = status.MPI_SOURCE;
	      if(rank_who_wants_to_steal == my_rank)
		throw std::runtime_error("Error: (rank_who_wants_to_steal == my_rank).");
	      HandleStealAttempt(rank_who_wants_to_steal, stealOnlyFromWithinOwnTasks, comm_worker_world);
	    } // END IF TAG_Steal_task_attempt
	    else if(status.MPI_TAG == TAG_Steal_failed) {
	      // OK, we have a TAG_Steal_failed message.
	      LockMutex_macro;
	      stealAttemptStatistics.add_time(cht::get_wall_seconds() - steal_attempt_start_time);
	      steal_attempt_in_progress = false;
	      UnlockMutex();
	    }
	    else if(status.MPI_TAG == TAG_Stolen_task_info) {
	      // Determine message length.
	      int msgSize;
	      MW._MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &msgSize);
	      // OK, we have a TAG_Stolen_task_info message.
	      HandleStolenTaskMessage(recvBufWorker, msgSize, status.MPI_SOURCE);
	    }
            else if(status.MPI_TAG == TAG_Stolen_task_finished) {
              // OK, we have a TAG_Stolen_task_finished message.
	      int msgSize;
	      MW._MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &msgSize);
	      HandleStolenTaskFinishedMessage(recvBufWorker, msgSize);
            }
	    else {
	      // This point should never be reached.
	      throw std::runtime_error("Error: unknown tag received from another worker.");
	    }
	    MW._MPI_Irecv(recvBufWorker, TASK_MESSAGE_MAX_SIZE, MPI_UNSIGNED_CHAR, 
			 MPI_ANY_SOURCE, MPI_ANY_TAG, *comm_worker_world, &requestList[1]);
	    LockMutex_macro;
	    RequestHandlerFuncStatisticsWorker.add_time(cht::get_wall_seconds() - startTime);
	    UnlockMutex();
	  } // end if something received from another worker
	} // END WHILE
      // Cleanup.
      for(int i = 0; i < nRequestsTot; i++) {
	MW._MPI_Cancel(&requestList[i]);
	MW._MPI_Request_free(&requestList[i]);
      }
    }

    bool Worker::getChildChunksRecursive(const Chunk & chunk, int level, std::vector<ChunkID> & childChunkIDVec) {
      std::list<ChunkID> childChunkIDs;
      chunk.getChildChunks(childChunkIDs);
      std::list<ChunkID>::iterator iter = childChunkIDs.begin();
      while(iter != childChunkIDs.end()) {
	if(*iter != CHUNK_ID_NULL) {
	  if(level == 0)
	    childChunkIDVec.push_back(*iter);
	  else {
	    // Level is not zero. Get the chunk and do a recursive call to get child chunks for next level.
	    cht::shared_ptr<Chunk const> childChunk;
	    if(!ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().getChunkIfLocal(*iter, childChunk))
	      return false;
	    if(!getChildChunksRecursive(*childChunk, level-1, childChunkIDVec))
	      return false;
	  }
	}
	iter++;
      }
      return true;
    }

    bool Worker::getChildChunksForLevel(const TaskInfoStruct & taskInfo, int level, std::vector<ChunkID> & childChunkIDVec) {
      std::vector<ChunkID> input_chunk_ids( taskInfo.inputChunks );
      for(int i = 0; i < input_chunk_ids.size(); i++) {
	ChunkID currChunkID = input_chunk_ids[i];
	if(currChunkID != CHUNK_ID_NULL) {
	  cht::shared_ptr<Chunk const> chunk;
	  if(!ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().getChunkIfLocal(currChunkID, chunk)) {
	    std::cerr << "Error in getChildChunksForLevel: getChunkIfLocal failed for currChunkID = " << currChunkID.str() << std::endl;
	    return false;
	  }
	  if(!getChildChunksRecursive(*chunk, level, childChunkIDVec)) {
	    std::cerr << "Error in getChildChunksForLevel: getChildChunksRecursive failed." << std::endl;
	    return false;
	  }
	}
      }
      return true;
    }

    void Worker::DataFetcherFunc() {
      LockMutex_macro;
      /* The number of tasks to fetch data for is determined by the number of worker threads plus n_prefetch_tasks,
	 since even without "prefetch" we still need to get data for one task per worker thread, otherwise the threads cannot execute tasks simultaneously. */
      int n_tasks_to_fetch_data_for = n_worker_threads + n_prefetch_tasks;
      while(1)
	{
	  // Check if we are done.
	  if(finishFlagForThreadCommunicationWorker)
	    {
	      // Done!
	      UnlockMutex();
	      break;
	    }
	  // Choose task from end of sorted list 
	  TaskListSortedMap::reverse_iterator task_iter = list_pending_for_transaction.rbegin();
	  int number_of_checked_tasks = 0;
	  bool did_nothing = true; 
	  while (task_iter != list_pending_for_transaction.rend() && 
		 number_of_checked_tasks < n_tasks_to_fetch_data_for) {
	    TaskInfoStruct* taskInfo = task_iter->second;	  
	    int nInpChunks = taskInfo->inputChunks.size();
	    if ( !taskInfo->ready_to_run ) {
	      did_nothing = false;
	      // We know already that the fallback execute will not be used
	      std::vector<bool> input_type_is_ChunkID(nInpChunks);
	      for(int ind = 0;ind < nInpChunks; ind ++) 
		input_type_is_ChunkID[ind] = 
		  cht::arg_manager<Task>::instance().check(taskInfo->get_taskTypeID(), "ChunkID", ind);
	      std::vector<ChunkID> input_chunk_ids( taskInfo->inputChunks );
	      std::vector<cht::shared_ptr<Chunk const> > chunk_vec(taskInfo->input_chunks);
	      TaskID tid = taskInfo->get_id();
	      UnlockMutex();
	      for(int ind = 0; ind < nInpChunks; ind++) {
		if ( chunk_vec[ind] == 0 && !input_type_is_ChunkID[ind] ) {
		  if(ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().getChunkIfExists(input_chunk_ids[ind], chunk_vec[ind]) == false) {
		    // getChunkIfExists returned false. This means the chunk has been deleted. 
		    // We accept this only if the task has been moved away from list_pending_for_transaction, then it could have been stolen which could lead to the chunk being deleted.
		    LockMutex_macro;
		    if( helper_list_pending_for_transaction.find(tid) != helper_list_pending_for_transaction.end() ) {
		      std::cerr << "Error in DataFetcherFunc: chunk does not exist for task still present in pending_for_transaction list. tid = " 
				<< tid.str() << ", input_chunk_ids[ind].str() = " << input_chunk_ids[ind].str() << std::endl;
		      throw std::runtime_error("Error in DataFetcherFunc: chunk does not exist for task still present in pending_for_transaction list.");
		    }
		    UnlockMutex();
		    // Apparently the task has been stolen, so there is no point in fetching any more data for it.
		    break;
		  }
		}
	      } // end for
	      LockMutex_macro;
	      // Check if task is still in list_pending_for_transaction.
	      // If it is not we do nothing with the chunks we just fetched.
	      TaskListMap::iterator it = helper_list_pending_for_transaction.find(tid);
	      if ( it != helper_list_pending_for_transaction.end() ) {
		// Task is still in list_pending_for_transaction!
		taskInfo = it->second;
		taskInfo->input_chunks = chunk_vec;
		taskInfo->ready_to_run = true;
	      }
	      break;
	    } // end if
	    task_iter++;
	    number_of_checked_tasks++;
	  } // end while
	  if ( did_nothing ) {
	    // At this point we know we have gone through all relevant tasks in list_pending_for_transaction but all of them had the data they need to run.
	    // In this case, we start fetching child chunks.
	    // We do this one child chunk at a time; after fetching one child chunk we go back to see if any more important chunks need to be fetched.
	    // The child_chunk_fetch_level value is used below.
	    // The meaning of child_chunk_fetch_level is as follows:
	    // -1 : child chunk fetching disabled (should give same behavior as before child chunk fetching was implemented)
	    //  0 : fetch chunks for ChunkID-type task input params, but still not for any "real" child chunks
	    //  1 : fetch child chunks (but no grand-children)
	    //  2 : fetch child chunks (including grand-children)
	    // ... and so on ...
	    //  n : fetch child chunks down to level n
	    //
	    if(child_chunk_fetch_level >= 0) {
	      // Fetch any chunks for ChunkID-type task input params
	      TaskListSortedMap::reverse_iterator task_iter = list_pending_for_transaction.rbegin();
	      int number_of_checked_tasks = 0;
	      while (task_iter != list_pending_for_transaction.rend() && 
		     number_of_checked_tasks < n_tasks_to_fetch_data_for) {
		TaskInfoStruct* taskInfo = task_iter->second;
		// First of all, make sure this TaskInfoStruct object has the correct size of the child_chunks and child_chunks_done_flags vectors.
		if(taskInfo->child_chunks.size() != child_chunk_fetch_level+1) {
		  taskInfo->child_chunks.resize(child_chunk_fetch_level+1);
		  taskInfo->child_chunks_done_flags.resize(child_chunk_fetch_level+1);
		  for(int k = 0; k < child_chunk_fetch_level+1; k++)
		    taskInfo->child_chunks_done_flags[k] = false;
		}
		// Check if the ChunkID-type task input params have already been taken care of for this task.
		if(taskInfo->child_chunks_done_flags[0]) {
		  task_iter++;
		  number_of_checked_tasks++;
		  continue;
		}
		int nInpChunks = taskInfo->inputChunks.size();
		if ( !taskInfo->ready_to_run ) {
		  std::cerr << "Error in DataFetcherFunc while trying to fetch child chunks: ( !taskInfo->ready_to_run )." << std::endl;
		  throw std::runtime_error("Error in DataFetcherFunc while trying to fetch child chunks: ( !taskInfo->ready_to_run ).");
		}
		std::vector<bool> input_type_is_ChunkID(nInpChunks);
		for(int ind = 0;ind < nInpChunks; ind ++) 
		  input_type_is_ChunkID[ind] = 
		    cht::arg_manager<Task>::instance().check(taskInfo->get_taskTypeID(), "ChunkID", ind);
		int count = 0;
		for(int ind = 0;ind < nInpChunks; ind ++) {
		  if(input_type_is_ChunkID[ind]) {
		    // OK, we found a ChunkID-type input chunk. Make sure it is fetched.
		    if(taskInfo->child_chunks[0].size() > count) {
		      // Chunk already fetched; just continue in this case.
		      count++;
		      continue;
		    }
		    // Chunk not already fetched. Consider fetching it now.
		    ChunkID ChunkID_to_fetch = taskInfo->inputChunks[ind];
		    if(ChunkID_to_fetch == CHUNK_ID_NULL)
		      continue;
		    did_nothing = false;
		    TaskID tid = taskInfo->get_id();
		    UnlockMutex();
		    cht::shared_ptr<Chunk const> fetchedChunk;
		    if(ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().getChunkIfExists(ChunkID_to_fetch, fetchedChunk) == false) {
		      // getChunkIfExists returned false. This means the chunk has been deleted. This can happen if the task was run while we were trying to fetch the chunk.
		      // FIXME: is this the right thing to do here?
		      LockMutex_macro;
		      break;
		    }
		    // OK, we have fetched a chunk. Add it to the list if possible.
		    LockMutex_macro;
		    // Check if task is still in list_pending_for_transaction.
		    // If it is not we do nothing with the chunk we just fetched.
		    TaskListMap::iterator it = helper_list_pending_for_transaction.find(tid);
		    if ( it != helper_list_pending_for_transaction.end() ) {
		      // Task is still in list_pending_for_transaction!
		      taskInfo = it->second;
		      taskInfo->child_chunks[0].push_back(fetchedChunk);
		    } // end if task still exists in list
		    count++;
		  } // end if input_type_is_ChunkID
		  if(!did_nothing)
		    break;
		} // end for go through list of input chunks for current task
		if(!did_nothing)
		  break;
		taskInfo->child_chunks_done_flags[0] = true; // we set this flag so that the work for this task can be skipped next time.
		task_iter++;
		number_of_checked_tasks++;
	      } // end while going through tasks in list_pending_for_transaction	      
	    } // end if (child_chunk_fetch_level >= 0)
	    if ( did_nothing ) {
	      for(int level = 0; level < child_chunk_fetch_level; level++) {
		TaskListSortedMap::reverse_iterator task_iter = list_pending_for_transaction.rbegin();
		int number_of_checked_tasks = 0;
		while (task_iter != list_pending_for_transaction.rend() && 
		       number_of_checked_tasks < n_tasks_to_fetch_data_for) {
		  TaskInfoStruct* taskInfo = task_iter->second;
		  if(taskInfo->child_chunks_done_flags[level+1]) {
		    task_iter++;
		    number_of_checked_tasks++;
		    continue;
		  }
		  int nInpChunks = taskInfo->inputChunks.size();
		  if ( !taskInfo->ready_to_run ) {
		    std::cerr << "Error in DataFetcherFunc while trying to fetch child chunks: ( !taskInfo->ready_to_run )." << std::endl;
		    throw std::runtime_error("Error in DataFetcherFunc while trying to fetch child chunks: ( !taskInfo->ready_to_run ).");
		  }
		  // Get ChunkIDs for child chunks at given level for this task.
		  std::vector<ChunkID> childChunkIDVec;
		  if(!getChildChunksForLevel(*taskInfo, level, childChunkIDVec)) {
		    std::cerr << "Error in DataFetcherFunc while trying to fetch child chunks: getChildChunksForLevel failed. Cache not enabled?" << std::endl;
		    std::cerr << "Note that chunk cache must be enabled when child chunk prefetching is used." << std::endl;
		    for(int i = 0; i < nInpChunks; i++) {
		      ChunkID ChunkID_tmp = taskInfo->inputChunks[i];
		      std::cerr << "taskInfo->inputChunks[" << i << "] <--> " << ChunkID_tmp.str() << std::endl;
		    }
		    throw std::runtime_error("Error in DataFetcherFunc while trying to fetch child chunks: getChildChunksForLevel failed.");
		  }
		  // Check if all those child chunks have already been fetched.
		  if(taskInfo->child_chunks[level+1].size() != childChunkIDVec.size()) {
		    // All child chunks have not been fetched. In this case, we fetch one more.
		    did_nothing = false;
		    ChunkID ChunkID_to_fetch = childChunkIDVec[taskInfo->child_chunks[level+1].size()];
		    TaskID tid = taskInfo->get_id();
		    UnlockMutex();
		    cht::shared_ptr<Chunk const> fetchedChildChunk;
		    if(ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().getChunkIfExists(ChunkID_to_fetch, fetchedChildChunk) == false) {
		      // getChunkIfExists returned false. This means the chunk has been deleted. This can happen if the task was run while we were trying to fetch the child chunk.
		      LockMutex_macro;
		      break;
		    }
		    // OK, we have fetched a child chunk. Add it to the list if possible.
		    LockMutex_macro;
		    // Check if task is still in list_pending_for_transaction.
		    // If it is not we do nothing with the chunk we just fetched.
		    TaskListMap::iterator it = helper_list_pending_for_transaction.find(tid);
		    if ( it != helper_list_pending_for_transaction.end() ) {
		      // Task is still in list_pending_for_transaction!
		      taskInfo = it->second;
		      taskInfo->child_chunks[level+1].push_back(fetchedChildChunk);
		    }
		  }
		  if(!did_nothing)
		    break;
		  taskInfo->child_chunks_done_flags[level+1] = true; // we set this flag so that the work for this task can be skipped next time.
		  task_iter++;
		  number_of_checked_tasks++;
		} // end while going through tasks in list_pending_for_transaction
		if(!did_nothing)
		  break;
	      } // end for child chunk fetch level
	    } // end if nothing done yet so consider fetching child chunks for different levels
	  } // end if nothing done yet so consider fetching child/ChunkID-input chunks
	  if ( did_nothing ) {
	    data_fetcher_cond.wait(mutex);
	    mutexIsLocked = true;
	    // After the wait, this thread has the mutex lock again. To help debbuging, we set the mutexLockDebugStr here.
	    mutexLockDebugStr = "DataFetcherFunc after data_fetcher_cond.wait(mutex).";
	  }
	} // end while
    } // end Worker::DataFetcherFunc()


    void Worker::MonitoringFunc() {
      double savedTimeForStatusReportFile = cht::get_wall_seconds();
      double savedTimeForStatusReportToOS = cht::get_wall_seconds();
      while(1)
	{
	  if(do_mutex_lock_checking) {
	    // We want to lock the mutex, but carefully using a
	    // separate thread, so that if the mutex seems locked by
	    // someone else for a long time, we can output an error
	    // message.
	    Threads::Thread* threadForMutexLock = new Threads::Thread("TS-mutexlock", global_thread_func_for_mutex_lock, this);
	    double savedTimeForMutexLockThread = cht::get_wall_seconds();
	    while(1) {
	      int one_milli_second = 1000;
	      usleep(1*one_milli_second);
	      if(threadForMutexLock->finished())
		break;
	      double secondsSinceMutexLockThreadStarted = cht::get_wall_seconds() - savedTimeForMutexLockThread;
	      if(secondsSinceMutexLockThreadStarted > 10.0) {
		std::cerr << "Error: Mutex lock takes too long time. mutexLockDebugStr = '" << mutexLockDebugStr 
			  << "'. Sleeping for 10 seconds before throwing exception, in case there are errors on other workers also." << std::endl;
		sleep(10);
		throw std::runtime_error("Error: Mutex lock takes too long time.");
	      }
	    }
	    // OK, thread finished. Now we know the mutex is locked.
	    delete threadForMutexLock;
	  }
	  else {
	    LockMutex_macro;
	  }
	  // Check if it is time to write a status report file.
	  if(status_report_file_interval_seconds != 0) {
	    double secondsSinceLastStatusReport = cht::get_wall_seconds() - savedTimeForStatusReportFile;
	    if(secondsSinceLastStatusReport > status_report_file_interval_seconds) {
	      savedTimeForStatusReportFile = cht::get_wall_seconds();
	      AccessKey key(this);
	      int my_rank = key.my_rank();
	      writeStatusReportFile(my_rank);
	    }
	  }
	  // Check if it is time to do a status report using the output service.
	  if(status_report_interval_milliseconds != 0) {
	    double secondsSinceLastStatusReport = cht::get_wall_seconds() - savedTimeForStatusReportToOS;
	    if(secondsSinceLastStatusReport*1000 > status_report_interval_milliseconds) {
	      savedTimeForStatusReportToOS = cht::get_wall_seconds();
	      doStatusReport();
	    }
	  }
	  bool doneFlag = finishFlagForThreadCommunicationWorker;
	  UnlockMutex();
	  if(doneFlag)
	    break;
	  // Do a rather long sleep call here to make sure the monitoring thread does not cause too much overhead.
	  int one_milli_second = 1000;
	  int sleepTime = 100*one_milli_second;
	  if(sleepTime > status_report_interval_milliseconds/10)
	    sleepTime = status_report_interval_milliseconds/10;
	  usleep(sleepTime);
	}
    }

    void Worker::MutexLockThreadFunc() {
      LockMutex_macro;
    }


    void Worker::setTaskInputChunkIDs( TaskInfoStruct & tis) {
      std::stringstream s;
      s << "Entered Worker::setTaskInputChunkIDs with taskTypeID = " 
	<< tis.get_taskTypeID();
      OS.outputDebug(s.str());
      tis.inputChunks.resize(tis.input_and_depend.size());
      std::vector<cht::ID>::iterator it = tis.input_and_depend.begin();
      size_t ind = 0;
      for(;it != tis.input_and_depend.end();it++,ind++) {
	if (it->is_taskID()) {
	  tis.inputChunks[ind] = getTaskResult(it->get_taskID());
	}
	else if (it->is_chunkID()) {
	  tis.inputChunks[ind] = it->get_chunkID();
	}
	else {
	  tis.inputChunks[ind] = CHUNK_ID_NULL;
	}
	// Check that the chunk with this chunk id is ready to use.
	if(ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().isRegisterChunkStillInProgress(tis.inputChunks[ind]))
	  throw std::runtime_error("Error in setTaskInputChunkIDs: isRegisterChunkStillInProgress returned true.");
	if(ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().isCopyChunkStillInProgress(tis.inputChunks[ind]))
	  throw std::runtime_error("Error in setTaskInputChunkIDs: isCopyChunkStillInProgress returned true.");
      } 
    }



    TaskID Worker::registerTaskStep1(std::string taskTypeID) {
      AccessKey key(this);
      int taskTypeInt = taskTypeIDToIntIDMap[taskTypeID];
      LockMutex_macro;
      // Determine new task ID
      TaskID taskID(key.my_rank(), ++idCounter, taskTypeInt);
      UnlockMutex();
      return taskID;
    }

    void Worker::registerTaskStep2(TaskID taskID,
				   TaskID creatorTaskID,
				   std::vector<cht::ID> const & input_and_dependencies)
    {
      AccessKey key(this);
      std::string taskTypeID = intIDToTaskTypeIDMap[taskID.taskTypeID];
      // Set ancestorDepth
      int ancestorDepth = 0; // creator == TASK_ID_NULL
      TaskInfoStruct const * creatorTask = NULL;
      if(creatorTaskID != TASK_ID_NULL)
	{
	  // Find creator.
	  LockMutex_macro;
	  TaskListMap::const_iterator it = list_doing_transaction_nonleaf_2.find(creatorTaskID);
	  if(it == list_doing_transaction_nonleaf_2.end())
	    throw std::runtime_error("Error! Task creator not found in list_doing_transaction_nonleaf_2.");
	  creatorTask = (*it).second;
	  ancestorDepth = creatorTask->get_ancestorDepth() + 1;
	  UnlockMutex();
	}
      
      // Create new TaskInfoStruct object.
      TaskInfoStruct* taskInfo = new TaskInfoStruct (taskID, 
						     creatorTaskID, 
						     creatorTask,
						     taskTypeID, 
						     input_and_dependencies,
						     taskID.ownerRank,
						     ancestorDepth);
      std::stringstream s;
      s << taskInfo->str() 
	<< " register_task" 
	<< " IDcreator = " << creatorTaskID.str() 
	<< " taskType = '" << taskTypeID << "'";
      OS.outputDebug(s.str());
      LockMutex_macro;
      assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
      // If some of the tasks in waitList are already finished, those should be removed from waitList.
      TaskIDList tasksToRemoveFromWaitList;
      TaskIDList::const_iterator it;
      for( it=taskInfo->waitList.begin() ; it != taskInfo->waitList.end(); it++ ) {
	if(ExistsInMap(list_finished, *it))
	  tasksToRemoveFromWaitList.push_back(*it);
      }
      for( it=tasksToRemoveFromWaitList.begin() ; it != tasksToRemoveFromWaitList.end(); it++ )
	taskInfo->waitList.remove(*it);
      // Now we need to add this task to either list_waiting_before_transaction or list_pending_for_transaction.
      // We need to place the task directly into list_pending_for_transaction if waitList is empty.
      if(taskInfo->waitList.empty()) {
	// Nothing in wait list. In this case we can add the new TaskInfoStruct object directly to list_pending_for_transaction.
	assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
	setTaskInputChunkIDs(*taskInfo);
	assert(list_pending_for_transaction.size() == helper_list_pending_for_transaction.size());
	AddTaskToPendingForTransactionLists(taskInfo);
      }
      else {
	// Add the new TaskInfoStruct object to list_waiting_before_execute.
	list_waiting_before_transaction.insert(TaskListMap::value_type(taskID, taskInfo));
      }    
      UnlockMutex();
    }


    bool Worker::checkIfTaskIsDoingFinalize(TaskID id) {
      AccessKey key(this);
      LockMutex_macro;
      bool result = ExistsInMap(list_running_finalize, id);
      UnlockMutex();
      return result;
    }


    // NOTE: Mutex should not be locked when calling this function
    ChunkID Worker::getTaskResult(TaskID const & calling_tid, 
						   TaskID const & tid) {
      AccessKey key(this);
      LockMutex_macro;
      if (!ExistsInMap(list_running_finalize, calling_tid))
	throw std::runtime_error("Calls to getTaskResult only allowed for "
				 "tasks running finalize()");
      // FIXME?: Check that calling_tid task is creator of tid task.
      ChunkID cid = getTaskResult(tid);
      UnlockMutex();
      return cid;
    }

    // NOTE: Mutex should be locked when calling this function
    ChunkID Worker::getTaskResult(TaskID const & tid) {
      AccessKey key(this);
      if(!ExistsInMap(list_finished, tid)) {
	int my_rank = key.my_rank();
	writeStatusReportFile(my_rank);
	std::cerr << "Error in Worker::getTaskResult: Task not in list_finished. tid.str() = " << tid.str() << std::endl;
	throw std::runtime_error("Error in Worker::getTaskResult: "
				 "Task not in list_finished");
      }
      ChunkID cid = list_finished[tid]->get_outputChunk();
      return cid;
    }


    Worker::Worker() 
      : OS( OutputService::Worker::instance() ),
	thisWorkerHasMotherTask(false),
	timeWhenTaskWasLastSetAsExecuting(0),
	accumulatedIdleTime(0),
	maxNoOfSimultaneousTasks(0),
 	pendingTaskCount_tot(0),
 	pendingTaskCount_readyToRun(0),
	n_worker_threads(0), // zero here, real value should be received from parent.
	n_prefetch_tasks(0), // zero here, real value should be received from parent.
	child_chunk_fetch_level(-1), // -1 here, real value should be received from parent.
	stealing_disabled_flag(false),
	steal_only_from_0_flag(false),
	status_report_file_interval_seconds(0),
	status_report_interval_milliseconds(0),
	do_mutex_lock_checking(false),
	mutexIsLocked(false),
	threadGroupWorkers(NULL),
	threadRequestHandler(NULL)
    {
      
    }

    void Worker::addTemporaryChunk(TaskID const & owner_task,
				   ID const & id) {
      AccessKey key(this);
      LockMutex_macro;
      if(!ExistsInMap(helper_list_pending_for_transaction, owner_task) && !ExistsInMap(list_stolen, owner_task))
	throw std::runtime_error("Error in Worker::addTemporaryChunk: "
				 "Task not in helper_list_pending_for_transaction");
      TaskInfoStruct* taskInfo = NULL;
      if(ExistsInMap(helper_list_pending_for_transaction, owner_task))
	taskInfo = helper_list_pending_for_transaction[owner_task];
      if(taskInfo == NULL) {
	// This may happen if the task was stolen
	if(!ExistsInMap(list_stolen, owner_task))
	  throw std::runtime_error("Error in addTemporaryChunk: (taskInfo == NULL) and not stolen.");
	// OK, it was stolen. In this case we do nothing.
	UnlockMutex();
	return;
      }
      taskInfo->addTemporaryChunk(id);
      UnlockMutex();
    }

  }; // end namespace
}; // end namespace

