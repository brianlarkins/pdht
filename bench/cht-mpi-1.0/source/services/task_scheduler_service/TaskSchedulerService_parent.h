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
#ifndef TASKSCHEDULERSERVICE_PARENT_HEADER
#define TASKSCHEDULERSERVICE_PARENT_HEADER
#include "services/Service_parent.h"
#include <list>
#include <fstream>
#include "utilities/Singleton.h"
#include "utilities/cht_utils.h"
#include "services/MPIWrapperInclude.h"
#include "services/output_service/OutputService_parent.h"
#include "services/task_scheduler_service/Task.h"
#include "services/task_scheduler_service/TaskSchedulerService.h"

namespace cht {
  namespace TaskSchedulerService {

    class Parent : 
    public Base<Service::Parent>,
      public Singleton<Parent> {
	friend class Singleton<Parent>;

    private:
	OutputService::Parent& OS;

	/* The following member variables exist for both parent and worker,
	 * but for the parent the vaiables are only used in order to send
	 * the values to the workers. */
	int n_worker_threads;  /* Number of worker threads */
	int n_prefetch_tasks; /* Max number of tasks for which to prefetch input chunks. */
	int child_chunk_fetch_level; /* Level of child-chunks to fetch (-1 means disabled). */
	bool stealing_disabled_flag; /* Mainly for debugging. */
	bool steal_only_from_0_flag; /* Mainly for debugging. */
	bool send_to_separate_workers_flag; /* Mainly for debugging. */
	int status_report_file_interval_seconds; /* Mainly for debugging. */
	int status_report_interval_milliseconds; /* Mainly for debugging. */
	bool do_mutex_lock_checking; /* Mainly for debugging. */

	// If send_to_separate_workers_flag is set, we use the following variable to keep track of which worker the next mother task should be sent to.
	int workerRankCounter;

    public:
	Parent();
	virtual void start_derived();
	virtual void stop_derived();
  
    public:
	ChunkID 
	  executeMotherTask(std::string taskTypeID,
			    std::vector<ChunkID> const & taskInputChunks);
	void setNoOfWorkerThreads(int n_worker_threads_);
	void setNoOfPrefetchTasks(int n_prefetch_tasks_);
	void setChildChunkFetchLevel(int child_chunk_fetch_level_);
	void setDebugParams(bool stealing_disabled_flag_,
			    bool steal_only_from_0_flag_,
			    bool send_to_separate_workers_flag_,
			    int status_report_file_interval_seconds_,
			    int status_report_interval_milliseconds_,
			    bool do_mutex_lock_checking);
  
      };


  }; // end namespace
}; // end namespace

#endif
