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
#include "services/ParentInterface.h"
#include <stdexcept>
#include <unistd.h>
#include "services/service_manager/spawn_workers_func.h"
#include "services/service_manager/ServiceManager.h"
#include "services/output_service/OutputService_parent.h"
#include "services/chunk_obj_service/ChunkObjService_parent.h"
#include "services/chunk_obj_cache_service/ChunkObjCacheService_parent.h"
#include "services/task_scheduler_service/TaskSchedulerService_parent.h"
#include "services/service_manager/determine_process_roles.h"
#include "services/service_manager/service_worker_func.h"
#include "services/services_registration_worker.h"
#include "services/services_registration_parent.h"
#include "license.h"

namespace cht {

#if BUILD_AS_SINGLE_PROGRAM
  static bool const services_registered_worker = register_internal_services_worker();
#endif
  static bool const services_registered_parent = register_internal_services_parent();

  static void WaitForFinalExitSignalFromParent() {
    MPI_Wrapper & MW = MPI_Wrapper::instance();
    // Note: we could just use MPI_Recv here, but that turned out to be
    // consuming a lot of cpu time for some MPI implementations,
    // e.g. openmpi/1.3.4 on isis. To be sure that "IDLE" processes
    // really do not use any cpu time, we use a loop with sleep() and
    // Iprobe() instead.
    while(1) {
      int flag;
      MW._MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
      if(flag)
	break;
      sleep(1);
    }
    MW._MPI_Recv(0, 0, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  ParentInterface::ParentInterface() 
    : service_worker_program_name("cht_worker"), n_workers(0), mainThread(0), access_granted(true) {    
  }
  ParentInterface::~ParentInterface() {
    // NOTE: "delete mainThread" cannot be called here because then
    // singleton objects (ParentInterface and Manager of threads) have
    // to be destructed in a particular order, which is difficult to
    // achieve.
    //
    // delete mainThread;
  }

  // Verify that access is allowed, meaning that the calling process
  // is the parent process and not a worker.
  void ParentInterface::verifyInstanceAccess() {
    if (!access_granted)
      throw std::runtime_error("Attempt to access cht function intended only for parent.");
  }     

  // Workers report that they are workers and are consequently not
  // allowed access to parent functions.
  void ParentInterface::i_am_worker() {
    access_granted = false;
  }

  void ParentInterface::start() {
    int required = MPI_WRAPPER_MPI_THREAD_LEVEL_REQUIRED;
    int provided = -1;
    MPI_Wrapper & MW = MPI_Wrapper::instance();
    MW._MPI_Init_thread(NULL, NULL, required, &provided);

    if (n_workers == 0) {
      int flag;
      int* universe_sizep;
      MW._MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_UNIVERSE_SIZE,  
			   &universe_sizep, &flag);
      if (!flag) 
        throw std::runtime_error("ParentInterface::start(): This MPI does not support MPI_UNIVERSE_SIZE. ( To aviod this problem, consider calling cht::extras::setNWorkers() before calling cht::start() )."); 
      n_workers = *universe_sizep;
    }

#if BUILD_AS_SINGLE_PROGRAM
    ProcessRole myRole;
    determine_process_roles(n_workers, myRole, intraComm, comm_to_workers);
    if(myRole == Idle) {
      WaitForFinalExitSignalFromParent();
      MW._MPI_Finalize();
      std::exit(0);
    }
    if(myRole == Worker) {
      service_worker_func(intraComm, comm_to_workers);
      WaitForFinalExitSignalFromParent();
      MW._MPI_Finalize();
      std::exit(0);
    }
    // At this point we know this is the parent.
#else
    spawnWorkers(service_worker_program_name, n_workers, comm_to_workers);
#endif
    mainThread = new Threads::Thread("main");

    ServiceManager::instance().Init(&comm_to_workers);    
    ServiceManager::instance().startService("OutputService");
    OutputService::Parent::instance().outputInfo("OutputService started.");
    OutputService::Parent::instance().outputInfo(CHT_LICENSE_TEXT);
    ServiceManager::instance().startService("ChunkObjService");
    ServiceManager::instance().startService("ChunkObjCacheService<ChunkObjService>");    
    ServiceManager::instance().startService("TaskSchedulerService");
  } // end start()

  // One could allow for starting and stopping the cht several times
  // within a single program. Then MPI_Finalize would be called
  // elsewhere, e.g. in the destructor. And MPI_Init would only be
  // called the first time the start() function is called.
  void ParentInterface::stop() {
    ServiceManager::instance().stopService("TaskSchedulerService");
    ServiceManager::instance().stopService("ChunkObjCacheService<ChunkObjService>");
    ServiceManager::instance().stopService("ChunkObjService");
    ServiceManager::instance().stopService("OutputService");
    ServiceManager::instance().Finalize();
    MPI_Wrapper & MW = MPI_Wrapper::instance();
#ifdef BUILD_AS_SINGLE_PROGRAM
    int nProcsTot;
    MW._MPI_Comm_size(MPI_COMM_WORLD, &nProcsTot);
    int myRank;
    MW._MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    for(int i = 0; i < nProcsTot; i++) {
      if(i != myRank)
	MW._MPI_Send(0, 0, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
    }
#endif
    MW._MPI_Finalize();    
    delete mainThread; // see comment in ParentInterface destructor
  }

  void ParentInterface::setNWorkers(unsigned int n_workers_) {
    n_workers = n_workers_;
  }    
  void ParentInterface::setOutputMode(Output::Mode mode) {
    OutputService::Parent::instance().setMode(mode);
  }
  void ParentInterface::setOutputLevel(Output::Priority prio) {
    OutputService::Parent::instance().setLevel(prio);
  }

  void ParentInterface::setCacheMode(extras::Cache::Mode mode) {
    ChunkObjCacheService::Parent<ChunkObjService::Parent>::instance().setCacheMode(mode);    
  }
  void ParentInterface::setCacheSize(size_t cache_memory_usage_limit) {
    ChunkObjCacheService::Parent<ChunkObjService::Parent>::instance().setCacheSize(cache_memory_usage_limit);    
  }
  void ParentInterface::setNoOfWorkerThreads(unsigned int n_worker_threads) {
    TaskSchedulerService::Parent::instance().setNoOfWorkerThreads(n_worker_threads);
  }
  void ParentInterface::setNoOfPrefetchTasks(unsigned int n_prefetch_tasks) {
    TaskSchedulerService::Parent::instance().setNoOfPrefetchTasks(n_prefetch_tasks);
  }
  void ParentInterface::setChildChunkFetchLevel(int child_chunk_fetch_level) {
    TaskSchedulerService::Parent::instance().setChildChunkFetchLevel(child_chunk_fetch_level);
  }
  void ParentInterface::setDebugParams(bool stealing_disabled_flag,
				       bool steal_only_from_0_flag,
				       bool send_to_separate_workers_flag,
				       int status_report_file_interval_seconds,
				       int status_report_interval_milliseconds,
				       bool do_mutex_lock_checking,
				       double probability_to_create_chunks_locally) {
    TaskSchedulerService::Parent::instance().setDebugParams(stealing_disabled_flag,
							    steal_only_from_0_flag,
							    send_to_separate_workers_flag,
							    status_report_file_interval_seconds,
							    status_report_interval_milliseconds,
							    do_mutex_lock_checking);
    ChunkObjService::Parent::instance().setDebugParams(probability_to_create_chunks_locally);
  }


} // end namespace cht
