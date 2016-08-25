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
#include "ParentInterface.h"
#include "chunks_and_tasks_parent.h"
#include "services/task_scheduler_service/TaskSchedulerService_parent.h"
#include "services/service_manager/ServiceManager.h"
#include "services/chunk_obj_service/ChunkObjService_parent.h"
#include "services/chunk_obj_cache_service/ChunkObjCacheService_parent.h"

namespace cht {
  void start() {
    cht::ParentInterface::instance().start();
  }
  void stop() {
    cht::ParentInterface::instance().stop();
  }
  namespace extras {
    void setNWorkers(unsigned int n_workers) {
      cht::ParentInterface::instance().setNWorkers(n_workers);
    }
    void setCacheMode(Cache::Mode mode) {
      cht::ParentInterface::instance().setCacheMode(mode);
    }
    void setCacheSize(size_t cache_memory_usage_limit) {
      cht::ParentInterface::instance().setCacheSize(cache_memory_usage_limit);
    }
    void setNoOfWorkerThreads(unsigned int n_worker_threads) {
      cht::ParentInterface::instance().setNoOfWorkerThreads(n_worker_threads);
    }
    void setNoOfPrefetchTasks(unsigned int n_prefetch_tasks) {
      cht::ParentInterface::instance().setNoOfPrefetchTasks(n_prefetch_tasks);
    }
    void setChildChunkFetchLevel(int child_chunk_fetch_level) {
      cht::ParentInterface::instance().setChildChunkFetchLevel(child_chunk_fetch_level);
    }
    void setDebugParams(bool stealing_disabled_flag,
			bool steal_only_from_0_flag,
			bool send_to_separate_workers_flag,
			int status_report_file_interval_seconds,
			int status_report_interval_milliseconds,
			bool do_mutex_lock_checking,
			double probability_to_register_chunks_locally) {
      cht::ParentInterface::instance().setDebugParams(stealing_disabled_flag, 
						      steal_only_from_0_flag, 
						      send_to_separate_workers_flag,
						      status_report_file_interval_seconds,
						      status_report_interval_milliseconds,
						      do_mutex_lock_checking, 
						      probability_to_register_chunks_locally);
    }
  } // end namespace extras
  void setOutputMode(Output::Mode mode) {
    cht::ParentInterface::instance().setOutputMode(mode);
  }
  void setOutputLevel(Output::Priority prio) {
    cht::ParentInterface::instance().setOutputLevel(prio);
  }

  void deleteChunk(ChunkID id) {
    // Make sure that access is allowed, i.e. the parent is calling.
    cht::ParentInterface::instance(); 
    ChunkObjCacheService::Parent<ChunkObjService::Parent>::instance().
      deleteChunk(id);
  }

  void output(std::string message, Output::Priority prio) {
    // Make sure that access is allowed, i.e. the parent is calling.
    cht::ParentInterface::instance();
    OutputService::Parent::instance().output(prio, message);
  }

  void resetStatistics() {
    // Make sure that access is allowed, i.e. the parent is calling.
    cht::ParentInterface::instance();
    cht::ServiceManager::instance().resetStatistics("TaskSchedulerService");
    cht::ServiceManager::instance().resetStatistics("ChunkObjService");
    cht::ServiceManager::instance().resetStatistics("ChunkObjCacheService<ChunkObjService>");
  }
  void reportStatistics() {
    // Make sure that access is allowed, i.e. the parent is calling.
    cht::ParentInterface::instance();
    cht::ServiceManager::instance().reportStatistics("TaskSchedulerService");
    cht::ServiceManager::instance().reportStatistics("ChunkObjService");
    cht::ServiceManager::instance().reportStatistics("ChunkObjCacheService<ChunkObjService>");

  }

  namespace internal {
    // Help functions needed to hide internal stuff
    ChunkID executeMotherTask(std::string taskTypeID, 
			      std::vector<ChunkID> const & task_input_chunks) {
      // Make sure that access is allowed, i.e. the parent is calling.
      cht::ParentInterface::instance();
      return TaskSchedulerService::Parent::instance().executeMotherTask(taskTypeID, task_input_chunks);      
    } 

    ChunkID registerChunk(Chunk const * tmp_ptr, 
			  std::string class_id_str) {
      // Make sure that access is allowed, i.e. the parent is calling.
      cht::ParentInterface::instance();
      return ChunkObjCacheService::Parent<ChunkObjService::Parent>::instance().
	registerChunkDirectly(tmp_ptr, class_id_str);
    }    
    

    void getChunk(ChunkID id, 
		  cht::shared_ptr<Chunk const> & objPtr,
		  std::string class_id_str) {
      // Make sure that access is allowed, i.e. the parent is calling.
      cht::ParentInterface::instance();
      int idInt = ChunkObjService::Parent::instance().getChunkTypeIDInt(class_id_str);
      if (idInt != id.chunkTypeID)
	throw std::runtime_error("getChunk(...) (parent): Chunk types do not match.");
      ChunkObjCacheService::Parent<ChunkObjService::Parent>::instance().
	getChunk(id, objPtr);
    }

  } // end namespace internal
  
} // end namespace cht
