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
#ifndef CHUNKS_AND_TASKS_PARENT_HEADER
#define CHUNKS_AND_TASKS_PARENT_HEADER
#include <cstdlib>
#include <memory>
#include "services/chunks_and_tasks_params.h"
#include "services/chunk_obj_service/ChunkObject.h"
#include "utilities/cht_shared_ptr.h"
#include "utilities/cht_static_check.h"

namespace cht {
  // CHT Start and Stop
  void start();
  void stop();

  namespace extras {  
    // CHT EXTRAS 
    void setNWorkers(unsigned int n_workers);
    void setCacheMode(Cache::Mode mode);
    void setCacheSize(size_t cache_memory_usage_limit);
    void setNoOfWorkerThreads(unsigned int n_worker_threads);
    void setNoOfPrefetchTasks(unsigned int n_prefetch_tasks);
    void setChildChunkFetchLevel(int child_chunk_fetch_level);
    void setDebugParams(bool, bool, bool, int, int, bool, double);
  } // end namespace extras
  void setOutputMode(Output::Mode mode);
  void setOutputLevel(Output::Priority prio);

  // Chunk Management Functions
  template<typename ChunkObjType>
    ChunkID registerChunk(ChunkObjType const * obj_ptr);
  template<typename ChunkObjType>
    void getChunk(ChunkID cid, cht::shared_ptr<ChunkObjType const> & obj_ptr);
  void deleteChunk(ChunkID id);

  // Function for task execution
  template<typename TaskType>
    ChunkID executeMotherTask(std::vector<ChunkID> const & task_input_chunks);

  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1);
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2);
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2, 
			      ChunkID const c3);  
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2, 
			      ChunkID const c3, 
			      ChunkID const c4);  
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2, 
			      ChunkID const c3, 
			      ChunkID const c4, 
			      ChunkID const c5);  
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2, 
			      ChunkID const c3, 
			      ChunkID const c4, 
			      ChunkID const c5, 
			      ChunkID const c6);  
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2, 
			      ChunkID const c3, 
			      ChunkID const c4, 
			      ChunkID const c5, 
			      ChunkID const c6, 
			      ChunkID const c7);  
  
  // Output functions  
  void output(std::string message, Output::Priority prio = Output::Info);
  
  // Statistics
  void resetStatistics();
  void reportStatistics();

  // Internal functions  
  namespace internal {
    ChunkID executeMotherTask(std::string taskTypeID, 
			      std::vector<ChunkID> const & task_input_chunks);
    ChunkID registerChunk(Chunk const * tmp_ptr, 
			  std::string class_id_str);    
    void getChunk(ChunkID id, 
		  cht::shared_ptr<Chunk const> & objPtr,
		  std::string class_id_str);
  } // end namespace internal


  // Definitions of function templates

  template<typename ChunkObjType>
    ChunkID registerChunk(ChunkObjType const * obj_ptr) {
    return internal::registerChunk(obj_ptr, ChunkObjType::get_class_id());    
  }
  template<typename ChunkObjType>
    void getChunk(ChunkID cid, cht::shared_ptr<ChunkObjType const> & obj_ptr) {
    std::string class_id_str = ChunkObjType::get_class_id();
    cht::shared_ptr<Chunk const> objCOPtr;
    internal::getChunk(cid, objCOPtr, class_id_str);
    // Ugly cast needed here!!
    obj_ptr = cht::shared_ptr<ChunkObjType const>( (ChunkObjType*)&(*objCOPtr), objCOPtr.getRefCountPtr() );
  }

  template<typename TaskType>
    ChunkID executeMotherTask(std::vector<ChunkID> const & task_input_chunks) {
    std::string taskTypeID = TaskType::get_class_id();
    return internal::executeMotherTask(taskTypeID, task_input_chunks);
  }
  
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1) {
    CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 1, 
			  WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
    std::vector<ChunkID> task_input_chunks(1);
    task_input_chunks[0] = c1;
    return executeMotherTask<TaskType>(task_input_chunks);
  }
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, ChunkID const c2) {
    CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 2, 
			  WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
    std::vector<ChunkID> task_input_chunks(2);
    task_input_chunks[0] = c1;
    task_input_chunks[1] = c2;
    return executeMotherTask<TaskType>(task_input_chunks);
  }
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2, 
			      ChunkID const c3) {
    CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 3, 
			  WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
    std::vector<ChunkID> task_input_chunks(3);
    task_input_chunks[0] = c1;
    task_input_chunks[1] = c2;
    task_input_chunks[2] = c3;
    return executeMotherTask<TaskType>(task_input_chunks);
  }
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2, 
			      ChunkID const c3,
			      ChunkID const c4) {
    CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 4, 
			  WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
    std::vector<ChunkID> task_input_chunks(4);
    task_input_chunks[0] = c1;
    task_input_chunks[1] = c2;
    task_input_chunks[2] = c3;
    task_input_chunks[3] = c4;
    return executeMotherTask<TaskType>(task_input_chunks);
  }
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2, 
			      ChunkID const c3,
			      ChunkID const c4,
			      ChunkID const c5) {
    CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 5, 
			  WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
    std::vector<ChunkID> task_input_chunks(5);
    task_input_chunks[0] = c1;
    task_input_chunks[1] = c2;
    task_input_chunks[2] = c3;
    task_input_chunks[3] = c4;
    task_input_chunks[4] = c5;
    return executeMotherTask<TaskType>(task_input_chunks);
  }
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2, 
			      ChunkID const c3,
			      ChunkID const c4,
			      ChunkID const c5,
			      ChunkID const c6) {
    CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 6, 
			  WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
    std::vector<ChunkID> task_input_chunks(6);
    task_input_chunks[0] = c1;
    task_input_chunks[1] = c2;
    task_input_chunks[2] = c3;
    task_input_chunks[3] = c4;
    task_input_chunks[4] = c5;
    task_input_chunks[5] = c6;
    return executeMotherTask<TaskType>(task_input_chunks);
  }
  template<typename TaskType>
    ChunkID executeMotherTask(ChunkID const c1, 
			      ChunkID const c2, 
			      ChunkID const c3,
			      ChunkID const c4,
			      ChunkID const c5,
			      ChunkID const c6,
			      ChunkID const c7) {
    CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 7, 
			  WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
    std::vector<ChunkID> task_input_chunks(7);
    task_input_chunks[0] = c1;
    task_input_chunks[1] = c2;
    task_input_chunks[2] = c3;
    task_input_chunks[3] = c4;
    task_input_chunks[4] = c5;
    task_input_chunks[5] = c6;
    task_input_chunks[6] = c7;
    return executeMotherTask<TaskType>(task_input_chunks);
  }
  


} // end namespace cht

#endif
