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
#include "services/chunk_obj_service/ChunkObjService_worker.h"
#include "services/chunk_obj_cache_service/ChunkObjCacheService_worker.h"
#include "TaskSchedulerService_worker.h"
#include "Task.h"

namespace cht {

  TaskID Task::registerTask(std::string taskTypeID,
			    std::vector<cht::ID> const & input_and_dependencies) {
    return addTaskToRegisterLater(taskTypeID, input_and_dependencies);
  }

  TaskID Task::addTaskToRegisterLater(std::string taskTypeID, std::vector<cht::ID> const & input_and_dependencies) {
    TaskID theTaskID = cht::TaskSchedulerService::Worker::instance().registerTaskStep1(taskTypeID);
    task_to_register new_task_to_register;
    new_task_to_register.id = theTaskID;
    new_task_to_register.taskTypeID = taskTypeID;
    std::vector<cht::ID>::const_iterator it = input_and_dependencies.begin();
    for(;it != input_and_dependencies.end();it++) {
      if (it->is_taskID()) {
	TaskID tid = it->get_taskID();
	cht::ID id = tid;
	new_task_to_register.input_and_dependencies.push_back(id);
      }
      else {
	ChunkID cid = it->get_chunkID();
	cht::ID id = cid;
	new_task_to_register.input_and_dependencies.push_back(id);
      }
    }
    tasksToRegister.push_back(new_task_to_register);
    return theTaskID;
  }

  void Task::internal_registerTasksInList() {
    // Go through list and register all tasks, using the saved information.
    std::list<task_to_register>::iterator it;
    for ( it=tasksToRegister.begin() ; it != tasksToRegister.end(); it++ )
      cht::TaskSchedulerService::Worker::instance().registerTaskStep2(it->id, myID, it->input_and_dependencies);
  }

  int Task::internal_getNoOfTasksToRegister() {
    return tasksToRegister.size();
  }

  void Task::internal_get_child_task_ids(std::list<TaskID> & list) const {
    list.clear();
    std::list<task_to_register>::const_iterator it;
    for ( it=tasksToRegister.begin() ; it != tasksToRegister.end(); it++ )
      list.push_back(it->id);
  }

  void Task::addTemporaryChunk(ID const & id) const {
    cht::TaskSchedulerService::Worker::instance().addTemporaryChunk(myID, id);
  }

  ChunkID Task::registerChunk(Chunk const * tmp_ptr, 
			      std::string class_id_str) {
    ChunkID cid = ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().getIdForRegisterChunk(tmp_ptr, class_id_str);
    // We need to save the new ChunkID so that we can check later that the RegisterChunk operation is finished, before we consider this task really done.
    chunk_ids_for_registerChunk_operations.push_back(cid);
    chunk_ptrs_for_registerChunk_operations.push_back(tmp_ptr);
    return cid;
  }

  ChunkID Task::copyChunk(ChunkID cid) {
    if ( cid == CHUNK_ID_NULL )
      return CHUNK_ID_NULL;
    ChunkID cid_new = ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().getIdForCopyChunk(cid);
    // We need to save the new ChunkID so that we can check later that the copy operation is finished, before we consider this task really done.
    chunk_ids_new_for_copyChunk_operations.push_back(cid_new);
    chunk_ids_old_for_copyChunk_operations.push_back(cid);
    return cid_new;
  }

  void Task::internal_startRegisterAndCopyChunkOps() {
    // Start registerChunk operations.
    {
      std::list<ChunkID>::iterator it_cid = chunk_ids_for_registerChunk_operations.begin();
      std::list<Chunk const *>::iterator it_ptr = chunk_ptrs_for_registerChunk_operations.begin();
      while(it_cid != chunk_ids_for_registerChunk_operations.end()) {
	ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().startRegisterChunk(*it_ptr, *it_cid);
	bool stillInProgress = ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().isRegisterChunkStillInProgress(*it_cid);
	it_cid++;
	it_ptr++;
      }
    }
    // Start copyChunk operations.
    {
      std::list<ChunkID>::iterator it_cid_new = chunk_ids_new_for_copyChunk_operations.begin();
      std::list<ChunkID>::iterator it_cid_old = chunk_ids_old_for_copyChunk_operations.begin();
      while(it_cid_new != chunk_ids_new_for_copyChunk_operations.end()) {
	ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().startCopyChunk(*it_cid_old, *it_cid_new);
	it_cid_new++;
	it_cid_old++;
      }
    }
  }

  bool Task::internal_allRegisterAndCopyChunkOpsHaveFinished() {
    bool allDone = true;
    std::list<ChunkID>::iterator it;
    for ( it=chunk_ids_for_registerChunk_operations.begin() ; it != chunk_ids_for_registerChunk_operations.end(); it++ ) {
      if(ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().isRegisterChunkStillInProgress(*it))
	allDone = false;
    }
    for ( it=chunk_ids_new_for_copyChunk_operations.begin() ; it != chunk_ids_new_for_copyChunk_operations.end(); it++ ) {
      if(ChunkObjCacheService::Worker<ChunkObjService::Worker>::instance().isCopyChunkStillInProgress(*it))
	allDone = false;
    }
    return allDone;
  }


  ChunkID Task::getTaskResult(TaskID const & tid) const {
    return cht::TaskSchedulerService::Worker::
      instance().getTaskResult(myID, tid);		      
  }

  void Task::deleteChunk(ChunkID cid) const {
    return ChunkObjCacheService::Worker<ChunkObjService::Worker>
      ::instance().deleteChunk(cid);
  }

  void Task::output(std::string message, Output::Priority prio) const {
    OutputService::Worker::instance().output(prio, message);    
  }

  cht::ChunkID Task::getInputChunkID(cht::Chunk const & chunk) const {
    for (unsigned int ind = 0;ind < inputChunks.size();ind++) 
      if (inputChunks[ind] == &chunk)
	return inputChunkIDs[ind];
    throw std::runtime_error("getInputChunkID(): Chunk not found in list of input chunks.");
  }


  cht::ID Task::internal_execute(TaskID myID_, bool fallback_) {
    // FIXME: Set ID earlier, for example when task is created.  Then
    //        this function would not be needed, one could call the
    //        0-parameter "execute()" directly.
    myID = myID_;
    fallback = fallback_;
    return internal_execute();
  }

  void Task::internal_setInputChunks(std::vector<cht::shared_ptr<Chunk const> > const & inpChunks) {
    inputChunks = inpChunks;
  }
  void Task::internal_setInputChunkIDs(std::vector<ChunkID> const & inpChunkIDs) {
    inputChunkIDs = inpChunkIDs;
  }
  void Task::internal_setinputChunksAndChunkIDs(std::vector<BaseObj const *> const & inpChunksAndChunkIDs) {
    inputChunksAndChunkIDs = inpChunksAndChunkIDs;
  }
  void Task::internal_clearInputVectors() {
    inputChunksAndChunkIDs.clear();
    inputChunks.clear();
    inputChunkIDs.clear();
  }

  bool internal::registerArgTypes(std::string & task_type_str, 
				  std::list<std::string> & task_input_arg_type_strs,
				  std::string & task_output_type_str) {
    return cht::arg_manager<Task>::instance().registerArgumentTypes(task_type_str, task_input_arg_type_strs, task_output_type_str);
  }

  cht::ID Task::execute(cht::ChunkID const & c1) {
    throw std::runtime_error("No execute(cht::ChunkID const &) fallback function implemented.");
  }   
  cht::ID Task::execute(cht::ChunkID const & c1, cht::ChunkID const & c2) {
    throw std::runtime_error("No execute(cht::ChunkID const &, cht::ChunkID const &) fallback function implemented.");
  }   
  cht::ID Task::execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3) {
    throw std::runtime_error("No execute(cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &) fallback function implemented.");
  }   
  cht::ID Task::execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3, 
			cht::ChunkID const & c4) {
    throw std::runtime_error("No execute(cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &) fallback function implemented.");
  }   
  cht::ID Task::execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3, 
			cht::ChunkID const & c4, cht::ChunkID const & c5) {
    throw std::runtime_error("No execute(cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &) fallback function implemented.");
  }   
  cht::ID Task::execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3, 
			cht::ChunkID const & c4, cht::ChunkID const & c5, cht::ChunkID const & c6) {
    throw std::runtime_error("No execute(cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &) fallback function implemented.");
  }   
  cht::ID Task::execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3, 
			cht::ChunkID const & c4, cht::ChunkID const & c5, cht::ChunkID const & c6, 
			cht::ChunkID const & c7) {
    throw std::runtime_error("No execute(cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &) fallback function implemented.");
  }   
  cht::ID Task::execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3, 
			cht::ChunkID const & c4, cht::ChunkID const & c5, cht::ChunkID const & c6, 
			cht::ChunkID const & c7, cht::ChunkID const & c8) {
    throw std::runtime_error("No execute(cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &) fallback function implemented.");
  }   

} // end namespace cht
