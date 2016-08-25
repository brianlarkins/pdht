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
#ifndef TASKSCHEDULERSERVICE_HEADER
#define TASKSCHEDULERSERVICE_HEADER
#include <list>
#include <fstream>
#include <cstring>
#include "utilities/Singleton.h"
#include "utilities/cht_utils.h"
#include "services/Service.h"
#include "services/MPIWrapperInclude.h"
#include "services/output_service/OutputService.h"
#include "services/task_scheduler_service/Task.h"

namespace cht {
  /** A service for management and execution of tasks. */
  namespace TaskSchedulerService {

    // Service_type is typically Service_worker or Service_parent
    template<typename Service_type> 
      class Base : public Service_type {
    public:
      static std::string object_type_id();
    protected:
      typedef std::map<std::string, int> strToIntMap;
      typedef std::map<int, std::string> intToStrMap;
      strToIntMap taskTypeIDToIntIDMap;
      intToStrMap intIDToTaskTypeIDMap;

      enum MPI_tags {
	TAG_Steal_task_attempt,
	TAG_Steal_failed,
	TAG_Stolen_task_info,
	TAG_Task_finished,
	TAG_Stolen_task_finished,
	TAG_Mother_task_info,
	TAG_Mother_task_finished,
	TAG_Thread_params,
	TAG_Workers_finished,
	TAG_TaskTypeID_map
      };

      int idCounter;
      Base() 
	: idCounter(0) {}
      
      typedef std::list<TaskID> TaskIDList;


      // Can this be implemented more smoothly? YES!
      // Elias note: this is tricky because of the way TaskInfoStruct
      // objects are copied. I got segfault crashes when I had waitList as
      // TaskIDList instead of TaskIDList*.
      class TaskInfoStruct 
      {
      private:	
	TaskID id;
	TaskID creator;
	TaskInfoStruct const * creator_tis; 
	std::string taskTypeID; 
	ID outputChunkOrTask; // local data (to be converted to ChunkID outputChunk)
	ChunkID outputChunk;  // local data
	Task* task;
	bool stolen;
	int ownerWorkerRank; // important if task is stolen, 
	                     // FIXME: This information is already 
	                     // available in the id member!
	int ancestorDepth; // Keeps track of how deep down in hierarchy the task is.
	std::list<ID> temporaryChunks; // (local data) List of temporary chunks owned by task 
      public:
	std::vector<ChunkID> inputChunks; // FIXME: change name to inputChunkIDs
	std::vector<cht::shared_ptr<Chunk const> > input_chunks; // local data
	std::vector< std::list< cht::shared_ptr<Chunk const> > > child_chunks; // local data
	std::vector<bool> child_chunks_done_flags; // local data
	bool fallback_execute_should_be_used; // local data (set when adding task to pending_for_execute_list)
	bool ready_to_run; // local data (true when either fallback execute should be used or all input chunks have been fetched)
	bool executing; // local data
	bool executed; // local data
	bool chunkOpsStarted; // local data
	int helperWorkerRank; // Local data. If the task is stolen, this is the rank of the worker that has stolen the task.
	std::list<ChunkID::MajorOwnerInfo> ownerList; // Local data. This list contains info about the biggest owners of input data for this task (recursively).
	size_t totInputDataSize; // Local data.
	size_t currWorkerOwnedSize; // Local data. This is simply the size in ownerList for this worker, or zero if current worker not present in ownerList.
	std::vector<ID> input_and_depend; // local data
	TaskIDList waitList;                   // local data
	TaskIDList childList; // Local data. The child list is set after execute(), then each TaskIDs is removed from childList when that child task is completely finished.
      public:
	bool operator==  ( TaskInfoStruct const & x ) const {
	  return id == x.id;
	}
	TaskID const & get_id() const {return id;}
	TaskID const & get_creator() const {return creator;}
	TaskInfoStruct const * get_creator_tis() const {return creator_tis;}
	std::string const & get_taskTypeID() const {return taskTypeID;} 
	bool outputChunkOrTask_is_null() const {return outputChunkOrTask.is_null();}
	ID const & get_outputChunkOrTask() const {return outputChunkOrTask;}
	void set_outputChunkOrTask(ID const & outputChunkOrTask_) {outputChunkOrTask = outputChunkOrTask_;}
	ChunkID const & get_outputChunk() const {return outputChunk;}
	void set_outputChunk(ChunkID const & outputChunk_) {outputChunk = outputChunk_;}
	Task * get_task() const {return task;}
	void set_task(Task* task_) {task = task_;}
	bool const & get_stolen() const {return stolen;}
	void set_stolen(bool const stolen_) {stolen = stolen_;}
	int const & get_ownerWorkerRank() const {return ownerWorkerRank;}
	int const & get_ancestorDepth() const {return ancestorDepth;}
	void addTemporaryChunk(ID const & id) {
	  temporaryChunks.push_back(id);
	}
	std::list<ID> const & getTemporaryChunks() const {
	  return temporaryChunks;
	}
	void clearTemporaryChunkList() {
	  temporaryChunks.clear();
	}

	std::string str() const {
	  std::stringstream ss;
	  ss << "TASKINFO ID = " << this->id.str();
	  return ss.str();
	}



      TaskInfoStruct(TaskID id_, 
		     TaskID creator_,
		     TaskInfoStruct const * creator_tis_, 
		     std::string taskTypeID_,
		     std::vector<cht::ID> input_and_depend_,
		     int ownerWorkerRank_,
		     int ancestorDepth_)
      :   id(id_), 
	  creator(creator_),
	  creator_tis(creator_tis_),
	  taskTypeID(taskTypeID_),
	  fallback_execute_should_be_used(false),
	  ready_to_run(false),
	  executing(false),
	  executed(false),
	  chunkOpsStarted(false),
	  helperWorkerRank(-1),
	  totInputDataSize(0),
	  currWorkerOwnedSize(0),
	  input_and_depend(input_and_depend_),
	  task(0),
	  stolen(false),
	  ownerWorkerRank(ownerWorkerRank_),
	  ancestorDepth(ancestorDepth_)
	{
	  std::vector<cht::ID>::iterator it = input_and_depend.begin();
	  for(;it != input_and_depend.end();it++) {
	    if (it->is_taskID())
	      waitList.push_back(it->get_taskID());
	  }
	}
      TaskInfoStruct(TaskID id_, 
		     TaskID creator_,
		     TaskInfoStruct const * creator_tis_, 
		     std::string taskTypeID_,
		     std::vector<ChunkID> inputChunks_,
		     int ownerWorkerRank_,
		     int ancestorDepth_)
	: id(id_), 
	  creator(creator_),
	  creator_tis(creator_tis_),
	  taskTypeID(taskTypeID_),
	  inputChunks(inputChunks_),
	  fallback_execute_should_be_used(false),
	  ready_to_run(false),
	  executing(false),
	  executed(false),
	  chunkOpsStarted(false),
	  helperWorkerRank(-1),
	  totInputDataSize(0),
	  currWorkerOwnedSize(0),
	  task(0),
	  stolen(false),
	  ownerWorkerRank(ownerWorkerRank_),
	  ancestorDepth(ancestorDepth_)
	  {}
      TaskInfoStruct() 
	: creator_tis(0), task(0), stolen(false), ownerWorkerRank(-1),
	  fallback_execute_should_be_used(false), ready_to_run(false),
	  executing(false), executed(false), chunkOpsStarted(false)
	  {
	  }
      ~TaskInfoStruct() { 
	delete task;
      }

        int pack_size() {
	  return
	    2*sizeof(TaskID) +          // id AND creator
	    (taskTypeID.length()+1) +   // taskTypeID
	    sizeof(size_t) +            // length of inputChunks
	    inputChunks.size()*sizeof(ChunkID) + // inputChunks
	    sizeof(bool) +                 // stolen
	    2*sizeof(int);              // ownerWorkerRank AND ancestorDepth
	}
	void pack(char* buffer) const {
	  char* p = buffer;
	  memcpy(p, &id,      sizeof(TaskID));
	  p += sizeof(TaskID);
	  memcpy(p, &creator, sizeof(TaskID));
	  p += sizeof(TaskID);
	  memcpy(p, taskTypeID.c_str(), taskTypeID.length()+1);
	  p += taskTypeID.length()+1;
	  size_t nInputChunks = inputChunks.size();
	  memcpy(p, &nInputChunks, sizeof(size_t));
	  p += sizeof(size_t);
	  memcpy(p, &inputChunks[0], 
		 inputChunks.size()*sizeof(ChunkID));
	  p += inputChunks.size()*sizeof(ChunkID);
	  memcpy(p, &stolen, sizeof(bool));
	  p += sizeof(bool);
	  memcpy(p, &ownerWorkerRank, sizeof(int));
	  p += sizeof(int);
	  memcpy(p, &ancestorDepth,   sizeof(int));
	}
	void unpack(char const * const buffer, int bufferSize) {
	  assert(bufferSize > 2*sizeof(TaskID) + sizeof(size_t) + sizeof(bool) + 2*sizeof(int));
	  char const * p = buffer;
	  memcpy(&id, p,      sizeof(TaskID));
	  p += sizeof(TaskID);
	  memcpy(&creator, p, sizeof(TaskID));
	  p += sizeof(TaskID);
	  // Check for end-of-string byte carefully, without going outside buffer.
	  int strLen = 0;
	  const char* endPtr = buffer + bufferSize;
	  const char* q = p;
	  while(q < endPtr) {
	    if(*q == '\0')
	      break;
	    strLen++;
	    q++;
	  }
	  int nBytesLeft = bufferSize - (p - buffer);
	  assert(nBytesLeft >= strLen + sizeof(size_t) + sizeof(bool) + 2*sizeof(int));
	  taskTypeID = std::string(p);
	  p += taskTypeID.length()+1;
	  size_t nInputChunks;
	  memcpy(&nInputChunks, p, sizeof(size_t));
	  p += sizeof(size_t);
	  nBytesLeft = bufferSize - (p - buffer);
	  assert(nBytesLeft == nInputChunks*sizeof(ChunkID) + sizeof(bool) + 2*sizeof(int));
	  inputChunks.resize(nInputChunks);
	  memcpy(&inputChunks[0], p, inputChunks.size()*sizeof(ChunkID));
	  p += inputChunks.size()*sizeof(ChunkID);
	  memcpy(&stolen, p, sizeof(bool));
	  p += sizeof(bool);
	  memcpy(&ownerWorkerRank, p, sizeof(int));
	  p += sizeof(int);
	  memcpy(&ancestorDepth, p, sizeof(int));
	  assert(task == NULL);
	  assert(creator_tis == NULL);
	  task        = NULL; // It should already be NULL.
	  creator_tis = NULL; // It should already be NULL.
	}

      }; // end TaskInfoStruct
    }; // end class template Base

    template<typename Service_type> 
      std::string Base<Service_type>::object_type_id() {
      return "TaskSchedulerService";
    }

  }; // end namespace
}; // end namespace

#endif
