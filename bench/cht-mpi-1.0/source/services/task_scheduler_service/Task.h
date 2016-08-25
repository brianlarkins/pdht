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
#ifndef TASK_HEADER
#define TASK_HEADER

#include <string>
#include <sstream>
#include "services/chunk_obj_service/ChunkObject.h"
#include "utilities/cht_typelist.h"
#include "utilities/cht_static_check.h"
#include "utilities/cht_obj_factory.h"

namespace cht {
  
  enum LifeSpan {persistent, temporary};

  // FIXME: move TaskID to its own source file?
  struct TaskID {
    int ownerRank;
    int localID;
    int taskTypeID;
  TaskID() : ownerRank(-1), localID(-1), taskTypeID(-1) {}
  TaskID(int ownerRank_, int localID_, int taskTypeID_) : 
    ownerRank(ownerRank_), localID(localID_), taskTypeID(taskTypeID_) {}
    bool operator<  ( const TaskID & x ) const {
      if(x.ownerRank == ownerRank)
	return (localID < x.localID);
      return (ownerRank < x.ownerRank);
    }
    bool operator>  ( const TaskID & x ) const {
      if(x.ownerRank == ownerRank)
	return (localID > x.localID);
      return (ownerRank > x.ownerRank);
    }
    bool operator==  ( const TaskID & x ) const {
      if(x.ownerRank == ownerRank && x.localID == localID)
	return true;
      return false;
    }
    bool operator!=  ( const TaskID & x ) const {
      if(x.ownerRank != ownerRank || x.localID != localID)
	return true;
      return false;
    }
    std::string str() const {
      std::stringstream ss;
      ss << "(" << ownerRank << "/" << localID << ")";
      return ss.str();
    }
  };

  static const TaskID TASK_ID_NULL;

  class ID {
  private:
    ChunkID chID;
    TaskID  taID;
  public:
    ID()
      : chID(CHUNK_ID_NULL), taID(TASK_ID_NULL) {}
    ID(ChunkID const & chID_)
      : chID(chID_), taID(TASK_ID_NULL) {}
    ID(TaskID const & taID_)
      : chID(CHUNK_ID_NULL), taID(taID_) {}
    ID & operator=(ChunkID const & chunk) {
      chID = chunk;
      taID = TASK_ID_NULL;
      return *this;
    }
    ID & operator=(TaskID const & task) {
      chID = CHUNK_ID_NULL;
      taID = task;
      return *this;
    }
    bool operator==(ID const & other) const {
      return (chID == other.chID) && (taID == other.taID);
    }
    bool is_taskID() const {return taID != TASK_ID_NULL;}
    TaskID get_taskID() const {return taID;}
    bool is_chunkID() const {return chID != CHUNK_ID_NULL;}
    ChunkID get_chunkID() const {return chID;}
    bool is_null() const {return !is_taskID() && !is_chunkID();}
  };

  class Task
  {
  private:
    TaskID myID;
  protected:
    template<typename TaskType>
      TaskID registerTask(std::vector<cht::ID> const & input_and_dependencies, 
			  LifeSpan lifespan = temporary);
    template<typename TaskType>
      TaskID registerTask(cht::ID const id1, 
			  LifeSpan lifespan = temporary);
    template<typename TaskType>
      TaskID registerTask(cht::ID const id1, 
			  cht::ID const id2, 
			  LifeSpan lifespan = temporary);
    template<typename TaskType>
      TaskID registerTask(cht::ID const id1, 
			  cht::ID const id2, 
			  cht::ID const id3, 
			  LifeSpan lifespan = temporary);
    template<typename TaskType>
      TaskID registerTask(cht::ID const id1, 
			  cht::ID const id2, 
			  cht::ID const id3, 
			  cht::ID const id4, 
			  LifeSpan lifespan = temporary);
    template<typename TaskType>
      TaskID registerTask(cht::ID const id1, 
			  cht::ID const id2, 
			  cht::ID const id3, 
			  cht::ID const id4, 
			  cht::ID const id5, 
			  LifeSpan lifespan = temporary);
    template<typename TaskType>
      TaskID registerTask(cht::ID const id1, 
			  cht::ID const id2, 
			  cht::ID const id3, 
			  cht::ID const id4, 
			  cht::ID const id5, 
			  cht::ID const id6, 
			  LifeSpan lifespan = temporary);
    template<typename TaskType>
      TaskID registerTask(cht::ID const id1, 
			  cht::ID const id2, 
			  cht::ID const id3, 
			  cht::ID const id4, 
			  cht::ID const id5, 
			  cht::ID const id6, 
			  cht::ID const id7, 
			  LifeSpan lifespan = temporary);
    template<typename TaskType>
      TaskID registerTask(cht::ID const id1, 
			  cht::ID const id2, 
			  cht::ID const id3, 
			  cht::ID const id4, 
			  cht::ID const id5, 
			  cht::ID const id6, 
			  cht::ID const id7, 
			  cht::ID const id8, 
			  LifeSpan lifespan = temporary);

    ChunkID getTaskResult(TaskID const & tid) const;
    template<typename Tobj>
      ChunkID registerChunk(Tobj const * objPtr, LifeSpan lifespan = temporary);
    ChunkID copyChunk(ChunkID cid);
    void deleteChunk(ChunkID cid) const;
    void output(std::string message, Output::Priority prio = Output::Info) const;

    cht::ChunkID getInputChunkID(cht::Chunk const &) const;
  private:
    TaskID registerTask(std::string taskTypeID,
    			std::vector<cht::ID> const & input_and_dependencies);
    struct task_to_register {
      TaskID id;
      std::string taskTypeID;
      std::vector<cht::ID> input_and_dependencies;
    };
    std::list<task_to_register> tasksToRegister;
    TaskID addTaskToRegisterLater(std::string taskTypeID, std::vector<cht::ID> const & input_and_dependencies);
    void addTemporaryChunk(ID const & id) const;
    ChunkID registerChunk(Chunk const * tmp_ptr, 
			  std::string class_id_str);
    std::vector<cht::shared_ptr<Chunk const> > inputChunks;
    std::vector<ChunkID> inputChunkIDs;    
    // Note: Use of the inputChunksAndChunkIDs vector relies on that
    // the inputChunks and inputChunkIDs vectors are populated. Its
    // elements are supposed to point to the elements in those
    // vectors.
    std::vector<BaseObj const * > inputChunksAndChunkIDs; 
    bool fallback; 
    /* Stuff saved for registerChunk operations. */
    std::list<ChunkID> chunk_ids_for_registerChunk_operations;
    std::list<Chunk const *> chunk_ptrs_for_registerChunk_operations;
    /* Stuff saved for copyChunk operations. */
    std::list<ChunkID> chunk_ids_new_for_copyChunk_operations;
    std::list<ChunkID> chunk_ids_old_for_copyChunk_operations;
  public:
    static std::string object_base_type_id() {
      return "Task";
    }
    virtual cht::ID internal_execute() = 0;
    // Inherited classes are required to define
    // static std::string get_class_id();

    // Note: a virtual destructor is needed here so that instances of 
    // task classes derived from this base class can be deleted using 
    // a pointer to the base class.
    virtual ~Task() { }
    Task() :fallback(false) {
    }
    virtual cht::ID internal_execute(TaskID myID, bool fallback_);
    /** Functions to be called in task scheduler. */
    void internal_setInputChunks(std::vector<cht::shared_ptr<Chunk const> > const & inpChunks);
    void internal_setInputChunkIDs(std::vector<ChunkID> const & inpChunkIDs);
    void internal_setinputChunksAndChunkIDs(std::vector<BaseObj const *> const & inpChunksAndChunkIDs);
    void internal_clearInputVectors();
    int internal_getNoOfTasksToRegister();
    void internal_startRegisterAndCopyChunkOps();
    bool internal_allRegisterAndCopyChunkOpsHaveFinished();
    void internal_registerTasksInList();
    void internal_get_child_task_ids(std::list<TaskID> & list) const;

  protected:
   virtual cht::ID execute(cht::ChunkID const & c1);
   virtual cht::ID execute(cht::ChunkID const & c1, cht::ChunkID const & c2);
   virtual cht::ID execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3);
   virtual cht::ID execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3, 
			   cht::ChunkID const & c4);
   virtual cht::ID execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3, 
			   cht::ChunkID const & c4, cht::ChunkID const & c5);
   virtual cht::ID execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3, 
			   cht::ChunkID const & c4, cht::ChunkID const & c5, cht::ChunkID const & c6);
   virtual cht::ID execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3, 
			   cht::ChunkID const & c4, cht::ChunkID const & c5, cht::ChunkID const & c6, 
			   cht::ChunkID const & c7);
   virtual cht::ID execute(cht::ChunkID const & c1, cht::ChunkID const & c2, cht::ChunkID const & c3, 
			   cht::ChunkID const & c4, cht::ChunkID const & c5, cht::ChunkID const & c6, 
			   cht::ChunkID const & c7, cht::ChunkID const & c8);

    template<int n_params>
      struct Exec;
    
  }; // end class Task

  template<typename TaskType>
    TaskID Task::registerTask(std::vector<cht::ID> const & input_and_dependencies, LifeSpan lifespan) {
    std::string taskTypeID = TaskType::get_class_id();
    TaskID tid = registerTask(taskTypeID, input_and_dependencies);
    if (lifespan == temporary)
      addTemporaryChunk(tid);
    return tid;
  }

    template<typename TaskType>
      TaskID Task::registerTask(cht::ID const id1, 
				LifeSpan lifespan) {
      CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 1, 
			    WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
      std::vector<cht::ID> input_and_dependencies(1);
      input_and_dependencies[0] = id1;
      return registerTask<TaskType>(input_and_dependencies, lifespan);
    }
    template<typename TaskType>
      TaskID Task::registerTask(cht::ID const id1, 
				cht::ID const id2, 
				LifeSpan lifespan) {
      CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 2, 
			    WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
      std::vector<cht::ID> input_and_dependencies(2);
      input_and_dependencies[0] = id1; 
      input_and_dependencies[1] = id2;
      return registerTask<TaskType>(input_and_dependencies, lifespan);      
    }
    template<typename TaskType>
      TaskID Task::registerTask(cht::ID const id1, 
				cht::ID const id2, 
				cht::ID const id3, 
				LifeSpan lifespan) {
      CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 3, 
			    WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
      std::vector<cht::ID> input_and_dependencies(3);
      input_and_dependencies[0] = id1; 
      input_and_dependencies[1] = id2;
      input_and_dependencies[2] = id3;
      return registerTask<TaskType>(input_and_dependencies, lifespan);      
    }
    template<typename TaskType>
      TaskID Task::registerTask(cht::ID const id1, 
				cht::ID const id2, 
				cht::ID const id3, 
				cht::ID const id4, 
				LifeSpan lifespan) {
      CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 4, 
			    WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
      std::vector<cht::ID> input_and_dependencies(4);
      input_and_dependencies[0] = id1; 
      input_and_dependencies[1] = id2;
      input_and_dependencies[2] = id3;
      input_and_dependencies[3] = id4;
      return registerTask<TaskType>(input_and_dependencies, lifespan);      
    }
    template<typename TaskType>
      TaskID Task::registerTask(cht::ID const id1, 
				cht::ID const id2, 
				cht::ID const id3, 
				cht::ID const id4, 
				cht::ID const id5, 
				LifeSpan lifespan) {
      CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 5, 
			    WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
      std::vector<cht::ID> input_and_dependencies(5);
      input_and_dependencies[0] = id1; 
      input_and_dependencies[1] = id2;
      input_and_dependencies[2] = id3;
      input_and_dependencies[3] = id4;
      input_and_dependencies[4] = id5;
      return registerTask<TaskType>(input_and_dependencies, lifespan);      
    }
    template<typename TaskType>
      TaskID Task::registerTask(cht::ID const id1, 
				cht::ID const id2, 
				cht::ID const id3, 
				cht::ID const id4, 
				cht::ID const id5, 
				cht::ID const id6, 
				LifeSpan lifespan) {
      CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 6, 
			    WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
      std::vector<cht::ID> input_and_dependencies(6);
      input_and_dependencies[0] = id1;
      input_and_dependencies[1] = id2;
      input_and_dependencies[2] = id3;
      input_and_dependencies[3] = id4;
      input_and_dependencies[4] = id5;
      input_and_dependencies[5] = id6;
      return registerTask<TaskType>(input_and_dependencies, lifespan);
    }
    template<typename TaskType>
      TaskID Task::registerTask(cht::ID const id1, 
				cht::ID const id2, 
				cht::ID const id3, 
				cht::ID const id4, 
				cht::ID const id5, 
				cht::ID const id6, 
				cht::ID const id7, 
				LifeSpan lifespan) {
      CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 7, 
			    WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
      std::vector<cht::ID> input_and_dependencies(7);
      input_and_dependencies[0] = id1;
      input_and_dependencies[1] = id2;
      input_and_dependencies[2] = id3;
      input_and_dependencies[3] = id4;
      input_and_dependencies[4] = id5;
      input_and_dependencies[5] = id6;
      input_and_dependencies[6] = id7;
      return registerTask<TaskType>(input_and_dependencies, lifespan);
    }
    template<typename TaskType>
      TaskID Task::registerTask(cht::ID const id1, 
				cht::ID const id2, 
				cht::ID const id3, 
				cht::ID const id4, 
				cht::ID const id5, 
				cht::ID const id6, 
				cht::ID const id7, 
				cht::ID const id8, 
				LifeSpan lifespan) {
      CHT_STATIC_CHECK_TRUE(TaskType::inputTypes::length == 8, 
			    WRONG_NUMBER_OF_ARGUMENTS_TO_TASK);
      std::vector<cht::ID> input_and_dependencies(8);
      input_and_dependencies[0] = id1;
      input_and_dependencies[1] = id2;
      input_and_dependencies[2] = id3;
      input_and_dependencies[3] = id4;
      input_and_dependencies[4] = id5;
      input_and_dependencies[5] = id6;
      input_and_dependencies[6] = id7;
      input_and_dependencies[7] = id8;
      return registerTask<TaskType>(input_and_dependencies, lifespan);
    }

 
  template<typename Tobj>
    ChunkID Task::registerChunk(Tobj const * objPtr, LifeSpan lifespan) {
    ChunkID cid = registerChunk(objPtr, Tobj::get_class_id());    
    if (lifespan == temporary)
      addTemporaryChunk(cid);
    return cid;
  }

  namespace internal {
    bool registerArgTypes(std::string & task_type_str, // FIXME should these params be const?
			  std::list<std::string> & task_args, // FIXME should these params be const?
			  std::string & task_output_type_str); // FIXME should these params be const?
  }

  template<typename TaskType>
    bool registerTaskType() {
    CHT_TYPECHECK_NOTEQUAL(typename TaskType::outputType, ChunkID, OUTPUT_TYPE_IS_CHUNKID); 
    bool const task_reg = 
      cht::obj_factory<Task>::instance().registerObjectType<TaskType>();
    std::string task_type_str = TaskType::get_class_id();
    std::list<std::string> task_args;
    TaskType::inputTypes::get_type_strings(task_args);
    std::string task_output_type_str = TaskType::outputType::get_class_id();
    bool const args_reg = internal::registerArgTypes(task_type_str, task_args, task_output_type_str);
    return task_reg && args_reg;
  }


  template<>
    struct Task::Exec<1> {
    template<typename TaskType>
      static cht::ID exec(TaskType * task_ptr) {
      Task* base_ptr = task_ptr;
      if (base_ptr->fallback)
	return base_ptr->execute( base_ptr->inputChunkIDs[0] );
      return task_ptr->execute( TaskType::inputTypes::template getElement<0>( base_ptr->inputChunksAndChunkIDs ) );
    }
  };
  template<>
    struct Task::Exec<2> {
    template<typename TaskType>
      static cht::ID exec(TaskType * task_ptr) {
      Task* base_ptr = task_ptr;
      if (base_ptr->fallback)
	return base_ptr->execute( base_ptr->inputChunkIDs[0],
				  base_ptr->inputChunkIDs[1]);
      return task_ptr->execute( TaskType::inputTypes::template getElement<0>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<1>( base_ptr->inputChunksAndChunkIDs ) );
    }
  };
  template<>
    struct Task::Exec<3> {
    template<typename TaskType>
      static cht::ID exec(TaskType * task_ptr) {
      Task* base_ptr = task_ptr;
      if (base_ptr->fallback)
	return base_ptr->execute( base_ptr->inputChunkIDs[0],
				  base_ptr->inputChunkIDs[1],
				  base_ptr->inputChunkIDs[2]);
      return task_ptr->execute( TaskType::inputTypes::template getElement<0>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<1>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<2>( base_ptr->inputChunksAndChunkIDs ));
    }
  };
  template<>
    struct Task::Exec<4> {
    template<typename TaskType>
      static cht::ID exec(TaskType * task_ptr) {
      Task* base_ptr = task_ptr;
      if (base_ptr->fallback)
	return base_ptr->execute( base_ptr->inputChunkIDs[0],
				  base_ptr->inputChunkIDs[1],
				  base_ptr->inputChunkIDs[2],
				  base_ptr->inputChunkIDs[3]);
      return task_ptr->execute( TaskType::inputTypes::template getElement<0>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<1>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<2>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<3>( base_ptr->inputChunksAndChunkIDs ));
    }
  };
  template<>
    struct Task::Exec<5> {
    template<typename TaskType>
      static cht::ID exec(TaskType * task_ptr) {
      Task* base_ptr = task_ptr;
      if (base_ptr->fallback)
	return base_ptr->execute( base_ptr->inputChunkIDs[0],
				  base_ptr->inputChunkIDs[1],
				  base_ptr->inputChunkIDs[2],
				  base_ptr->inputChunkIDs[3],
				  base_ptr->inputChunkIDs[4]);
      return task_ptr->execute( TaskType::inputTypes::template getElement<0>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<1>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<2>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<3>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<4>( base_ptr->inputChunksAndChunkIDs ));
    }
  };
  template<>
    struct Task::Exec<6> {
    template<typename TaskType>
      static cht::ID exec(TaskType * task_ptr) {
      Task* base_ptr = task_ptr;
      if (base_ptr->fallback)
	return base_ptr->execute( base_ptr->inputChunkIDs[0],
				  base_ptr->inputChunkIDs[1],
				  base_ptr->inputChunkIDs[2],
				  base_ptr->inputChunkIDs[3],
				  base_ptr->inputChunkIDs[4],
				  base_ptr->inputChunkIDs[5]);
      return task_ptr->execute( TaskType::inputTypes::template getElement<0>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<1>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<2>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<3>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<4>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<5>( base_ptr->inputChunksAndChunkIDs ));
    }
  };
  template<>
    struct Task::Exec<7> {
    template<typename TaskType>
      static cht::ID exec(TaskType * task_ptr) {
      Task* base_ptr = task_ptr;
      if (base_ptr->fallback)
	return base_ptr->execute( base_ptr->inputChunkIDs[0],
				  base_ptr->inputChunkIDs[1],
				  base_ptr->inputChunkIDs[2],
				  base_ptr->inputChunkIDs[3],
				  base_ptr->inputChunkIDs[4],
				  base_ptr->inputChunkIDs[5],
				  base_ptr->inputChunkIDs[6]);
      return task_ptr->execute( TaskType::inputTypes::template getElement<0>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<1>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<2>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<3>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<4>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<5>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<6>( base_ptr->inputChunksAndChunkIDs ));
    }
  };
  template<>
    struct Task::Exec<8> {
    template<typename TaskType>
      static cht::ID exec(TaskType * task_ptr) {
      Task* base_ptr = task_ptr;
      if (base_ptr->fallback)
	return base_ptr->execute( base_ptr->inputChunkIDs[0],
				  base_ptr->inputChunkIDs[1],
				  base_ptr->inputChunkIDs[2],
				  base_ptr->inputChunkIDs[3],
				  base_ptr->inputChunkIDs[4],
				  base_ptr->inputChunkIDs[5],
				  base_ptr->inputChunkIDs[6],
				  base_ptr->inputChunkIDs[7]);
      return task_ptr->execute( TaskType::inputTypes::template getElement<0>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<1>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<2>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<3>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<4>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<5>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<6>( base_ptr->inputChunksAndChunkIDs ),
				TaskType::inputTypes::template getElement<7>( base_ptr->inputChunksAndChunkIDs ));
    }
  };


} // end namespace cht
#endif
