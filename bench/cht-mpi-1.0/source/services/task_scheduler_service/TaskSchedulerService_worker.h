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
#ifndef TASKSCHEDULERSERVICE_WORKER_HEADER
#define TASKSCHEDULERSERVICE_WORKER_HEADER
#include "services/Service_worker.h"
#include <list>
#include <fstream>
#include "utilities/Singleton.h"
#include "utilities/cht_utils.h"
#include "services/MPIWrapperInclude.h"
#include "services/output_service/OutputService_worker.h"
#include "services/task_scheduler_service/Task.h"
#include "services/task_scheduler_service/TaskSchedulerService.h"

namespace cht {
  namespace TaskSchedulerService {

    static void* global_thread_func_main_worker(void*);
    static void* global_thread_func_request_handler(void*);
    static void* global_thread_func_data_fetcher(void*);
    static void* global_thread_func_monitoring(void*);
    static void* global_thread_func_for_mutex_lock(void*);
    
    class Worker : 
    public Base<Service::Worker>,
      public Singleton<Worker> {
	friend class Singleton<Worker>;
  

	struct TaskTypeStatistics {
	  int count_execute_leaf;  // number of times this task type has been executed (leaf)
	  int count_execute_nonleaf;  // number of times this task type has been executed (nonleaf)
	  int count_finalize; // number of times this task type has been finalized
	  Threads::ThreadGroup::State_stats stats_execute_leaf[Threads::ThreadGroup::numberOfPossibleStates];
	  Threads::ThreadGroup::State_stats stats_execute_nonleaf[Threads::ThreadGroup::numberOfPossibleStates];
	  Threads::ThreadGroup::State_stats stats_finalize[Threads::ThreadGroup::numberOfPossibleStates];
	  double execute_leaf_working_time_min;
	  double execute_leaf_working_time_max;
	  void clear() {
	    count_execute_leaf  = 0;
	    count_execute_nonleaf  = 0;
	    count_finalize = 0;
	    execute_leaf_working_time_min = 0; // min/max anyway reset first time
	    execute_leaf_working_time_max = 0; // min/max anyway reset first time
	    for (int i = 0; i < Threads::ThreadGroup::numberOfPossibleStates; ++i) {
	      stats_execute_leaf[i].clear();
	      stats_execute_nonleaf[i].clear();
	      stats_finalize[i].clear();
	    }
	  }
	  TaskTypeStatistics() {
	    clear();
	  }
	}; // end struct TaskTypeStatistics

	struct TaskInfoPtr {
	  TaskInfoPtr(TaskInfoStruct const & tis_) : tis(tis_) {}
	  TaskInfoStruct const & tis;
	  bool operator==  ( TaskInfoPtr const & x ) const {
	    return tis == x.tis;
	  }
	  bool operator!=  ( TaskInfoPtr const & x ) const {
	    return !(*this == x);
	  }

	  bool operator<  ( TaskInfoPtr const & x ) const {
	    /* FIXME: find out why compare_complicated() is not working properly, or try to make comparisons faster in some other way, if needed. */
	    // if(compare_complicated(x) != compare_simple(x))
            //   std::cout << "Warning! (compare_complicated(x) != compare_simple(x))" << std::endl;
	    return compare_simple(x);
	  }

	  void print_list(std::list<TaskInfoStruct const *> & list) const {
	    std::list<TaskInfoStruct const *>::iterator it  = list.begin();
	    while(it != list.end()) {
	      std::cout << "   " << (*it)->get_id().str();
	      it++;
	    }
	    std::cout << std::endl;
	  }

	  void get_ancestor_list(TaskInfoStruct const & taskInfo, std::list<TaskInfoStruct const *> & ancestors) const {
	    ancestors.push_front(&taskInfo);
	    TaskInfoStruct const * parentPtr = taskInfo.get_creator_tis();
	    if(!parentPtr)
	      return;
	    return get_ancestor_list(*parentPtr, ancestors);
	  }

	  bool compare_simple ( TaskInfoPtr const & x ) const {
	    if(tis.get_id() == x.tis.get_id())
	      return false;
	    // A stolen task should always be executed as soon as possible, so we check get_stolen() first of all here.
            if(tis.get_stolen() && !x.tis.get_stolen())
	      return false;
	    if(!tis.get_stolen() && x.tis.get_stolen())
	      return true;
	    std::list<TaskInfoStruct const *> ancestors_this;
	    get_ancestor_list(tis, ancestors_this);
	    std::list<TaskInfoStruct const *> ancestors_other;
	    get_ancestor_list(x.tis, ancestors_other);
	    // Now go through ancestor lists together to find the first point where they differ.
	    std::list<TaskInfoStruct const *>::iterator it_this  = ancestors_this .begin();
	    std::list<TaskInfoStruct const *>::iterator it_other = ancestors_other.begin();
	    // First check if already the topmost ancestors differ (in that case the two tasks have no (visible) common ancestor).
	    if((*it_this)->get_id() != (*it_other)->get_id()) {
	      // No (visible) common ancestor.
	      // A stolen task goes before others.
	      if((*it_this)->get_stolen() && !(*it_other)->get_stolen())
		return false;
	      if(!(*it_this)->get_stolen() && (*it_other)->get_stolen())
		return true;
	      // Both stolen, or none of them stolen. In this case, go by ancestorDepth.
	      if((*it_this)->get_ancestorDepth() > (*it_other)->get_ancestorDepth())
		return false;
              if((*it_this)->get_ancestorDepth() < (*it_other)->get_ancestorDepth())
                return true;
	      // Same ancestorDepth. In this case, just compare the ids.
              return ( (*it_this)->get_id() < (*it_other)->get_id() );
	    }
	    // If we get to this point, it means the two tasks have a (visible) common ancestor. Wee loop through the ancestor lists until we find the first point where they differ.
	    while(1) {
	      TaskInfoStruct const * taskInfo_this = *it_this;
	      TaskInfoStruct const * taskInfo_other = *it_other;
	      if(taskInfo_this->get_id() != taskInfo_other->get_id()) {
		// OK, we have found a point where the ancestor lists differ.
		int ownedPercentage_this  = 100;
		if(taskInfo_this ->totInputDataSize > 0)
		  ownedPercentage_this  = (int)((100*taskInfo_this ->currWorkerOwnedSize)/taskInfo_this ->totInputDataSize);
		int ownedPercentage_other = 100;
		if(taskInfo_other->totInputDataSize > 0)
		  ownedPercentage_other = (int)((100*taskInfo_other->currWorkerOwnedSize)/taskInfo_other->totInputDataSize);
		if(ownedPercentage_this != ownedPercentage_other)
		  return ( ownedPercentage_this < ownedPercentage_other );
		return ( taskInfo_this->get_id() < taskInfo_other->get_id() );
	      }
	      it_this++;
	      it_other++;
	      if(it_this == ancestors_this.end() && it_other == ancestors_other.end()) {
		// This should not happen.
		std::cout << "Error, both ancestor lists reached end."
			  << " ancestors_this.size() = " << ancestors_this.size() 
			  << " ancestors_other.size() = " << ancestors_other.size() 
			  << std::endl;
		std::cout << "ancestors_this: ";
		print_list(ancestors_this);
		std::cout << "ancestors_other: ";
		print_list(ancestors_other);
		throw std::runtime_error("Error, both ancestor lists reached end.");
	      }
	      if(it_this == ancestors_this.end())
		return true;
	      if(it_other == ancestors_other.end())
		return false;
	    }
	  }

	  /* FIXME: this routine is wrong and nobody knows how to fix it. */
	  bool compare_complicated ( TaskInfoPtr const & x ) const {
	    std::list<TaskInfoStruct const *> ancestors_this;
	    ancestors_this.push_back(&tis);
	    ancestors_this.push_back(tis.get_creator_tis());
	    std::list<TaskInfoStruct const *> ancestors_other;
	    ancestors_other.push_back(&x.tis);
	    ancestors_other.push_back(x.tis.get_creator_tis());
	    while (ancestors_this.back() != ancestors_other.back()) {
	      if (ancestors_this.back() == NULL || ancestors_other.back() == NULL) {
		if ( ancestors_this.back() == NULL )
		  ancestors_this.pop_back();
		if ( ancestors_other.back() == NULL )
		  ancestors_other.pop_back();
		if ( ancestors_this.back()->get_id() == ancestors_other.back()->get_id() )
		  // This means that one of the tasks was an ancestor to the other!
		  return tis.get_ancestorDepth() < x.tis.get_ancestorDepth();
		return ( ancestors_this.back()->get_id() < ancestors_other.back()->get_id() );
	      }		
	      if ( ancestors_this.back()->get_ancestorDepth() > ancestors_other.back()->get_ancestorDepth() ) 
		ancestors_this.push_back(ancestors_this.back()->get_creator_tis());		
	      else
		ancestors_other.push_back(ancestors_other.back()->get_creator_tis());		
	    } // end while
	    ancestors_this.pop_back();
	    ancestors_other.pop_back();
	    if ( ancestors_this.back()->get_id() == ancestors_other.back()->get_id() )
	      // This means that one of the tasks was an ancestor to the other!
	      return tis.get_ancestorDepth() < x.tis.get_ancestorDepth();
	    return ( ancestors_this.back()->get_id() < ancestors_other.back()->get_id() );
	  }
	  
	}; // end struct TaskInfoPtr

	typedef std::map<std::string, TaskTypeStatistics> TaskStatsMap;
	TaskStatsMap taskStatistics;
	cht::timer timerForEverything;
	bool thisWorkerHasMotherTask;
	double timeWhenTaskWasLastSetAsExecuting;
	double accumulatedIdleTime;
	void sleepAndAddToAccumulatedIdleTime(int microsecondsToSleep);
	void updateStatistics(const TaskInfoStruct & taskInfo,
			      cht::vector<Threads::ThreadGroup::State_stats> const & threadStats_start,
			      cht::vector<Threads::ThreadGroup::State_stats> const & threadStats_stop,
			      std::string execute_or_finalize);
    private:
	OutputService::Worker& OS;

	/* The following member variables exist for both parent and worker,
	 * but for the parent the vaiables are only used in order to send
	 * the values to the workers. */
	int n_worker_threads;  /* Number of worker threads */
	int n_prefetch_tasks; /* Max number of tasks for which to prefetch input chunks. */
	int child_chunk_fetch_level; /* Level of child-chunks to fetch (-1 means disabled). */
	bool stealing_disabled_flag; /* Mainly for debugging. */
	bool steal_only_from_0_flag; /* Mainly for debugging. */
	int status_report_file_interval_seconds; /* Mainly for debugging. */
	int status_report_interval_milliseconds; /* Mainly for debugging. */
	bool do_mutex_lock_checking; /* Mainly for debugging. */

	friend void* global_thread_func_main_worker(void*);
	friend void* global_thread_func_request_handler(void*);
	friend void* global_thread_func_data_fetcher(void*);
	friend void* global_thread_func_monitoring(void*);
	friend void* global_thread_func_for_mutex_lock(void*);

	virtual void start_derived();
	virtual void stop_derived();

	Worker();

	int maxNoOfSimultaneousTasks;
	// The pendingTaskCount_tot and pendingTaskCount_readyToRun variables are used to keep statistics of how often a task is ready to run directly when the ChunkIDs for it are known.
	int pendingTaskCount_tot;
	int pendingTaskCount_readyToRun;	
	// A more optimized map could possibly speed up lookups if needed
	typedef std::map<TaskID, TaskInfoStruct*> TaskListMap;
	typedef std::map<TaskInfoPtr, TaskInfoStruct*> TaskListSortedMap;
	void outputListContents(std::ostream & stream, const char* s, TaskListMap & theList);
	void writeStatusReportFile(int my_rank);
	void doStatusReport();

	// Statistics about stealing operations
	int transactionCount; /* Counts how many transactions have been performed. Gives an indication of how much work has been done since previous steal. */
	int noOfStealsFromThisWorker; /* Number of times other workers have successfully stolen from this worker. */
	int noOfStealsByThisWorker; /* Number of times this worker has successfully stolen from other workers. */
	struct StealInfoStruct {
	  int victimRank;
	  int transactionCountWhenStealOccurred;
	  int ancestorDepthOfStolenTask;
	  double wallTimeWhenStealOccurred;
	  bool stealOnlyOwnStuff;
	  std::string taskTypeString;
	StealInfoStruct() : victimRank(-1), transactionCountWhenStealOccurred(-1), ancestorDepthOfStolenTask(-1), wallTimeWhenStealOccurred(-1), stealOnlyOwnStuff(false) { }
	};
	std::list<StealInfoStruct> stealStatisticsList;

	/**************************************************************
	  Lists of tasks in different states.
	***************************************************************/

	/* "Waiting before transaction". This is the state of a task
	   that has been registered but is depending on the output
	   chunks of some other tasks, so it must wait until those
	   other tasks have produced their output chunks. Another way
	   of saying this is that all input chunk ids for the task are
	   not known yet. */
	TaskListMap list_waiting_before_transaction;

	/* - "Pending for transaction". This state means that all
	   input chunk ids for the task are known. While a task is in
	   this state, data may be fetched for it since the input
	   chunk ids are known. Also, it may be speculatively
	   executed. After a task has been speculatively executed, we
	   see if that was a leaf task or not. If it was a leaf task,
	   create and copy operations are started and the task is then
	   moved directly to list_waiting_before_finalize. A task in
	   this state can also be stolen. */
	TaskListSortedMap list_pending_for_transaction;
	TaskListMap helper_list_pending_for_transaction;

	/* "Doing transaction non-leaf, part 1". This state means that
	   createChunk and copyChunk operations are performed being
	   performed, for a non-leaf task. */
	TaskListMap list_doing_transaction_nonleaf_1;

	/* "Doing transaction non-leaf, part 2". This state means that
	   the registerTask operations are performed, so that any
	   child tasks are really registered. */
	TaskListMap list_doing_transaction_nonleaf_2;

	/* "Doing transaction leaf". This state means that createChunk
	   and copyChunk operations are performed being performed, for
	   a leaf task. */
	TaskListMap list_doing_transaction_leaf;

	/* "Waiting before finalize". In this state, we are waiting
	   for any child tasks to complete and also (for leaf tasks)
	   for createChunk and copyChunk operations to be
	   completed. */
	TaskListMap list_waiting_before_finalize;
	TaskListMap list_pending_for_finalize;
	TaskListMap list_running_finalize;
	TaskListMap list_finished;
	TaskListMap list_stolen;
	// Mutex to avoid conflicting simultaneous access.
	Threads::Mutex mutex;
	// The LockMutex function takes a string argument used for debugging.
	bool mutexIsLocked;
	std::string mutexLockDebugStr;
	void LockMutex(std::string debugInfoStr);
	void UnlockMutex();
	// Flag set when time to finish.
	bool finishFlagForThreadCommunicationWorker;
	bool steal_attempt_in_progress;
	bool currStealAttemptIsOnlyForOwnStuff; // Relevant only when steal_attempt_in_progress=true
	// Steal attempt statistics.
	double steal_attempt_start_time;
	cht::simple_statistics stealAttemptStatistics;
	cht::simple_statistics stealAttemptSendStatistics;
	cht::simple_statistics RequestHandlerFuncStatisticsParent;
	cht::simple_statistics RequestHandlerFuncStatisticsWorker;
  
	// Function counting all tasks in all lists.
	int GetTotNoOfTasks();
	// Function to check if a task exists in any of the lists.
	bool ExistsInMap(const TaskListMap & map, TaskID id);
	int CountChildTasks(const TaskListMap & map, TaskID id);
	int CountChildTasks(TaskID id);
	void AddTaskToPendingForTransactionLists(TaskInfoStruct* taskInfo);
	void RemoveTaskFromPendingForTransactionLists(TaskInfoStruct* taskInfo);
	void ReinsertTaskToPendingForTransactionLists(TaskInfoStruct* taskInfo);
	void RemoveTaskFromWaitLists(TaskID tid);
	void RemoveTaskFromChildListOfCreator(TaskInfoStruct* taskInfo);
	bool UpdateLists();
	void GetHelperRanks(const TaskListMap & taskListMap, std::list<int> & resultList);
	Threads::ThreadGroup* threadGroupWorkers; /**< Worker thread group pointer. @warning Not protected by mutex. */
	Threads::Thread* threadRequestHandler; /**< Request handler thread pointer. @warning Not protected by mutex. */
	Threads::Thread* threadDataFetcher; /**< Data fetcher thread pointer. @warning Not protected by mutex. */
	Threads::Thread* threadMonitoring; /**< Monitoring thread pointer. @warning Not protected by mutex. */
	void HandleMessageFromParent(int tag, int msgSize, char* recvBufParent, MPI_Comm* comm_to_parent, MPI_Request & request);
	void HandleStealAttempt(int rank_who_wants_to_steal, bool stealOnlyFromWithinOwnTasks, MPI_Comm* comm_worker_world);
	void HandleStolenTaskMessage(const char* recvBufWorker, int msgSize, int victimRank);
	void HandleStolenTaskFinishedMessage(const char* recvBufWorker, int msgSize);
	void WorkerFunc();
	void RequestHandlerFunc();
	void DataFetcherFunc();
	void MonitoringFunc();
	void MutexLockThreadFunc();
	Threads::Cond data_fetcher_cond;

	void setTaskInputChunkIDs( TaskInfoStruct & tis);

	bool getChildChunksRecursive(const Chunk & chunk, int level, std::vector<ChunkID> & childChunkIDVec);
	bool getChildChunksForLevel(const TaskInfoStruct & taskInfo, int level, std::vector<ChunkID> & childChunkIDVec);

    public:
        void call_mpi_abort();
  
	// Parent uses default version
	void resetStatistics();
	void reportStatistics(std::string messageHeader);

	TaskID registerTaskStep1(std::string taskTypeID);
	void registerTaskStep2(TaskID taskID,
			       TaskID creatorTaskID,
			       std::vector<cht::ID> const & input_and_dependencies);
	
	ChunkID getTaskResult(TaskID const & calling_tid, 
			      TaskID const & tid);
	ChunkID getTaskResult(TaskID const & tid);

	bool checkIfTaskIsDoingFinalize(TaskID id);

	void addTemporaryChunk(TaskID const & owner_task,
			       ID const & id);
  
      };

  }; // end namespace
}; // end namespace

#endif
