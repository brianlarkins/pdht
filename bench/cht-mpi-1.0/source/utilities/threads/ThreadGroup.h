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
#ifndef THREADGROUP_HEADER
#define THREADGROUP_HEADER
#include <map>
#include <queue>
#include "../cht_vector.h"
#include "Cond.h"
#include "Mutex.h"
#include "Thread.h"
namespace cht {
  namespace Threads {
    
    /** Class used to start, synchronize, and stop a group of threads
	that (supposedly) work together to solve a task.  A key
	feature is that the threads can and should report their
	'thread state' to the thread group so that an educated
	decision can be made when to activate and deactivate different
	threads.
    */
    class ThreadGroup {
    public:
      struct State_stats {
	int counter;
	double time_acc;
	void clear() {
	  counter = 0;
	  time_acc = 0;
	}
        State_stats() : counter(0), time_acc(0) {}
      };

      static const int beforeAttach = 0;
      static const int afterDetach = 1;
      static const int doingBlockingOperation = 2;
      static const int readyAfterBlockingOperation = 3;
      static const int working = 4;
      static const int waiting = 5;
      static const int numberOfPossibleStates = 6;

      static std::string get_thread_state_name(int ts);      

    private:      
      class ThreadInfo {
      private:
	int myState; 
	double timer_checkpoint;
	State_stats stats[numberOfPossibleStates];
      public:
	bool highPriority;
	int blockingOperationCounter;
	Cond cond;
	
	void changeState(int state); /**< Change thread state. @warning thread group mutex must be locked before call!! */
	void updateStats(); 	/**< Update statistics. @warning thread group mutex must be locked before call!!     */
	void getStats(cht::vector<State_stats> & threadStats); 	/**< Aquire thread group statistics
								   @warning thread group mutex must be locked before call!! */
	int getState() const { return myState; }
	ThreadInfo();
      }; // end class ThreadInfo


      typedef std::map< pthread_t, ThreadInfo* > ThreadInfoMap;
      typedef std::queue<pthread_t> ThreadQueue;
      int const n_total;
      int n_max_working;
      cht::vector<Thread*> threads;
      typedef void* (*threadFun)(void*); 
      threadFun tfun;
      void* arg;
      ThreadInfoMap infoMap;
      ThreadQueue readyAfterBlockingQueue;
      Mutex mutex_group;
      static void* thread_func(void* arg);
      void attach();
      void detach();
      int n_working();               /**< Number of working threads. @warning thread group mutex must be locked before call! */
      void tryToActivateOneWorker(); /**< Activate another worker. @warning thread group mutex must be locked before call!!   */
    public:
      ThreadGroup(std::string identifier, threadFun tfun, void* arg, int n_total, int n_max_working);
      ~ThreadGroup();
      void yield();       
      void startingBlockingOperation();
      void blockingOperationFinished();
      void setPriorityHigh();
      void setPriorityNormal();
      void getThreadStats(cht::vector<State_stats> & threadStats);
      void increaseNMaxWorking();
      bool decreaseNMaxWorking();
      int getNMaxWorking();
    };
    
  } // end namespace Threads
} // end namespace cht
#endif
