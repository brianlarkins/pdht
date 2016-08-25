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
#include <iostream>
#include "../cht_time.h"
#include "ThreadGroup.h"
#include "Manager.h"

namespace cht {
  namespace Threads {

    std::string ThreadGroup::get_thread_state_name(int ts) {
      switch (ts) {
      case beforeAttach:
	return "before-attach             ";
      case afterDetach:
	return "after-detach              ";
      case doingBlockingOperation:
	return "doing-blocking-op         ";
      case readyAfterBlockingOperation:
	return "ready-after-blocking-op   ";
      case working:
	return "working                   ";
      case waiting:
	return "waiting                   ";
      };
      throw std::runtime_error("Invalid thread state");
    }
    
    void* ThreadGroup::thread_func(void* arg) {
      ThreadGroup* p = (ThreadGroup*) arg;
      try {
	p->attach();
	Manager::instance().registerThreadGroup(p);
	(p->tfun)(p->arg);
	Manager::instance().unregisterThreadGroup();
	p->detach();
      }
      catch (std::runtime_error e) {
	std::cerr << "Error! Exception std::runtime_error caught in "
	  "cht::threads::thread_func()." << std::endl 
		  << "what() = '" << e.what() << "'" << std::endl;
      }
      catch ( ... ) {
	std::cerr << "Error! Exception caught in cht::threads::thread_func()." 
		  << std::endl;
      }
      return NULL;
    }

    void ThreadGroup::ThreadInfo::changeState(int state) {
      updateStats(); // Note: important to update before changing state.
      stats[state].counter++;
      myState = state;
    }
    void ThreadGroup::ThreadInfo::updateStats() {
      stats[myState].time_acc += (cht::get_wall_seconds() - timer_checkpoint);
      timer_checkpoint = cht::get_wall_seconds();
    }
    void ThreadGroup::ThreadInfo::getStats
    (cht::vector<State_stats> & threadStats) {
      // thread group mutex must be locked when calling this function!!
      updateStats();
      threadStats.setsize(numberOfPossibleStates);
      for (int i = 0; i < numberOfPossibleStates; ++i) 
	threadStats[i] = stats[i];
    }
    ThreadGroup::ThreadInfo::ThreadInfo() 
      : myState(beforeAttach), blockingOperationCounter(0),
	highPriority(false) {
      for (int i = 0; i < numberOfPossibleStates; ++i) {
	stats[i].clear();
      }
      stats[beforeAttach].counter++;
      timer_checkpoint = cht::get_wall_seconds();
    }
    

    ThreadGroup::ThreadGroup(std::string identifier, 
			     threadFun tfun, void* arg, 
			     int n_total, int n_max_working) 
      : n_total(n_total), n_max_working(n_max_working), 
	tfun(tfun), arg(arg), threads(n_total) {
      // We need to lock here since 'this' is passed to the Thread constructor 
      mutex_group.lock(); 
      for (int i = 0; i < n_total; ++i)
	threads[i] = new Thread(identifier, thread_func, this);
      mutex_group.unlock();
    }

    ThreadGroup::~ThreadGroup() {
      mutex_group.lock();
      n_max_working = n_total;
      pthread_t threadToActivate;
      // Try to activate thread in waiting state
      ThreadInfoMap::iterator it;
      for ( it=infoMap.begin() ; it != infoMap.end(); it++ ) {
	if ( (*it).second->getState() == waiting ) {
	  threadToActivate = (*it).first;
	  infoMap[threadToActivate]->changeState(working);
	  infoMap[threadToActivate]->cond.signal();
	}
      } // end for
      mutex_group.unlock();
      for (int i = 0; i < n_total; ++i) 
	delete threads[i];      
    }
    
    void ThreadGroup::attach() {
      pthread_t myID = pthread_self();
      mutex_group.lock();
      if (infoMap.find(myID) != infoMap.end())
	throw std::runtime_error("attach() called for thread already "
				 "in thread group");      
      infoMap[myID] = new ThreadInfo();
      if ( n_working() < n_max_working ) {
	infoMap[myID]->changeState(working);
      }
      else {
	infoMap[myID]->changeState(waiting);
	while(infoMap[myID]->getState() == waiting) {
	  infoMap[myID]->cond.wait(mutex_group);
	}
      }
      mutex_group.unlock();
      // Should be in working state here
      mutex_group.lock();
      if ( infoMap[myID]->getState() != working) 
	throw std::runtime_error("Thread not in working state when "
				 "leaving attach function");
      mutex_group.unlock();      
    }

    void ThreadGroup::detach() {
      pthread_t myID = pthread_self();
      mutex_group.lock();      
      if ( infoMap[myID]->getState() != working)
	throw std::runtime_error("Thread group detach called for thread in"
				 " nonworking state");
      infoMap[myID]->changeState(afterDetach);
      delete infoMap[myID];
      if ( infoMap.erase(myID) != 1 )
	throw std::runtime_error("Thread group detach: error in erase");
      if (n_working() >= n_max_working)
	throw std::runtime_error("Thread group detach called with "
				 "n_working() >= n_working_max after erase");
      tryToActivateOneWorker();
      mutex_group.unlock();      
    }

    int ThreadGroup::n_working() {
      ThreadInfoMap::iterator it;
      int count = 0;
      for ( it=infoMap.begin() ; it != infoMap.end(); it++ ) {
	if ( (*it).second->getState() == working)
	  ++count;
      }
      return count;
    }
    
    void ThreadGroup::tryToActivateOneWorker() {
      // thread group mutex already locked here
      pthread_t threadToActivate;
      // First, try to activate thread in ready after blocking state
      if ( !readyAfterBlockingQueue.empty() ) {
	threadToActivate = readyAfterBlockingQueue.front();
	readyAfterBlockingQueue.pop();
      }
      else {
	// No threads in ready after blocking state, try
	// to activate thread in waiting state
	ThreadInfoMap::iterator it;
	for ( it=infoMap.begin() ; it != infoMap.end(); it++ ) {
	  if ( (*it).second->getState() == waiting ) {
	    threadToActivate = (*it).first;
	    break;
	  }
	} // end for
	if (it == infoMap.end()) 
	  // No threads to activate
	  return; 
      }
      // There is a thread to activate
      infoMap[threadToActivate]->changeState(working);
      infoMap[threadToActivate]->cond.signal();
    }

    void ThreadGroup::yield() {
      pthread_t myID = pthread_self();
      mutex_group.lock();
      if (infoMap.find(myID) == infoMap.end())
	throw std::runtime_error("Thread group yield() called"
				 " for thread not in thread group");      
      if ( infoMap[myID]->getState() != working )
	throw std::runtime_error("Thread group yield() called "
				 "for thread in nonworking state");    
      if ( !readyAfterBlockingQueue.empty() ) {
	// Activate another thread
	pthread_t threadToActivate = readyAfterBlockingQueue.front();
	readyAfterBlockingQueue.pop();
	infoMap[threadToActivate]->changeState(working);
	infoMap[threadToActivate]->cond.signal();
	// Go to waiting state
	infoMap[myID]->changeState(waiting);
	while(infoMap[myID]->getState() == waiting) {
	  infoMap[myID]->cond.wait(mutex_group);
	}
      }
      // If queue is empty, there is no thread to activate, 
      // just continue "working"
      mutex_group.unlock();
    }

    void ThreadGroup::startingBlockingOperation() {
      pthread_t myID = pthread_self();
      mutex_group.lock();
      if (infoMap.find(myID) == infoMap.end())
	throw std::runtime_error("Thread group startingBlockingOperation called"
				 " for thread not in thread group");
      if (infoMap[myID]->highPriority) {
	/* Do nothing if thread has high priority. */
	mutex_group.unlock();
	return;
      }
      infoMap[myID]->blockingOperationCounter++;
      /* Note: only change state if blockingOperationCounter is 1 here. */
      if(infoMap[myID]->blockingOperationCounter == 1) {
	infoMap[myID]->changeState(doingBlockingOperation);
	tryToActivateOneWorker();
      }
      mutex_group.unlock();
    }
    
    void ThreadGroup::blockingOperationFinished() {
      pthread_t myID = pthread_self();
      mutex_group.lock();
      if (infoMap.find(myID) == infoMap.end())
	throw std::runtime_error("Thread group blockingOperationFinished called"
				 " for thread not in thread group.");
      if (infoMap[myID]->highPriority) {
	/* Do nothing if thread has high priority. */
	mutex_group.unlock();
	return;
      }
      if ( infoMap[myID]->getState() != doingBlockingOperation) {
	/* Do nothing in this case. We allow this because such "extra"
	   calls to blockingOperationFinished() can happen when
	   getting several objects at once. */
	mutex_group.unlock();
	return;
      }
      infoMap[myID]->blockingOperationCounter = 0;
      if (n_working() == n_max_working) {
	readyAfterBlockingQueue.push(myID);
	infoMap[myID]->changeState(readyAfterBlockingOperation);
	while(infoMap[myID]->getState() == readyAfterBlockingOperation)
	  infoMap[myID]->cond.wait(mutex_group);
      }
      else {
	// go to working state
	infoMap[myID]->changeState(working);
      }
      /* Check that thread is now in working state. */
      if ( infoMap[myID]->getState() != working)
	throw std::runtime_error("Thread not in working state when released"
				 " from blockingOperationFinished ");
      mutex_group.unlock();
    }

    void ThreadGroup::setPriorityHigh() {
      pthread_t myID = pthread_self();
      mutex_group.lock();
      if (infoMap.find(myID) == infoMap.end())
	throw std::runtime_error("Thread group setPriorityHigh called"
				 " for thread not in thread group.");
      if ( infoMap[myID]->getState() != working) 
	throw std::runtime_error("Thread group setPriorityHigh called"
				 " for thread not working state.");
      if(infoMap[myID]->highPriority)
	throw std::runtime_error("Thread group setPriorityHigh called"
				 " for thread that already has high priority.");
      infoMap[myID]->highPriority = true;
      mutex_group.unlock();
    }
    
    void ThreadGroup::setPriorityNormal() {
      pthread_t myID = pthread_self();
      mutex_group.lock();
      if (infoMap.find(myID) == infoMap.end())
	throw std::runtime_error("Thread group setPriorityNormal called"
				 " for thread not in thread group.");
      if ( infoMap[myID]->getState() != working) 
	throw std::runtime_error("Thread group setPriorityNormal called"
				 " for thread not working state.");
      if(!infoMap[myID]->highPriority)
	throw std::runtime_error("Thread group setPriorityNormal called"
				 " for thread that does not have high priority.");
      infoMap[myID]->highPriority = false;
      mutex_group.unlock();
    }
    
    void ThreadGroup::getThreadStats(cht::vector<State_stats> & threadStats) {
      pthread_t myID = pthread_self();
      mutex_group.lock();
      infoMap[myID]->getStats(threadStats);
      mutex_group.unlock();
    }

    void ThreadGroup::increaseNMaxWorking() {
      pthread_t myID = pthread_self();
      mutex_group.lock();
      if (infoMap.find(myID) == infoMap.end())
	throw std::runtime_error("Thread group increaseNMaxWorking() called"
				 " for thread not in thread group");      
      if ( infoMap[myID]->getState() != working )
	throw std::runtime_error("Thread group increaseNMaxWorking() called "
				 "for thread in nonworking state");    
      tryToActivateOneWorker();
      n_max_working++;
      mutex_group.unlock();
    }
    bool ThreadGroup::decreaseNMaxWorking() {
      pthread_t myID = pthread_self();
      mutex_group.lock();
      if (infoMap.find(myID) == infoMap.end())
	throw std::runtime_error("Thread group yield() called"
				 " for thread not in thread group");      
      if ( infoMap[myID]->getState() != working )
	throw std::runtime_error("Thread group yield() called "
				 "for thread in nonworking state");    
      if (n_max_working <= 1) {
	mutex_group.unlock();      
	return false; // Only decrease if more than 1 thread
      }
      n_max_working--;
      if ( n_working() > n_max_working ) {
	// Go to waiting state
	infoMap[myID]->changeState(waiting);
	while(infoMap[myID]->getState() == waiting) {
	  infoMap[myID]->cond.wait(mutex_group);
	}
      }
      mutex_group.unlock();      
      return true;
    }

    int ThreadGroup::getNMaxWorking() {
      mutex_group.lock();
      int tmp = n_max_working;
      mutex_group.unlock();      
      return tmp;
    }
    
  } // end namespace Threads
} // end namespace cht
