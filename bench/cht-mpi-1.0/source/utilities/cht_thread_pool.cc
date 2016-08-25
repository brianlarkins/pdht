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
#include <stdexcept>
#include <iostream>
#include "cht_thread_pool.h"

namespace cht {

  void thread_pool::set_n_working_max(int const n_max) {
    if ( !infoMap.empty() ) 
      throw std::runtime_error("set_n_working_max called when when "
			       "list is not empty");
    n_working_max = n_max;
  }

  void thread_pool::join() {
    pthread_t myID = pthread_self();
    if (infoMap.find(myID) != infoMap.end())
      throw std::runtime_error("Join called for thread already "
			       "in thread pool");      
    pthread_mutex_lock(&mutex_pool);
    size_t pool_tid = ++threadCounter;
    infoMap[myID] = new ThreadInfo( pool_tid );
    if ( n_working() < n_working_max ) {
      infoMap[myID]->changeState(working);
    }
    else {
      infoMap[myID]->changeState(waiting);
      while(infoMap[myID]->getState() == waiting) {
	pthread_cond_wait(&infoMap[myID]->cond, &mutex_pool);
      }
    }
    pthread_mutex_unlock(&mutex_pool);
    // Should be in working state here
    pthread_mutex_lock(&mutex_pool); 
    if ( infoMap[myID]->getState() != working) 
      throw std::runtime_error("Thread not in working state when "
			       "leaving join function");
    pthread_mutex_unlock(&mutex_pool); 
  }

  void thread_pool::leave() {
    pthread_t myID = pthread_self();
    pthread_mutex_lock(&mutex_pool); 
    if ( infoMap[myID]->getState() != working)
      throw std::runtime_error("Thread pool leave called for thread in"
			       " nonworking state");
    infoMap[myID]->changeState(afterLeave);
    delete infoMap[myID];
    if ( infoMap.erase(myID) != 1 )
      throw std::runtime_error("Thread pool leave: error in erase");
    if (n_working() >= n_working_max)
      throw std::runtime_error("Thread pool leave called with "
			       "n_working() >= n_working_max after erase");
    tryToActivateOneWorker();
    pthread_mutex_unlock(&mutex_pool); 
  }

  void thread_pool::yield() {
    pthread_t myID = pthread_self();
    pthread_mutex_lock(&mutex_pool); 
    if ( infoMap[myID]->getState() != working)
      throw std::runtime_error("Thread pool yield() called"
			       " for thread in nonworking state");    
    if ( !readyAfterBlockingQueue.empty() ) {
      // Activate another thread
      pthread_t threadToActivate = readyAfterBlockingQueue.front();
      readyAfterBlockingQueue.pop();
      infoMap[threadToActivate]->changeState(working);
      pthread_cond_signal(&infoMap[threadToActivate]->cond);
      // Go to waiting state
      infoMap[myID]->changeState(waiting);
      while(infoMap[myID]->getState() == waiting) {
	pthread_cond_wait(&infoMap[myID]->cond, &mutex_pool);
      }
      pthread_mutex_unlock(&mutex_pool);
    }
    else
      // No thread to activate, just continue "working"
      pthread_mutex_unlock(&mutex_pool); 
    pthread_mutex_lock(&mutex_pool); 
    pthread_mutex_unlock(&mutex_pool);
  }


  void thread_pool::startingBlockingOperation() {
    pthread_t myID = pthread_self();
    pthread_mutex_lock(&mutex_pool); 
    if (infoMap.find(myID) == infoMap.end())
      throw std::runtime_error("Thread pool startingBlockingOperation called"
			       " for thread not in thread pool");      
    if ( infoMap[myID]->getState() != working)
      throw std::runtime_error("Thread pool startingBlockingOperation called"
			       " for thread in nonworking state");
    infoMap[myID]->changeState(doingBlockingOperation);
    tryToActivateOneWorker();
    pthread_mutex_unlock(&mutex_pool);
  }

  void thread_pool::blockingOperationFinished() {
    pthread_t myID = pthread_self();
    pthread_mutex_lock(&mutex_pool); 
    if (infoMap.find(myID) == infoMap.end())
      throw std::runtime_error("Thread pool blockingOperationFinished called"
			       " for thread not in thread pool");      
    if ( infoMap[myID]->getState() != doingBlockingOperation)
      throw std::runtime_error("Thread pool blockingOperationFinished called"
			       " for thread not in blocking operation state");
    if (n_working() == n_working_max) {
      // lock mutex - go to readyAfterBlockingOperation state
      readyAfterBlockingQueue.push(myID);
      infoMap[myID]->changeState(readyAfterBlockingOperation);
      while(infoMap[myID]->getState() == readyAfterBlockingOperation)
	pthread_cond_wait(&infoMap[myID]->cond, &mutex_pool);
      pthread_mutex_unlock(&mutex_pool);
    }
    else {
      // go to working state
      infoMap[myID]->changeState(working);
      pthread_mutex_unlock(&mutex_pool);     
    }    
    pthread_mutex_lock(&mutex_pool); 
    if ( infoMap[myID]->getState() != working)
      throw std::runtime_error("Thread not in working state when released from"
			       " blockingOperationFinished ");
    pthread_mutex_unlock(&mutex_pool);
  }
  
  int thread_pool::n_working() {
    // thread pool mutex already locked here
    ThreadInfoMap::iterator it;
    int count = 0;
    for ( it=infoMap.begin() ; it != infoMap.end(); it++ ) {
      if ( (*it).second->getState() == working)
	++count;
    }
    return count;
  }

  void thread_pool::tryToActivateOneWorker() {
    // thread pool mutex already locked here
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
	if ( (*it).second->getState() == waiting) {
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
    pthread_cond_signal(&infoMap[threadToActivate]->cond);
  }

  thread_pool::thread_pool() 
    : threadCounter(0), n_working_max(1)
  {
    pthread_mutex_init(&mutex_pool, NULL);
  }

  
};
