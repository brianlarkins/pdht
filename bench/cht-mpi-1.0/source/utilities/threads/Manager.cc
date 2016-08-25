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
#include "Manager.h"
namespace cht {
  namespace Threads {

    void Manager::registerThread(Thread* thr) {
      pthread_t myID = pthread_self();
      mutex.lock();
      threadMap_all_threads[myID] = thr;
      mutex.unlock();
    }
    void Manager::unregisterThread() {
      pthread_t myID = pthread_self();
      mutex.lock();
      if ( threadMap_all_threads.erase(myID) != 1 )
	throw std::runtime_error("Threads::Manager unregisterThread: "
				 "error in erase");
      mutex.unlock();
    }
    void Manager::registerThreadGroup(ThreadGroup* thrGrp) {
      pthread_t myID = pthread_self();
      mutex.lock();
      threadGroupMap_all_threads[myID] = thrGrp;
      mutex.unlock();
    }
    void Manager::unregisterThreadGroup() {
      pthread_t myID = pthread_self();
      mutex.lock();
      if ( threadGroupMap_all_threads.erase(myID) != 1 )
	throw std::runtime_error("Threads::Manager unregisterThreadGroup: "
				 "error in erase");
      mutex.unlock();
    }

    Thread& Manager::getThread() {
      pthread_t myID = pthread_self();
      mutex.lock();
      if ( threadMap_all_threads.find(myID) == threadMap_all_threads.end() )
	throw std::runtime_error("Thread not registered."); 
      Thread* tmp = threadMap_all_threads[myID];
      mutex.unlock();
      return *tmp;
    }
    ThreadGroup& Manager::getThreadGroup() {
      pthread_t myID = pthread_self();
      mutex.lock();
      if ( threadGroupMap_all_threads.find(myID) == 
	   threadGroupMap_all_threads.end() )
	throw std::runtime_error("Thread not in any thread group."); 
      ThreadGroup* tmp = threadGroupMap_all_threads[myID];
      mutex.unlock();
      return *tmp;
    }
    
    
  } // end namespace Threads
} // end namespace cht
