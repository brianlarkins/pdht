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
#include "Thread.h"
#include "Manager.h"
namespace cht {
  namespace Threads {

    void* Thread::thread_func( void* arg ) {
      Thread* thrPtr = (Thread*) arg;
      Manager::instance().registerThread(thrPtr);
      (thrPtr->tfun)(thrPtr->arg);
      thrPtr->setFinishedFlag();
      Manager::instance().unregisterThread();
      return NULL;
    }
    Thread::Thread( std::string identifier ) 
      : identifier( identifier ), tfun( NULL ), arg( NULL ), 
	cht_ID( Manager::instance().acquireThreadID() ),
	pthr_ID_creator( pthread_self() ), finishedFlag(false) {
      pthr_ID = pthread_self();
      Manager::instance().registerThread(this);
    }
    Thread::Thread( std::string identifier, threadFun tfun, void* arg ) 
      : identifier( identifier ), tfun( tfun ), arg( arg ), 
	cht_ID( Manager::instance().acquireThreadID() ),
	pthr_ID_creator( pthread_self() ), finishedFlag(false) {
      // Note that the Thread object is completely initialized at this point.
      pthread_create(&pthr_ID, NULL, thread_func, this);	
    }
    Thread::~Thread() {
      pthread_t myID = pthread_self();
      if ( !pthread_equal(pthr_ID_creator, myID) )
	throw std::runtime_error("Thread destructor not called by the "
				 "same thread that called the constructor.");
      if ( pthread_equal(pthr_ID_creator, pthr_ID) )
	Manager::instance().unregisterThread();	
      else
	pthread_join(pthr_ID, NULL);
    }
    size_t Thread::get_ID_number() const { return cht_ID; }
    std::string Thread::get_ID_string() const { return identifier; }
    void Thread::setFinishedFlag() {
      mutex.lock();
      finishedFlag = true;
      mutex.unlock();
    }
    bool Thread::finished() {
      bool result;
      mutex.lock();
      result = finishedFlag;
      mutex.unlock();
      return result;
    }

  } // end namespace Threads
} // end namespace cht
