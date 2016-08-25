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
#include "MPIWrapperFunneled.h"
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <pthread.h>

// Proxy base classes

struct Proxy_MPI_Base {
  pthread_mutex_t proxy_mutex;
  void LockMutex() {
    pthread_mutex_lock(&proxy_mutex);
  }
  void UnlockMutex() {
    pthread_mutex_unlock(&proxy_mutex);
  }
  Proxy_MPI_Base() {
    pthread_mutex_init(&proxy_mutex, NULL);
    LockMutex();
  }
};

struct Proxy_MPI : public Proxy_MPI_Base {
  virtual void callBack() = 0;
  virtual int id_for_debug() = 0;
};

struct Proxy_MPI_Wait_calls : public Proxy_MPI_Base {
  virtual int callBackReturn() = 0;
};

#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
struct Proxy_MPI_Non_blocking_send_calls : public Proxy_MPI_Base {
  virtual void createGRequest() = 0;
  virtual void callNonBlockingSend() = 0;
  virtual bool tryToFinish() = 0;
};
int non_blocking_send_query_fn(void *extra_state, MPI_Status *status);
int non_blocking_send_free_fn(void *extra_state);
int non_blocking_send_cancel_fn(void *extra_state, int complete);
#endif

void* global_thread_func(void* arg)
{
  MPIWrapperFunneled* p = (MPIWrapperFunneled*) arg;
  p->executor_thread_func();
  return NULL;
}

MPIWrapperFunneled::MPIWrapperFunneled() 
  : init_proxy(0), finalize_proxy(0)
{
  for (int i = 0; i < number_of_mutexes; ++i)
    pthread_mutex_init(&mutexes[i], NULL);
  pthread_cond_init(&executor_wakeup_cond, NULL);
  // Start up thread
  pthread_create(&thread, NULL, global_thread_func, this);
}

MPIWrapperFunneled::~MPIWrapperFunneled() {
#if 0
  // Signal exit to executor thread
  LockMutex();
  wake_up_executor();
  UnlockMutex();
#endif
  // Join thread here!!
  pthread_join(thread, NULL);
}

void MPIWrapperFunneled::LockMutex(int const mutex_id)
{
  pthread_mutex_lock(&mutexes[mutex_id]);
}

void MPIWrapperFunneled::UnlockMutex(int const mutex_id)
{
  pthread_mutex_unlock(&mutexes[mutex_id]);
}


/* 
   This function is called by the thread that is supposed to make all
   MPI calls. No other thread is supposed to make any MPI calls. 
*/
void MPIWrapperFunneled::executor_thread_func() {
  std::list<Proxy_MPI_Wait_calls*> pendingForCompletionWaitCalls;
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  std::list<Proxy_MPI_Non_blocking_send_calls*> waitingToCallSendNonBlockingSendCalls[2];
  std::list<Proxy_MPI_Non_blocking_send_calls*> pendingForCompletionNonBlockingSendCalls[2];
#endif

#if 0
  std::cout<<"Beginning of MPIWrapperFunneled::executor_thread_func()" << std::endl;
#endif

  // Waiting for initialization loop
  bool mpi_initialized_in_wrapper = false;
  while(1) {
    // Check MPI INIT -> break
    int flag;    
    ::MPI_Initialized(&flag);
    if (flag) {
      break;
    }
    // Check init_proxy -> init, break
    LockMutex(mutex_id_init);
    if (init_proxy) {
      init_proxy->callBack();
      init_proxy->UnlockMutex(); 
      init_proxy = 0;
      UnlockMutex(mutex_id_init);  
      mpi_initialized_in_wrapper = true;
      break;
    }
    UnlockMutex(mutex_id_init);  
  }
  
  while(1) {

    LockMutex(mutex_id_init);
    if (init_proxy) {
      throw std::runtime_error("Attempt to initialize MPI a second time");
    }
    UnlockMutex(mutex_id_init);  

#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
    int max_no_simultaneous_sends[2] = {1, 1};
    for (int ind = 0; ind < 2; ++ind) {
      {
	// Loop over send calls pending for start of generalized request
	LockMutex(mutex_id_non_blocking_send_list);
	Proxy_MPI_Non_blocking_send_calls* proxy;
	while ( !pendingNonBlockingSendCalls[ind].empty() ) {
	  proxy = pendingNonBlockingSendCalls[ind].front();
	  // Create Grequest
	  proxy->createGRequest();
	  // Call is ready, unlock mutex
	  proxy->UnlockMutex(); 
	  pendingNonBlockingSendCalls[ind].pop_front();
	  waitingToCallSendNonBlockingSendCalls[ind].push_back(proxy);
	}
	UnlockMutex(mutex_id_non_blocking_send_list);  
      }
      if (pendingForCompletionNonBlockingSendCalls[ind].size() < max_no_simultaneous_sends[ind] && 
	  !waitingToCallSendNonBlockingSendCalls[ind].empty()) {
	Proxy_MPI_Non_blocking_send_calls* proxy;
	proxy = waitingToCallSendNonBlockingSendCalls[ind].front();
	proxy->callNonBlockingSend();
	waitingToCallSendNonBlockingSendCalls[ind].pop_front();
	pendingForCompletionNonBlockingSendCalls[ind].push_back(proxy);      
      }

      {
	Proxy_MPI_Non_blocking_send_calls* proxy;
	std::list<Proxy_MPI_Non_blocking_send_calls*>::iterator it
	  = pendingForCompletionNonBlockingSendCalls[ind].begin();
	while ( it != pendingForCompletionNonBlockingSendCalls[ind].end() ) {
	  proxy = (*it);
	  if ( proxy->tryToFinish() ) {
	    it = pendingForCompletionNonBlockingSendCalls[ind].erase(it);
	  }
	  else {
	    it++;
	  }
	}
      }
    }
    // non-blocking send calls stuff done 
#endif

    {
      // Loop over non blocking calls pending for execution
      LockMutex(mutex_id_non_blocking_list);
      Proxy_MPI* proxy;
      while ( !pendingNonBlocking.empty() ) {
	proxy = pendingNonBlocking.front();
	// NOTE: It can be bad to unlock the list mutex here since this makes
	// it possible to add more proxies to the list. Then it can take
	// a long time to get out of the while loop to process the other
	// lists.
	proxy->callBack();
	// Call is ready, unlock mutex
	proxy->UnlockMutex(); 
	pendingNonBlocking.pop_front();
      }
      UnlockMutex(mutex_id_non_blocking_list);  
    }
    
    // Move pending wait calls to local list (without mutex)
    LockMutex(mutex_id_waitcalls_list);
    while ( !pendingWaitCalls.empty() ) {
      pendingForCompletionWaitCalls.push_back(pendingWaitCalls.front());
      pendingWaitCalls.pop_front();
    }
    UnlockMutex(mutex_id_waitcalls_list);  
    
    {
      // Local waitcalls list processed here, (no list mutex locked)
      Proxy_MPI_Wait_calls* proxy_wait_calls;
      std::list<Proxy_MPI_Wait_calls*>::iterator it
	= pendingForCompletionWaitCalls.begin();
      while ( it != pendingForCompletionWaitCalls.end() ) {
	proxy_wait_calls = (*it);
	int flag = proxy_wait_calls->callBackReturn();
	if (flag) {
	  // Call is ready, unlock mutex
	  proxy_wait_calls->UnlockMutex(); 
	  it = pendingForCompletionWaitCalls.erase(it);
	}
	else {
	  it++;
	}
      }
    }
    
    bool allListsEmpty = pendingForCompletionWaitCalls.empty();
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
    for (int ind = 0; ind < 2; ++ind) {
      if (allListsEmpty) 
	allListsEmpty = waitingToCallSendNonBlockingSendCalls[ind].empty() && pendingForCompletionNonBlockingSendCalls[ind].empty();
    }
    if (allListsEmpty) {
      LockMutex(mutex_id_non_blocking_send_list);
      for (int ind = 0; ind < 2; ++ind) 
	allListsEmpty = allListsEmpty && pendingNonBlockingSendCalls[ind].empty();      
      UnlockMutex(mutex_id_non_blocking_send_list);
    }
#endif
    if (allListsEmpty) {
      LockMutex(mutex_id_non_blocking_list);
      allListsEmpty = pendingNonBlocking.empty();
      UnlockMutex(mutex_id_non_blocking_list);
    }
    if (allListsEmpty) {
      LockMutex(mutex_id_waitcalls_list);
      allListsEmpty = pendingWaitCalls.empty();
      UnlockMutex(mutex_id_waitcalls_list);
    }
    
    if ( allListsEmpty ) {
      // Check for finalize
      int flag;
      ::MPI_Finalized(&flag);
      if (flag) {
	if (mpi_initialized_in_wrapper) {
	  std::cout<<"WARNING (MPIWrapperFunneled): Since MPI was initialized via MPIWrapperFunneled MPI_Finalize should be called via wrapper as well."<<std::endl;
	}
	break;
      }
      // Check finalize proxy
      LockMutex(mutex_id_finalize);
      if (finalize_proxy != 0) {
	if (!mpi_initialized_in_wrapper) {
	  std::cout<<"WARNING (MPIWrapperFunneled): Since MPI was not initialized via MPIWrapperFunneled MPI_Finalize should be called outside wrapper as well."<<std::endl;
	}
	finalize_proxy->callBack();
	finalize_proxy->UnlockMutex(); 
	finalize_proxy = 0;
	UnlockMutex(mutex_id_finalize);  
	break;
      }
      UnlockMutex(mutex_id_finalize);
#if 1
      // Check carefully that lists are still empty before setting
      // executor_wakeup_flag to false.
      for(int i = 0; i < number_of_mutexes; i++)
	LockMutex(i);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
      for (int ind = 0; ind < 2; ++ind) 
	allListsEmpty = allListsEmpty && pendingNonBlockingSendCalls[ind].empty();
#endif
      allListsEmpty = allListsEmpty && pendingNonBlocking.empty();
      allListsEmpty = allListsEmpty && pendingWaitCalls.empty();
      if(allListsEmpty && (!init_proxy) )
	executor_wakeup_flag = false;
      for(int i = 0; i < number_of_mutexes; i++)
	UnlockMutex(i);
      // Wait until someone calls wake_up_executor().
      LockMutex(mutex_id_executor_wakeup_flag);
      while(executor_wakeup_flag == false) {
	pthread_cond_wait(&executor_wakeup_cond, &mutexes[mutex_id_executor_wakeup_flag]);
      }
      UnlockMutex(mutex_id_executor_wakeup_flag);
#endif
    } // end if all lists empty
  } // end while
  
  
  // Elias note: what we wanted to do here was to create a MPI
  // "generalized request" and wait for it. We thought is would then
  // be possible for another thread to "wake up" this thread without
  // making any MPI call, by using some kind of pthread call. But now
  // after reading in the MPI 2.2 specification about how "generalized
  // requests" work, I think this will not be possible. It seems that
  // the only way to "wake up" a generalized request is by using the
  // MPI call MPI_Grequest_complete() but that means the other thread
  // would have to make an MPI call as well. So now I think it is
  // impossible to make only one thread handle MPI calls. What we need
  // is an MPI implementation that really supports
  // MPI_THREAD_MULTIPLE.
  //  throw std::runtime_error("Error in MPIWrapperFunneled::executor_thread_func(): "
  //			   "damn, this doesn't work!");
#if 0
  if (!mpi_was_already_initialized) {
    ::MPI_Finalize();
  }
#endif
}


MPIWrapperFunneled& MPIWrapperFunneled::instance()
{
  // FIXME: This is not thread safe!!
  // Move theInstance to create() function and use the
  // double checked locking pattern
  // Where should thread be created? 
  // Probably inside create() or in constructor? 
  static MPIWrapperFunneled theInstance; // local static variable
  return theInstance;
}  

// New code

// MPI calls locking mutex with id mutex_id_init

struct Proxy_MPI_Init_thread : public Proxy_MPI {
  int  return_value;
  int * argc;
  char *** argv;
  int  required;
  int * provided;
  Proxy_MPI_Init_thread(int * argc, char *** argv, int  required, int * provided)
  : argc(argc), argv(argv), required(required), provided(provided)
  {}
  void callBack() {
    return_value = ::MPI_Init_thread(argc, argv, required, provided);
  }
  int id_for_debug() { return 1; }
};
int  MPIWrapperFunneled::_MPI_Init_thread(int * argc, char *** argv, int  required, int * provided) {
  Proxy_MPI_Init_thread proxyObj(argc, argv, required, provided);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_init);
  init_proxy = proxy;
  UnlockMutex(mutex_id_init);
  wake_up_executor();
  proxy->LockMutex();
  // Checking provided thread level
  if (*provided < MPI_THREAD_FUNNELED) {
    std::cout << "WARNING (MPIWrapperFunneled) : MPIWrapperFunneled requires MPI thread level MPI_THREAD_FUNNELED which"
      " is not provided" << std::endl;
  }
  return proxyObj.return_value;
}

// MPI calls locking mutex with id mutex_id_finalize

struct Proxy_MPI_Finalize : public Proxy_MPI {
  int  return_value;
  Proxy_MPI_Finalize()
  {}
  void callBack() {
    return_value = ::MPI_Finalize();
  }
  int id_for_debug() { return 2; }
};
int  MPIWrapperFunneled::_MPI_Finalize() {
  Proxy_MPI_Finalize proxyObj;
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_finalize);
  finalize_proxy = proxy;
  UnlockMutex(mutex_id_finalize);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



// Blocking calls that are handled specifically by calling corresponding internal nonblocking routines

int MPIWrapperFunneled::_MPI_Recv(void * buf, int  count, MPI_Datatype  datatype, int  source, int  tag, MPI_Comm  comm, MPI_Status * status) {
  MPI_Request request;
  this->_MPI_Irecv(buf, count, datatype, source, tag, comm, &request); 
  return this->_MPI_Wait(&request, status);
}

int MPIWrapperFunneled::_MPI_Bsend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm) {
  MPI_Request request;
  this->_MPI_Ibsend(buf, count, datatype, dest, tag, comm, &request); 
  return this->_MPI_Wait(&request, MPI_STATUS_IGNORE);
}

int MPIWrapperFunneled::_MPI_Rsend(void * ibuf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm) {
  MPI_Request request;
  this->_MPI_Irsend(ibuf, count, datatype, dest, tag, comm, &request); 
  return this->_MPI_Wait(&request, MPI_STATUS_IGNORE);
}

int MPIWrapperFunneled::_MPI_Send(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm) {
  MPI_Request request;
  this->_MPI_Isend(buf, count, datatype, dest, tag, comm, &request); 
  return this->_MPI_Wait(&request, MPI_STATUS_IGNORE);
}

int MPIWrapperFunneled::_MPI_Ssend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm) {
  MPI_Request request;
  this->_MPI_Issend(buf, count, datatype, dest, tag, comm, &request); 
  return this->_MPI_Wait(&request, MPI_STATUS_IGNORE);  
}

// MPI calls locking mutex with id mutex_id_waitcalls_list
// Particularly handled wait functions


struct Proxy_MPI_Waitall : public Proxy_MPI_Wait_calls {
  int  return_value;
  int  count;
  MPI_Request * array_of_requests;
  MPI_Status * array_of_statuses;
  Proxy_MPI_Waitall(int  count, MPI_Request * array_of_requests, MPI_Status * array_of_statuses)
  : count(count), array_of_requests(array_of_requests), array_of_statuses(array_of_statuses)
  {}
  int callBackReturn() {
    int flag;
    return_value = ::MPI_Testall(count, array_of_requests, 
				 &flag, array_of_statuses);
    return flag;
  }
};
int  MPIWrapperFunneled::_MPI_Waitall(int  count, MPI_Request * array_of_requests, MPI_Status * array_of_statuses) {
  Proxy_MPI_Waitall proxyObj(count, array_of_requests, array_of_statuses);
  Proxy_MPI_Wait_calls* proxy = &proxyObj;
  LockMutex(mutex_id_waitcalls_list);
  pendingWaitCalls.push_back(proxy);
  UnlockMutex(mutex_id_waitcalls_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Waitany : public Proxy_MPI_Wait_calls {
  int  return_value;
  int  count;
  MPI_Request * array_of_requests;
  int * index;
  MPI_Status * status;
  Proxy_MPI_Waitany(int  count, MPI_Request * array_of_requests, int * index, MPI_Status * status)
  : count(count), array_of_requests(array_of_requests), index(index), status(status)
  {}
  int callBackReturn() {
    int flag;
    return_value = ::MPI_Testany(count, array_of_requests, index, &flag, status);
    return flag;
  }
};
int  MPIWrapperFunneled::_MPI_Waitany(int  count, MPI_Request * array_of_requests, int * index, MPI_Status * status) {
  Proxy_MPI_Waitany proxyObj(count, array_of_requests, index, status);
  Proxy_MPI_Wait_calls* proxy = &proxyObj;
  LockMutex(mutex_id_waitcalls_list);
  pendingWaitCalls.push_back(proxy);
  UnlockMutex(mutex_id_waitcalls_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}


struct Proxy_MPI_Wait : public Proxy_MPI_Wait_calls {
  int  return_value;
  MPI_Request * request;
  MPI_Status * status;
  Proxy_MPI_Wait(MPI_Request * request, MPI_Status * status)
  : request(request), status(status)
  {}
  int callBackReturn() {
    int flag;
    return_value = ::MPI_Test(request, &flag, status);
    return flag;
  }
};
int  MPIWrapperFunneled::_MPI_Wait(MPI_Request * request, MPI_Status * status) {
  Proxy_MPI_Wait proxyObj(request, status);
  Proxy_MPI_Wait_calls* proxy = &proxyObj;
  LockMutex(mutex_id_waitcalls_list);
  pendingWaitCalls.push_back(proxy);
  UnlockMutex(mutex_id_waitcalls_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

struct Proxy_MPI_Waitsome : public Proxy_MPI_Wait_calls {
  int  return_value;
  int  incount;
  MPI_Request * array_of_requests;
  int * outcount;
  int * array_of_indices;
  MPI_Status * array_of_statuses;
  Proxy_MPI_Waitsome(int  incount, MPI_Request * array_of_requests, int * outcount, int * array_of_indices, MPI_Status * array_of_statuses)
    : incount(incount), array_of_requests(array_of_requests), outcount(outcount), array_of_indices(array_of_indices), array_of_statuses(array_of_statuses)
  {}
  int callBackReturn() {
    return_value = ::MPI_Testsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
    return *outcount;
  }
};
int  MPIWrapperFunneled::_MPI_Waitsome(int  incount, MPI_Request * array_of_requests, int * outcount, int * array_of_indices, MPI_Status * array_of_statuses) {
  Proxy_MPI_Waitsome proxyObj(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  Proxy_MPI_Wait_calls* proxy = &proxyObj;
  LockMutex(mutex_id_waitcalls_list);
  pendingWaitCalls.push_back(proxy);
  UnlockMutex(mutex_id_waitcalls_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

struct Proxy_MPI_Probe : public Proxy_MPI_Wait_calls {
  int  return_value;
  int  source;
  int  tag;
  MPI_Comm  comm;
  MPI_Status * status;
  Proxy_MPI_Probe(int  source, int  tag, MPI_Comm  comm, MPI_Status * status)
  : source(source), tag(tag), comm(comm), status(status)
  {}
  int callBackReturn() {
    int flag;
    return_value = ::MPI_Iprobe(source, tag, comm, &flag, status);
    return flag;
  }
};
int MPIWrapperFunneled::_MPI_Probe(int  source, int  tag, MPI_Comm  comm, MPI_Status * status) {
  Proxy_MPI_Probe proxyObj(source, tag, comm, status);
  Proxy_MPI_Wait_calls* proxy = &proxyObj;
  LockMutex(mutex_id_waitcalls_list);
  pendingWaitCalls.push_back(proxy);
  UnlockMutex(mutex_id_waitcalls_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}





// MPI calls locking mutex with id mutex_id_non_blocking_list
// Blocking and nonblocking calls that are only funneled to the executor thread,
// not handled in a special manner even if they block the executor.


struct Proxy_MPI_Abort : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  errorcode;
  Proxy_MPI_Abort(MPI_Comm  comm, int  errorcode)
  : comm(comm), errorcode(errorcode)
  {}
  void callBack() {
    return_value = ::MPI_Abort(comm, errorcode);
  }
  int id_for_debug() { return 3; }
};
int  MPIWrapperFunneled::_MPI_Abort(MPI_Comm  comm, int  errorcode) {
  Proxy_MPI_Abort proxyObj(comm, errorcode);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}




#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Accumulate : public Proxy_MPI {
  int  return_value;
  void * origin_addr;
  int  origin_count;
  MPI_Datatype  origin_datatype;
  int  target_rank;
  MPI_Aint  target_disp;
  int  target_count;
  MPI_Datatype  target_datatype;
  MPI_Op  op;
  MPI_Win  win;
  Proxy_MPI_Accumulate(void * origin_addr, int  origin_count, MPI_Datatype  origin_datatype, int  target_rank, MPI_Aint  target_disp, int  target_count, MPI_Datatype  target_datatype, MPI_Op  op, MPI_Win  win)
  : origin_addr(origin_addr), origin_count(origin_count), origin_datatype(origin_datatype), target_rank(target_rank), target_disp(target_disp), target_count(target_count), target_datatype(target_datatype), op(op), win(win)
  {}
  void callBack() {
    return_value = ::MPI_Accumulate(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, op, win);
  }
  int id_for_debug() { return 4; }
};
int  MPIWrapperFunneled::_MPI_Accumulate(void * origin_addr, int  origin_count, MPI_Datatype  origin_datatype, int  target_rank, MPI_Aint  target_disp, int  target_count, MPI_Datatype  target_datatype, MPI_Op  op, MPI_Win  win) {
  Proxy_MPI_Accumulate proxyObj(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, op, win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Add_error_class : public Proxy_MPI {
  int  return_value;
  int * errorclass;
  Proxy_MPI_Add_error_class(int * errorclass)
  : errorclass(errorclass)
  {}
  void callBack() {
    return_value = ::MPI_Add_error_class(errorclass);
  }
  int id_for_debug() { return 5; }
};
int  MPIWrapperFunneled::_MPI_Add_error_class(int * errorclass) {
  Proxy_MPI_Add_error_class proxyObj(errorclass);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Add_error_code : public Proxy_MPI {
  int  return_value;
  int  errorclass;
  int * errorcode;
  Proxy_MPI_Add_error_code(int  errorclass, int * errorcode)
  : errorclass(errorclass), errorcode(errorcode)
  {}
  void callBack() {
    return_value = ::MPI_Add_error_code(errorclass, errorcode);
  }
  int id_for_debug() { return 6; }
};
int  MPIWrapperFunneled::_MPI_Add_error_code(int  errorclass, int * errorcode) {
  Proxy_MPI_Add_error_code proxyObj(errorclass, errorcode);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Add_error_string : public Proxy_MPI {
  int  return_value;
  int  errorcode;
  char * string;
  Proxy_MPI_Add_error_string(int  errorcode, char * string)
  : errorcode(errorcode), string(string)
  {}
  void callBack() {
    return_value = ::MPI_Add_error_string(errorcode, string);
  }
  int id_for_debug() { return 7; }
};
int  MPIWrapperFunneled::_MPI_Add_error_string(int  errorcode, char * string) {
  Proxy_MPI_Add_error_string proxyObj(errorcode, string);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Allgather : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  int  sendcount;
  MPI_Datatype  sendtype;
  void * recvbuf;
  int  recvcount;
  MPI_Datatype  recvtype;
  MPI_Comm  comm;
  Proxy_MPI_Allgather(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, MPI_Comm  comm)
  : sendbuf(sendbuf), sendcount(sendcount), sendtype(sendtype), recvbuf(recvbuf), recvcount(recvcount), recvtype(recvtype), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }
  int id_for_debug() { return 9; }
};
int  MPIWrapperFunneled::_MPI_Allgather(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, MPI_Comm  comm) {
  Proxy_MPI_Allgather proxyObj(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Allgatherv : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  int  sendcount;
  MPI_Datatype  sendtype;
  void * recvbuf;
  int * recvcounts;
  int * displs;
  MPI_Datatype  recvtype;
  MPI_Comm  comm;
  Proxy_MPI_Allgatherv(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int * recvcounts, int * displs, MPI_Datatype  recvtype, MPI_Comm  comm)
  : sendbuf(sendbuf), sendcount(sendcount), sendtype(sendtype), recvbuf(recvbuf), recvcounts(recvcounts), displs(displs), recvtype(recvtype), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
  }
  int id_for_debug() { return 10; }
};
int  MPIWrapperFunneled::_MPI_Allgatherv(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int * recvcounts, int * displs, MPI_Datatype  recvtype, MPI_Comm  comm) {
  Proxy_MPI_Allgatherv proxyObj(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Alloc_mem : public Proxy_MPI {
  int  return_value;
  MPI_Aint  size;
  MPI_Info  info;
  void * baseptr;
  Proxy_MPI_Alloc_mem(MPI_Aint  size, MPI_Info  info, void * baseptr)
  : size(size), info(info), baseptr(baseptr)
  {}
  void callBack() {
    return_value = ::MPI_Alloc_mem(size, info, baseptr);
  }
  int id_for_debug() { return 11; }
};
int  MPIWrapperFunneled::_MPI_Alloc_mem(MPI_Aint  size, MPI_Info  info, void * baseptr) {
  Proxy_MPI_Alloc_mem proxyObj(size, info, baseptr);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Allreduce : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  void * recvbuf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Op  op;
  MPI_Comm  comm;
  Proxy_MPI_Allreduce(void * sendbuf, void * recvbuf, int  count, MPI_Datatype  datatype, MPI_Op  op, MPI_Comm  comm)
  : sendbuf(sendbuf), recvbuf(recvbuf), count(count), datatype(datatype), op(op), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  }
  int id_for_debug() { return 12; }
};
int  MPIWrapperFunneled::_MPI_Allreduce(void * sendbuf, void * recvbuf, int  count, MPI_Datatype  datatype, MPI_Op  op, MPI_Comm  comm) {
  Proxy_MPI_Allreduce proxyObj(sendbuf, recvbuf, count, datatype, op, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Alltoall : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  int  sendcount;
  MPI_Datatype  sendtype;
  void * recvbuf;
  int  recvcount;
  MPI_Datatype  recvtype;
  MPI_Comm  comm;
  Proxy_MPI_Alltoall(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, MPI_Comm  comm)
  : sendbuf(sendbuf), sendcount(sendcount), sendtype(sendtype), recvbuf(recvbuf), recvcount(recvcount), recvtype(recvtype), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }
  int id_for_debug() { return 13; }
};
int  MPIWrapperFunneled::_MPI_Alltoall(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, MPI_Comm  comm) {
  Proxy_MPI_Alltoall proxyObj(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Alltoallv : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  int * sendcounts;
  int * sdispls;
  MPI_Datatype  sendtype;
  void * recvbuf;
  int * recvcounts;
  int * rdispls;
  MPI_Datatype  recvtype;
  MPI_Comm  comm;
  Proxy_MPI_Alltoallv(void * sendbuf, int * sendcounts, int * sdispls, MPI_Datatype  sendtype, void * recvbuf, int * recvcounts, int * rdispls, MPI_Datatype  recvtype, MPI_Comm  comm)
  : sendbuf(sendbuf), sendcounts(sendcounts), sdispls(sdispls), sendtype(sendtype), recvbuf(recvbuf), recvcounts(recvcounts), rdispls(rdispls), recvtype(recvtype), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
  }
  int id_for_debug() { return 14; }
};
int  MPIWrapperFunneled::_MPI_Alltoallv(void * sendbuf, int * sendcounts, int * sdispls, MPI_Datatype  sendtype, void * recvbuf, int * recvcounts, int * rdispls, MPI_Datatype  recvtype, MPI_Comm  comm) {
  Proxy_MPI_Alltoallv proxyObj(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Alltoallw : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  int * sendcounts;
  int * sdispls;
  MPI_Datatype * sendtypes;
  void * recvbuf;
  int * recvcounts;
  int * rdispls;
  MPI_Datatype * recvtypes;
  MPI_Comm  comm;
  Proxy_MPI_Alltoallw(void * sendbuf, int * sendcounts, int * sdispls, MPI_Datatype * sendtypes, void * recvbuf, int * recvcounts, int * rdispls, MPI_Datatype * recvtypes, MPI_Comm  comm)
  : sendbuf(sendbuf), sendcounts(sendcounts), sdispls(sdispls), sendtypes(sendtypes), recvbuf(recvbuf), recvcounts(recvcounts), rdispls(rdispls), recvtypes(recvtypes), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm);
  }
  int id_for_debug() { return 15; }
};
int  MPIWrapperFunneled::_MPI_Alltoallw(void * sendbuf, int * sendcounts, int * sdispls, MPI_Datatype * sendtypes, void * recvbuf, int * recvcounts, int * rdispls, MPI_Datatype * recvtypes, MPI_Comm  comm) {
  Proxy_MPI_Alltoallw proxyObj(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Barrier : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  Proxy_MPI_Barrier(MPI_Comm  comm)
  : comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Barrier(comm);
  }
  int id_for_debug() { return 19; }
};
int  MPIWrapperFunneled::_MPI_Barrier(MPI_Comm  comm) {
  Proxy_MPI_Barrier proxyObj(comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Bcast : public Proxy_MPI {
  int  return_value;
  void * buffer;
  int  count;
  MPI_Datatype  datatype;
  int  root;
  MPI_Comm  comm;
  Proxy_MPI_Bcast(void * buffer, int  count, MPI_Datatype  datatype, int  root, MPI_Comm  comm)
  : buffer(buffer), count(count), datatype(datatype), root(root), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Bcast(buffer, count, datatype, root, comm);
  }
  int id_for_debug() { return 20; }
};
int  MPIWrapperFunneled::_MPI_Bcast(void * buffer, int  count, MPI_Datatype  datatype, int  root, MPI_Comm  comm) {
  Proxy_MPI_Bcast proxyObj(buffer, count, datatype, root, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Bsend_init : public Proxy_MPI {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  dest;
  int  tag;
  MPI_Comm  comm;
  MPI_Request * request;
  Proxy_MPI_Bsend_init(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request)
  : buf(buf), count(count), datatype(datatype), dest(dest), tag(tag), comm(comm), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Bsend_init(buf, count, datatype, dest, tag, comm, request);
  }
  int id_for_debug() { return 21; }
};
int  MPIWrapperFunneled::_MPI_Bsend_init(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Bsend_init proxyObj(buf, count, datatype, dest, tag, comm, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Buffer_attach : public Proxy_MPI {
  int  return_value;
  void * buffer;
  int  size;
  Proxy_MPI_Buffer_attach(void * buffer, int  size)
  : buffer(buffer), size(size)
  {}
  void callBack() {
    return_value = ::MPI_Buffer_attach(buffer, size);
  }
  int id_for_debug() { return 22; }
};
int  MPIWrapperFunneled::_MPI_Buffer_attach(void * buffer, int  size) {
  Proxy_MPI_Buffer_attach proxyObj(buffer, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Buffer_detach : public Proxy_MPI {
  int  return_value;
  void * buffer;
  int * size;
  Proxy_MPI_Buffer_detach(void * buffer, int * size)
  : buffer(buffer), size(size)
  {}
  void callBack() {
    return_value = ::MPI_Buffer_detach(buffer, size);
  }
  int id_for_debug() { return 23; }
};
int  MPIWrapperFunneled::_MPI_Buffer_detach(void * buffer, int * size) {
  Proxy_MPI_Buffer_detach proxyObj(buffer, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Cancel : public Proxy_MPI {
  int  return_value;
  MPI_Request * request;
  Proxy_MPI_Cancel(MPI_Request * request)
  : request(request)
  {}
  void callBack() {
    return_value = ::MPI_Cancel(request);
  }
  int id_for_debug() { return 24; }
};
int  MPIWrapperFunneled::_MPI_Cancel(MPI_Request * request) {
  Proxy_MPI_Cancel proxyObj(request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Cart_coords : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  rank;
  int  maxdims;
  int * coords;
  Proxy_MPI_Cart_coords(MPI_Comm  comm, int  rank, int  maxdims, int * coords)
  : comm(comm), rank(rank), maxdims(maxdims), coords(coords)
  {}
  void callBack() {
    return_value = ::MPI_Cart_coords(comm, rank, maxdims, coords);
  }
  int id_for_debug() { return 25; }
};
int  MPIWrapperFunneled::_MPI_Cart_coords(MPI_Comm  comm, int  rank, int  maxdims, int * coords) {
  Proxy_MPI_Cart_coords proxyObj(comm, rank, maxdims, coords);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Cart_create : public Proxy_MPI {
  int  return_value;
  MPI_Comm  old_comm;
  int  ndims;
  int * dims;
  int * periods;
  int  reorder;
  MPI_Comm * comm_cart;
  Proxy_MPI_Cart_create(MPI_Comm  old_comm, int  ndims, int * dims, int * periods, int  reorder, MPI_Comm * comm_cart)
  : old_comm(old_comm), ndims(ndims), dims(dims), periods(periods), reorder(reorder), comm_cart(comm_cart)
  {}
  void callBack() {
    return_value = ::MPI_Cart_create(old_comm, ndims, dims, periods, reorder, comm_cart);
  }
  int id_for_debug() { return 26; }
};
int  MPIWrapperFunneled::_MPI_Cart_create(MPI_Comm  old_comm, int  ndims, int * dims, int * periods, int  reorder, MPI_Comm * comm_cart) {
  Proxy_MPI_Cart_create proxyObj(old_comm, ndims, dims, periods, reorder, comm_cart);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Cart_get : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  maxdims;
  int * dims;
  int * periods;
  int * coords;
  Proxy_MPI_Cart_get(MPI_Comm  comm, int  maxdims, int * dims, int * periods, int * coords)
  : comm(comm), maxdims(maxdims), dims(dims), periods(periods), coords(coords)
  {}
  void callBack() {
    return_value = ::MPI_Cart_get(comm, maxdims, dims, periods, coords);
  }
  int id_for_debug() { return 27; }
};
int  MPIWrapperFunneled::_MPI_Cart_get(MPI_Comm  comm, int  maxdims, int * dims, int * periods, int * coords) {
  Proxy_MPI_Cart_get proxyObj(comm, maxdims, dims, periods, coords);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Cart_map : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  ndims;
  int * dims;
  int * periods;
  int * newrank;
  Proxy_MPI_Cart_map(MPI_Comm  comm, int  ndims, int * dims, int * periods, int * newrank)
  : comm(comm), ndims(ndims), dims(dims), periods(periods), newrank(newrank)
  {}
  void callBack() {
    return_value = ::MPI_Cart_map(comm, ndims, dims, periods, newrank);
  }
  int id_for_debug() { return 28; }
};
int  MPIWrapperFunneled::_MPI_Cart_map(MPI_Comm  comm, int  ndims, int * dims, int * periods, int * newrank) {
  Proxy_MPI_Cart_map proxyObj(comm, ndims, dims, periods, newrank);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Cart_rank : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int * coords;
  int * rank;
  Proxy_MPI_Cart_rank(MPI_Comm  comm, int * coords, int * rank)
  : comm(comm), coords(coords), rank(rank)
  {}
  void callBack() {
    return_value = ::MPI_Cart_rank(comm, coords, rank);
  }
  int id_for_debug() { return 29; }
};
int  MPIWrapperFunneled::_MPI_Cart_rank(MPI_Comm  comm, int * coords, int * rank) {
  Proxy_MPI_Cart_rank proxyObj(comm, coords, rank);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Cart_shift : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  direction;
  int  disp;
  int * rank_source;
  int * rank_dest;
  Proxy_MPI_Cart_shift(MPI_Comm  comm, int  direction, int  disp, int * rank_source, int * rank_dest)
  : comm(comm), direction(direction), disp(disp), rank_source(rank_source), rank_dest(rank_dest)
  {}
  void callBack() {
    return_value = ::MPI_Cart_shift(comm, direction, disp, rank_source, rank_dest);
  }
  int id_for_debug() { return 30; }
};
int  MPIWrapperFunneled::_MPI_Cart_shift(MPI_Comm  comm, int  direction, int  disp, int * rank_source, int * rank_dest) {
  Proxy_MPI_Cart_shift proxyObj(comm, direction, disp, rank_source, rank_dest);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Cart_sub : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int * remain_dims;
  MPI_Comm * new_comm;
  Proxy_MPI_Cart_sub(MPI_Comm  comm, int * remain_dims, MPI_Comm * new_comm)
  : comm(comm), remain_dims(remain_dims), new_comm(new_comm)
  {}
  void callBack() {
    return_value = ::MPI_Cart_sub(comm, remain_dims, new_comm);
  }
  int id_for_debug() { return 31; }
};
int  MPIWrapperFunneled::_MPI_Cart_sub(MPI_Comm  comm, int * remain_dims, MPI_Comm * new_comm) {
  Proxy_MPI_Cart_sub proxyObj(comm, remain_dims, new_comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Cartdim_get : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int * ndims;
  Proxy_MPI_Cartdim_get(MPI_Comm  comm, int * ndims)
  : comm(comm), ndims(ndims)
  {}
  void callBack() {
    return_value = ::MPI_Cartdim_get(comm, ndims);
  }
  int id_for_debug() { return 32; }
};
int  MPIWrapperFunneled::_MPI_Cartdim_get(MPI_Comm  comm, int * ndims) {
  Proxy_MPI_Cartdim_get proxyObj(comm, ndims);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Close_port : public Proxy_MPI {
  int  return_value;
  char * port_name;
  Proxy_MPI_Close_port(char * port_name)
  : port_name(port_name)
  {}
  void callBack() {
    return_value = ::MPI_Close_port(port_name);
  }
  int id_for_debug() { return 33; }
};
int  MPIWrapperFunneled::_MPI_Close_port(char * port_name) {
  Proxy_MPI_Close_port proxyObj(port_name);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Comm_accept : public Proxy_MPI {
  int  return_value;
  char * port_name;
  MPI_Info  info;
  int  root;
  MPI_Comm  comm;
  MPI_Comm * newcomm;
  Proxy_MPI_Comm_accept(char * port_name, MPI_Info  info, int  root, MPI_Comm  comm, MPI_Comm * newcomm)
  : port_name(port_name), info(info), root(root), comm(comm), newcomm(newcomm)
  {}
  void callBack() {
    return_value = ::MPI_Comm_accept(port_name, info, root, comm, newcomm);
  }
  int id_for_debug() { return 34; }
};
int  MPIWrapperFunneled::_MPI_Comm_accept(char * port_name, MPI_Info  info, int  root, MPI_Comm  comm, MPI_Comm * newcomm) {
  Proxy_MPI_Comm_accept proxyObj(port_name, info, root, comm, newcomm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Comm_call_errhandler : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  errorcode;
  Proxy_MPI_Comm_call_errhandler(MPI_Comm  comm, int  errorcode)
  : comm(comm), errorcode(errorcode)
  {}
  void callBack() {
    return_value = ::MPI_Comm_call_errhandler(comm, errorcode);
  }
  int id_for_debug() { return 36; }
};
int  MPIWrapperFunneled::_MPI_Comm_call_errhandler(MPI_Comm  comm, int  errorcode) {
  Proxy_MPI_Comm_call_errhandler proxyObj(comm, errorcode);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Comm_compare : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm1;
  MPI_Comm  comm2;
  int * result;
  Proxy_MPI_Comm_compare(MPI_Comm  comm1, MPI_Comm  comm2, int * result)
  : comm1(comm1), comm2(comm2), result(result)
  {}
  void callBack() {
    return_value = ::MPI_Comm_compare(comm1, comm2, result);
  }
  int id_for_debug() { return 37; }
};
int  MPIWrapperFunneled::_MPI_Comm_compare(MPI_Comm  comm1, MPI_Comm  comm2, int * result) {
  Proxy_MPI_Comm_compare proxyObj(comm1, comm2, result);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Comm_connect : public Proxy_MPI {
  int  return_value;
  char * port_name;
  MPI_Info  info;
  int  root;
  MPI_Comm  comm;
  MPI_Comm * newcomm;
  Proxy_MPI_Comm_connect(char * port_name, MPI_Info  info, int  root, MPI_Comm  comm, MPI_Comm * newcomm)
  : port_name(port_name), info(info), root(root), comm(comm), newcomm(newcomm)
  {}
  void callBack() {
    return_value = ::MPI_Comm_connect(port_name, info, root, comm, newcomm);
  }
  int id_for_debug() { return 38; }
};
int  MPIWrapperFunneled::_MPI_Comm_connect(char * port_name, MPI_Info  info, int  root, MPI_Comm  comm, MPI_Comm * newcomm) {
  Proxy_MPI_Comm_connect proxyObj(port_name, info, root, comm, newcomm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Comm_create_keyval : public Proxy_MPI {
  int  return_value;
  MPI_Comm_copy_attr_function * comm_copy_attr_fn;
  MPI_Comm_delete_attr_function * comm_delete_attr_fn;
  int * comm_keyval;
  void * extra_state;
  Proxy_MPI_Comm_create_keyval(MPI_Comm_copy_attr_function * comm_copy_attr_fn, MPI_Comm_delete_attr_function * comm_delete_attr_fn, int * comm_keyval, void * extra_state)
  : comm_copy_attr_fn(comm_copy_attr_fn), comm_delete_attr_fn(comm_delete_attr_fn), comm_keyval(comm_keyval), extra_state(extra_state)
  {}
  void callBack() {
    return_value = ::MPI_Comm_create_keyval(comm_copy_attr_fn, comm_delete_attr_fn, comm_keyval, extra_state);
  }
  int id_for_debug() { return 40; }
};
int  MPIWrapperFunneled::_MPI_Comm_create_keyval(MPI_Comm_copy_attr_function * comm_copy_attr_fn, MPI_Comm_delete_attr_function * comm_delete_attr_fn, int * comm_keyval, void * extra_state) {
  Proxy_MPI_Comm_create_keyval proxyObj(comm_copy_attr_fn, comm_delete_attr_fn, comm_keyval, extra_state);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Comm_create : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  MPI_Group  group;
  MPI_Comm * newcomm;
  Proxy_MPI_Comm_create(MPI_Comm  comm, MPI_Group  group, MPI_Comm * newcomm)
  : comm(comm), group(group), newcomm(newcomm)
  {}
  void callBack() {
    return_value = ::MPI_Comm_create(comm, group, newcomm);
  }
  int id_for_debug() { return 41; }
};
int  MPIWrapperFunneled::_MPI_Comm_create(MPI_Comm  comm, MPI_Group  group, MPI_Comm * newcomm) {
  Proxy_MPI_Comm_create proxyObj(comm, group, newcomm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Comm_delete_attr : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  comm_keyval;
  Proxy_MPI_Comm_delete_attr(MPI_Comm  comm, int  comm_keyval)
  : comm(comm), comm_keyval(comm_keyval)
  {}
  void callBack() {
    return_value = ::MPI_Comm_delete_attr(comm, comm_keyval);
  }
  int id_for_debug() { return 42; }
};
int  MPIWrapperFunneled::_MPI_Comm_delete_attr(MPI_Comm  comm, int  comm_keyval) {
  Proxy_MPI_Comm_delete_attr proxyObj(comm, comm_keyval);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Comm_disconnect : public Proxy_MPI {
  int  return_value;
  MPI_Comm * comm;
  Proxy_MPI_Comm_disconnect(MPI_Comm * comm)
  : comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Comm_disconnect(comm);
  }
  int id_for_debug() { return 43; }
};
int  MPIWrapperFunneled::_MPI_Comm_disconnect(MPI_Comm * comm) {
  Proxy_MPI_Comm_disconnect proxyObj(comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Comm_dup : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  MPI_Comm * newcomm;
  Proxy_MPI_Comm_dup(MPI_Comm  comm, MPI_Comm * newcomm)
  : comm(comm), newcomm(newcomm)
  {}
  void callBack() {
    return_value = ::MPI_Comm_dup(comm, newcomm);
  }
  int id_for_debug() { return 44; }
};
int  MPIWrapperFunneled::_MPI_Comm_dup(MPI_Comm  comm, MPI_Comm * newcomm) {
  Proxy_MPI_Comm_dup proxyObj(comm, newcomm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1


struct Proxy_MPI_Comm_free_keyval : public Proxy_MPI {
  int  return_value;
  int * comm_keyval;
  Proxy_MPI_Comm_free_keyval(int * comm_keyval)
  : comm_keyval(comm_keyval)
  {}
  void callBack() {
    return_value = ::MPI_Comm_free_keyval(comm_keyval);
  }
  int id_for_debug() { return 46; }
};
int  MPIWrapperFunneled::_MPI_Comm_free_keyval(int * comm_keyval) {
  Proxy_MPI_Comm_free_keyval proxyObj(comm_keyval);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Comm_free : public Proxy_MPI {
  int  return_value;
  MPI_Comm * comm;
  Proxy_MPI_Comm_free(MPI_Comm * comm)
  : comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Comm_free(comm);
  }
  int id_for_debug() { return 47; }
};
int  MPIWrapperFunneled::_MPI_Comm_free(MPI_Comm * comm) {
  Proxy_MPI_Comm_free proxyObj(comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Comm_get_attr : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  comm_keyval;
  void * attribute_val;
  int * flag;
  Proxy_MPI_Comm_get_attr(MPI_Comm  comm, int  comm_keyval, void * attribute_val, int * flag)
  : comm(comm), comm_keyval(comm_keyval), attribute_val(attribute_val), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Comm_get_attr(comm, comm_keyval, attribute_val, flag);
  }
  int id_for_debug() { return 48; }
};
int  MPIWrapperFunneled::_MPI_Comm_get_attr(MPI_Comm  comm, int  comm_keyval, void * attribute_val, int * flag) {
  Proxy_MPI_Comm_get_attr proxyObj(comm, comm_keyval, attribute_val, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Comm_get_errhandler : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  MPI_Errhandler * erhandler;
  Proxy_MPI_Comm_get_errhandler(MPI_Comm  comm, MPI_Errhandler * erhandler)
  : comm(comm), erhandler(erhandler)
  {}
  void callBack() {
    return_value = ::MPI_Comm_get_errhandler(comm, erhandler);
  }
  int id_for_debug() { return 49; }
};
int  MPIWrapperFunneled::_MPI_Comm_get_errhandler(MPI_Comm  comm, MPI_Errhandler * erhandler) {
  Proxy_MPI_Comm_get_errhandler proxyObj(comm, erhandler);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Comm_get_name : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  char * comm_name;
  int * resultlen;
  Proxy_MPI_Comm_get_name(MPI_Comm  comm, char * comm_name, int * resultlen)
  : comm(comm), comm_name(comm_name), resultlen(resultlen)
  {}
  void callBack() {
    return_value = ::MPI_Comm_get_name(comm, comm_name, resultlen);
  }
  int id_for_debug() { return 50; }
};
int  MPIWrapperFunneled::_MPI_Comm_get_name(MPI_Comm  comm, char * comm_name, int * resultlen) {
  Proxy_MPI_Comm_get_name proxyObj(comm, comm_name, resultlen);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Comm_get_parent : public Proxy_MPI {
  int  return_value;
  MPI_Comm * parent;
  Proxy_MPI_Comm_get_parent(MPI_Comm * parent)
  : parent(parent)
  {}
  void callBack() {
    return_value = ::MPI_Comm_get_parent(parent);
  }
  int id_for_debug() { return 51; }
};
int  MPIWrapperFunneled::_MPI_Comm_get_parent(MPI_Comm * parent) {
  Proxy_MPI_Comm_get_parent proxyObj(parent);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Comm_group : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  MPI_Group * group;
  Proxy_MPI_Comm_group(MPI_Comm  comm, MPI_Group * group)
  : comm(comm), group(group)
  {}
  void callBack() {
    return_value = ::MPI_Comm_group(comm, group);
  }
  int id_for_debug() { return 52; }
};
int  MPIWrapperFunneled::_MPI_Comm_group(MPI_Comm  comm, MPI_Group * group) {
  Proxy_MPI_Comm_group proxyObj(comm, group);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Comm_join : public Proxy_MPI {
  int  return_value;
  int  fd;
  MPI_Comm * intercomm;
  Proxy_MPI_Comm_join(int  fd, MPI_Comm * intercomm)
  : fd(fd), intercomm(intercomm)
  {}
  void callBack() {
    return_value = ::MPI_Comm_join(fd, intercomm);
  }
  int id_for_debug() { return 53; }
};
int  MPIWrapperFunneled::_MPI_Comm_join(int  fd, MPI_Comm * intercomm) {
  Proxy_MPI_Comm_join proxyObj(fd, intercomm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Comm_rank : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int * rank;
  Proxy_MPI_Comm_rank(MPI_Comm  comm, int * rank)
  : comm(comm), rank(rank)
  {}
  void callBack() {
    return_value = ::MPI_Comm_rank(comm, rank);
  }
  int id_for_debug() { return 54; }
};
int  MPIWrapperFunneled::_MPI_Comm_rank(MPI_Comm  comm, int * rank) {
  Proxy_MPI_Comm_rank proxyObj(comm, rank);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Comm_remote_group : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  MPI_Group * group;
  Proxy_MPI_Comm_remote_group(MPI_Comm  comm, MPI_Group * group)
  : comm(comm), group(group)
  {}
  void callBack() {
    return_value = ::MPI_Comm_remote_group(comm, group);
  }
  int id_for_debug() { return 55; }
};
int  MPIWrapperFunneled::_MPI_Comm_remote_group(MPI_Comm  comm, MPI_Group * group) {
  Proxy_MPI_Comm_remote_group proxyObj(comm, group);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Comm_remote_size : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int * size;
  Proxy_MPI_Comm_remote_size(MPI_Comm  comm, int * size)
  : comm(comm), size(size)
  {}
  void callBack() {
    return_value = ::MPI_Comm_remote_size(comm, size);
  }
  int id_for_debug() { return 56; }
};
int  MPIWrapperFunneled::_MPI_Comm_remote_size(MPI_Comm  comm, int * size) {
  Proxy_MPI_Comm_remote_size proxyObj(comm, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Comm_set_attr : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  comm_keyval;
  void * attribute_val;
  Proxy_MPI_Comm_set_attr(MPI_Comm  comm, int  comm_keyval, void * attribute_val)
  : comm(comm), comm_keyval(comm_keyval), attribute_val(attribute_val)
  {}
  void callBack() {
    return_value = ::MPI_Comm_set_attr(comm, comm_keyval, attribute_val);
  }
  int id_for_debug() { return 57; }
};
int  MPIWrapperFunneled::_MPI_Comm_set_attr(MPI_Comm  comm, int  comm_keyval, void * attribute_val) {
  Proxy_MPI_Comm_set_attr proxyObj(comm, comm_keyval, attribute_val);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Comm_set_errhandler : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  MPI_Errhandler  errhandler;
  Proxy_MPI_Comm_set_errhandler(MPI_Comm  comm, MPI_Errhandler  errhandler)
  : comm(comm), errhandler(errhandler)
  {}
  void callBack() {
    return_value = ::MPI_Comm_set_errhandler(comm, errhandler);
  }
  int id_for_debug() { return 58; }
};
int  MPIWrapperFunneled::_MPI_Comm_set_errhandler(MPI_Comm  comm, MPI_Errhandler  errhandler) {
  Proxy_MPI_Comm_set_errhandler proxyObj(comm, errhandler);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Comm_set_name : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  char * comm_name;
  Proxy_MPI_Comm_set_name(MPI_Comm  comm, char * comm_name)
  : comm(comm), comm_name(comm_name)
  {}
  void callBack() {
    return_value = ::MPI_Comm_set_name(comm, comm_name);
  }
  int id_for_debug() { return 59; }
};
int  MPIWrapperFunneled::_MPI_Comm_set_name(MPI_Comm  comm, char * comm_name) {
  Proxy_MPI_Comm_set_name proxyObj(comm, comm_name);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Comm_size : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int * size;
  Proxy_MPI_Comm_size(MPI_Comm  comm, int * size)
  : comm(comm), size(size)
  {}
  void callBack() {
    return_value = ::MPI_Comm_size(comm, size);
  }
  int id_for_debug() { return 60; }
};
int  MPIWrapperFunneled::_MPI_Comm_size(MPI_Comm  comm, int * size) {
  Proxy_MPI_Comm_size proxyObj(comm, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Comm_spawn : public Proxy_MPI {
  int  return_value;
  char * command;
  char ** argv;
  int  maxprocs;
  MPI_Info  info;
  int  root;
  MPI_Comm  comm;
  MPI_Comm * intercomm;
  int * array_of_errcodes;
  Proxy_MPI_Comm_spawn(char * command, char ** argv, int  maxprocs, MPI_Info  info, int  root, MPI_Comm  comm, MPI_Comm * intercomm, int * array_of_errcodes)
  : command(command), argv(argv), maxprocs(maxprocs), info(info), root(root), comm(comm), intercomm(intercomm), array_of_errcodes(array_of_errcodes)
  {}
  void callBack() {
    if (comm == MPI_COMM_NULL)
      throw std::runtime_error("spawn: callBack: comm == MPI_COMM_NULL");
    return_value = ::MPI_Comm_spawn(command, argv, maxprocs, info, root, comm, intercomm, array_of_errcodes);
  }
  int id_for_debug() { return 61; }
};
int  MPIWrapperFunneled::_MPI_Comm_spawn(char * command, char ** argv, int  maxprocs, MPI_Info  info, int  root, MPI_Comm  comm, MPI_Comm * intercomm, int * array_of_errcodes) {
  if (comm == MPI_COMM_NULL)
    throw std::runtime_error("spawn: initially: comm == MPI_COMM_NULL");
  Proxy_MPI_Comm_spawn proxyObj(command, argv, maxprocs, info, root, comm, intercomm, array_of_errcodes);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Comm_spawn_multiple : public Proxy_MPI {
  int  return_value;
  int  count;
  char ** array_of_commands;
  char *** array_of_argv;
  int * array_of_maxprocs;
  MPI_Info * array_of_info;
  int  root;
  MPI_Comm  comm;
  MPI_Comm * intercomm;
  int * array_of_errcodes;
  Proxy_MPI_Comm_spawn_multiple(int  count, char ** array_of_commands, char *** array_of_argv, int * array_of_maxprocs, MPI_Info * array_of_info, int  root, MPI_Comm  comm, MPI_Comm * intercomm, int * array_of_errcodes)
  : count(count), array_of_commands(array_of_commands), array_of_argv(array_of_argv), array_of_maxprocs(array_of_maxprocs), array_of_info(array_of_info), root(root), comm(comm), intercomm(intercomm), array_of_errcodes(array_of_errcodes)
  {}
  void callBack() {
    return_value = ::MPI_Comm_spawn_multiple(count, array_of_commands, array_of_argv, array_of_maxprocs, array_of_info, root, comm, intercomm, array_of_errcodes);
  }
  int id_for_debug() { return 62; }
};
int  MPIWrapperFunneled::_MPI_Comm_spawn_multiple(int  count, char ** array_of_commands, char *** array_of_argv, int * array_of_maxprocs, MPI_Info * array_of_info, int  root, MPI_Comm  comm, MPI_Comm * intercomm, int * array_of_errcodes) {
  Proxy_MPI_Comm_spawn_multiple proxyObj(count, array_of_commands, array_of_argv, array_of_maxprocs, array_of_info, root, comm, intercomm, array_of_errcodes);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Comm_split : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  color;
  int  key;
  MPI_Comm * newcomm;
  Proxy_MPI_Comm_split(MPI_Comm  comm, int  color, int  key, MPI_Comm * newcomm)
  : comm(comm), color(color), key(key), newcomm(newcomm)
  {}
  void callBack() {
    return_value = ::MPI_Comm_split(comm, color, key, newcomm);
  }
  int id_for_debug() { return 63; }
};
int  MPIWrapperFunneled::_MPI_Comm_split(MPI_Comm  comm, int  color, int  key, MPI_Comm * newcomm) {
  Proxy_MPI_Comm_split proxyObj(comm, color, key, newcomm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Comm_test_inter : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int * flag;
  Proxy_MPI_Comm_test_inter(MPI_Comm  comm, int * flag)
  : comm(comm), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Comm_test_inter(comm, flag);
  }
  int id_for_debug() { return 64; }
};
int  MPIWrapperFunneled::_MPI_Comm_test_inter(MPI_Comm  comm, int * flag) {
  Proxy_MPI_Comm_test_inter proxyObj(comm, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Dims_create : public Proxy_MPI {
  int  return_value;
  int  nnodes;
  int  ndims;
  int * dims;
  Proxy_MPI_Dims_create(int  nnodes, int  ndims, int * dims)
  : nnodes(nnodes), ndims(ndims), dims(dims)
  {}
  void callBack() {
    return_value = ::MPI_Dims_create(nnodes, ndims, dims);
  }
  int id_for_debug() { return 65; }
};
int  MPIWrapperFunneled::_MPI_Dims_create(int  nnodes, int  ndims, int * dims) {
  Proxy_MPI_Dims_create proxyObj(nnodes, ndims, dims);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Errhandler_free : public Proxy_MPI {
  int  return_value;
  MPI_Errhandler * errhandler;
  Proxy_MPI_Errhandler_free(MPI_Errhandler * errhandler)
  : errhandler(errhandler)
  {}
  void callBack() {
    return_value = ::MPI_Errhandler_free(errhandler);
  }
  int id_for_debug() { return 69; }
};
int  MPIWrapperFunneled::_MPI_Errhandler_free(MPI_Errhandler * errhandler) {
  Proxy_MPI_Errhandler_free proxyObj(errhandler);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Error_class : public Proxy_MPI {
  int  return_value;
  int  errorcode;
  int * errorclass;
  Proxy_MPI_Error_class(int  errorcode, int * errorclass)
  : errorcode(errorcode), errorclass(errorclass)
  {}
  void callBack() {
    return_value = ::MPI_Error_class(errorcode, errorclass);
  }
  int id_for_debug() { return 72; }
};
int  MPIWrapperFunneled::_MPI_Error_class(int  errorcode, int * errorclass) {
  Proxy_MPI_Error_class proxyObj(errorcode, errorclass);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Error_string : public Proxy_MPI {
  int  return_value;
  int  errorcode;
  char * string;
  int * resultlen;
  Proxy_MPI_Error_string(int  errorcode, char * string, int * resultlen)
  : errorcode(errorcode), string(string), resultlen(resultlen)
  {}
  void callBack() {
    return_value = ::MPI_Error_string(errorcode, string, resultlen);
  }
  int id_for_debug() { return 73; }
};
int  MPIWrapperFunneled::_MPI_Error_string(int  errorcode, char * string, int * resultlen) {
  Proxy_MPI_Error_string proxyObj(errorcode, string, resultlen);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Exscan : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  void * recvbuf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Op  op;
  MPI_Comm  comm;
  Proxy_MPI_Exscan(void * sendbuf, void * recvbuf, int  count, MPI_Datatype  datatype, MPI_Op  op, MPI_Comm  comm)
  : sendbuf(sendbuf), recvbuf(recvbuf), count(count), datatype(datatype), op(op), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Exscan(sendbuf, recvbuf, count, datatype, op, comm);
  }
  int id_for_debug() { return 74; }
};
int  MPIWrapperFunneled::_MPI_Exscan(void * sendbuf, void * recvbuf, int  count, MPI_Datatype  datatype, MPI_Op  op, MPI_Comm  comm) {
  Proxy_MPI_Exscan proxyObj(sendbuf, recvbuf, count, datatype, op, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_c2f : public Proxy_MPI {
  MPI_Fint  return_value;
  MPI_File  file;
  Proxy_MPI_File_c2f(MPI_File  file)
  : file(file)
  {}
  void callBack() {
    return_value = ::MPI_File_c2f(file);
  }
  int id_for_debug() { return 75; }
};
MPI_Fint  MPIWrapperFunneled::_MPI_File_c2f(MPI_File  file) {
  Proxy_MPI_File_c2f proxyObj(file);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_f2c : public Proxy_MPI {
  MPI_File  return_value;
  MPI_Fint  file;
  Proxy_MPI_File_f2c(MPI_Fint  file)
  : file(file)
  {}
  void callBack() {
    return_value = ::MPI_File_f2c(file);
  }
  int id_for_debug() { return 76; }
};
MPI_File  MPIWrapperFunneled::_MPI_File_f2c(MPI_Fint  file) {
  Proxy_MPI_File_f2c proxyObj(file);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_call_errhandler : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  int  errorcode;
  Proxy_MPI_File_call_errhandler(MPI_File  fh, int  errorcode)
  : fh(fh), errorcode(errorcode)
  {}
  void callBack() {
    return_value = ::MPI_File_call_errhandler(fh, errorcode);
  }
  int id_for_debug() { return 77; }
};
int  MPIWrapperFunneled::_MPI_File_call_errhandler(MPI_File  fh, int  errorcode) {
  Proxy_MPI_File_call_errhandler proxyObj(fh, errorcode);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_set_errhandler : public Proxy_MPI {
  int  return_value;
  MPI_File  file;
  MPI_Errhandler  errhandler;
  Proxy_MPI_File_set_errhandler(MPI_File  file, MPI_Errhandler  errhandler)
  : file(file), errhandler(errhandler)
  {}
  void callBack() {
    return_value = ::MPI_File_set_errhandler(file, errhandler);
  }
  int id_for_debug() { return 79; }
};
int  MPIWrapperFunneled::_MPI_File_set_errhandler(MPI_File  file, MPI_Errhandler  errhandler) {
  Proxy_MPI_File_set_errhandler proxyObj(file, errhandler);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_errhandler : public Proxy_MPI {
  int  return_value;
  MPI_File  file;
  MPI_Errhandler * errhandler;
  Proxy_MPI_File_get_errhandler(MPI_File  file, MPI_Errhandler * errhandler)
  : file(file), errhandler(errhandler)
  {}
  void callBack() {
    return_value = ::MPI_File_get_errhandler(file, errhandler);
  }
  int id_for_debug() { return 80; }
};
int  MPIWrapperFunneled::_MPI_File_get_errhandler(MPI_File  file, MPI_Errhandler * errhandler) {
  Proxy_MPI_File_get_errhandler proxyObj(file, errhandler);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_open : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  char * filename;
  int  amode;
  MPI_Info  info;
  MPI_File * fh;
  Proxy_MPI_File_open(MPI_Comm  comm, char * filename, int  amode, MPI_Info  info, MPI_File * fh)
  : comm(comm), filename(filename), amode(amode), info(info), fh(fh)
  {}
  void callBack() {
    return_value = ::MPI_File_open(comm, filename, amode, info, fh);
  }
  int id_for_debug() { return 81; }
};
int  MPIWrapperFunneled::_MPI_File_open(MPI_Comm  comm, char * filename, int  amode, MPI_Info  info, MPI_File * fh) {
  Proxy_MPI_File_open proxyObj(comm, filename, amode, info, fh);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_close : public Proxy_MPI {
  int  return_value;
  MPI_File * fh;
  Proxy_MPI_File_close(MPI_File * fh)
  : fh(fh)
  {}
  void callBack() {
    return_value = ::MPI_File_close(fh);
  }
  int id_for_debug() { return 82; }
};
int  MPIWrapperFunneled::_MPI_File_close(MPI_File * fh) {
  Proxy_MPI_File_close proxyObj(fh);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_delete : public Proxy_MPI {
  int  return_value;
  char * filename;
  MPI_Info  info;
  Proxy_MPI_File_delete(char * filename, MPI_Info  info)
  : filename(filename), info(info)
  {}
  void callBack() {
    return_value = ::MPI_File_delete(filename, info);
  }
  int id_for_debug() { return 83; }
};
int  MPIWrapperFunneled::_MPI_File_delete(char * filename, MPI_Info  info) {
  Proxy_MPI_File_delete proxyObj(filename, info);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_set_size : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  size;
  Proxy_MPI_File_set_size(MPI_File  fh, MPI_Offset  size)
  : fh(fh), size(size)
  {}
  void callBack() {
    return_value = ::MPI_File_set_size(fh, size);
  }
  int id_for_debug() { return 84; }
};
int  MPIWrapperFunneled::_MPI_File_set_size(MPI_File  fh, MPI_Offset  size) {
  Proxy_MPI_File_set_size proxyObj(fh, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_preallocate : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  size;
  Proxy_MPI_File_preallocate(MPI_File  fh, MPI_Offset  size)
  : fh(fh), size(size)
  {}
  void callBack() {
    return_value = ::MPI_File_preallocate(fh, size);
  }
  int id_for_debug() { return 85; }
};
int  MPIWrapperFunneled::_MPI_File_preallocate(MPI_File  fh, MPI_Offset  size) {
  Proxy_MPI_File_preallocate proxyObj(fh, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_size : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset * size;
  Proxy_MPI_File_get_size(MPI_File  fh, MPI_Offset * size)
  : fh(fh), size(size)
  {}
  void callBack() {
    return_value = ::MPI_File_get_size(fh, size);
  }
  int id_for_debug() { return 86; }
};
int  MPIWrapperFunneled::_MPI_File_get_size(MPI_File  fh, MPI_Offset * size) {
  Proxy_MPI_File_get_size proxyObj(fh, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_group : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Group * group;
  Proxy_MPI_File_get_group(MPI_File  fh, MPI_Group * group)
  : fh(fh), group(group)
  {}
  void callBack() {
    return_value = ::MPI_File_get_group(fh, group);
  }
  int id_for_debug() { return 87; }
};
int  MPIWrapperFunneled::_MPI_File_get_group(MPI_File  fh, MPI_Group * group) {
  Proxy_MPI_File_get_group proxyObj(fh, group);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_amode : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  int * amode;
  Proxy_MPI_File_get_amode(MPI_File  fh, int * amode)
  : fh(fh), amode(amode)
  {}
  void callBack() {
    return_value = ::MPI_File_get_amode(fh, amode);
  }
  int id_for_debug() { return 88; }
};
int  MPIWrapperFunneled::_MPI_File_get_amode(MPI_File  fh, int * amode) {
  Proxy_MPI_File_get_amode proxyObj(fh, amode);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_set_info : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Info  info;
  Proxy_MPI_File_set_info(MPI_File  fh, MPI_Info  info)
  : fh(fh), info(info)
  {}
  void callBack() {
    return_value = ::MPI_File_set_info(fh, info);
  }
  int id_for_debug() { return 89; }
};
int  MPIWrapperFunneled::_MPI_File_set_info(MPI_File  fh, MPI_Info  info) {
  Proxy_MPI_File_set_info proxyObj(fh, info);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_info : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Info * info_used;
  Proxy_MPI_File_get_info(MPI_File  fh, MPI_Info * info_used)
  : fh(fh), info_used(info_used)
  {}
  void callBack() {
    return_value = ::MPI_File_get_info(fh, info_used);
  }
  int id_for_debug() { return 90; }
};
int  MPIWrapperFunneled::_MPI_File_get_info(MPI_File  fh, MPI_Info * info_used) {
  Proxy_MPI_File_get_info proxyObj(fh, info_used);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_set_view : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  disp;
  MPI_Datatype  etype;
  MPI_Datatype  filetype;
  char * datarep;
  MPI_Info  info;
  Proxy_MPI_File_set_view(MPI_File  fh, MPI_Offset  disp, MPI_Datatype  etype, MPI_Datatype  filetype, char * datarep, MPI_Info  info)
  : fh(fh), disp(disp), etype(etype), filetype(filetype), datarep(datarep), info(info)
  {}
  void callBack() {
    return_value = ::MPI_File_set_view(fh, disp, etype, filetype, datarep, info);
  }
  int id_for_debug() { return 91; }
};
int  MPIWrapperFunneled::_MPI_File_set_view(MPI_File  fh, MPI_Offset  disp, MPI_Datatype  etype, MPI_Datatype  filetype, char * datarep, MPI_Info  info) {
  Proxy_MPI_File_set_view proxyObj(fh, disp, etype, filetype, datarep, info);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_view : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset * disp;
  MPI_Datatype * etype;
  MPI_Datatype * filetype;
  char * datarep;
  Proxy_MPI_File_get_view(MPI_File  fh, MPI_Offset * disp, MPI_Datatype * etype, MPI_Datatype * filetype, char * datarep)
  : fh(fh), disp(disp), etype(etype), filetype(filetype), datarep(datarep)
  {}
  void callBack() {
    return_value = ::MPI_File_get_view(fh, disp, etype, filetype, datarep);
  }
  int id_for_debug() { return 92; }
};
int  MPIWrapperFunneled::_MPI_File_get_view(MPI_File  fh, MPI_Offset * disp, MPI_Datatype * etype, MPI_Datatype * filetype, char * datarep) {
  Proxy_MPI_File_get_view proxyObj(fh, disp, etype, filetype, datarep);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_at : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_read_at(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), offset(offset), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_read_at(fh, offset, buf, count, datatype, status);
  }
  int id_for_debug() { return 93; }
};
int  MPIWrapperFunneled::_MPI_File_read_at(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_read_at proxyObj(fh, offset, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_at_all : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_read_at_all(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), offset(offset), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_read_at_all(fh, offset, buf, count, datatype, status);
  }
  int id_for_debug() { return 94; }
};
int  MPIWrapperFunneled::_MPI_File_read_at_all(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_read_at_all proxyObj(fh, offset, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_at : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_write_at(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), offset(offset), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_write_at(fh, offset, buf, count, datatype, status);
  }
  int id_for_debug() { return 95; }
};
int  MPIWrapperFunneled::_MPI_File_write_at(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_write_at proxyObj(fh, offset, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_at_all : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_write_at_all(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), offset(offset), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_write_at_all(fh, offset, buf, count, datatype, status);
  }
  int id_for_debug() { return 96; }
};
int  MPIWrapperFunneled::_MPI_File_write_at_all(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_write_at_all proxyObj(fh, offset, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_iread_at : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Request * request;
  Proxy_MPI_File_iread_at(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request)
  : fh(fh), offset(offset), buf(buf), count(count), datatype(datatype), request(request)
  {}
  void callBack() {
    return_value = ::MPI_File_iread_at(fh, offset, buf, count, datatype, request);
  }
  int id_for_debug() { return 97; }
};
int  MPIWrapperFunneled::_MPI_File_iread_at(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request) {
  Proxy_MPI_File_iread_at proxyObj(fh, offset, buf, count, datatype, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_iwrite_at : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Request * request;
  Proxy_MPI_File_iwrite_at(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request)
  : fh(fh), offset(offset), buf(buf), count(count), datatype(datatype), request(request)
  {}
  void callBack() {
    return_value = ::MPI_File_iwrite_at(fh, offset, buf, count, datatype, request);
  }
  int id_for_debug() { return 98; }
};
int  MPIWrapperFunneled::_MPI_File_iwrite_at(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request) {
  Proxy_MPI_File_iwrite_at proxyObj(fh, offset, buf, count, datatype, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_read(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_read(fh, buf, count, datatype, status);
  }
  int id_for_debug() { return 99; }
};
int  MPIWrapperFunneled::_MPI_File_read(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_read proxyObj(fh, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_all : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_read_all(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_read_all(fh, buf, count, datatype, status);
  }
  int id_for_debug() { return 100; }
};
int  MPIWrapperFunneled::_MPI_File_read_all(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_read_all proxyObj(fh, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_write(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_write(fh, buf, count, datatype, status);
  }
  int id_for_debug() { return 101; }
};
int  MPIWrapperFunneled::_MPI_File_write(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_write proxyObj(fh, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_all : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_write_all(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_write_all(fh, buf, count, datatype, status);
  }
  int id_for_debug() { return 102; }
};
int  MPIWrapperFunneled::_MPI_File_write_all(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_write_all proxyObj(fh, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_iread : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Request * request;
  Proxy_MPI_File_iread(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request)
  : fh(fh), buf(buf), count(count), datatype(datatype), request(request)
  {}
  void callBack() {
    return_value = ::MPI_File_iread(fh, buf, count, datatype, request);
  }
  int id_for_debug() { return 103; }
};
int  MPIWrapperFunneled::_MPI_File_iread(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request) {
  Proxy_MPI_File_iread proxyObj(fh, buf, count, datatype, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_iwrite : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Request * request;
  Proxy_MPI_File_iwrite(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request)
  : fh(fh), buf(buf), count(count), datatype(datatype), request(request)
  {}
  void callBack() {
    return_value = ::MPI_File_iwrite(fh, buf, count, datatype, request);
  }
  int id_for_debug() { return 104; }
};
int  MPIWrapperFunneled::_MPI_File_iwrite(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request) {
  Proxy_MPI_File_iwrite proxyObj(fh, buf, count, datatype, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_seek : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  int  whence;
  Proxy_MPI_File_seek(MPI_File  fh, MPI_Offset  offset, int  whence)
  : fh(fh), offset(offset), whence(whence)
  {}
  void callBack() {
    return_value = ::MPI_File_seek(fh, offset, whence);
  }
  int id_for_debug() { return 105; }
};
int  MPIWrapperFunneled::_MPI_File_seek(MPI_File  fh, MPI_Offset  offset, int  whence) {
  Proxy_MPI_File_seek proxyObj(fh, offset, whence);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_position : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset * offset;
  Proxy_MPI_File_get_position(MPI_File  fh, MPI_Offset * offset)
  : fh(fh), offset(offset)
  {}
  void callBack() {
    return_value = ::MPI_File_get_position(fh, offset);
  }
  int id_for_debug() { return 106; }
};
int  MPIWrapperFunneled::_MPI_File_get_position(MPI_File  fh, MPI_Offset * offset) {
  Proxy_MPI_File_get_position proxyObj(fh, offset);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_byte_offset : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  MPI_Offset * disp;
  Proxy_MPI_File_get_byte_offset(MPI_File  fh, MPI_Offset  offset, MPI_Offset * disp)
  : fh(fh), offset(offset), disp(disp)
  {}
  void callBack() {
    return_value = ::MPI_File_get_byte_offset(fh, offset, disp);
  }
  int id_for_debug() { return 107; }
};
int  MPIWrapperFunneled::_MPI_File_get_byte_offset(MPI_File  fh, MPI_Offset  offset, MPI_Offset * disp) {
  Proxy_MPI_File_get_byte_offset proxyObj(fh, offset, disp);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_shared : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_read_shared(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_read_shared(fh, buf, count, datatype, status);
  }
  int id_for_debug() { return 108; }
};
int  MPIWrapperFunneled::_MPI_File_read_shared(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_read_shared proxyObj(fh, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_shared : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_write_shared(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_write_shared(fh, buf, count, datatype, status);
  }
  int id_for_debug() { return 109; }
};
int  MPIWrapperFunneled::_MPI_File_write_shared(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_write_shared proxyObj(fh, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_iread_shared : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Request * request;
  Proxy_MPI_File_iread_shared(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request)
  : fh(fh), buf(buf), count(count), datatype(datatype), request(request)
  {}
  void callBack() {
    return_value = ::MPI_File_iread_shared(fh, buf, count, datatype, request);
  }
  int id_for_debug() { return 110; }
};
int  MPIWrapperFunneled::_MPI_File_iread_shared(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request) {
  Proxy_MPI_File_iread_shared proxyObj(fh, buf, count, datatype, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_iwrite_shared : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Request * request;
  Proxy_MPI_File_iwrite_shared(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request)
  : fh(fh), buf(buf), count(count), datatype(datatype), request(request)
  {}
  void callBack() {
    return_value = ::MPI_File_iwrite_shared(fh, buf, count, datatype, request);
  }
  int id_for_debug() { return 111; }
};
int  MPIWrapperFunneled::_MPI_File_iwrite_shared(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Request * request) {
  Proxy_MPI_File_iwrite_shared proxyObj(fh, buf, count, datatype, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_ordered : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_read_ordered(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_read_ordered(fh, buf, count, datatype, status);
  }
  int id_for_debug() { return 112; }
};
int  MPIWrapperFunneled::_MPI_File_read_ordered(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_read_ordered proxyObj(fh, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_ordered : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Status * status;
  Proxy_MPI_File_write_ordered(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status)
  : fh(fh), buf(buf), count(count), datatype(datatype), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_write_ordered(fh, buf, count, datatype, status);
  }
  int id_for_debug() { return 113; }
};
int  MPIWrapperFunneled::_MPI_File_write_ordered(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype, MPI_Status * status) {
  Proxy_MPI_File_write_ordered proxyObj(fh, buf, count, datatype, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_seek_shared : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  int  whence;
  Proxy_MPI_File_seek_shared(MPI_File  fh, MPI_Offset  offset, int  whence)
  : fh(fh), offset(offset), whence(whence)
  {}
  void callBack() {
    return_value = ::MPI_File_seek_shared(fh, offset, whence);
  }
  int id_for_debug() { return 114; }
};
int  MPIWrapperFunneled::_MPI_File_seek_shared(MPI_File  fh, MPI_Offset  offset, int  whence) {
  Proxy_MPI_File_seek_shared proxyObj(fh, offset, whence);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_position_shared : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset * offset;
  Proxy_MPI_File_get_position_shared(MPI_File  fh, MPI_Offset * offset)
  : fh(fh), offset(offset)
  {}
  void callBack() {
    return_value = ::MPI_File_get_position_shared(fh, offset);
  }
  int id_for_debug() { return 115; }
};
int  MPIWrapperFunneled::_MPI_File_get_position_shared(MPI_File  fh, MPI_Offset * offset) {
  Proxy_MPI_File_get_position_shared proxyObj(fh, offset);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_at_all_begin : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  Proxy_MPI_File_read_at_all_begin(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype)
  : fh(fh), offset(offset), buf(buf), count(count), datatype(datatype)
  {}
  void callBack() {
    return_value = ::MPI_File_read_at_all_begin(fh, offset, buf, count, datatype);
  }
  int id_for_debug() { return 116; }
};
int  MPIWrapperFunneled::_MPI_File_read_at_all_begin(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype) {
  Proxy_MPI_File_read_at_all_begin proxyObj(fh, offset, buf, count, datatype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_at_all_end : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  MPI_Status * status;
  Proxy_MPI_File_read_at_all_end(MPI_File  fh, void * buf, MPI_Status * status)
  : fh(fh), buf(buf), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_read_at_all_end(fh, buf, status);
  }
  int id_for_debug() { return 117; }
};
int  MPIWrapperFunneled::_MPI_File_read_at_all_end(MPI_File  fh, void * buf, MPI_Status * status) {
  Proxy_MPI_File_read_at_all_end proxyObj(fh, buf, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_at_all_begin : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Offset  offset;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  Proxy_MPI_File_write_at_all_begin(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype)
  : fh(fh), offset(offset), buf(buf), count(count), datatype(datatype)
  {}
  void callBack() {
    return_value = ::MPI_File_write_at_all_begin(fh, offset, buf, count, datatype);
  }
  int id_for_debug() { return 118; }
};
int  MPIWrapperFunneled::_MPI_File_write_at_all_begin(MPI_File  fh, MPI_Offset  offset, void * buf, int  count, MPI_Datatype  datatype) {
  Proxy_MPI_File_write_at_all_begin proxyObj(fh, offset, buf, count, datatype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_at_all_end : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  MPI_Status * status;
  Proxy_MPI_File_write_at_all_end(MPI_File  fh, void * buf, MPI_Status * status)
  : fh(fh), buf(buf), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_write_at_all_end(fh, buf, status);
  }
  int id_for_debug() { return 119; }
};
int  MPIWrapperFunneled::_MPI_File_write_at_all_end(MPI_File  fh, void * buf, MPI_Status * status) {
  Proxy_MPI_File_write_at_all_end proxyObj(fh, buf, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_all_begin : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  Proxy_MPI_File_read_all_begin(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype)
  : fh(fh), buf(buf), count(count), datatype(datatype)
  {}
  void callBack() {
    return_value = ::MPI_File_read_all_begin(fh, buf, count, datatype);
  }
  int id_for_debug() { return 120; }
};
int  MPIWrapperFunneled::_MPI_File_read_all_begin(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype) {
  Proxy_MPI_File_read_all_begin proxyObj(fh, buf, count, datatype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_all_end : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  MPI_Status * status;
  Proxy_MPI_File_read_all_end(MPI_File  fh, void * buf, MPI_Status * status)
  : fh(fh), buf(buf), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_read_all_end(fh, buf, status);
  }
  int id_for_debug() { return 121; }
};
int  MPIWrapperFunneled::_MPI_File_read_all_end(MPI_File  fh, void * buf, MPI_Status * status) {
  Proxy_MPI_File_read_all_end proxyObj(fh, buf, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_all_begin : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  Proxy_MPI_File_write_all_begin(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype)
  : fh(fh), buf(buf), count(count), datatype(datatype)
  {}
  void callBack() {
    return_value = ::MPI_File_write_all_begin(fh, buf, count, datatype);
  }
  int id_for_debug() { return 122; }
};
int  MPIWrapperFunneled::_MPI_File_write_all_begin(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype) {
  Proxy_MPI_File_write_all_begin proxyObj(fh, buf, count, datatype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_all_end : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  MPI_Status * status;
  Proxy_MPI_File_write_all_end(MPI_File  fh, void * buf, MPI_Status * status)
  : fh(fh), buf(buf), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_write_all_end(fh, buf, status);
  }
  int id_for_debug() { return 123; }
};
int  MPIWrapperFunneled::_MPI_File_write_all_end(MPI_File  fh, void * buf, MPI_Status * status) {
  Proxy_MPI_File_write_all_end proxyObj(fh, buf, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_ordered_begin : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  Proxy_MPI_File_read_ordered_begin(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype)
  : fh(fh), buf(buf), count(count), datatype(datatype)
  {}
  void callBack() {
    return_value = ::MPI_File_read_ordered_begin(fh, buf, count, datatype);
  }
  int id_for_debug() { return 124; }
};
int  MPIWrapperFunneled::_MPI_File_read_ordered_begin(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype) {
  Proxy_MPI_File_read_ordered_begin proxyObj(fh, buf, count, datatype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_read_ordered_end : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  MPI_Status * status;
  Proxy_MPI_File_read_ordered_end(MPI_File  fh, void * buf, MPI_Status * status)
  : fh(fh), buf(buf), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_read_ordered_end(fh, buf, status);
  }
  int id_for_debug() { return 125; }
};
int  MPIWrapperFunneled::_MPI_File_read_ordered_end(MPI_File  fh, void * buf, MPI_Status * status) {
  Proxy_MPI_File_read_ordered_end proxyObj(fh, buf, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_ordered_begin : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  Proxy_MPI_File_write_ordered_begin(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype)
  : fh(fh), buf(buf), count(count), datatype(datatype)
  {}
  void callBack() {
    return_value = ::MPI_File_write_ordered_begin(fh, buf, count, datatype);
  }
  int id_for_debug() { return 126; }
};
int  MPIWrapperFunneled::_MPI_File_write_ordered_begin(MPI_File  fh, void * buf, int  count, MPI_Datatype  datatype) {
  Proxy_MPI_File_write_ordered_begin proxyObj(fh, buf, count, datatype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_write_ordered_end : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  void * buf;
  MPI_Status * status;
  Proxy_MPI_File_write_ordered_end(MPI_File  fh, void * buf, MPI_Status * status)
  : fh(fh), buf(buf), status(status)
  {}
  void callBack() {
    return_value = ::MPI_File_write_ordered_end(fh, buf, status);
  }
  int id_for_debug() { return 127; }
};
int  MPIWrapperFunneled::_MPI_File_write_ordered_end(MPI_File  fh, void * buf, MPI_Status * status) {
  Proxy_MPI_File_write_ordered_end proxyObj(fh, buf, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_type_extent : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  MPI_Datatype  datatype;
  MPI_Aint * extent;
  Proxy_MPI_File_get_type_extent(MPI_File  fh, MPI_Datatype  datatype, MPI_Aint * extent)
  : fh(fh), datatype(datatype), extent(extent)
  {}
  void callBack() {
    return_value = ::MPI_File_get_type_extent(fh, datatype, extent);
  }
  int id_for_debug() { return 128; }
};
int  MPIWrapperFunneled::_MPI_File_get_type_extent(MPI_File  fh, MPI_Datatype  datatype, MPI_Aint * extent) {
  Proxy_MPI_File_get_type_extent proxyObj(fh, datatype, extent);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_set_atomicity : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  int  flag;
  Proxy_MPI_File_set_atomicity(MPI_File  fh, int  flag)
  : fh(fh), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_File_set_atomicity(fh, flag);
  }
  int id_for_debug() { return 129; }
};
int  MPIWrapperFunneled::_MPI_File_set_atomicity(MPI_File  fh, int  flag) {
  Proxy_MPI_File_set_atomicity proxyObj(fh, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_get_atomicity : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  int * flag;
  Proxy_MPI_File_get_atomicity(MPI_File  fh, int * flag)
  : fh(fh), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_File_get_atomicity(fh, flag);
  }
  int id_for_debug() { return 130; }
};
int  MPIWrapperFunneled::_MPI_File_get_atomicity(MPI_File  fh, int * flag) {
  Proxy_MPI_File_get_atomicity proxyObj(fh, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_File_sync : public Proxy_MPI {
  int  return_value;
  MPI_File  fh;
  Proxy_MPI_File_sync(MPI_File  fh)
  : fh(fh)
  {}
  void callBack() {
    return_value = ::MPI_File_sync(fh);
  }
  int id_for_debug() { return 131; }
};
int  MPIWrapperFunneled::_MPI_File_sync(MPI_File  fh) {
  Proxy_MPI_File_sync proxyObj(fh);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif




struct Proxy_MPI_Finalized : public Proxy_MPI {
  int  return_value;
  int * flag;
  Proxy_MPI_Finalized(int * flag)
  : flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Finalized(flag);
  }
  int id_for_debug() { return 132; }
};
int  MPIWrapperFunneled::_MPI_Finalized(int * flag) {
  Proxy_MPI_Finalized proxyObj(flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Free_mem : public Proxy_MPI {
  int  return_value;
  void * base;
  Proxy_MPI_Free_mem(void * base)
  : base(base)
  {}
  void callBack() {
    return_value = ::MPI_Free_mem(base);
  }
  int id_for_debug() { return 133; }
};
int  MPIWrapperFunneled::_MPI_Free_mem(void * base) {
  Proxy_MPI_Free_mem proxyObj(base);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Gather : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  int  sendcount;
  MPI_Datatype  sendtype;
  void * recvbuf;
  int  recvcount;
  MPI_Datatype  recvtype;
  int  root;
  MPI_Comm  comm;
  Proxy_MPI_Gather(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, int  root, MPI_Comm  comm)
  : sendbuf(sendbuf), sendcount(sendcount), sendtype(sendtype), recvbuf(recvbuf), recvcount(recvcount), recvtype(recvtype), root(root), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
  int id_for_debug() { return 134; }
};
int  MPIWrapperFunneled::_MPI_Gather(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, int  root, MPI_Comm  comm) {
  Proxy_MPI_Gather proxyObj(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Gatherv : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  int  sendcount;
  MPI_Datatype  sendtype;
  void * recvbuf;
  int * recvcounts;
  int * displs;
  MPI_Datatype  recvtype;
  int  root;
  MPI_Comm  comm;
  Proxy_MPI_Gatherv(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int * recvcounts, int * displs, MPI_Datatype  recvtype, int  root, MPI_Comm  comm)
  : sendbuf(sendbuf), sendcount(sendcount), sendtype(sendtype), recvbuf(recvbuf), recvcounts(recvcounts), displs(displs), recvtype(recvtype), root(root), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
  }
  int id_for_debug() { return 135; }
};
int  MPIWrapperFunneled::_MPI_Gatherv(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int * recvcounts, int * displs, MPI_Datatype  recvtype, int  root, MPI_Comm  comm) {
  Proxy_MPI_Gatherv proxyObj(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Get_address : public Proxy_MPI {
  int  return_value;
  void * location;
  MPI_Aint * address;
  Proxy_MPI_Get_address(void * location, MPI_Aint * address)
  : location(location), address(address)
  {}
  void callBack() {
    return_value = ::MPI_Get_address(location, address);
  }
  int id_for_debug() { return 136; }
};
int  MPIWrapperFunneled::_MPI_Get_address(void * location, MPI_Aint * address) {
  Proxy_MPI_Get_address proxyObj(location, address);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Get_count : public Proxy_MPI {
  int  return_value;
  MPI_Status * status;
  MPI_Datatype  datatype;
  int * count;
  Proxy_MPI_Get_count(MPI_Status * status, MPI_Datatype  datatype, int * count)
  : status(status), datatype(datatype), count(count)
  {}
  void callBack() {
    return_value = ::MPI_Get_count(status, datatype, count);
  }
  int id_for_debug() { return 137; }
};
int  MPIWrapperFunneled::_MPI_Get_count(MPI_Status * status, MPI_Datatype  datatype, int * count) {
  Proxy_MPI_Get_count proxyObj(status, datatype, count);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Get_elements : public Proxy_MPI {
  int  return_value;
  MPI_Status * status;
  MPI_Datatype  datatype;
  int * count;
  Proxy_MPI_Get_elements(MPI_Status * status, MPI_Datatype  datatype, int * count)
  : status(status), datatype(datatype), count(count)
  {}
  void callBack() {
    return_value = ::MPI_Get_elements(status, datatype, count);
  }
  int id_for_debug() { return 138; }
};
int  MPIWrapperFunneled::_MPI_Get_elements(MPI_Status * status, MPI_Datatype  datatype, int * count) {
  Proxy_MPI_Get_elements proxyObj(status, datatype, count);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Get : public Proxy_MPI {
  int  return_value;
  void * origin_addr;
  int  origin_count;
  MPI_Datatype  origin_datatype;
  int  target_rank;
  MPI_Aint  target_disp;
  int  target_count;
  MPI_Datatype  target_datatype;
  MPI_Win  win;
  Proxy_MPI_Get(void * origin_addr, int  origin_count, MPI_Datatype  origin_datatype, int  target_rank, MPI_Aint  target_disp, int  target_count, MPI_Datatype  target_datatype, MPI_Win  win)
  : origin_addr(origin_addr), origin_count(origin_count), origin_datatype(origin_datatype), target_rank(target_rank), target_disp(target_disp), target_count(target_count), target_datatype(target_datatype), win(win)
  {}
  void callBack() {
    return_value = ::MPI_Get(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, win);
  }
  int id_for_debug() { return 139; }
};
int  MPIWrapperFunneled::_MPI_Get(void * origin_addr, int  origin_count, MPI_Datatype  origin_datatype, int  target_rank, MPI_Aint  target_disp, int  target_count, MPI_Datatype  target_datatype, MPI_Win  win) {
  Proxy_MPI_Get proxyObj(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Get_processor_name : public Proxy_MPI {
  int  return_value;
  char * name;
  int * resultlen;
  Proxy_MPI_Get_processor_name(char * name, int * resultlen)
  : name(name), resultlen(resultlen)
  {}
  void callBack() {
    return_value = ::MPI_Get_processor_name(name, resultlen);
  }
  int id_for_debug() { return 140; }
};
int  MPIWrapperFunneled::_MPI_Get_processor_name(char * name, int * resultlen) {
  Proxy_MPI_Get_processor_name proxyObj(name, resultlen);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Get_version : public Proxy_MPI {
  int  return_value;
  int * version;
  int * subversion;
  Proxy_MPI_Get_version(int * version, int * subversion)
  : version(version), subversion(subversion)
  {}
  void callBack() {
    return_value = ::MPI_Get_version(version, subversion);
  }
  int id_for_debug() { return 141; }
};
int  MPIWrapperFunneled::_MPI_Get_version(int * version, int * subversion) {
  Proxy_MPI_Get_version proxyObj(version, subversion);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Graph_create : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm_old;
  int  nnodes;
  int * index;
  int * edges;
  int  reorder;
  MPI_Comm * comm_graph;
  Proxy_MPI_Graph_create(MPI_Comm  comm_old, int  nnodes, int * index, int * edges, int  reorder, MPI_Comm * comm_graph)
  : comm_old(comm_old), nnodes(nnodes), index(index), edges(edges), reorder(reorder), comm_graph(comm_graph)
  {}
  void callBack() {
    return_value = ::MPI_Graph_create(comm_old, nnodes, index, edges, reorder, comm_graph);
  }
  int id_for_debug() { return 142; }
};
int  MPIWrapperFunneled::_MPI_Graph_create(MPI_Comm  comm_old, int  nnodes, int * index, int * edges, int  reorder, MPI_Comm * comm_graph) {
  Proxy_MPI_Graph_create proxyObj(comm_old, nnodes, index, edges, reorder, comm_graph);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Graph_get : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  maxindex;
  int  maxedges;
  int * index;
  int * edges;
  Proxy_MPI_Graph_get(MPI_Comm  comm, int  maxindex, int  maxedges, int * index, int * edges)
  : comm(comm), maxindex(maxindex), maxedges(maxedges), index(index), edges(edges)
  {}
  void callBack() {
    return_value = ::MPI_Graph_get(comm, maxindex, maxedges, index, edges);
  }
  int id_for_debug() { return 143; }
};
int  MPIWrapperFunneled::_MPI_Graph_get(MPI_Comm  comm, int  maxindex, int  maxedges, int * index, int * edges) {
  Proxy_MPI_Graph_get proxyObj(comm, maxindex, maxedges, index, edges);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Graph_map : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  nnodes;
  int * index;
  int * edges;
  int * newrank;
  Proxy_MPI_Graph_map(MPI_Comm  comm, int  nnodes, int * index, int * edges, int * newrank)
  : comm(comm), nnodes(nnodes), index(index), edges(edges), newrank(newrank)
  {}
  void callBack() {
    return_value = ::MPI_Graph_map(comm, nnodes, index, edges, newrank);
  }
  int id_for_debug() { return 144; }
};
int  MPIWrapperFunneled::_MPI_Graph_map(MPI_Comm  comm, int  nnodes, int * index, int * edges, int * newrank) {
  Proxy_MPI_Graph_map proxyObj(comm, nnodes, index, edges, newrank);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Graph_neighbors_count : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  rank;
  int * nneighbors;
  Proxy_MPI_Graph_neighbors_count(MPI_Comm  comm, int  rank, int * nneighbors)
  : comm(comm), rank(rank), nneighbors(nneighbors)
  {}
  void callBack() {
    return_value = ::MPI_Graph_neighbors_count(comm, rank, nneighbors);
  }
  int id_for_debug() { return 145; }
};
int  MPIWrapperFunneled::_MPI_Graph_neighbors_count(MPI_Comm  comm, int  rank, int * nneighbors) {
  Proxy_MPI_Graph_neighbors_count proxyObj(comm, rank, nneighbors);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Graph_neighbors : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int  rank;
  int  maxneighbors;
  int * neighbors;
  Proxy_MPI_Graph_neighbors(MPI_Comm  comm, int  rank, int  maxneighbors, int * neighbors)
  : comm(comm), rank(rank), maxneighbors(maxneighbors), neighbors(neighbors)
  {}
  void callBack() {
    return_value = ::MPI_Graph_neighbors(comm, rank, maxneighbors, neighbors);
  }
  int id_for_debug() { return 146; }
};
int  MPIWrapperFunneled::_MPI_Graph_neighbors(MPI_Comm  comm, int  rank, int  maxneighbors, int * neighbors) {
  Proxy_MPI_Graph_neighbors proxyObj(comm, rank, maxneighbors, neighbors);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Graphdims_get : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int * nnodes;
  int * nedges;
  Proxy_MPI_Graphdims_get(MPI_Comm  comm, int * nnodes, int * nedges)
  : comm(comm), nnodes(nnodes), nedges(nedges)
  {}
  void callBack() {
    return_value = ::MPI_Graphdims_get(comm, nnodes, nedges);
  }
  int id_for_debug() { return 147; }
};
int  MPIWrapperFunneled::_MPI_Graphdims_get(MPI_Comm  comm, int * nnodes, int * nedges) {
  Proxy_MPI_Graphdims_get proxyObj(comm, nnodes, nedges);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Grequest_complete : public Proxy_MPI {
  int  return_value;
  MPI_Request  request;
  Proxy_MPI_Grequest_complete(MPI_Request  request)
  : request(request)
  {}
  void callBack() {
    return_value = ::MPI_Grequest_complete(request);
  }
  int id_for_debug() { return 148; }
};
int  MPIWrapperFunneled::_MPI_Grequest_complete(MPI_Request  request) {
  Proxy_MPI_Grequest_complete proxyObj(request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Grequest_start : public Proxy_MPI {
  int  return_value;
  MPI_Grequest_query_function * query_fn;
  MPI_Grequest_free_function * free_fn;
  MPI_Grequest_cancel_function * cancel_fn;
  void * extra_state;
  MPI_Request * request;
  Proxy_MPI_Grequest_start(MPI_Grequest_query_function * query_fn, MPI_Grequest_free_function * free_fn, MPI_Grequest_cancel_function * cancel_fn, void * extra_state, MPI_Request * request)
  : query_fn(query_fn), free_fn(free_fn), cancel_fn(cancel_fn), extra_state(extra_state), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Grequest_start(query_fn, free_fn, cancel_fn, extra_state, request);
  }
  int id_for_debug() { return 149; }
};
int  MPIWrapperFunneled::_MPI_Grequest_start(MPI_Grequest_query_function * query_fn, MPI_Grequest_free_function * free_fn, MPI_Grequest_cancel_function * cancel_fn, void * extra_state, MPI_Request * request) {
  Proxy_MPI_Grequest_start proxyObj(query_fn, free_fn, cancel_fn, extra_state, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}


#endif



struct Proxy_MPI_Group_compare : public Proxy_MPI {
  int  return_value;
  MPI_Group  group1;
  MPI_Group  group2;
  int * result;
  Proxy_MPI_Group_compare(MPI_Group  group1, MPI_Group  group2, int * result)
  : group1(group1), group2(group2), result(result)
  {}
  void callBack() {
    return_value = ::MPI_Group_compare(group1, group2, result);
  }
  int id_for_debug() { return 151; }
};
int  MPIWrapperFunneled::_MPI_Group_compare(MPI_Group  group1, MPI_Group  group2, int * result) {
  Proxy_MPI_Group_compare proxyObj(group1, group2, result);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Group_difference : public Proxy_MPI {
  int  return_value;
  MPI_Group  group1;
  MPI_Group  group2;
  MPI_Group * newgroup;
  Proxy_MPI_Group_difference(MPI_Group  group1, MPI_Group  group2, MPI_Group * newgroup)
  : group1(group1), group2(group2), newgroup(newgroup)
  {}
  void callBack() {
    return_value = ::MPI_Group_difference(group1, group2, newgroup);
  }
  int id_for_debug() { return 152; }
};
int  MPIWrapperFunneled::_MPI_Group_difference(MPI_Group  group1, MPI_Group  group2, MPI_Group * newgroup) {
  Proxy_MPI_Group_difference proxyObj(group1, group2, newgroup);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Group_excl : public Proxy_MPI {
  int  return_value;
  MPI_Group  group;
  int  n;
  int * ranks;
  MPI_Group * newgroup;
  Proxy_MPI_Group_excl(MPI_Group  group, int  n, int * ranks, MPI_Group * newgroup)
  : group(group), n(n), ranks(ranks), newgroup(newgroup)
  {}
  void callBack() {
    return_value = ::MPI_Group_excl(group, n, ranks, newgroup);
  }
  int id_for_debug() { return 153; }
};
int  MPIWrapperFunneled::_MPI_Group_excl(MPI_Group  group, int  n, int * ranks, MPI_Group * newgroup) {
  Proxy_MPI_Group_excl proxyObj(group, n, ranks, newgroup);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Group_free : public Proxy_MPI {
  int  return_value;
  MPI_Group * group;
  Proxy_MPI_Group_free(MPI_Group * group)
  : group(group)
  {}
  void callBack() {
    return_value = ::MPI_Group_free(group);
  }
  int id_for_debug() { return 155; }
};
int  MPIWrapperFunneled::_MPI_Group_free(MPI_Group * group) {
  Proxy_MPI_Group_free proxyObj(group);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Group_incl : public Proxy_MPI {
  int  return_value;
  MPI_Group  group;
  int  n;
  int * ranks;
  MPI_Group * newgroup;
  Proxy_MPI_Group_incl(MPI_Group  group, int  n, int * ranks, MPI_Group * newgroup)
  : group(group), n(n), ranks(ranks), newgroup(newgroup)
  {}
  void callBack() {
    return_value = ::MPI_Group_incl(group, n, ranks, newgroup);
  }
  int id_for_debug() { return 156; }
};
int  MPIWrapperFunneled::_MPI_Group_incl(MPI_Group  group, int  n, int * ranks, MPI_Group * newgroup) {
  Proxy_MPI_Group_incl proxyObj(group, n, ranks, newgroup);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Group_intersection : public Proxy_MPI {
  int  return_value;
  MPI_Group  group1;
  MPI_Group  group2;
  MPI_Group * newgroup;
  Proxy_MPI_Group_intersection(MPI_Group  group1, MPI_Group  group2, MPI_Group * newgroup)
  : group1(group1), group2(group2), newgroup(newgroup)
  {}
  void callBack() {
    return_value = ::MPI_Group_intersection(group1, group2, newgroup);
  }
  int id_for_debug() { return 157; }
};
int  MPIWrapperFunneled::_MPI_Group_intersection(MPI_Group  group1, MPI_Group  group2, MPI_Group * newgroup) {
  Proxy_MPI_Group_intersection proxyObj(group1, group2, newgroup);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Group_rank : public Proxy_MPI {
  int  return_value;
  MPI_Group  group;
  int * rank;
  Proxy_MPI_Group_rank(MPI_Group  group, int * rank)
  : group(group), rank(rank)
  {}
  void callBack() {
    return_value = ::MPI_Group_rank(group, rank);
  }
  int id_for_debug() { return 158; }
};
int  MPIWrapperFunneled::_MPI_Group_rank(MPI_Group  group, int * rank) {
  Proxy_MPI_Group_rank proxyObj(group, rank);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Group_size : public Proxy_MPI {
  int  return_value;
  MPI_Group  group;
  int * size;
  Proxy_MPI_Group_size(MPI_Group  group, int * size)
  : group(group), size(size)
  {}
  void callBack() {
    return_value = ::MPI_Group_size(group, size);
  }
  int id_for_debug() { return 159; }
};
int  MPIWrapperFunneled::_MPI_Group_size(MPI_Group  group, int * size) {
  Proxy_MPI_Group_size proxyObj(group, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Group_translate_ranks : public Proxy_MPI {
  int  return_value;
  MPI_Group  group1;
  int  n;
  int * ranks1;
  MPI_Group  group2;
  int * ranks2;
  Proxy_MPI_Group_translate_ranks(MPI_Group  group1, int  n, int * ranks1, MPI_Group  group2, int * ranks2)
  : group1(group1), n(n), ranks1(ranks1), group2(group2), ranks2(ranks2)
  {}
  void callBack() {
    return_value = ::MPI_Group_translate_ranks(group1, n, ranks1, group2, ranks2);
  }
  int id_for_debug() { return 160; }
};
int  MPIWrapperFunneled::_MPI_Group_translate_ranks(MPI_Group  group1, int  n, int * ranks1, MPI_Group  group2, int * ranks2) {
  Proxy_MPI_Group_translate_ranks proxyObj(group1, n, ranks1, group2, ranks2);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Group_union : public Proxy_MPI {
  int  return_value;
  MPI_Group  group1;
  MPI_Group  group2;
  MPI_Group * newgroup;
  Proxy_MPI_Group_union(MPI_Group  group1, MPI_Group  group2, MPI_Group * newgroup)
  : group1(group1), group2(group2), newgroup(newgroup)
  {}
  void callBack() {
    return_value = ::MPI_Group_union(group1, group2, newgroup);
  }
  int id_for_debug() { return 161; }
};
int  MPIWrapperFunneled::_MPI_Group_union(MPI_Group  group1, MPI_Group  group2, MPI_Group * newgroup) {
  Proxy_MPI_Group_union proxyObj(group1, group2, newgroup);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Ibsend : public Proxy_MPI {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  dest;
  int  tag;
  MPI_Comm  comm;
  MPI_Request * request;
  Proxy_MPI_Ibsend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request)
  : buf(buf), count(count), datatype(datatype), dest(dest), tag(tag), comm(comm), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Ibsend(buf, count, datatype, dest, tag, comm, request);
  }
  int id_for_debug() { return 162; }
};
int  MPIWrapperFunneled::_MPI_Ibsend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Ibsend proxyObj(buf, count, datatype, dest, tag, comm, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Info_create : public Proxy_MPI {
  int  return_value;
  MPI_Info * info;
  Proxy_MPI_Info_create(MPI_Info * info)
  : info(info)
  {}
  void callBack() {
    return_value = ::MPI_Info_create(info);
  }
  int id_for_debug() { return 164; }
};
int  MPIWrapperFunneled::_MPI_Info_create(MPI_Info * info) {
  Proxy_MPI_Info_create proxyObj(info);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Info_delete : public Proxy_MPI {
  int  return_value;
  MPI_Info  info;
  char * key;
  Proxy_MPI_Info_delete(MPI_Info  info, char * key)
  : info(info), key(key)
  {}
  void callBack() {
    return_value = ::MPI_Info_delete(info, key);
  }
  int id_for_debug() { return 165; }
};
int  MPIWrapperFunneled::_MPI_Info_delete(MPI_Info  info, char * key) {
  Proxy_MPI_Info_delete proxyObj(info, key);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Info_dup : public Proxy_MPI {
  int  return_value;
  MPI_Info  info;
  MPI_Info * newinfo;
  Proxy_MPI_Info_dup(MPI_Info  info, MPI_Info * newinfo)
  : info(info), newinfo(newinfo)
  {}
  void callBack() {
    return_value = ::MPI_Info_dup(info, newinfo);
  }
  int id_for_debug() { return 166; }
};
int  MPIWrapperFunneled::_MPI_Info_dup(MPI_Info  info, MPI_Info * newinfo) {
  Proxy_MPI_Info_dup proxyObj(info, newinfo);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}




struct Proxy_MPI_Info_free : public Proxy_MPI {
  int  return_value;
  MPI_Info * info;
  Proxy_MPI_Info_free(MPI_Info * info)
  : info(info)
  {}
  void callBack() {
    return_value = ::MPI_Info_free(info);
  }
  int id_for_debug() { return 168; }
};
int  MPIWrapperFunneled::_MPI_Info_free(MPI_Info * info) {
  Proxy_MPI_Info_free proxyObj(info);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Info_get : public Proxy_MPI {
  int  return_value;
  MPI_Info  info;
  char * key;
  int  valuelen;
  char * value;
  int * flag;
  Proxy_MPI_Info_get(MPI_Info  info, char * key, int  valuelen, char * value, int * flag)
  : info(info), key(key), valuelen(valuelen), value(value), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Info_get(info, key, valuelen, value, flag);
  }
  int id_for_debug() { return 169; }
};
int  MPIWrapperFunneled::_MPI_Info_get(MPI_Info  info, char * key, int  valuelen, char * value, int * flag) {
  Proxy_MPI_Info_get proxyObj(info, key, valuelen, value, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Info_get_nkeys : public Proxy_MPI {
  int  return_value;
  MPI_Info  info;
  int * nkeys;
  Proxy_MPI_Info_get_nkeys(MPI_Info  info, int * nkeys)
  : info(info), nkeys(nkeys)
  {}
  void callBack() {
    return_value = ::MPI_Info_get_nkeys(info, nkeys);
  }
  int id_for_debug() { return 170; }
};
int  MPIWrapperFunneled::_MPI_Info_get_nkeys(MPI_Info  info, int * nkeys) {
  Proxy_MPI_Info_get_nkeys proxyObj(info, nkeys);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Info_get_nthkey : public Proxy_MPI {
  int  return_value;
  MPI_Info  info;
  int  n;
  char * key;
  Proxy_MPI_Info_get_nthkey(MPI_Info  info, int  n, char * key)
  : info(info), n(n), key(key)
  {}
  void callBack() {
    return_value = ::MPI_Info_get_nthkey(info, n, key);
  }
  int id_for_debug() { return 171; }
};
int  MPIWrapperFunneled::_MPI_Info_get_nthkey(MPI_Info  info, int  n, char * key) {
  Proxy_MPI_Info_get_nthkey proxyObj(info, n, key);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Info_get_valuelen : public Proxy_MPI {
  int  return_value;
  MPI_Info  info;
  char * key;
  int * valuelen;
  int * flag;
  Proxy_MPI_Info_get_valuelen(MPI_Info  info, char * key, int * valuelen, int * flag)
  : info(info), key(key), valuelen(valuelen), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Info_get_valuelen(info, key, valuelen, flag);
  }
  int id_for_debug() { return 172; }
};
int  MPIWrapperFunneled::_MPI_Info_get_valuelen(MPI_Info  info, char * key, int * valuelen, int * flag) {
  Proxy_MPI_Info_get_valuelen proxyObj(info, key, valuelen, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Info_set : public Proxy_MPI {
  int  return_value;
  MPI_Info  info;
  char * key;
  char * value;
  Proxy_MPI_Info_set(MPI_Info  info, char * key, char * value)
  : info(info), key(key), value(value)
  {}
  void callBack() {
    return_value = ::MPI_Info_set(info, key, value);
  }
  int id_for_debug() { return 173; }
};
int  MPIWrapperFunneled::_MPI_Info_set(MPI_Info  info, char * key, char * value) {
  Proxy_MPI_Info_set proxyObj(info, key, value);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}


int  MPIWrapperFunneled::_MPI_Init(int * argc, char *** argv) {
  int required = MPI_THREAD_FUNNELED;
  int provided;
  return this->_MPI_Init_thread(argc, argv, required, &provided);
}



struct Proxy_MPI_Initialized : public Proxy_MPI {
  int  return_value;
  int * flag;
  Proxy_MPI_Initialized(int * flag)
  : flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Initialized(flag);
  }
  int id_for_debug() { return 174; }
};
int  MPIWrapperFunneled::_MPI_Initialized(int * flag) {
  Proxy_MPI_Initialized proxyObj(flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Intercomm_create : public Proxy_MPI {
  int  return_value;
  MPI_Comm  local_comm;
  int  local_leader;
  MPI_Comm  bridge_comm;
  int  remote_leader;
  int  tag;
  MPI_Comm * newintercomm;
  Proxy_MPI_Intercomm_create(MPI_Comm  local_comm, int  local_leader, MPI_Comm  bridge_comm, int  remote_leader, int  tag, MPI_Comm * newintercomm)
  : local_comm(local_comm), local_leader(local_leader), bridge_comm(bridge_comm), remote_leader(remote_leader), tag(tag), newintercomm(newintercomm)
  {}
  void callBack() {
    return_value = ::MPI_Intercomm_create(local_comm, local_leader, bridge_comm, remote_leader, tag, newintercomm);
  }
  int id_for_debug() { return 175; }
};
int  MPIWrapperFunneled::_MPI_Intercomm_create(MPI_Comm  local_comm, int  local_leader, MPI_Comm  bridge_comm, int  remote_leader, int  tag, MPI_Comm * newintercomm) {
  Proxy_MPI_Intercomm_create proxyObj(local_comm, local_leader, bridge_comm, remote_leader, tag, newintercomm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Intercomm_merge : public Proxy_MPI {
  int  return_value;
  MPI_Comm  intercomm;
  int  high;
  MPI_Comm * newintercomm;
  Proxy_MPI_Intercomm_merge(MPI_Comm  intercomm, int  high, MPI_Comm * newintercomm)
  : intercomm(intercomm), high(high), newintercomm(newintercomm)
  {}
  void callBack() {
    return_value = ::MPI_Intercomm_merge(intercomm, high, newintercomm);
  }
  int id_for_debug() { return 176; }
};
int  MPIWrapperFunneled::_MPI_Intercomm_merge(MPI_Comm  intercomm, int  high, MPI_Comm * newintercomm) {
  Proxy_MPI_Intercomm_merge proxyObj(intercomm, high, newintercomm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Iprobe : public Proxy_MPI {
  int  return_value;
  int  source;
  int  tag;
  MPI_Comm  comm;
  int * flag;
  MPI_Status * status;
  Proxy_MPI_Iprobe(int  source, int  tag, MPI_Comm  comm, int * flag, MPI_Status * status)
  : source(source), tag(tag), comm(comm), flag(flag), status(status)
  {}
  void callBack() {
    return_value = ::MPI_Iprobe(source, tag, comm, flag, status);
  }
  int id_for_debug() { return 177; }
};
int  MPIWrapperFunneled::_MPI_Iprobe(int  source, int  tag, MPI_Comm  comm, int * flag, MPI_Status * status) {
  Proxy_MPI_Iprobe proxyObj(source, tag, comm, flag, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Irecv : public Proxy_MPI {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  source;
  int  tag;
  MPI_Comm  comm;
  MPI_Request * request;
  Proxy_MPI_Irecv(void * buf, int  count, MPI_Datatype  datatype, int  source, int  tag, MPI_Comm  comm, MPI_Request * request)
  : buf(buf), count(count), datatype(datatype), source(source), tag(tag), comm(comm), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Irecv(buf, count, datatype, source, tag, comm, request);
  }
  int id_for_debug() { return 178; }
};
int  MPIWrapperFunneled::_MPI_Irecv(void * buf, int  count, MPI_Datatype  datatype, int  source, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Irecv proxyObj(buf, count, datatype, source, tag, comm, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Irsend : public Proxy_MPI {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  dest;
  int  tag;
  MPI_Comm  comm;
  MPI_Request * request;
  Proxy_MPI_Irsend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request)
  : buf(buf), count(count), datatype(datatype), dest(dest), tag(tag), comm(comm), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Irsend(buf, count, datatype, dest, tag, comm, request);
  }
  int id_for_debug() { return 179; }
};
int  MPIWrapperFunneled::_MPI_Irsend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Irsend proxyObj(buf, count, datatype, dest, tag, comm, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

 


#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Isend : public Proxy_MPI_Non_blocking_send_calls {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  dest;
  int  tag;
  MPI_Comm  comm;
  MPI_Request request_internal;
  MPI_Request * request_generalized;
  int cancelled;
  Proxy_MPI_Isend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, 
		  MPI_Request * request_generalized)
    : buf(buf), count(count), datatype(datatype), dest(dest), tag(tag), comm(comm), 
      request_generalized(request_generalized)
  {}
  void createGRequest() {
    ::MPI_Grequest_start(non_blocking_send_query_fn, 
			 non_blocking_send_free_fn, 
			 non_blocking_send_cancel_fn, 
			 this, this->request_generalized);
  }
  void callNonBlockingSend() {
    return_value = ::MPI_Isend(buf, count, datatype, dest, tag, comm, &request_internal);
  }
  bool tryToFinish() {
    int flag;
    ::MPI_Test(&request_internal, &flag, MPI_STATUS_IGNORE);
    if ( flag ) {
      ::MPI_Grequest_complete(*request_generalized);
      return true;
    }
    return false;
  }
};
int  MPIWrapperFunneled::_MPI_Isend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Isend* proxy = new Proxy_MPI_Isend(buf, count, datatype, dest, tag, comm, request);
  LockMutex(mutex_id_non_blocking_send_list);
  if (count > 1000)
    pendingNonBlockingSendCalls[0].push_back(proxy);
  else
    pendingNonBlockingSendCalls[1].push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_send_list);
  wake_up_executor();
  proxy->LockMutex();
  return MPI_SUCCESS; 
}


int non_blocking_send_query_fn(void *extra_state, MPI_Status *status) {
  Proxy_MPI_Isend* proxy = (Proxy_MPI_Isend*)extra_state;
  MPI_Status_set_elements(status, proxy->datatype, proxy->count);
  MPI_Status_set_cancelled(status, proxy->cancelled); 
  status->MPI_SOURCE = proxy->dest; 
  status->MPI_TAG    = proxy->tag; 
  // status->MPI_ERROR  = ;  Since noone knows what this is we do not bother either 
  return MPI_SUCCESS; 
}
int non_blocking_send_free_fn(void *extra_state) {
  Proxy_MPI_Isend* proxy = (Proxy_MPI_Isend*)extra_state;
  delete proxy;
  return MPI_SUCCESS;   
}
int non_blocking_send_cancel_fn(void *extra_state, int complete) {
  Proxy_MPI_Isend* proxy = (Proxy_MPI_Isend*)extra_state;
  if (!complete) { 
    fprintf(stderr, "Cannot cancel generalized request - aborting program\n"); 
    MPI_Abort(MPI_COMM_WORLD, 99); 
  } 
  proxy->cancelled = 1;
  return MPI_SUCCESS;     
}



#else


struct Proxy_MPI_Isend : public Proxy_MPI {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  dest;
  int  tag;
  MPI_Comm  comm;
  MPI_Request * request;
  Proxy_MPI_Isend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request)
  : buf(buf), count(count), datatype(datatype), dest(dest), tag(tag), comm(comm), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Isend(buf, count, datatype, dest, tag, comm, request);
  }
  int id_for_debug() { return 180; }
};
int  MPIWrapperFunneled::_MPI_Isend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Isend proxyObj(buf, count, datatype, dest, tag, comm, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}


#endif








struct Proxy_MPI_Issend : public Proxy_MPI {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  dest;
  int  tag;
  MPI_Comm  comm;
  MPI_Request * request;
  Proxy_MPI_Issend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request)
  : buf(buf), count(count), datatype(datatype), dest(dest), tag(tag), comm(comm), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Issend(buf, count, datatype, dest, tag, comm, request);
  }
  int id_for_debug() { return 181; }
};
int  MPIWrapperFunneled::_MPI_Issend(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Issend proxyObj(buf, count, datatype, dest, tag, comm, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Is_thread_main : public Proxy_MPI {
  int  return_value;
  int * flag;
  Proxy_MPI_Is_thread_main(int * flag)
  : flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Is_thread_main(flag);
  }
  int id_for_debug() { return 182; }
};
int  MPIWrapperFunneled::_MPI_Is_thread_main(int * flag) {
  Proxy_MPI_Is_thread_main proxyObj(flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Lookup_name : public Proxy_MPI {
  int  return_value;
  char * service_name;
  MPI_Info  info;
  char * port_name;
  Proxy_MPI_Lookup_name(char * service_name, MPI_Info  info, char * port_name)
  : service_name(service_name), info(info), port_name(port_name)
  {}
  void callBack() {
    return_value = ::MPI_Lookup_name(service_name, info, port_name);
  }
  int id_for_debug() { return 185; }
};
int  MPIWrapperFunneled::_MPI_Lookup_name(char * service_name, MPI_Info  info, char * port_name) {
  Proxy_MPI_Lookup_name proxyObj(service_name, info, port_name);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}


#endif



struct Proxy_MPI_Op_create : public Proxy_MPI {
  int  return_value;
  MPI_User_function * function;
  int  commute;
  MPI_Op * op;
  Proxy_MPI_Op_create(MPI_User_function * function, int  commute, MPI_Op * op)
  : function(function), commute(commute), op(op)
  {}
  void callBack() {
    return_value = ::MPI_Op_create(function, commute, op);
  }
  int id_for_debug() { return 187; }
};
int  MPIWrapperFunneled::_MPI_Op_create(MPI_User_function * function, int  commute, MPI_Op * op) {
  Proxy_MPI_Op_create proxyObj(function, commute, op);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Open_port : public Proxy_MPI {
  int  return_value;
  MPI_Info  info;
  char * port_name;
  Proxy_MPI_Open_port(MPI_Info  info, char * port_name)
  : info(info), port_name(port_name)
  {}
  void callBack() {
    return_value = ::MPI_Open_port(info, port_name);
  }
  int id_for_debug() { return 188; }
};
int  MPIWrapperFunneled::_MPI_Open_port(MPI_Info  info, char * port_name) {
  Proxy_MPI_Open_port proxyObj(info, port_name);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}


#endif



struct Proxy_MPI_Op_free : public Proxy_MPI {
  int  return_value;
  MPI_Op * op;
  Proxy_MPI_Op_free(MPI_Op * op)
  : op(op)
  {}
  void callBack() {
    return_value = ::MPI_Op_free(op);
  }
  int id_for_debug() { return 190; }
};
int  MPIWrapperFunneled::_MPI_Op_free(MPI_Op * op) {
  Proxy_MPI_Op_free proxyObj(op);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Pack_external : public Proxy_MPI {
  int  return_value;
  char * datarep;
  void * inbuf;
  int  incount;
  MPI_Datatype  datatype;
  void * outbuf;
  MPI_Aint  outsize;
  MPI_Aint * position;
  Proxy_MPI_Pack_external(char * datarep, void * inbuf, int  incount, MPI_Datatype  datatype, void * outbuf, MPI_Aint  outsize, MPI_Aint * position)
  : datarep(datarep), inbuf(inbuf), incount(incount), datatype(datatype), outbuf(outbuf), outsize(outsize), position(position)
  {}
  void callBack() {
    return_value = ::MPI_Pack_external(datarep, inbuf, incount, datatype, outbuf, outsize, position);
  }
  int id_for_debug() { return 191; }
};
int  MPIWrapperFunneled::_MPI_Pack_external(char * datarep, void * inbuf, int  incount, MPI_Datatype  datatype, void * outbuf, MPI_Aint  outsize, MPI_Aint * position) {
  Proxy_MPI_Pack_external proxyObj(datarep, inbuf, incount, datatype, outbuf, outsize, position);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Pack_external_size : public Proxy_MPI {
  int  return_value;
  char * datarep;
  int  incount;
  MPI_Datatype  datatype;
  MPI_Aint * size;
  Proxy_MPI_Pack_external_size(char * datarep, int  incount, MPI_Datatype  datatype, MPI_Aint * size)
  : datarep(datarep), incount(incount), datatype(datatype), size(size)
  {}
  void callBack() {
    return_value = ::MPI_Pack_external_size(datarep, incount, datatype, size);
  }
  int id_for_debug() { return 192; }
};
int  MPIWrapperFunneled::_MPI_Pack_external_size(char * datarep, int  incount, MPI_Datatype  datatype, MPI_Aint * size) {
  Proxy_MPI_Pack_external_size proxyObj(datarep, incount, datatype, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Pack : public Proxy_MPI {
  int  return_value;
  void * inbuf;
  int  incount;
  MPI_Datatype  datatype;
  void * outbuf;
  int  outsize;
  int * position;
  MPI_Comm  comm;
  Proxy_MPI_Pack(void * inbuf, int  incount, MPI_Datatype  datatype, void * outbuf, int  outsize, int * position, MPI_Comm  comm)
  : inbuf(inbuf), incount(incount), datatype(datatype), outbuf(outbuf), outsize(outsize), position(position), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Pack(inbuf, incount, datatype, outbuf, outsize, position, comm);
  }
  int id_for_debug() { return 193; }
};
int  MPIWrapperFunneled::_MPI_Pack(void * inbuf, int  incount, MPI_Datatype  datatype, void * outbuf, int  outsize, int * position, MPI_Comm  comm) {
  Proxy_MPI_Pack proxyObj(inbuf, incount, datatype, outbuf, outsize, position, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Pack_size : public Proxy_MPI {
  int  return_value;
  int  incount;
  MPI_Datatype  datatype;
  MPI_Comm  comm;
  int * size;
  Proxy_MPI_Pack_size(int  incount, MPI_Datatype  datatype, MPI_Comm  comm, int * size)
  : incount(incount), datatype(datatype), comm(comm), size(size)
  {}
  void callBack() {
    return_value = ::MPI_Pack_size(incount, datatype, comm, size);
  }
  int id_for_debug() { return 194; }
};
int  MPIWrapperFunneled::_MPI_Pack_size(int  incount, MPI_Datatype  datatype, MPI_Comm  comm, int * size) {
  Proxy_MPI_Pack_size proxyObj(incount, datatype, comm, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Publish_name : public Proxy_MPI {
  int  return_value;
  char * service_name;
  MPI_Info  info;
  char * port_name;
  Proxy_MPI_Publish_name(char * service_name, MPI_Info  info, char * port_name)
  : service_name(service_name), info(info), port_name(port_name)
  {}
  void callBack() {
    return_value = ::MPI_Publish_name(service_name, info, port_name);
  }
  int id_for_debug() { return 195; }
};
int  MPIWrapperFunneled::_MPI_Publish_name(char * service_name, MPI_Info  info, char * port_name) {
  Proxy_MPI_Publish_name proxyObj(service_name, info, port_name);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Put : public Proxy_MPI {
  int  return_value;
  void * origin_addr;
  int  origin_count;
  MPI_Datatype  origin_datatype;
  int  target_rank;
  MPI_Aint  target_disp;
  int  target_count;
  MPI_Datatype  target_datatype;
  MPI_Win  win;
  Proxy_MPI_Put(void * origin_addr, int  origin_count, MPI_Datatype  origin_datatype, int  target_rank, MPI_Aint  target_disp, int  target_count, MPI_Datatype  target_datatype, MPI_Win  win)
  : origin_addr(origin_addr), origin_count(origin_count), origin_datatype(origin_datatype), target_rank(target_rank), target_disp(target_disp), target_count(target_count), target_datatype(target_datatype), win(win)
  {}
  void callBack() {
    return_value = ::MPI_Put(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, win);
  }
  int id_for_debug() { return 196; }
};
int  MPIWrapperFunneled::_MPI_Put(void * origin_addr, int  origin_count, MPI_Datatype  origin_datatype, int  target_rank, MPI_Aint  target_disp, int  target_count, MPI_Datatype  target_datatype, MPI_Win  win) {
  Proxy_MPI_Put proxyObj(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Query_thread : public Proxy_MPI {
  int  return_value;
  int * provided;
  Proxy_MPI_Query_thread(int * provided)
  : provided(provided)
  {}
  void callBack() {
    return_value = ::MPI_Query_thread(provided);
  }
  int id_for_debug() { return 197; }
};
int  MPIWrapperFunneled::_MPI_Query_thread(int * provided) {
  Proxy_MPI_Query_thread proxyObj(provided);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Recv_init : public Proxy_MPI {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  source;
  int  tag;
  MPI_Comm  comm;
  MPI_Request * request;
  Proxy_MPI_Recv_init(void * buf, int  count, MPI_Datatype  datatype, int  source, int  tag, MPI_Comm  comm, MPI_Request * request)
  : buf(buf), count(count), datatype(datatype), source(source), tag(tag), comm(comm), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Recv_init(buf, count, datatype, source, tag, comm, request);
  }
  int id_for_debug() { return 198; }
};
int  MPIWrapperFunneled::_MPI_Recv_init(void * buf, int  count, MPI_Datatype  datatype, int  source, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Recv_init proxyObj(buf, count, datatype, source, tag, comm, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Reduce : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  void * recvbuf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Op  op;
  int  root;
  MPI_Comm  comm;
  Proxy_MPI_Reduce(void * sendbuf, void * recvbuf, int  count, MPI_Datatype  datatype, MPI_Op  op, int  root, MPI_Comm  comm)
  : sendbuf(sendbuf), recvbuf(recvbuf), count(count), datatype(datatype), op(op), root(root), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  }
  int id_for_debug() { return 199; }
};
int  MPIWrapperFunneled::_MPI_Reduce(void * sendbuf, void * recvbuf, int  count, MPI_Datatype  datatype, MPI_Op  op, int  root, MPI_Comm  comm) {
  Proxy_MPI_Reduce proxyObj(sendbuf, recvbuf, count, datatype, op, root, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Reduce_scatter : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  void * recvbuf;
  int * recvcounts;
  MPI_Datatype  datatype;
  MPI_Op  op;
  MPI_Comm  comm;
  Proxy_MPI_Reduce_scatter(void * sendbuf, void * recvbuf, int * recvcounts, MPI_Datatype  datatype, MPI_Op  op, MPI_Comm  comm)
  : sendbuf(sendbuf), recvbuf(recvbuf), recvcounts(recvcounts), datatype(datatype), op(op), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
  }
  int id_for_debug() { return 200; }
};
int  MPIWrapperFunneled::_MPI_Reduce_scatter(void * sendbuf, void * recvbuf, int * recvcounts, MPI_Datatype  datatype, MPI_Op  op, MPI_Comm  comm) {
  Proxy_MPI_Reduce_scatter proxyObj(sendbuf, recvbuf, recvcounts, datatype, op, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Register_datarep : public Proxy_MPI {
  int  return_value;
  char * datarep;
  MPI_Datarep_conversion_function * read_conversion_fn;
  MPI_Datarep_conversion_function * write_conversion_fn;
  MPI_Datarep_extent_function * dtype_file_extent_fn;
  void * extra_state;
  Proxy_MPI_Register_datarep(char * datarep, MPI_Datarep_conversion_function * read_conversion_fn, MPI_Datarep_conversion_function * write_conversion_fn, MPI_Datarep_extent_function * dtype_file_extent_fn, void * extra_state)
  : datarep(datarep), read_conversion_fn(read_conversion_fn), write_conversion_fn(write_conversion_fn), dtype_file_extent_fn(dtype_file_extent_fn), extra_state(extra_state)
  {}
  void callBack() {
    return_value = ::MPI_Register_datarep(datarep, read_conversion_fn, write_conversion_fn, dtype_file_extent_fn, extra_state);
  }
  int id_for_debug() { return 201; }
};
int  MPIWrapperFunneled::_MPI_Register_datarep(char * datarep, MPI_Datarep_conversion_function * read_conversion_fn, MPI_Datarep_conversion_function * write_conversion_fn, MPI_Datarep_extent_function * dtype_file_extent_fn, void * extra_state) {
  Proxy_MPI_Register_datarep proxyObj(datarep, read_conversion_fn, write_conversion_fn, dtype_file_extent_fn, extra_state);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#endif



struct Proxy_MPI_Request_free : public Proxy_MPI {
  int  return_value;
  MPI_Request * request;
  Proxy_MPI_Request_free(MPI_Request * request)
  : request(request)
  {}
  void callBack() {
    return_value = ::MPI_Request_free(request);
  }
  int id_for_debug() { return 204; }
};
int  MPIWrapperFunneled::_MPI_Request_free(MPI_Request * request) {
  Proxy_MPI_Request_free proxyObj(request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Request_get_status : public Proxy_MPI {
  int  return_value;
  MPI_Request  request;
  int * flag;
  MPI_Status * status;
  Proxy_MPI_Request_get_status(MPI_Request  request, int * flag, MPI_Status * status)
  : request(request), flag(flag), status(status)
  {}
  void callBack() {
    return_value = ::MPI_Request_get_status(request, flag, status);
  }
  int id_for_debug() { return 205; }
};
int  MPIWrapperFunneled::_MPI_Request_get_status(MPI_Request  request, int * flag, MPI_Status * status) {
  Proxy_MPI_Request_get_status proxyObj(request, flag, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Rsend_init : public Proxy_MPI {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  dest;
  int  tag;
  MPI_Comm  comm;
  MPI_Request * request;
  Proxy_MPI_Rsend_init(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request)
  : buf(buf), count(count), datatype(datatype), dest(dest), tag(tag), comm(comm), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Rsend_init(buf, count, datatype, dest, tag, comm, request);
  }
  int id_for_debug() { return 206; }
};
int  MPIWrapperFunneled::_MPI_Rsend_init(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Rsend_init proxyObj(buf, count, datatype, dest, tag, comm, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Scan : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  void * recvbuf;
  int  count;
  MPI_Datatype  datatype;
  MPI_Op  op;
  MPI_Comm  comm;
  Proxy_MPI_Scan(void * sendbuf, void * recvbuf, int  count, MPI_Datatype  datatype, MPI_Op  op, MPI_Comm  comm)
  : sendbuf(sendbuf), recvbuf(recvbuf), count(count), datatype(datatype), op(op), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
  }
  int id_for_debug() { return 207; }
};
int  MPIWrapperFunneled::_MPI_Scan(void * sendbuf, void * recvbuf, int  count, MPI_Datatype  datatype, MPI_Op  op, MPI_Comm  comm) {
  Proxy_MPI_Scan proxyObj(sendbuf, recvbuf, count, datatype, op, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Scatter : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  int  sendcount;
  MPI_Datatype  sendtype;
  void * recvbuf;
  int  recvcount;
  MPI_Datatype  recvtype;
  int  root;
  MPI_Comm  comm;
  Proxy_MPI_Scatter(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, int  root, MPI_Comm  comm)
  : sendbuf(sendbuf), sendcount(sendcount), sendtype(sendtype), recvbuf(recvbuf), recvcount(recvcount), recvtype(recvtype), root(root), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
  int id_for_debug() { return 208; }
};
int  MPIWrapperFunneled::_MPI_Scatter(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, int  root, MPI_Comm  comm) {
  Proxy_MPI_Scatter proxyObj(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Scatterv : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  int * sendcounts;
  int * displs;
  MPI_Datatype  sendtype;
  void * recvbuf;
  int  recvcount;
  MPI_Datatype  recvtype;
  int  root;
  MPI_Comm  comm;
  Proxy_MPI_Scatterv(void * sendbuf, int * sendcounts, int * displs, MPI_Datatype  sendtype, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, int  root, MPI_Comm  comm)
  : sendbuf(sendbuf), sendcounts(sendcounts), displs(displs), sendtype(sendtype), recvbuf(recvbuf), recvcount(recvcount), recvtype(recvtype), root(root), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
  int id_for_debug() { return 209; }
};
int  MPIWrapperFunneled::_MPI_Scatterv(void * sendbuf, int * sendcounts, int * displs, MPI_Datatype  sendtype, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, int  root, MPI_Comm  comm) {
  Proxy_MPI_Scatterv proxyObj(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Send_init : public Proxy_MPI {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  dest;
  int  tag;
  MPI_Comm  comm;
  MPI_Request * request;
  Proxy_MPI_Send_init(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request)
  : buf(buf), count(count), datatype(datatype), dest(dest), tag(tag), comm(comm), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Send_init(buf, count, datatype, dest, tag, comm, request);
  }
  int id_for_debug() { return 210; }
};
int  MPIWrapperFunneled::_MPI_Send_init(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Send_init proxyObj(buf, count, datatype, dest, tag, comm, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Sendrecv : public Proxy_MPI {
  int  return_value;
  void * sendbuf;
  int  sendcount;
  MPI_Datatype  sendtype;
  int  dest;
  int  sendtag;
  void * recvbuf;
  int  recvcount;
  MPI_Datatype  recvtype;
  int  source;
  int  recvtag;
  MPI_Comm  comm;
  MPI_Status * status;
  Proxy_MPI_Sendrecv(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, int  dest, int  sendtag, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, int  source, int  recvtag, MPI_Comm  comm, MPI_Status * status)
  : sendbuf(sendbuf), sendcount(sendcount), sendtype(sendtype), dest(dest), sendtag(sendtag), recvbuf(recvbuf), recvcount(recvcount), recvtype(recvtype), source(source), recvtag(recvtag), comm(comm), status(status)
  {}
  void callBack() {
    return_value = ::MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
  }
  int id_for_debug() { return 211; }
};
int  MPIWrapperFunneled::_MPI_Sendrecv(void * sendbuf, int  sendcount, MPI_Datatype  sendtype, int  dest, int  sendtag, void * recvbuf, int  recvcount, MPI_Datatype  recvtype, int  source, int  recvtag, MPI_Comm  comm, MPI_Status * status) {
  Proxy_MPI_Sendrecv proxyObj(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Sendrecv_replace : public Proxy_MPI {
  int  return_value;
  void *  buf;
  int  count;
  MPI_Datatype  datatype;
  int  dest;
  int  sendtag;
  int  source;
  int  recvtag;
  MPI_Comm  comm;
  MPI_Status * status;
  Proxy_MPI_Sendrecv_replace(void *  buf, int  count, MPI_Datatype  datatype, int  dest, int  sendtag, int  source, int  recvtag, MPI_Comm  comm, MPI_Status * status)
  : buf(buf), count(count), datatype(datatype), dest(dest), sendtag(sendtag), source(source), recvtag(recvtag), comm(comm), status(status)
  {}
  void callBack() {
    return_value = ::MPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
  }
  int id_for_debug() { return 212; }
};
int  MPIWrapperFunneled::_MPI_Sendrecv_replace(void *  buf, int  count, MPI_Datatype  datatype, int  dest, int  sendtag, int  source, int  recvtag, MPI_Comm  comm, MPI_Status * status) {
  Proxy_MPI_Sendrecv_replace proxyObj(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Ssend_init : public Proxy_MPI {
  int  return_value;
  void * buf;
  int  count;
  MPI_Datatype  datatype;
  int  dest;
  int  tag;
  MPI_Comm  comm;
  MPI_Request * request;
  Proxy_MPI_Ssend_init(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request)
  : buf(buf), count(count), datatype(datatype), dest(dest), tag(tag), comm(comm), request(request)
  {}
  void callBack() {
    return_value = ::MPI_Ssend_init(buf, count, datatype, dest, tag, comm, request);
  }
  int id_for_debug() { return 213; }
};
int  MPIWrapperFunneled::_MPI_Ssend_init(void * buf, int  count, MPI_Datatype  datatype, int  dest, int  tag, MPI_Comm  comm, MPI_Request * request) {
  Proxy_MPI_Ssend_init proxyObj(buf, count, datatype, dest, tag, comm, request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Start : public Proxy_MPI {
  int  return_value;
  MPI_Request * request;
  Proxy_MPI_Start(MPI_Request * request)
  : request(request)
  {}
  void callBack() {
    return_value = ::MPI_Start(request);
  }
  int id_for_debug() { return 214; }
};
int  MPIWrapperFunneled::_MPI_Start(MPI_Request * request) {
  Proxy_MPI_Start proxyObj(request);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Startall : public Proxy_MPI {
  int  return_value;
  int  count;
  MPI_Request * array_of_requests;
  Proxy_MPI_Startall(int  count, MPI_Request * array_of_requests)
  : count(count), array_of_requests(array_of_requests)
  {}
  void callBack() {
    return_value = ::MPI_Startall(count, array_of_requests);
  }
  int id_for_debug() { return 215; }
};
int  MPIWrapperFunneled::_MPI_Startall(int  count, MPI_Request * array_of_requests) {
  Proxy_MPI_Startall proxyObj(count, array_of_requests);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Status_c2f : public Proxy_MPI {
  int  return_value;
  MPI_Status * c_status;
  MPI_Fint * f_status;
  Proxy_MPI_Status_c2f(MPI_Status * c_status, MPI_Fint * f_status)
  : c_status(c_status), f_status(f_status)
  {}
  void callBack() {
    return_value = ::MPI_Status_c2f(c_status, f_status);
  }
  int id_for_debug() { return 216; }
};
int  MPIWrapperFunneled::_MPI_Status_c2f(MPI_Status * c_status, MPI_Fint * f_status) {
  Proxy_MPI_Status_c2f proxyObj(c_status, f_status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Status_f2c : public Proxy_MPI {
  int  return_value;
  MPI_Fint * f_status;
  MPI_Status * c_status;
  Proxy_MPI_Status_f2c(MPI_Fint * f_status, MPI_Status * c_status)
  : f_status(f_status), c_status(c_status)
  {}
  void callBack() {
    return_value = ::MPI_Status_f2c(f_status, c_status);
  }
  int id_for_debug() { return 217; }
};
int  MPIWrapperFunneled::_MPI_Status_f2c(MPI_Fint * f_status, MPI_Status * c_status) {
  Proxy_MPI_Status_f2c proxyObj(f_status, c_status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Status_set_cancelled : public Proxy_MPI {
  int  return_value;
  MPI_Status * status;
  int  flag;
  Proxy_MPI_Status_set_cancelled(MPI_Status * status, int  flag)
  : status(status), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Status_set_cancelled(status, flag);
  }
  int id_for_debug() { return 218; }
};
int  MPIWrapperFunneled::_MPI_Status_set_cancelled(MPI_Status * status, int  flag) {
  Proxy_MPI_Status_set_cancelled proxyObj(status, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Status_set_elements : public Proxy_MPI {
  int  return_value;
  MPI_Status * status;
  MPI_Datatype  datatype;
  int  count;
  Proxy_MPI_Status_set_elements(MPI_Status * status, MPI_Datatype  datatype, int  count)
  : status(status), datatype(datatype), count(count)
  {}
  void callBack() {
    return_value = ::MPI_Status_set_elements(status, datatype, count);
  }
  int id_for_debug() { return 219; }
};
int  MPIWrapperFunneled::_MPI_Status_set_elements(MPI_Status * status, MPI_Datatype  datatype, int  count) {
  Proxy_MPI_Status_set_elements proxyObj(status, datatype, count);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Testall : public Proxy_MPI {
  int  return_value;
  int  count;
  MPI_Request * array_of_requests;
  int * flag;
  MPI_Status * array_of_statuses;
  Proxy_MPI_Testall(int  count, MPI_Request * array_of_requests, int * flag, MPI_Status * array_of_statuses)
  : count(count), array_of_requests(array_of_requests), flag(flag), array_of_statuses(array_of_statuses)
  {}
  void callBack() {
    return_value = ::MPI_Testall(count, array_of_requests, flag, array_of_statuses);
  }
  int id_for_debug() { return 220; }
};
int  MPIWrapperFunneled::_MPI_Testall(int  count, MPI_Request * array_of_requests, int * flag, MPI_Status * array_of_statuses) {
  Proxy_MPI_Testall proxyObj(count, array_of_requests, flag, array_of_statuses);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Testany : public Proxy_MPI {
  int  return_value;
  int  count;
  MPI_Request * array_of_requests;
  int * index;
  int * flag;
  MPI_Status * status;
  Proxy_MPI_Testany(int  count, MPI_Request * array_of_requests, int * index, int * flag, MPI_Status * status)
  : count(count), array_of_requests(array_of_requests), index(index), flag(flag), status(status)
  {}
  void callBack() {
    return_value = ::MPI_Testany(count, array_of_requests, index, flag, status);
  }
  int id_for_debug() { return 221; }
};
int  MPIWrapperFunneled::_MPI_Testany(int  count, MPI_Request * array_of_requests, int * index, int * flag, MPI_Status * status) {
  Proxy_MPI_Testany proxyObj(count, array_of_requests, index, flag, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Test : public Proxy_MPI {
  int  return_value;
  MPI_Request * request;
  int * flag;
  MPI_Status * status;
  Proxy_MPI_Test(MPI_Request * request, int * flag, MPI_Status * status)
  : request(request), flag(flag), status(status)
  {}
  void callBack() {
    return_value = ::MPI_Test(request, flag, status);
  }
  int id_for_debug() { return 222; }
};
int  MPIWrapperFunneled::_MPI_Test(MPI_Request * request, int * flag, MPI_Status * status) {
  Proxy_MPI_Test proxyObj(request, flag, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Test_cancelled : public Proxy_MPI {
  int  return_value;
  MPI_Status * status;
  int * flag;
  Proxy_MPI_Test_cancelled(MPI_Status * status, int * flag)
  : status(status), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Test_cancelled(status, flag);
  }
  int id_for_debug() { return 223; }
};
int  MPIWrapperFunneled::_MPI_Test_cancelled(MPI_Status * status, int * flag) {
  Proxy_MPI_Test_cancelled proxyObj(status, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Testsome : public Proxy_MPI {
  int  return_value;
  int  incount;
  MPI_Request * array_of_requests;
  int * outcount;
  int * array_of_indices;
  MPI_Status * array_of_statuses;
  Proxy_MPI_Testsome(int  incount, MPI_Request * array_of_requests, int * outcount, int * array_of_indices, MPI_Status * array_of_statuses)
  : incount(incount), array_of_requests(array_of_requests), outcount(outcount), array_of_indices(array_of_indices), array_of_statuses(array_of_statuses)
  {}
  void callBack() {
    return_value = ::MPI_Testsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  }
  int id_for_debug() { return 224; }
};
int  MPIWrapperFunneled::_MPI_Testsome(int  incount, MPI_Request * array_of_requests, int * outcount, int * array_of_indices, MPI_Status * array_of_statuses) {
  Proxy_MPI_Testsome proxyObj(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Topo_test : public Proxy_MPI {
  int  return_value;
  MPI_Comm  comm;
  int * status;
  Proxy_MPI_Topo_test(MPI_Comm  comm, int * status)
  : comm(comm), status(status)
  {}
  void callBack() {
    return_value = ::MPI_Topo_test(comm, status);
  }
  int id_for_debug() { return 225; }
};
int  MPIWrapperFunneled::_MPI_Topo_test(MPI_Comm  comm, int * status) {
  Proxy_MPI_Topo_test proxyObj(comm, status);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}




struct Proxy_MPI_Type_commit : public Proxy_MPI {
  int  return_value;
  MPI_Datatype * type;
  Proxy_MPI_Type_commit(MPI_Datatype * type)
  : type(type)
  {}
  void callBack() {
    return_value = ::MPI_Type_commit(type);
  }
  int id_for_debug() { return 227; }
};
int  MPIWrapperFunneled::_MPI_Type_commit(MPI_Datatype * type) {
  Proxy_MPI_Type_commit proxyObj(type);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_contiguous : public Proxy_MPI {
  int  return_value;
  int  count;
  MPI_Datatype  oldtype;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_contiguous(int  count, MPI_Datatype  oldtype, MPI_Datatype * newtype)
  : count(count), oldtype(oldtype), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_contiguous(count, oldtype, newtype);
  }
  int id_for_debug() { return 228; }
};
int  MPIWrapperFunneled::_MPI_Type_contiguous(int  count, MPI_Datatype  oldtype, MPI_Datatype * newtype) {
  Proxy_MPI_Type_contiguous proxyObj(count, oldtype, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_create_darray : public Proxy_MPI {
  int  return_value;
  int  size;
  int  rank;
  int  ndims;
  int * gsize_array;
  int * distrib_array;
  int * darg_array;
  int * psize_array;
  int  order;
  MPI_Datatype  oldtype;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_create_darray(int  size, int  rank, int  ndims, int * gsize_array, int * distrib_array, int * darg_array, int * psize_array, int  order, MPI_Datatype  oldtype, MPI_Datatype * newtype)
  : size(size), rank(rank), ndims(ndims), gsize_array(gsize_array), distrib_array(distrib_array), darg_array(darg_array), psize_array(psize_array), order(order), oldtype(oldtype), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_darray(size, rank, ndims, gsize_array, distrib_array, darg_array, psize_array, order, oldtype, newtype);
  }
  int id_for_debug() { return 229; }
};
int  MPIWrapperFunneled::_MPI_Type_create_darray(int  size, int  rank, int  ndims, int * gsize_array, int * distrib_array, int * darg_array, int * psize_array, int  order, MPI_Datatype  oldtype, MPI_Datatype * newtype) {
  Proxy_MPI_Type_create_darray proxyObj(size, rank, ndims, gsize_array, distrib_array, darg_array, psize_array, order, oldtype, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Type_create_f90_complex : public Proxy_MPI {
  int  return_value;
  int  p;
  int  r;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_create_f90_complex(int  p, int  r, MPI_Datatype * newtype)
  : p(p), r(r), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_f90_complex(p, r, newtype);
  }
  int id_for_debug() { return 230; }
};
int  MPIWrapperFunneled::_MPI_Type_create_f90_complex(int  p, int  r, MPI_Datatype * newtype) {
  Proxy_MPI_Type_create_f90_complex proxyObj(p, r, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_create_f90_integer : public Proxy_MPI {
  int  return_value;
  int  r;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_create_f90_integer(int  r, MPI_Datatype * newtype)
  : r(r), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_f90_integer(r, newtype);
  }
  int id_for_debug() { return 231; }
};
int  MPIWrapperFunneled::_MPI_Type_create_f90_integer(int  r, MPI_Datatype * newtype) {
  Proxy_MPI_Type_create_f90_integer proxyObj(r, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_create_f90_real : public Proxy_MPI {
  int  return_value;
  int  p;
  int  r;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_create_f90_real(int  p, int  r, MPI_Datatype * newtype)
  : p(p), r(r), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_f90_real(p, r, newtype);
  }
  int id_for_debug() { return 232; }
};
int  MPIWrapperFunneled::_MPI_Type_create_f90_real(int  p, int  r, MPI_Datatype * newtype) {
  Proxy_MPI_Type_create_f90_real proxyObj(p, r, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Type_create_hindexed : public Proxy_MPI {
  int  return_value;
  int  count;
  int * array_of_blocklengths;
  MPI_Aint * array_of_displacements;
  MPI_Datatype  oldtype;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_create_hindexed(int  count, int * array_of_blocklengths, MPI_Aint * array_of_displacements, MPI_Datatype  oldtype, MPI_Datatype * newtype)
  : count(count), array_of_blocklengths(array_of_blocklengths), array_of_displacements(array_of_displacements), oldtype(oldtype), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
  }
  int id_for_debug() { return 233; }
};
int  MPIWrapperFunneled::_MPI_Type_create_hindexed(int  count, int * array_of_blocklengths, MPI_Aint * array_of_displacements, MPI_Datatype  oldtype, MPI_Datatype * newtype) {
  Proxy_MPI_Type_create_hindexed proxyObj(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_create_hvector : public Proxy_MPI {
  int  return_value;
  int  count;
  int  blocklength;
  MPI_Aint  stride;
  MPI_Datatype  oldtype;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_create_hvector(int  count, int  blocklength, MPI_Aint  stride, MPI_Datatype  oldtype, MPI_Datatype * newtype)
  : count(count), blocklength(blocklength), stride(stride), oldtype(oldtype), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_hvector(count, blocklength, stride, oldtype, newtype);
  }
  int id_for_debug() { return 234; }
};
int  MPIWrapperFunneled::_MPI_Type_create_hvector(int  count, int  blocklength, MPI_Aint  stride, MPI_Datatype  oldtype, MPI_Datatype * newtype) {
  Proxy_MPI_Type_create_hvector proxyObj(count, blocklength, stride, oldtype, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}




#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Type_create_keyval : public Proxy_MPI {
  int  return_value;
  MPI_Type_copy_attr_function * type_copy_attr_fn;
  MPI_Type_delete_attr_function * type_delete_attr_fn;
  int * type_keyval;
  void * extra_state;
  Proxy_MPI_Type_create_keyval(MPI_Type_copy_attr_function * type_copy_attr_fn, MPI_Type_delete_attr_function * type_delete_attr_fn, int * type_keyval, void * extra_state)
  : type_copy_attr_fn(type_copy_attr_fn), type_delete_attr_fn(type_delete_attr_fn), type_keyval(type_keyval), extra_state(extra_state)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_keyval(type_copy_attr_fn, type_delete_attr_fn, type_keyval, extra_state);
  }
  int id_for_debug() { return 235; }
};
int  MPIWrapperFunneled::_MPI_Type_create_keyval(MPI_Type_copy_attr_function * type_copy_attr_fn, MPI_Type_delete_attr_function * type_delete_attr_fn, int * type_keyval, void * extra_state) {
  Proxy_MPI_Type_create_keyval proxyObj(type_copy_attr_fn, type_delete_attr_fn, type_keyval, extra_state);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_create_indexed_block : public Proxy_MPI {
  int  return_value;
  int  count;
  int  blocklength;
  int * array_of_displacements;
  MPI_Datatype  oldtype;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_create_indexed_block(int  count, int  blocklength, int * array_of_displacements, MPI_Datatype  oldtype, MPI_Datatype * newtype)
  : count(count), blocklength(blocklength), array_of_displacements(array_of_displacements), oldtype(oldtype), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_indexed_block(count, blocklength, array_of_displacements, oldtype, newtype);
  }
  int id_for_debug() { return 236; }
};
int  MPIWrapperFunneled::_MPI_Type_create_indexed_block(int  count, int  blocklength, int * array_of_displacements, MPI_Datatype  oldtype, MPI_Datatype * newtype) {
  Proxy_MPI_Type_create_indexed_block proxyObj(count, blocklength, array_of_displacements, oldtype, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Type_create_struct : public Proxy_MPI {
  int  return_value;
  int  count;
  int * array_of_block_lengths;
  MPI_Aint * array_of_displacements;
  MPI_Datatype * array_of_types;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_create_struct(int  count, int * array_of_block_lengths, MPI_Aint * array_of_displacements, MPI_Datatype * array_of_types, MPI_Datatype * newtype)
  : count(count), array_of_block_lengths(array_of_block_lengths), array_of_displacements(array_of_displacements), array_of_types(array_of_types), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_struct(count, array_of_block_lengths, array_of_displacements, array_of_types, newtype);
  }
  int id_for_debug() { return 237; }
};
int  MPIWrapperFunneled::_MPI_Type_create_struct(int  count, int * array_of_block_lengths, MPI_Aint * array_of_displacements, MPI_Datatype * array_of_types, MPI_Datatype * newtype) {
  Proxy_MPI_Type_create_struct proxyObj(count, array_of_block_lengths, array_of_displacements, array_of_types, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_create_subarray : public Proxy_MPI {
  int  return_value;
  int  ndims;
  int * size_array;
  int * subsize_array;
  int * start_array;
  int  order;
  MPI_Datatype  oldtype;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_create_subarray(int  ndims, int * size_array, int * subsize_array, int * start_array, int  order, MPI_Datatype  oldtype, MPI_Datatype * newtype)
  : ndims(ndims), size_array(size_array), subsize_array(subsize_array), start_array(start_array), order(order), oldtype(oldtype), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_subarray(ndims, size_array, subsize_array, start_array, order, oldtype, newtype);
  }
  int id_for_debug() { return 238; }
};
int  MPIWrapperFunneled::_MPI_Type_create_subarray(int  ndims, int * size_array, int * subsize_array, int * start_array, int  order, MPI_Datatype  oldtype, MPI_Datatype * newtype) {
  Proxy_MPI_Type_create_subarray proxyObj(ndims, size_array, subsize_array, start_array, order, oldtype, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_create_resized : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  oldtype;
  MPI_Aint  lb;
  MPI_Aint  extent;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_create_resized(MPI_Datatype  oldtype, MPI_Aint  lb, MPI_Aint  extent, MPI_Datatype * newtype)
  : oldtype(oldtype), lb(lb), extent(extent), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_create_resized(oldtype, lb, extent, newtype);
  }
  int id_for_debug() { return 239; }
};
int  MPIWrapperFunneled::_MPI_Type_create_resized(MPI_Datatype  oldtype, MPI_Aint  lb, MPI_Aint  extent, MPI_Datatype * newtype) {
  Proxy_MPI_Type_create_resized proxyObj(oldtype, lb, extent, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Type_delete_attr : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  type;
  int  type_keyval;
  Proxy_MPI_Type_delete_attr(MPI_Datatype  type, int  type_keyval)
  : type(type), type_keyval(type_keyval)
  {}
  void callBack() {
    return_value = ::MPI_Type_delete_attr(type, type_keyval);
  }
  int id_for_debug() { return 240; }
};
int  MPIWrapperFunneled::_MPI_Type_delete_attr(MPI_Datatype  type, int  type_keyval) {
  Proxy_MPI_Type_delete_attr proxyObj(type, type_keyval);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Type_dup : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  type;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_dup(MPI_Datatype  type, MPI_Datatype * newtype)
  : type(type), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_dup(type, newtype);
  }
  int id_for_debug() { return 241; }
};
int  MPIWrapperFunneled::_MPI_Type_dup(MPI_Datatype  type, MPI_Datatype * newtype) {
  Proxy_MPI_Type_dup proxyObj(type, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_free : public Proxy_MPI {
  int  return_value;
  MPI_Datatype * type;
  Proxy_MPI_Type_free(MPI_Datatype * type)
  : type(type)
  {}
  void callBack() {
    return_value = ::MPI_Type_free(type);
  }
  int id_for_debug() { return 243; }
};
int  MPIWrapperFunneled::_MPI_Type_free(MPI_Datatype * type) {
  Proxy_MPI_Type_free proxyObj(type);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Type_free_keyval : public Proxy_MPI {
  int  return_value;
  int * type_keyval;
  Proxy_MPI_Type_free_keyval(int * type_keyval)
  : type_keyval(type_keyval)
  {}
  void callBack() {
    return_value = ::MPI_Type_free_keyval(type_keyval);
  }
  int id_for_debug() { return 244; }
};
int  MPIWrapperFunneled::_MPI_Type_free_keyval(int * type_keyval) {
  Proxy_MPI_Type_free_keyval proxyObj(type_keyval);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_get_attr : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  type;
  int  type_keyval;
  void * attribute_val;
  int * flag;
  Proxy_MPI_Type_get_attr(MPI_Datatype  type, int  type_keyval, void * attribute_val, int * flag)
  : type(type), type_keyval(type_keyval), attribute_val(attribute_val), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Type_get_attr(type, type_keyval, attribute_val, flag);
  }
  int id_for_debug() { return 246; }
};
int  MPIWrapperFunneled::_MPI_Type_get_attr(MPI_Datatype  type, int  type_keyval, void * attribute_val, int * flag) {
  Proxy_MPI_Type_get_attr proxyObj(type, type_keyval, attribute_val, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Type_get_contents : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  mtype;
  int  max_integers;
  int  max_addresses;
  int  max_datatypes;
  int * array_of_integers;
  MPI_Aint * array_of_addresses;
  MPI_Datatype * array_of_datatypes;
  Proxy_MPI_Type_get_contents(MPI_Datatype  mtype, int  max_integers, int  max_addresses, int  max_datatypes, int * array_of_integers, MPI_Aint * array_of_addresses, MPI_Datatype * array_of_datatypes)
  : mtype(mtype), max_integers(max_integers), max_addresses(max_addresses), max_datatypes(max_datatypes), array_of_integers(array_of_integers), array_of_addresses(array_of_addresses), array_of_datatypes(array_of_datatypes)
  {}
  void callBack() {
    return_value = ::MPI_Type_get_contents(mtype, max_integers, max_addresses, max_datatypes, array_of_integers, array_of_addresses, array_of_datatypes);
  }
  int id_for_debug() { return 247; }
};
int  MPIWrapperFunneled::_MPI_Type_get_contents(MPI_Datatype  mtype, int  max_integers, int  max_addresses, int  max_datatypes, int * array_of_integers, MPI_Aint * array_of_addresses, MPI_Datatype * array_of_datatypes) {
  Proxy_MPI_Type_get_contents proxyObj(mtype, max_integers, max_addresses, max_datatypes, array_of_integers, array_of_addresses, array_of_datatypes);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_get_envelope : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  type;
  int * num_integers;
  int * num_addresses;
  int * num_datatypes;
  int * combiner;
  Proxy_MPI_Type_get_envelope(MPI_Datatype  type, int * num_integers, int * num_addresses, int * num_datatypes, int * combiner)
  : type(type), num_integers(num_integers), num_addresses(num_addresses), num_datatypes(num_datatypes), combiner(combiner)
  {}
  void callBack() {
    return_value = ::MPI_Type_get_envelope(type, num_integers, num_addresses, num_datatypes, combiner);
  }
  int id_for_debug() { return 248; }
};
int  MPIWrapperFunneled::_MPI_Type_get_envelope(MPI_Datatype  type, int * num_integers, int * num_addresses, int * num_datatypes, int * combiner) {
  Proxy_MPI_Type_get_envelope proxyObj(type, num_integers, num_addresses, num_datatypes, combiner);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_get_extent : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  type;
  MPI_Aint * lb;
  MPI_Aint * extent;
  Proxy_MPI_Type_get_extent(MPI_Datatype  type, MPI_Aint * lb, MPI_Aint * extent)
  : type(type), lb(lb), extent(extent)
  {}
  void callBack() {
    return_value = ::MPI_Type_get_extent(type, lb, extent);
  }
  int id_for_debug() { return 249; }
};
int  MPIWrapperFunneled::_MPI_Type_get_extent(MPI_Datatype  type, MPI_Aint * lb, MPI_Aint * extent) {
  Proxy_MPI_Type_get_extent proxyObj(type, lb, extent);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_get_name : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  type;
  char * type_name;
  int * resultlen;
  Proxy_MPI_Type_get_name(MPI_Datatype  type, char * type_name, int * resultlen)
  : type(type), type_name(type_name), resultlen(resultlen)
  {}
  void callBack() {
    return_value = ::MPI_Type_get_name(type, type_name, resultlen);
  }
  int id_for_debug() { return 250; }
};
int  MPIWrapperFunneled::_MPI_Type_get_name(MPI_Datatype  type, char * type_name, int * resultlen) {
  Proxy_MPI_Type_get_name proxyObj(type, type_name, resultlen);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Type_get_true_extent : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  datatype;
  MPI_Aint * true_lb;
  MPI_Aint * true_extent;
  Proxy_MPI_Type_get_true_extent(MPI_Datatype  datatype, MPI_Aint * true_lb, MPI_Aint * true_extent)
  : datatype(datatype), true_lb(true_lb), true_extent(true_extent)
  {}
  void callBack() {
    return_value = ::MPI_Type_get_true_extent(datatype, true_lb, true_extent);
  }
  int id_for_debug() { return 251; }
};
int  MPIWrapperFunneled::_MPI_Type_get_true_extent(MPI_Datatype  datatype, MPI_Aint * true_lb, MPI_Aint * true_extent) {
  Proxy_MPI_Type_get_true_extent proxyObj(datatype, true_lb, true_extent);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Type_indexed : public Proxy_MPI {
  int  return_value;
  int  count;
  int * array_of_blocklengths;
  int * array_of_displacements;
  MPI_Datatype  oldtype;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_indexed(int  count, int * array_of_blocklengths, int * array_of_displacements, MPI_Datatype  oldtype, MPI_Datatype * newtype)
  : count(count), array_of_blocklengths(array_of_blocklengths), array_of_displacements(array_of_displacements), oldtype(oldtype), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_indexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
  }
  int id_for_debug() { return 254; }
};
int  MPIWrapperFunneled::_MPI_Type_indexed(int  count, int * array_of_blocklengths, int * array_of_displacements, MPI_Datatype  oldtype, MPI_Datatype * newtype) {
  Proxy_MPI_Type_indexed proxyObj(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Type_match_size : public Proxy_MPI {
  int  return_value;
  int  typeclass;
  int  size;
  MPI_Datatype * type;
  Proxy_MPI_Type_match_size(int  typeclass, int  size, MPI_Datatype * type)
  : typeclass(typeclass), size(size), type(type)
  {}
  void callBack() {
    return_value = ::MPI_Type_match_size(typeclass, size, type);
  }
  int id_for_debug() { return 256; }
};
int  MPIWrapperFunneled::_MPI_Type_match_size(int  typeclass, int  size, MPI_Datatype * type) {
  Proxy_MPI_Type_match_size proxyObj(typeclass, size, type);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_set_attr : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  type;
  int  type_keyval;
  void * attr_val;
  Proxy_MPI_Type_set_attr(MPI_Datatype  type, int  type_keyval, void * attr_val)
  : type(type), type_keyval(type_keyval), attr_val(attr_val)
  {}
  void callBack() {
    return_value = ::MPI_Type_set_attr(type, type_keyval, attr_val);
  }
  int id_for_debug() { return 257; }
};
int  MPIWrapperFunneled::_MPI_Type_set_attr(MPI_Datatype  type, int  type_keyval, void * attr_val) {
  Proxy_MPI_Type_set_attr proxyObj(type, type_keyval, attr_val);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

#endif



struct Proxy_MPI_Type_set_name : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  type;
  char * type_name;
  Proxy_MPI_Type_set_name(MPI_Datatype  type, char * type_name)
  : type(type), type_name(type_name)
  {}
  void callBack() {
    return_value = ::MPI_Type_set_name(type, type_name);
  }
  int id_for_debug() { return 258; }
};
int  MPIWrapperFunneled::_MPI_Type_set_name(MPI_Datatype  type, char * type_name) {
  Proxy_MPI_Type_set_name proxyObj(type, type_name);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_size : public Proxy_MPI {
  int  return_value;
  MPI_Datatype  type;
  int * size;
  Proxy_MPI_Type_size(MPI_Datatype  type, int * size)
  : type(type), size(size)
  {}
  void callBack() {
    return_value = ::MPI_Type_size(type, size);
  }
  int id_for_debug() { return 259; }
};
int  MPIWrapperFunneled::_MPI_Type_size(MPI_Datatype  type, int * size) {
  Proxy_MPI_Type_size proxyObj(type, size);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Type_vector : public Proxy_MPI {
  int  return_value;
  int  count;
  int  blocklength;
  int  stride;
  MPI_Datatype  oldtype;
  MPI_Datatype * newtype;
  Proxy_MPI_Type_vector(int  count, int  blocklength, int  stride, MPI_Datatype  oldtype, MPI_Datatype * newtype)
  : count(count), blocklength(blocklength), stride(stride), oldtype(oldtype), newtype(newtype)
  {}
  void callBack() {
    return_value = ::MPI_Type_vector(count, blocklength, stride, oldtype, newtype);
  }
  int id_for_debug() { return 262; }
};
int  MPIWrapperFunneled::_MPI_Type_vector(int  count, int  blocklength, int  stride, MPI_Datatype  oldtype, MPI_Datatype * newtype) {
  Proxy_MPI_Type_vector proxyObj(count, blocklength, stride, oldtype, newtype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Unpack : public Proxy_MPI {
  int  return_value;
  void * inbuf;
  int  insize;
  int * position;
  void * outbuf;
  int  outcount;
  MPI_Datatype  datatype;
  MPI_Comm  comm;
  Proxy_MPI_Unpack(void * inbuf, int  insize, int * position, void * outbuf, int  outcount, MPI_Datatype  datatype, MPI_Comm  comm)
  : inbuf(inbuf), insize(insize), position(position), outbuf(outbuf), outcount(outcount), datatype(datatype), comm(comm)
  {}
  void callBack() {
    return_value = ::MPI_Unpack(inbuf, insize, position, outbuf, outcount, datatype, comm);
  }
  int id_for_debug() { return 263; }
};
int  MPIWrapperFunneled::_MPI_Unpack(void * inbuf, int  insize, int * position, void * outbuf, int  outcount, MPI_Datatype  datatype, MPI_Comm  comm) {
  Proxy_MPI_Unpack proxyObj(inbuf, insize, position, outbuf, outcount, datatype, comm);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1

struct Proxy_MPI_Unpublish_name : public Proxy_MPI {
  int  return_value;
  char * service_name;
  MPI_Info  info;
  char * port_name;
  Proxy_MPI_Unpublish_name(char * service_name, MPI_Info  info, char * port_name)
  : service_name(service_name), info(info), port_name(port_name)
  {}
  void callBack() {
    return_value = ::MPI_Unpublish_name(service_name, info, port_name);
  }
  int id_for_debug() { return 264; }
};
int  MPIWrapperFunneled::_MPI_Unpublish_name(char * service_name, MPI_Info  info, char * port_name) {
  Proxy_MPI_Unpublish_name proxyObj(service_name, info, port_name);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Unpack_external  : public Proxy_MPI {
  int  return_value;
  char * datarep;
  void * inbuf;
  MPI_Aint  insize;
  MPI_Aint * position;
  void * outbuf;
  int  outcount;
  MPI_Datatype  datatype;
  Proxy_MPI_Unpack_external (char * datarep, void * inbuf, MPI_Aint  insize, MPI_Aint * position, void * outbuf, int  outcount, MPI_Datatype  datatype)
  : datarep(datarep), inbuf(inbuf), insize(insize), position(position), outbuf(outbuf), outcount(outcount), datatype(datatype)
  {}
  void callBack() {
    return_value = ::MPI_Unpack_external (datarep, inbuf, insize, position, outbuf, outcount, datatype);
  }
  int id_for_debug() { return 265; }
};
int  MPIWrapperFunneled::_MPI_Unpack_external (char * datarep, void * inbuf, MPI_Aint  insize, MPI_Aint * position, void * outbuf, int  outcount, MPI_Datatype  datatype) {
  Proxy_MPI_Unpack_external  proxyObj(datarep, inbuf, insize, position, outbuf, outcount, datatype);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}


#endif



#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1



struct Proxy_MPI_Win_call_errhandler : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  int  errorcode;
  Proxy_MPI_Win_call_errhandler(MPI_Win  win, int  errorcode)
  : win(win), errorcode(errorcode)
  {}
  void callBack() {
    return_value = ::MPI_Win_call_errhandler(win, errorcode);
  }
  int id_for_debug() { return 267; }
};
int  MPIWrapperFunneled::_MPI_Win_call_errhandler(MPI_Win  win, int  errorcode) {
  Proxy_MPI_Win_call_errhandler proxyObj(win, errorcode);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_complete : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  Proxy_MPI_Win_complete(MPI_Win  win)
  : win(win)
  {}
  void callBack() {
    return_value = ::MPI_Win_complete(win);
  }
  int id_for_debug() { return 268; }
};
int  MPIWrapperFunneled::_MPI_Win_complete(MPI_Win  win) {
  Proxy_MPI_Win_complete proxyObj(win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_create : public Proxy_MPI {
  int  return_value;
  void * base;
  MPI_Aint  size;
  int  disp_unit;
  MPI_Info  info;
  MPI_Comm  comm;
  MPI_Win * win;
  Proxy_MPI_Win_create(void * base, MPI_Aint  size, int  disp_unit, MPI_Info  info, MPI_Comm  comm, MPI_Win * win)
  : base(base), size(size), disp_unit(disp_unit), info(info), comm(comm), win(win)
  {}
  void callBack() {
    return_value = ::MPI_Win_create(base, size, disp_unit, info, comm, win);
  }
  int id_for_debug() { return 269; }
};
int  MPIWrapperFunneled::_MPI_Win_create(void * base, MPI_Aint  size, int  disp_unit, MPI_Info  info, MPI_Comm  comm, MPI_Win * win) {
  Proxy_MPI_Win_create proxyObj(base, size, disp_unit, info, comm, win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_create_keyval : public Proxy_MPI {
  int  return_value;
  MPI_Win_copy_attr_function * win_copy_attr_fn;
  MPI_Win_delete_attr_function * win_delete_attr_fn;
  int * win_keyval;
  void * extra_state;
  Proxy_MPI_Win_create_keyval(MPI_Win_copy_attr_function * win_copy_attr_fn, MPI_Win_delete_attr_function * win_delete_attr_fn, int * win_keyval, void * extra_state)
  : win_copy_attr_fn(win_copy_attr_fn), win_delete_attr_fn(win_delete_attr_fn), win_keyval(win_keyval), extra_state(extra_state)
  {}
  void callBack() {
    return_value = ::MPI_Win_create_keyval(win_copy_attr_fn, win_delete_attr_fn, win_keyval, extra_state);
  }
  int id_for_debug() { return 271; }
};
int  MPIWrapperFunneled::_MPI_Win_create_keyval(MPI_Win_copy_attr_function * win_copy_attr_fn, MPI_Win_delete_attr_function * win_delete_attr_fn, int * win_keyval, void * extra_state) {
  Proxy_MPI_Win_create_keyval proxyObj(win_copy_attr_fn, win_delete_attr_fn, win_keyval, extra_state);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_delete_attr : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  int  win_keyval;
  Proxy_MPI_Win_delete_attr(MPI_Win  win, int  win_keyval)
  : win(win), win_keyval(win_keyval)
  {}
  void callBack() {
    return_value = ::MPI_Win_delete_attr(win, win_keyval);
  }
  int id_for_debug() { return 272; }
};
int  MPIWrapperFunneled::_MPI_Win_delete_attr(MPI_Win  win, int  win_keyval) {
  Proxy_MPI_Win_delete_attr proxyObj(win, win_keyval);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_fence : public Proxy_MPI {
  int  return_value;
  int  assert_;
  MPI_Win  win;
  Proxy_MPI_Win_fence(int  assert_, MPI_Win  win)
  : assert_(assert_), win(win)
  {}
  void callBack() {
    return_value = ::MPI_Win_fence(assert_, win);
  }
  int id_for_debug() { return 274; }
};
int  MPIWrapperFunneled::_MPI_Win_fence(int  assert_, MPI_Win  win) {
  Proxy_MPI_Win_fence proxyObj(assert_, win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_free : public Proxy_MPI {
  int  return_value;
  MPI_Win * win;
  Proxy_MPI_Win_free(MPI_Win * win)
  : win(win)
  {}
  void callBack() {
    return_value = ::MPI_Win_free(win);
  }
  int id_for_debug() { return 275; }
};
int  MPIWrapperFunneled::_MPI_Win_free(MPI_Win * win) {
  Proxy_MPI_Win_free proxyObj(win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_free_keyval : public Proxy_MPI {
  int  return_value;
  int * win_keyval;
  Proxy_MPI_Win_free_keyval(int * win_keyval)
  : win_keyval(win_keyval)
  {}
  void callBack() {
    return_value = ::MPI_Win_free_keyval(win_keyval);
  }
  int id_for_debug() { return 276; }
};
int  MPIWrapperFunneled::_MPI_Win_free_keyval(int * win_keyval) {
  Proxy_MPI_Win_free_keyval proxyObj(win_keyval);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_get_attr : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  int  win_keyval;
  void * attribute_val;
  int * flag;
  Proxy_MPI_Win_get_attr(MPI_Win  win, int  win_keyval, void * attribute_val, int * flag)
  : win(win), win_keyval(win_keyval), attribute_val(attribute_val), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Win_get_attr(win, win_keyval, attribute_val, flag);
  }
  int id_for_debug() { return 277; }
};
int  MPIWrapperFunneled::_MPI_Win_get_attr(MPI_Win  win, int  win_keyval, void * attribute_val, int * flag) {
  Proxy_MPI_Win_get_attr proxyObj(win, win_keyval, attribute_val, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_get_errhandler : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  MPI_Errhandler * errhandler;
  Proxy_MPI_Win_get_errhandler(MPI_Win  win, MPI_Errhandler * errhandler)
  : win(win), errhandler(errhandler)
  {}
  void callBack() {
    return_value = ::MPI_Win_get_errhandler(win, errhandler);
  }
  int id_for_debug() { return 278; }
};
int  MPIWrapperFunneled::_MPI_Win_get_errhandler(MPI_Win  win, MPI_Errhandler * errhandler) {
  Proxy_MPI_Win_get_errhandler proxyObj(win, errhandler);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_get_group : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  MPI_Group * group;
  Proxy_MPI_Win_get_group(MPI_Win  win, MPI_Group * group)
  : win(win), group(group)
  {}
  void callBack() {
    return_value = ::MPI_Win_get_group(win, group);
  }
  int id_for_debug() { return 279; }
};
int  MPIWrapperFunneled::_MPI_Win_get_group(MPI_Win  win, MPI_Group * group) {
  Proxy_MPI_Win_get_group proxyObj(win, group);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_get_name : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  char * win_name;
  int * resultlen;
  Proxy_MPI_Win_get_name(MPI_Win  win, char * win_name, int * resultlen)
  : win(win), win_name(win_name), resultlen(resultlen)
  {}
  void callBack() {
    return_value = ::MPI_Win_get_name(win, win_name, resultlen);
  }
  int id_for_debug() { return 280; }
};
int  MPIWrapperFunneled::_MPI_Win_get_name(MPI_Win  win, char * win_name, int * resultlen) {
  Proxy_MPI_Win_get_name proxyObj(win, win_name, resultlen);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_lock : public Proxy_MPI {
  int  return_value;
  int  lock_type;
  int  rank;
  int  assert_;
  MPI_Win  win;
  Proxy_MPI_Win_lock(int  lock_type, int  rank, int  assert_, MPI_Win  win)
  : lock_type(lock_type), rank(rank), assert_(assert_), win(win)
  {}
  void callBack() {
    return_value = ::MPI_Win_lock(lock_type, rank, assert_, win);
  }
  int id_for_debug() { return 281; }
};
int  MPIWrapperFunneled::_MPI_Win_lock(int  lock_type, int  rank, int  assert_, MPI_Win  win) {
  Proxy_MPI_Win_lock proxyObj(lock_type, rank, assert_, win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_post : public Proxy_MPI {
  int  return_value;
  MPI_Group  group;
  int  assert_;
  MPI_Win  win;
  Proxy_MPI_Win_post(MPI_Group  group, int  assert_, MPI_Win  win)
  : group(group), assert_(assert_), win(win)
  {}
  void callBack() {
    return_value = ::MPI_Win_post(group, assert_, win);
  }
  int id_for_debug() { return 282; }
};
int  MPIWrapperFunneled::_MPI_Win_post(MPI_Group  group, int  assert_, MPI_Win  win) {
  Proxy_MPI_Win_post proxyObj(group, assert_, win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_set_attr : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  int  win_keyval;
  void * attribute_val;
  Proxy_MPI_Win_set_attr(MPI_Win  win, int  win_keyval, void * attribute_val)
  : win(win), win_keyval(win_keyval), attribute_val(attribute_val)
  {}
  void callBack() {
    return_value = ::MPI_Win_set_attr(win, win_keyval, attribute_val);
  }
  int id_for_debug() { return 283; }
};
int  MPIWrapperFunneled::_MPI_Win_set_attr(MPI_Win  win, int  win_keyval, void * attribute_val) {
  Proxy_MPI_Win_set_attr proxyObj(win, win_keyval, attribute_val);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_set_errhandler : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  MPI_Errhandler  errhandler;
  Proxy_MPI_Win_set_errhandler(MPI_Win  win, MPI_Errhandler  errhandler)
  : win(win), errhandler(errhandler)
  {}
  void callBack() {
    return_value = ::MPI_Win_set_errhandler(win, errhandler);
  }
  int id_for_debug() { return 284; }
};
int  MPIWrapperFunneled::_MPI_Win_set_errhandler(MPI_Win  win, MPI_Errhandler  errhandler) {
  Proxy_MPI_Win_set_errhandler proxyObj(win, errhandler);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_set_name : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  char * win_name;
  Proxy_MPI_Win_set_name(MPI_Win  win, char * win_name)
  : win(win), win_name(win_name)
  {}
  void callBack() {
    return_value = ::MPI_Win_set_name(win, win_name);
  }
  int id_for_debug() { return 285; }
};
int  MPIWrapperFunneled::_MPI_Win_set_name(MPI_Win  win, char * win_name) {
  Proxy_MPI_Win_set_name proxyObj(win, win_name);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_start : public Proxy_MPI {
  int  return_value;
  MPI_Group  group;
  int  assert_;
  MPI_Win  win;
  Proxy_MPI_Win_start(MPI_Group  group, int  assert_, MPI_Win  win)
  : group(group), assert_(assert_), win(win)
  {}
  void callBack() {
    return_value = ::MPI_Win_start(group, assert_, win);
  }
  int id_for_debug() { return 286; }
};
int  MPIWrapperFunneled::_MPI_Win_start(MPI_Group  group, int  assert_, MPI_Win  win) {
  Proxy_MPI_Win_start proxyObj(group, assert_, win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_test : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  int * flag;
  Proxy_MPI_Win_test(MPI_Win  win, int * flag)
  : win(win), flag(flag)
  {}
  void callBack() {
    return_value = ::MPI_Win_test(win, flag);
  }
  int id_for_debug() { return 287; }
};
int  MPIWrapperFunneled::_MPI_Win_test(MPI_Win  win, int * flag) {
  Proxy_MPI_Win_test proxyObj(win, flag);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_unlock : public Proxy_MPI {
  int  return_value;
  int  rank;
  MPI_Win  win;
  Proxy_MPI_Win_unlock(int  rank, MPI_Win  win)
  : rank(rank), win(win)
  {}
  void callBack() {
    return_value = ::MPI_Win_unlock(rank, win);
  }
  int id_for_debug() { return 288; }
};
int  MPIWrapperFunneled::_MPI_Win_unlock(int  rank, MPI_Win  win) {
  Proxy_MPI_Win_unlock proxyObj(rank, win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Win_wait : public Proxy_MPI {
  int  return_value;
  MPI_Win  win;
  Proxy_MPI_Win_wait(MPI_Win  win)
  : win(win)
  {}
  void callBack() {
    return_value = ::MPI_Win_wait(win);
  }
  int id_for_debug() { return 289; }
};
int  MPIWrapperFunneled::_MPI_Win_wait(MPI_Win  win) {
  Proxy_MPI_Win_wait proxyObj(win);
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

// End of MPI_Win functions.

#endif



struct Proxy_MPI_Wtick : public Proxy_MPI {
  double  return_value;
  Proxy_MPI_Wtick()
  {}
  void callBack() {
    return_value = ::MPI_Wtick();
  }
  int id_for_debug() { return 290; }
};
double  MPIWrapperFunneled::_MPI_Wtick() {
  Proxy_MPI_Wtick proxyObj;
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}



struct Proxy_MPI_Wtime : public Proxy_MPI {
  double  return_value;
  Proxy_MPI_Wtime()
  {}
  void callBack() {
    return_value = ::MPI_Wtime();
  }
  int id_for_debug() { return 291; }
};
double  MPIWrapperFunneled::_MPI_Wtime() {
  Proxy_MPI_Wtime proxyObj;
  Proxy_MPI* proxy = &proxyObj;
  LockMutex(mutex_id_non_blocking_list);
  pendingNonBlocking.push_back(proxy);
  UnlockMutex(mutex_id_non_blocking_list);
  wake_up_executor();
  proxy->LockMutex();
  return proxyObj.return_value;
}

// Extra non-MPI routine: probe for any of the tags in given list. 
int MPIWrapperFunneled::NonMPI_Iprobe_multiple(int src, const int tagList[], 
				       int nTags, MPI_Comm comm, 
				       int *flag, MPI_Status *stat)
{
  int rc;
  for(int i = 0; i < nTags; i++)
    {
      rc = MPI_Iprobe(src, tagList[i], comm, flag, stat);
      if(*flag) return rc;
    }
  return rc;
}





