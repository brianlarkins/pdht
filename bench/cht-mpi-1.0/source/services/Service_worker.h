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
#ifndef SERVICE_WORKER_HEADER
#define SERVICE_WORKER_HEADER

#include "services/MPIWrapperInclude.h"
#include "services/Service.h"
#include "utilities/cht_threads.h"

namespace cht {
  namespace Service {

    class Worker : public Base
    {
    public:
      Worker();
      static std::string object_base_type_id();
      void start(MPI_Comm* comm_to_parent, 
		 MPI_Comm* comm_worker_world);
      void stop();
      virtual void resetStatistics();
      virtual void reportStatistics(std::string messageHeader);  
    private:
      bool service_running;
      bool key_access_allowed;
      int number_of_access_keys;
      Threads::Mutex mutex;   
      void LockMutex() {
	mutex.lock();
      }
      void UnlockMutex() {
	mutex.unlock();
      }
      MPI_Comm* comm_to_parent;
      MPI_Comm* comm_worker_world;
      int n_workers;
      int my_rank;

    protected:
      class AccessKey {
	Worker* ptr;
      private:
	// no copying of keys
	AccessKey( AccessKey const & ); 
	AccessKey const & operator=( AccessKey const & );
      public:
      AccessKey(Worker* ptr_)
	:ptr(ptr_) {
	  ptr->LockMutex();
	  if ( !ptr->key_access_allowed )
	    throw std::runtime_error("Attempt to acquire access key when not allowed.");
	  ++ptr->number_of_access_keys;
	  ptr->UnlockMutex();
	}
	~AccessKey() {
	  ptr->LockMutex();
	  --ptr->number_of_access_keys;
	  ptr->UnlockMutex();
	}
	MPI_Comm* comm_to_parent() const { return ptr->comm_to_parent; }
	MPI_Comm* comm_worker_world() const { return ptr->comm_worker_world; }
	int n_workers() const { return ptr->n_workers; }
	int my_rank() const { return ptr->my_rank; }
      };

      MPI_Wrapper& MW;
      virtual void start_derived() = 0;
      virtual void stop_derived() = 0;
      bool serviceIsRunning() {
	bool tmp;
	LockMutex();
	tmp = service_running;
	UnlockMutex();
	return tmp;
      }
    };

  }; // end namespace
}; // end namespace

#endif 
