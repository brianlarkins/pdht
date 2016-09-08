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
#include "Service_worker.h"
#include <string>
#include <stdexcept>

namespace cht {
  namespace Service {

    Worker::Worker() 
      :  MW( MPI_Wrapper::instance() ),
	 service_running(false), 
	 key_access_allowed(false),
	 number_of_access_keys(0),
	 comm_to_parent(NULL), comm_worker_world(NULL),
	 n_workers(0),
	 my_rank(0) { 
    }



    std::string Worker::object_base_type_id() {
      return "cht::Service::Worker";
    }
    void Worker::start(MPI_Comm* comm_to_parent_, 
		       MPI_Comm* comm_worker_world_) { 
      LockMutex();
      if ( service_running )
	throw std::runtime_error("Service start attempt on already running service.");
      if ( number_of_access_keys != 0 )
	throw std::runtime_error("Service start attempt when number_of_access_keys != 1");
      if (comm_to_parent_ == NULL)
	throw std::runtime_error("error in start(): comm_to_parent_ == NULL !!!");
      if (comm_worker_world_ == NULL)
	throw std::runtime_error("error in start(): comm_worker_world_ == NULL !!!");
      if (*comm_to_parent_ == MPI_COMM_NULL)
	throw std::runtime_error("error in start(): *comm_to_parent_ == MPI_COMM_NULL !!!");
      if (*comm_worker_world_ == MPI_COMM_NULL)
	throw std::runtime_error("error in start(): *comm_worker_world_ == MPI_COMM_NULL !!!");
      comm_to_parent = comm_to_parent_; 
      comm_worker_world = comm_worker_world_;
      MW._MPI_Comm_size(*comm_worker_world, &n_workers);
      MW._MPI_Comm_rank(*comm_worker_world, &my_rank);
      key_access_allowed = true;
      UnlockMutex();
      // Access granted to derived class via start/stop key
      this->start_derived();
      LockMutex();
      // Service is not really running until here
      service_running = true;
      UnlockMutex();
    }
    void Worker::stop() { 
      LockMutex();
      if ( !service_running )
	throw std::runtime_error("Service stop attempt on service not running.");
      service_running = false;
      UnlockMutex();
      // Access granted to derived class via start/stop key
      this->stop_derived();
      LockMutex();
      if ( number_of_access_keys != 0 )
	throw std::runtime_error("Service stop attempt when number_of_access_keys != 1");
      key_access_allowed = false;
      comm_to_parent = NULL; 
      comm_worker_world = NULL;
      n_workers = 0;
      my_rank = 0;
      UnlockMutex();
    }
    void Worker::resetStatistics() { }
    void Worker::reportStatistics(std::string messageHeader) { }
#if 0
    void Worker::start_derived() { }
    void Worker::stop_derived() { }
#endif
  }; // end namespace
}; // end namespace