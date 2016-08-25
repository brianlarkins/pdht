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
#include "services/Service_parent.h"
#include <cassert>
#include <stdexcept>
#include <cstring>
#include "utilities/cht_utils.h"
#include "services/service_manager/ServiceManager.h"

namespace cht {

  void ServiceManager::Init(MPI_Comm* comm_to_workers_) {
    comm_to_workers = comm_to_workers_;
    MW._MPI_Comm_remote_size(*comm_to_workers, &n_workers);
    if(n_workers < 1)
      throw std::runtime_error("Error in ServiceManager::Init(): (n_workers < 1).");
  }

  void ServiceManager::Finalize() {
    for (int dest = 0; dest < n_workers; ++dest) {
      MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, dest, 
		  Tag_worker_shutdown, *comm_to_workers);
    }
    //  MW._MPI_Comm_free(&comm_to_workers);
  }

  void ServiceManager::startService(std::string serviceID) {
    if ( comm_to_workers == NULL ) 
      throw std::runtime_error("ServiceManager not initialized");
    ServiceID_Comm_Map::iterator i = runningServices.find( serviceID );
    if ( i != runningServices.end() ) 
      throw std::runtime_error("ServiceManager::startService - Service already started");
    copyServiceIDToBuf(serviceID);
    for (int dest = 0; dest < n_workers; ++dest)
      MW._MPI_Send(serviceIDBuf, Service::max_length_serviceID_string, 
		  MPI_UNSIGNED_CHAR, dest, 
		  Tag_start_service, *comm_to_workers);
    for (int src = 0; src < n_workers; ++src)
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, src, 
		  Tag_start_service_ack1, *comm_to_workers, 
		  MPI_STATUS_IGNORE);
    MPI_Comm* comm_to_workers_copy = new MPI_Comm;
    // NOTE: workers must call this function as well!!  ("in the correct order"??)
    MW._MPI_Comm_dup(*comm_to_workers, comm_to_workers_copy);
    runningServices.insert
      ( ServiceID_Comm_Map::value_type(serviceID, comm_to_workers_copy) );
    Service::Parent* new_service = cht::obj_factory<Service::Parent>::instance().createObject(serviceID);
    new_service->start(comm_to_workers_copy);
    // Wait for second acknowledgement from each worker. This is to make
    // sure parent does not return until workers have really started the
    // service.
    for (int src = 0; src < n_workers; ++src)
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, src, 
		  Tag_start_service_ack2, *comm_to_workers, 
		  MPI_STATUS_IGNORE);
  }

  void ServiceManager::stopService(std::string serviceID) {
    if ( comm_to_workers == NULL ) 
      throw std::runtime_error("ServiceManager not initialized");
    ServiceID_Comm_Map::iterator i = runningServices.find( serviceID );
    if ( i == runningServices.end() ) 
      throw std::runtime_error("Service not started");
    copyServiceIDToBuf(serviceID);
    for (int dest = 0; dest < n_workers; ++dest) 
      MW._MPI_Send(serviceIDBuf, Service::max_length_serviceID_string, 
		  MPI_UNSIGNED_CHAR, dest, 
		  Tag_stop_service, *comm_to_workers);
    for (int src = 0; src < n_workers; ++src) 
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, src, 
		  Tag_stop_service_ack1, *comm_to_workers, 
		  MPI_STATUS_IGNORE);
    Service::Parent* service = cht::obj_factory<Service::Parent>::instance().createObject(serviceID);
    service->stop();
    MW._MPI_Comm_free( (*i).second );
    delete (*i).second;
    runningServices.erase(i);
    // Wait for second acknowledgement from each worker. This is to make
    // sure parent does not return until workers have really stopped the
    // service.
    for (int src = 0; src < n_workers; ++src) 
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, src, 
		  Tag_stop_service_ack2, *comm_to_workers, 
		  MPI_STATUS_IGNORE);
  }

  void ServiceManager::resetStatistics(std::string serviceID) {
    if ( comm_to_workers == NULL ) 
      throw std::runtime_error("ServiceManager not initialized");
    ServiceID_Comm_Map::iterator i = runningServices.find( serviceID );
    if ( i == runningServices.end() ) 
      throw std::runtime_error("Service not started");
    copyServiceIDToBuf(serviceID);
    for (int dest = 0; dest < n_workers; ++dest) 
      MW._MPI_Send(serviceIDBuf, Service::max_length_serviceID_string, 
		  MPI_UNSIGNED_CHAR, dest, 
		  Tag_reset_statistics, *comm_to_workers);
    for (int src = 0; src < n_workers; ++src) 
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, src, 
		  Tag_reset_statistics_ack, *comm_to_workers, 
		  MPI_STATUS_IGNORE);
    Service::Parent* service = cht::obj_factory<Service::Parent>::instance().createObject(serviceID);
    service->resetStatistics();
  }
  void ServiceManager::reportStatistics(std::string serviceID) {
    if ( comm_to_workers == NULL ) 
      throw std::runtime_error("ServiceManager not initialized");
    ServiceID_Comm_Map::iterator i = runningServices.find( serviceID );
    if ( i == runningServices.end() ) 
      throw std::runtime_error("Service not started");
    copyServiceIDToBuf(serviceID);
    for (int dest = 0; dest < n_workers; ++dest) 
      MW._MPI_Send(serviceIDBuf, Service::max_length_serviceID_string, 
		  MPI_UNSIGNED_CHAR, dest, 
		  Tag_report_statistics, *comm_to_workers);
    for (int src = 0; src < n_workers; ++src) 
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, src, 
		  Tag_report_statistics_ack, *comm_to_workers, 
		  MPI_STATUS_IGNORE);
    Service::Parent* service = cht::obj_factory<Service::Parent>::instance().createObject(serviceID);
    std::string messageHeader;
    messageHeader = "=====  " + serviceID + "::reportStatistics()  =====";
    service->reportStatistics(messageHeader);
  }


  ServiceManager::ServiceManager() 
    : n_workers(0), comm_to_workers(NULL), MW( MPI_Wrapper::instance() )
  {
  
  }


  void ServiceManager::copyServiceIDToBuf(std::string serviceID) {
    assert(serviceID.length() <= Service::max_length_serviceID_string);
    strcpy ( serviceIDBuf, serviceID.c_str() );
  }

}; // end namespace
