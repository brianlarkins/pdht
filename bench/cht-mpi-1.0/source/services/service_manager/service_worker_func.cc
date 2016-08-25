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
#include "services/Service_worker.h"
#include <stdexcept>
#include <iostream>
#include <unistd.h>
#include "utilities/cht_utils.h"
#include "services/service_manager/ServiceManager.h"
#include "services/service_manager/service_worker_func.h"
#include "services/ParentInterface.h"

namespace cht {
  
  void service_worker_func(MPI_Comm comm_worker_world,
			   MPI_Comm comm_parent) {
    Threads::Thread workerMainThread("main");
    char serviceIDBuf[Service::max_length_serviceID_string];
    if (comm_parent == MPI_COMM_NULL) throw std::runtime_error("Error: (comm_parent == MPI_COMM_NULL).");
    MPI_Wrapper& MW = MPI_Wrapper::instance();
    ServiceManager::ServiceID_Comm_Map runningServicesParentComms;
    ServiceManager::ServiceID_Comm_Map runningServicesWorkerComms;

    // Prevent access to parent functions from workers
    ParentInterface::instance().i_am_worker();

    try {

      int myRank;
      MW._MPI_Comm_rank(comm_parent, &myRank);
      int myRank2;
      MW._MPI_Comm_rank(comm_worker_world, &myRank2);
      if(myRank != myRank2)
	throw std::runtime_error("Error in service_worker_func: (myRank != myRank2)");    

      while(1) {
	MPI_Status status;
	if (comm_parent == MPI_COMM_NULL)
	  throw std::runtime_error("service_worker: comm_parent == MPI_COMM_NULL !!!");
	// Note: we could just use MPI_Recv here, but that turned out to
	// be consuming a lot of cpu time for some MPI implementations,
	// e.g. openmpi/1.3.4 on isis. To be sure that the
	// service_worker_func thread really does not use any cpu time,
	// we use a loop with sleep() and Iprobe() instead.
	while(1) {
	  int flag;
	  MW._MPI_Iprobe(0, MPI_ANY_TAG, comm_parent, &flag, MPI_STATUS_IGNORE);
	  if(flag)
	    break;
	  int one_milli_second = 1000;
	  usleep(100*one_milli_second);
	}
	MW._MPI_Recv(serviceIDBuf, Service::max_length_serviceID_string, 
		    MPI_UNSIGNED_CHAR, 0, MPI_ANY_TAG, 
		    comm_parent, &status);
	if (status.MPI_TAG == ServiceManager::Tag_worker_shutdown) {
	  if ( !runningServicesParentComms.empty() || 
	       !runningServicesWorkerComms.empty()    )
	    throw std::runtime_error("Services still running, cannot shutdown");
	  {
	    break;
	  }
	}
	else if (status.MPI_TAG == ServiceManager::Tag_start_service) {
	  MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, 
		      ServiceManager::Tag_start_service_ack1, comm_parent);
	  std::string serviceID(serviceIDBuf);
	  ServiceManager::ServiceID_Comm_Map::iterator i_p = 
	    runningServicesParentComms.find( serviceID );
	  if ( i_p != runningServicesParentComms.end() ) 
	    throw std::runtime_error("service_worker - Service already started - p");
	  ServiceManager::ServiceID_Comm_Map::iterator i_w = 
	    runningServicesWorkerComms.find( serviceID );
	  if ( i_w != runningServicesWorkerComms.end() ) 
	    throw std::runtime_error("service_worker - Service already started - w");
      
	  MPI_Comm* comm_parent_copy = new MPI_Comm;
	  MW._MPI_Comm_dup(comm_parent, comm_parent_copy);
	  runningServicesParentComms.insert
	    ( ServiceManager::ServiceID_Comm_Map::value_type(serviceID, comm_parent_copy) );
	  MPI_Comm* comm_worker_world_copy = new MPI_Comm;
	  MW._MPI_Comm_dup(comm_worker_world, comm_worker_world_copy);
	  runningServicesWorkerComms.insert
	    ( ServiceManager::ServiceID_Comm_Map::value_type(serviceID, comm_worker_world_copy) );

	  Service::Worker* service = cht::obj_factory<Service::Worker>::instance().createObject(serviceID);
	  service->start(comm_parent_copy, comm_worker_world_copy);
	  MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, 
		      ServiceManager::Tag_start_service_ack2, comm_parent);
	}
	else if (status.MPI_TAG == ServiceManager::Tag_stop_service) {
	  MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, 
		      ServiceManager::Tag_stop_service_ack1, comm_parent);
	  std::string serviceID(serviceIDBuf);
	  ServiceManager::ServiceID_Comm_Map::iterator i_p = 
	    runningServicesParentComms.find( serviceID );
	  if ( i_p == runningServicesParentComms.end() ) 
	    throw std::runtime_error("Service not started");
	  ServiceManager::ServiceID_Comm_Map::iterator i_w = 
	    runningServicesWorkerComms.find( serviceID );
	  if ( i_w == runningServicesWorkerComms.end() ) 
	    throw std::runtime_error("Service not started");
      
	  Service::Worker* service = cht::obj_factory<Service::Worker>::instance().createObject(serviceID);
	  service->stop();

	  MW._MPI_Comm_free( (*i_w).second );
	  delete (*i_w).second;
	  runningServicesWorkerComms.erase(i_w);
	  MW._MPI_Comm_free( (*i_p).second );
	  delete (*i_p).second;
	  runningServicesParentComms.erase(i_p);
	  MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, 
		      ServiceManager::Tag_stop_service_ack2, comm_parent);
	}
	else if (status.MPI_TAG == ServiceManager::Tag_reset_statistics) {
	  MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, 
		      ServiceManager::Tag_reset_statistics_ack, comm_parent);
	  std::string serviceID(serviceIDBuf);
	  if ( runningServicesParentComms.find( serviceID ) == runningServicesParentComms.end() )
	    throw std::runtime_error("Service not active");
	  if ( runningServicesWorkerComms.find( serviceID ) == runningServicesWorkerComms.end() ) 
	    throw std::runtime_error("Service not active");
	  Service::Worker* service = cht::obj_factory<Service::Worker>::instance().createObject(serviceID);
	  service->resetStatistics();
	}
	else if (status.MPI_TAG == ServiceManager::Tag_report_statistics) {
	  MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, 
		      ServiceManager::Tag_report_statistics_ack, comm_parent);
	  std::string serviceID(serviceIDBuf);
	  if ( runningServicesParentComms.find( serviceID ) == runningServicesParentComms.end() )
	    throw std::runtime_error("Service not active");
	  if ( runningServicesWorkerComms.find( serviceID ) == runningServicesWorkerComms.end() ) 
	    throw std::runtime_error("Service not active");
	  Service::Worker* service = cht::obj_factory<Service::Worker>::instance().createObject(serviceID);
	  std::string messageHeader;
	  messageHeader = "=====  " + serviceID + "::reportStatistics()  =====";
	  service->reportStatistics(messageHeader);
	}
	else {
	  throw std::runtime_error("Unknown tag from parent");
	}
    
      } // end while

    } catch ( std::runtime_error e ) {
      std::cerr << "======= ERROR ERROR ERROR ===============" << std::endl;
      std::cerr << "Error! Exception std::runtime_error caught in service_worker_func!" << std::endl;
      std::cerr << "what(): " << e.what() << std::endl;
      std::cerr << "=========================================" << std::endl;
    } catch ( std::exception e ) {
      std::cerr << "======= ERROR ERROR ERROR ===============" << std::endl;
      std::cerr << "Error! Exception std::exception caught in service_worker_func!" << std::endl;
      std::cerr << "what(): " << e.what() << std::endl;
      std::cerr << "=========================================" << std::endl;
    } catch ( ... ) {
      std::cerr << "======= ERROR ERROR ERROR ===============" << std::endl;
      std::cerr << "Error! Exception caught in service_worker_func!" << std::endl;
      std::cerr << "=========================================" << std::endl;
    }

  }

}; // end namespace
