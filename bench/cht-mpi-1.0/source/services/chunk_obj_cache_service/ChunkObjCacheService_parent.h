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
#ifndef CHUNK_OBJ_CACHE_PARENT_HEADER
#define CHUNK_OBJ_CACHE_PARENT_HEADER

#include <cassert>
#include <iomanip>
#include <memory>
#include <set>
#include "utilities/Singleton.h"
#include "utilities/cht_shared_ptr.h"
#include "services/Service_parent.h"
#include "services/chunk_obj_cache_service/ChunkObjCacheService.h"
#include "services/chunk_obj_service/ChunkObject.h"

namespace cht {
  namespace ChunkObjCacheService {

    template<typename T_CO_Service>
      class Parent : 
    public Base<Service::Parent, T_CO_Service>,
      public Singleton<Parent<T_CO_Service> > {
	friend class Singleton<Parent<T_CO_Service> >;
    public:

    protected:
	typedef typename Base<Service::Parent, T_CO_Service>::MPI_Tags MPI_Tags;
	// Stuff to make this a service
	virtual void start_derived();
	virtual void stop_derived();
	// End, service stuff
    public:
  

	//  static void onRefCountChangeCallback( size_t ptrId, int ref_count, int previous_ref_count);
  
      };

    template<typename T_CO_Service> 
      void Parent<T_CO_Service>::start_derived() {
      Service::Parent::AccessKey key(this);
      // Send memoryUsageDanglingLimit to workers
      for(int i = 0; i < key.n_workers(); i++) {
	size_t valueToSend = this->memoryUsageDanglingLimit;
	this->MW._MPI_Send(&valueToSend, sizeof(size_t), 
			  MPI_UNSIGNED_CHAR, i, 
			  Parent<T_CO_Service>::Tag_params, *key.comm_to_workers());
	this->MW._MPI_Send(&this->mode, sizeof(extras::Cache::Mode), 
			  MPI_UNSIGNED_CHAR, i, 
			  Parent<T_CO_Service>::Tag_params, *key.comm_to_workers());
      }
    }


    template<typename T_CO_Service> 
      void Parent<T_CO_Service>::stop_derived() {
      // Empty cache here?
    }

  }; // end namespace
}; // end namespace

#endif
