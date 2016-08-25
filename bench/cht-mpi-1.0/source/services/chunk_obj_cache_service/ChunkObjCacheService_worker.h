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
#ifndef CHUNK_OBJ_CACHE_WORKER_HEADER
#define CHUNK_OBJ_CACHE_WORKER_HEADER
#include <cassert>
#include <iomanip>
#include <memory>
#include <set>
#include "utilities/Singleton.h"
#include "services/chunk_obj_service/ChunkObject.h"
#include "utilities/cht_shared_ptr.h"
#include "services/chunk_obj_cache_service/ChunkObjCacheService.h"
#include "services/Service_parent.h"

namespace cht {
  namespace ChunkObjCacheService {

    template<typename T_CO_Service>
      class Worker : 
    public Base<Service::Worker, T_CO_Service>,
      public Singleton<Worker<T_CO_Service> > {
	friend class Singleton<Worker<T_CO_Service> >;
    public:
    protected:
	typedef typename ChunkObjCacheService::Base<Service::Worker, T_CO_Service>::MPI_Tags MPI_Tags;
	typedef Service::Worker::AccessKey AccessKey; 
	// Stuff to make this a service
	virtual void start_derived();
	virtual void stop_derived();
	// End, service stuff
    public:
  
	// Parent uses default functions (empty ones)
	void resetStatistics();
	void reportStatistics(std::string messageHeader);

      };


    template<typename T_CO_Service> 
      void Worker<T_CO_Service>::start_derived() {
      AccessKey key(this);
      // Receive memoryUsageDanglingLimit from parent
      this->MW._MPI_Recv(&this->memoryUsageDanglingLimit, sizeof(size_t), MPI_UNSIGNED_CHAR, 0, 
			Worker<T_CO_Service>::Tag_params, 
			*key.comm_to_parent(), MPI_STATUS_IGNORE);
      // Receive mode from parent
      this->MW._MPI_Recv(&this->mode, sizeof(extras::Cache::Mode), MPI_UNSIGNED_CHAR, 
			0, Worker<T_CO_Service>::Tag_params, 
			*key.comm_to_parent(), 
			MPI_STATUS_IGNORE);
    }

    template<typename T_CO_Service> 
      void Worker<T_CO_Service>::stop_derived() {
      // Empty cache here?
    }


    template<typename T_CO_Service>
      void Worker<T_CO_Service>::resetStatistics() {
      this->LockMutex();
      this->setOfChunkTypes.clear();
      this->map_stat_being_fetched.clear();
      this->map_stat_used.clear();
      this->map_stat_dangling.clear();
      this->map_stat_not_locally.clear();
      this->number_of_chunks_found_in_being_fetched = 0;
      this->number_of_chunks_found_in_used  = 0;
      this->number_of_chunks_found_in_dangling = 0;
      this->number_of_chunks_not_found_locally = 0;
      this->memoryUsageDanglingMaxUsed = this->totalMemoryUsageDangling;
      this->UnlockMutex();
    }

    template<typename T_CO_Service>
      void Worker<T_CO_Service>::reportStatistics(std::string messageHeader) {
      // Output some statistics.
      std::stringstream s;
      s << messageHeader << "\n";
      this->LockMutex();
      size_t cache_hits = 
	this->number_of_chunks_found_in_used + 
	this->number_of_chunks_found_in_dangling;
      size_t total_number_of_accesses = cache_hits + 
	this->number_of_chunks_found_in_being_fetched + 
	this->number_of_chunks_not_found_locally;
      s<< "=== Cache hit statistics, all chunk types ===\n";
      s << "  Cache hit rate ( used + dangling )     : " 
	<< std::fixed << std::setw(6) << std::setprecision(2)    
	<< 100*(cache_hits / (double)total_number_of_accesses)
	<< "% \n";
      s << "  No. of chunks found in map_used        :  " 
	<< this->number_of_chunks_found_in_used << "\n";
      s << "  No. of chunks found in map_dangling    :  " 
	<< this->number_of_chunks_found_in_dangling << "\n";
      s << "  No. of chunks being_fetched            :  " 
	<< this->number_of_chunks_found_in_being_fetched << "\n";
      s << "  No. of chunks not found locally        :  " 
	<< this->number_of_chunks_not_found_locally << "\n";
      s << "--------------\n";
  
      // Output stats per chunk type:
      s<< "=== Cache hit statistics per chunk type ===\n";
      size_t n_used, n_dangling, n_being_fetched, n_not_locally;
      std::set<std::string>::iterator it = this->setOfChunkTypes.begin();
      for (;it != this->setOfChunkTypes.end();++it) {
	n_used          = this->map_stat_used[*it];
	n_dangling      = this->map_stat_dangling[*it];
	n_being_fetched = this->map_stat_being_fetched[*it];
	n_not_locally   = this->map_stat_not_locally[*it];
	cache_hits      = n_used + n_dangling;
	total_number_of_accesses = cache_hits + n_being_fetched + n_not_locally;
	s << *it << ":\n";
	s << "  Cache hit rate ( used + dangling )     : " 
	  << std::fixed << std::setw(6) << std::setprecision(2)    
	  << 100*(cache_hits / (double)total_number_of_accesses)
	  << "% \n";
	s << "  No. of chunks found in map_used        :  " << n_used << "\n";
	s << "  No. of chunks found in map_dangling    :  " << n_dangling << "\n";
	s << "  No. of chunks being_fetched            :  " << n_being_fetched << "\n";
	s << "  No. of chunks not found locally        :  " << n_not_locally << "\n";
	s << "--------------\n";
      }
      // Output max memory usage for dangling objects
      s << "=== Cache memory usage statistics (dangling) ===\n";
      s 
	<< " max   " 
	<< std::setw(8) << std::setprecision(6)    
	<< this->memoryUsageDanglingMaxUsed/(double)1000000 << " MBytes\n"
	<< " limit " 
	<< std::setw(8) << std::setprecision(6)    
	<< this->memoryUsageDanglingLimit/(double)1000000 << " MBytes\n";
      s << "--------------\n";  
      OutputService::Worker::instance().outputInfo(s.str());
      this->UnlockMutex();
    }

  }; // end namespace
}; // end namespace

#endif
