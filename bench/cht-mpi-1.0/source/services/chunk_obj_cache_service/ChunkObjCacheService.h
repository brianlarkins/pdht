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
#ifndef CHUNK_OBJ_CACHE_HEADER
#define CHUNK_OBJ_CACHE_HEADER

#include <iostream>
#include <cassert>
#include <iomanip>
#include <memory>
#include <set>
#include "utilities/Singleton.h"
#include "services/chunk_obj_service/ChunkObject.h"
#include "utilities/cht_shared_ptr.h"
#include "utilities/cht_threads.h"

namespace cht {
  /** A service for caching of Chunk objects. Requires a service for
      management of Chunk object (template argument). One example of
      such a service is ChunkObjService. */
  namespace ChunkObjCacheService {
 
    // The template parameters are as follows:
    // T_CO_Service : ChunkObjCache is a layer on top of the Chunk object 
    //                service T_CO_Service
    // Service_type : Determines whether this is the parent or worker part 
    //                of the service
    template<typename Service_type, typename T_CO_Service>
      class Base : public Service_type {
    public:
      static std::string object_type_id() {
	return "ChunkObjCacheService<" + T_CO_Service::object_type_id() + ">";
      }


    protected:
      enum MPI_Tags {
	Tag_params
      };
      // End, service stuff
    public:


    Base()
      : mode(extras::Cache::Enabled), accessCounter(0), totalMemoryUsageDangling(0), memoryUsageDanglingLimit(0), memoryUsageDanglingMaxUsed(0), ptrIdCounter(0), number_of_chunks_found_in_being_fetched(0), number_of_chunks_found_in_used(0), number_of_chunks_found_in_dangling(0), number_of_chunks_not_found_locally(0) {
      };

      void deleteChunk(ChunkID id);
  
      bool getChunkIfExists(ChunkID id, cht::shared_ptr<Chunk const> & objPtr);
      template<typename Tobj>
	bool getChunkIfExists(ChunkID id, cht::shared_ptr<Tobj const> & objPtr);
      void getChunk(ChunkID id, 
		    cht::shared_ptr<Chunk const> & objPtr);
      template<typename Tobj>
	void getChunk(ChunkID id, cht::shared_ptr<Tobj const> & objPtr);
      template<typename Tobj>
	void getChunks(ChunkID cid1, cht::shared_ptr<Tobj const> & objPtr1,
		       ChunkID cid2, cht::shared_ptr<Tobj const> & objPtr2);

      bool getChunkIfInCache(ChunkID cid, 
			     cht::shared_ptr<Chunk const> & objPtr);
      bool getChunkIfLocal(ChunkID id, 
			   cht::shared_ptr<Chunk const> & objPtr);
      template<typename Tobj>
	bool getChunkIfLocal(ChunkID id, cht::shared_ptr<Tobj const> & objPtr);

      void setCacheMode(extras::Cache::Mode newMode);
      void setCacheSize(size_t memoryUsageLimit_);

      ChunkID getIdForRegisterChunk(Chunk const * objPtr,
				  std::string class_id_str);
      void startRegisterChunk(Chunk const * objPtr, ChunkID cid);
      bool isRegisterChunkStillInProgress(ChunkID cid);
      ChunkID registerChunkDirectly(Chunk const * objPtr,
				  std::string class_id_str);

      ChunkID getIdForCopyChunk(ChunkID cid_old);
      void startCopyChunk(ChunkID cid_old, ChunkID cid_new);
      bool isCopyChunkStillInProgress(ChunkID cid);

    protected:
      struct CachedChunk {
	ChunkID const id;
	cht::shared_ptr<Chunk const> objPtr;
	std::string const class_id_str;
	size_t latestAccess;
	Threads::Mutex mutex;    
	void LockMutex() {
	  mutex.lock();
	}
	void UnlockMutex() {
	  mutex.unlock();
	}
      CachedChunk(ChunkID id, 
		  std::string class_id_str,
		  size_t latestAccess) 
        : id(id), objPtr(0), class_id_str(class_id_str), 
	  latestAccess(latestAccess) { 
	}
      };
  
      // Mutex must not be locked when calling this function
      void checkMapSizes(int i) {
	LockMutex();
	assert( map_dangling.size() == map_dangling_latestAccess.size() );
	if(map_used.size() + map_dangling.size() != map_all_cached.size())
	  std::cerr << "   map_used.size() = " << map_used.size() 
		    << "   map_dangling.size() = " << map_dangling.size() 
		    << "   map_all_cached.size() = " << map_all_cached.size()
		    << "   i = " << i 
		    << std::endl;
	assert( map_used.size() + map_dangling.size() == map_all_cached.size() );
	UnlockMutex();
      }

      extras::Cache::Mode mode;

      typedef std::map<ChunkID, CachedChunk*> AccessMap;
      typedef std::map<size_t, CachedChunk*> LatestAccessMap;
      typedef std::map<size_t, CachedChunk*> PtrMap;
      AccessMap  map_used;   /* Used objects are not accounted for in
				total memory usage for cache (they do not
				"belong" to the cache). */
      AccessMap  map_dangling;
      LatestAccessMap map_dangling_latestAccess;
      PtrMap map_all_cached;
      size_t accessCounter;
      size_t totalMemoryUsageDangling;
      size_t memoryUsageDanglingLimit;
      size_t memoryUsageDanglingMaxUsed;
      Threads::Mutex mutex;
      size_t ptrIdCounter;

      std::set<ChunkID> chunksBeingFetched;
      Threads::Cond chunksBeingFetched_cond;

      typedef std::map<std::string, size_t> StatisticsMap;
      StatisticsMap map_stat_being_fetched;
      StatisticsMap map_stat_used;
      StatisticsMap map_stat_dangling;
      StatisticsMap map_stat_not_locally;
      std::set<std::string> setOfChunkTypes;

      size_t number_of_chunks_found_in_being_fetched;
      size_t number_of_chunks_found_in_used;
      size_t number_of_chunks_found_in_dangling;
      size_t number_of_chunks_not_found_locally;

      void LockMutex();
      void UnlockMutex();
      void freeMemoryIfNeeded(); 
      template<typename Tobj>
	void assignObjPtrViaUglyCast(CachedChunk const * const theCachedChunk,
				     cht::shared_ptr<Tobj const> & objPtr) const;

      static void onRefCountChangeCallback( void* contextPtr, size_t ptrId, int ref_count, int previous_ref_count);
  
    }; // end class Base


    template<typename Service_type, typename T_CO_Service>
      void Base<Service_type, T_CO_Service>::LockMutex() {
      mutex.lock();
    }

    template<typename Service_type, typename T_CO_Service>
      void Base<Service_type, T_CO_Service>::UnlockMutex() {
      mutex.unlock();
    }



    // Mutex must not be locked when calling this function
    template<typename Service_type, typename T_CO_Service>
      void Base<Service_type, T_CO_Service>::setCacheMode(extras::Cache::Mode newMode) {
      if ( Service_type::serviceIsRunning() )
	throw std::runtime_error("Call to cht::ChunkObjCacheService::Base<...>::setCacheMode while service is running.");
      LockMutex();
      mode = newMode;
      UnlockMutex();
    }
    // Mutex must not be locked when calling this function
    template<typename Service_type, typename T_CO_Service>
      void Base<Service_type, T_CO_Service>::setCacheSize(size_t memoryUsageLimit_) {
      if ( Service_type::serviceIsRunning() )
	throw std::runtime_error("Call to cht::ChunkObjCacheService::Base<...>::setCacheSize while service is running.");
      LockMutex();
      memoryUsageDanglingLimit = memoryUsageLimit_;
      UnlockMutex();
    }

    // Mutex must not be locked when calling this function
    template<typename Service_type, typename T_CO_Service>
      void Base<Service_type, T_CO_Service>::freeMemoryIfNeeded() {
      checkMapSizes(1);
      typename LatestAccessMap::iterator it;
      LockMutex();
      while( totalMemoryUsageDangling > memoryUsageDanglingLimit ) {
	if ( map_dangling_latestAccess.empty() )
	  throw std::runtime_error("freeMemoryIfNeeded(): map_dangling_latestAccess "
				   "is empty but more memory needs to be freed!!");
	// Remove oldest unused object from cache. 
	// (One could if there are no unused objects remove one of the used objects? )
	it = map_dangling_latestAccess.begin(); // begin() here gives oldest object.
	CachedChunk * theCachedChunk = it->second;
	theCachedChunk->LockMutex();
	if(theCachedChunk->objPtr == 0)
	  throw std::runtime_error("Error in freeMemoryIfNeeded(): (theCachedChunk->objPtr == 0).");
	totalMemoryUsageDangling -= theCachedChunk->objPtr->memoryUsage();
	if(map_all_cached.erase( theCachedChunk->objPtr.get_callback_id() ) != 1)
	  throw std::runtime_error("Error: failed to erase from map_all_cached.");
	if(map_dangling.erase(theCachedChunk->id) != 1)
	  throw std::runtime_error("Error: failed to erase from map_dangling.");
	map_dangling_latestAccess.erase(it);
	UnlockMutex();
	// We do not unlock theCachedChunk->mutex here
	delete theCachedChunk; // Note that ref count changes for objPtr in it->second object
	LockMutex();
      }
      UnlockMutex();
      checkMapSizes(2);
    }

    // Mutex must not be locked when calling this function
    // theCachedChunk->mutex must be locked
    template<typename Service_type, typename T_CO_Service>
      template<typename Tobj>
      void Base<Service_type, T_CO_Service>::assignObjPtrViaUglyCast(CachedChunk const * const theCachedChunk,
								     cht::shared_ptr<Tobj const> & objPtr) const {
      // Check that correct mutexes are locked or not here in debug mode
      std::string class_id_str = Tobj::get_class_id();
      if (theCachedChunk->class_id_str != class_id_str)
	throw std::runtime_error("Type given to assignObjPtrViaUglyCast does not "
				 "match type stored in cache.");
      // Ugly cast needed here!!
      objPtr = cht::shared_ptr<Tobj const>( (Tobj*)&(*(theCachedChunk->objPtr)), 
					    theCachedChunk->objPtr.getRefCountPtr() );
    }

    template<typename Service_type, typename T_CO_Service>
      void Base<Service_type, T_CO_Service>::deleteChunk(ChunkID id) {
      if (mode == extras::Cache::Disabled) 
	return T_CO_Service::instance().deleteChunk(id);
      // If the object is cache resident we remove it
      LockMutex();
      typename AccessMap::iterator it = map_used.find( id );
      if ( it != map_used.end() ) {
	CachedChunk* theCachedChunk = it->second;
	theCachedChunk->LockMutex();
	map_used.erase(it);
	if(map_all_cached.erase( theCachedChunk->objPtr.get_callback_id() ) != 1)
	  throw std::runtime_error("Error: failed to erase from map_all_cached.");
	// Now the cached chunk has been removed from the two lists it was in 
	UnlockMutex();
	// We do not unlock theCachedChunk->mutex here
	delete theCachedChunk; // Note that ref count changes for objPtr in it->second object
      }
      else {
	it = map_dangling.find( id );  
	if ( it != map_dangling.end() ) {
	  CachedChunk* theCachedChunk = it->second;
	  theCachedChunk->LockMutex();
	  map_dangling.erase(it);
	  if(map_dangling_latestAccess.erase(theCachedChunk->latestAccess) != 1)
	    throw std::runtime_error("Error: failed to erase from map_dangling_latestAccess.");
	  if(map_all_cached.erase( theCachedChunk->objPtr.get_callback_id() ) != 1)
	    throw std::runtime_error("Error: failed to erase from map_all_cached.");
	  // Now the cached chunk has been removed from the three lists it was in 
	  UnlockMutex();
	  // We do not unlock theCachedChunk->mutex here
	  delete theCachedChunk; // Note that ref count changes for objPtr in it->second object
	}
	else {
	  UnlockMutex();
	}
      }
      T_CO_Service::instance().deleteChunk(id);
    }



    template<typename Service_type, typename T_CO_Service>
      ChunkID Base<Service_type, T_CO_Service>::getIdForRegisterChunk(Chunk const * objPtr,
								     std::string class_id_str) {
      return T_CO_Service::instance().getIdForRegisterChunk( objPtr, class_id_str );
    }

    template<typename Service_type, typename T_CO_Service>
      void Base<Service_type, T_CO_Service>::startRegisterChunk(Chunk const * objPtr, ChunkID cid) {
      if (mode == extras::Cache::Disabled) 
	return T_CO_Service::instance().startRegisterChunk(objPtr, cid);
      // Check now if the chunk will be created locally, then we do not need to cache it
      if ( T_CO_Service::instance().chunkResidesLocally(cid) )
	return T_CO_Service::instance().startRegisterChunk(objPtr, cid);
      // If we reach this point it means the chunk will be created locally, so we will need to cache it.
      T_CO_Service::instance().startRegisterChunk(objPtr, cid);
      // Now make sure we keep a copy of the chunk in cache.
      checkMapSizes(3);
      CachedChunk* theCachedChunk;
      LockMutex();
      std::string class_id_str = T_CO_Service::instance().getChunkTypeIDStr(cid.chunkTypeID);
      theCachedChunk = new CachedChunk(cid, class_id_str, ++accessCounter);
      UnlockMutex();
      theCachedChunk->objPtr = objPtr;
      LockMutex();
      /* ref count for objPtr == 1 here (this is important for the
	 callBack function).  */      
      // The ref count for theCachedChunk->objPtr must now be 1. 
      assert(theCachedChunk->objPtr.getRefCountPtr()->get_count() == 1);
      // We have the only ref to the object. In this case we treat it as
      // "dangling".
      map_dangling[cid] = theCachedChunk;
      map_dangling_latestAccess[theCachedChunk->latestAccess] = theCachedChunk;
      if(theCachedChunk->objPtr == 0)
	throw std::runtime_error("Error: (theCachedChunk->objPtr == 0)");
      totalMemoryUsageDangling += theCachedChunk->objPtr->memoryUsage();
      memoryUsageDanglingMaxUsed = memoryUsageDanglingMaxUsed > totalMemoryUsageDangling ? memoryUsageDanglingMaxUsed : totalMemoryUsageDangling;
      size_t ptrId = ++ptrIdCounter;
      theCachedChunk->objPtr.setRefCountChangeCallback(onRefCountChangeCallback, this, ptrId);
      map_all_cached[ptrId] = theCachedChunk;
      UnlockMutex();
      // Only now can theCachedChunk be accessed by other threads
      freeMemoryIfNeeded();
      checkMapSizes(4);
    }

    template<typename Service_type, typename T_CO_Service>
      bool Base<Service_type, T_CO_Service>::isRegisterChunkStillInProgress(ChunkID cid) {
      return T_CO_Service::instance().isRegisterChunkStillInProgress(cid);
    }

    template<typename Service_type, typename T_CO_Service>
      ChunkID Base<Service_type, T_CO_Service>::registerChunkDirectly(Chunk const * objPtr,
								    std::string class_id_str) {
      return T_CO_Service::instance().registerChunkDirectly(objPtr, class_id_str);
    }

    template<typename Service_type, typename T_CO_Service>
      ChunkID Base<Service_type, T_CO_Service>::getIdForCopyChunk(ChunkID cid_old) {
      return T_CO_Service::instance().getIdForCopyChunk(cid_old);
    }

    template<typename Service_type, typename T_CO_Service>
      void Base<Service_type, T_CO_Service>::startCopyChunk(ChunkID cid_old, ChunkID cid_new) {
      return T_CO_Service::instance().startCopyChunk(cid_old, cid_new);
    }

    template<typename Service_type, typename T_CO_Service>
      bool Base<Service_type, T_CO_Service>::isCopyChunkStillInProgress(ChunkID cid) {
      return T_CO_Service::instance().isCopyChunkStillInProgress(cid);
    }


    template<typename Service_type, typename T_CO_Service>
      template<typename Tobj>
      void Base<Service_type, T_CO_Service>::getChunks(ChunkID cid1, cht::shared_ptr<Tobj const> & objPtr1,
						       ChunkID cid2, cht::shared_ptr<Tobj const> & objPtr2) {
      getChunk(cid1, objPtr1);
      getChunk(cid2, objPtr2);
    }
    
    template<typename Service_type, typename T_CO_Service>
      bool Base<Service_type, T_CO_Service>::getChunkIfExists(ChunkID id, 
							      cht::shared_ptr<Chunk const> & objPtr) {
      if (mode == extras::Cache::Disabled) 
	return T_CO_Service::instance().getChunkIfExists( id, objPtr );
      // Check now if the chunk resides locally, then it should not be cached and we just call the T_CO_Service
      if ( T_CO_Service::instance().chunkResidesLocally(id) )
	return T_CO_Service::instance().getChunkIfExists( id, objPtr );
      checkMapSizes(5);
      std::string cid_class_id_str = T_CO_Service::instance().getChunkTypeIDStr(id.chunkTypeID);
      if ( getChunkIfLocal( id, objPtr) )
	return true;
      LockMutex();
      if ( getChunkIfInCache( id, objPtr) ) 
	return true;
      // Object is not in cache.
      // Check if some other thread is already getting this chunk.
      if(chunksBeingFetched.find(id) != chunksBeingFetched.end()) {
	// Chunk found in chunksBeingFetched
	number_of_chunks_found_in_being_fetched++;
	map_stat_being_fetched[cid_class_id_str]++;
	// Some other thread is getting the chunk. We just wait.
	while(chunksBeingFetched.find(id) != chunksBeingFetched.end()) {
	  // Some other thread is getting the chunk. We just wait.
	  chunksBeingFetched_cond.wait(mutex);
	}
	UnlockMutex();
	return getChunkIfExists(id, objPtr);
      }
      // Chunk not found in locally
      number_of_chunks_not_found_locally++;
      map_stat_not_locally[cid_class_id_str]++;
      // No other thread is currently getting the chunk. Add id to
      // chunksBeingFetched set to be sure noone else starts to get it while we
      // are getting it.
      chunksBeingFetched.insert(id);
      UnlockMutex();
      // Note that now that the mutex is unlocked, another thread may
      // modify the maps. For example the chunk we are about to get here may be
      // fetched by another thread and added to map_used.
      if(T_CO_Service::instance().getChunkIfExists(id, objPtr) == false) {
	// Chunk does not exist.
	LockMutex();
	chunksBeingFetched.erase(id);
	chunksBeingFetched_cond.broadcast();
	UnlockMutex();
	return false;
      }
      LockMutex();
      if(map_used.find(id) != map_used.end() || map_dangling.find(id) != map_dangling.end())
	throw std::runtime_error("Error: Some other thread has already added this chunk. This should not happen.");
      chunksBeingFetched.erase(id);
      chunksBeingFetched_cond.broadcast();
      CachedChunk* theCachedChunk = new CachedChunk(id, cid_class_id_str, ++accessCounter);
      if(map_used.find(id) != map_used.end())
	throw std::runtime_error("Error: id already found in map_used (3).");
      if(map_dangling.find(id) != map_dangling.end())
	throw std::runtime_error("Error: id to be added to map_used found in map_dangling (3).");
      map_used[id] = theCachedChunk;
      size_t ptrId = ++ptrIdCounter;
      objPtr.setRefCountChangeCallback(onRefCountChangeCallback, this, ptrId);
      map_all_cached[ptrId] = theCachedChunk;
      theCachedChunk->LockMutex();
      UnlockMutex();
      // Note that theCachedChunk->mutex is locked here
      theCachedChunk->objPtr = objPtr;
      theCachedChunk->UnlockMutex();
      checkMapSizes(9);
      return true;
    }

    template<typename Service_type, typename T_CO_Service>
      template<typename Tobj>
      bool Base<Service_type, T_CO_Service>::getChunkIfExists(ChunkID id, 
								 cht::shared_ptr<Tobj const> & objPtr) {
      cht::shared_ptr<Chunk const> objCOPtr;
      if(getChunkIfExists(id, objCOPtr) == false)
	return false;
      // Ugly cast needed here!!
      objPtr = cht::shared_ptr<Tobj const>( (Tobj*)&(*objCOPtr), objCOPtr.getRefCountPtr() );
      return true;
    }

    template<typename Service_type, typename T_CO_Service>
      void Base<Service_type, T_CO_Service>::getChunk(ChunkID id, 
						      cht::shared_ptr<Chunk const> & objPtr) {
      bool success = getChunkIfExists(id, objPtr);
      if(success == false)
	throw std::runtime_error("Error in ChunkObjCacheService getChunk: chunk does not seem to exist.");
    }

    template<typename Service_type, typename T_CO_Service>
      template<typename Tobj>
      void Base<Service_type, T_CO_Service>::getChunk(ChunkID id, 
						      cht::shared_ptr<Tobj const> & objPtr) {
      cht::shared_ptr<Chunk const> objCOPtr;
      getChunk(id, objCOPtr);
      // Ugly cast needed here!!
      objPtr = cht::shared_ptr<Tobj const>( (Tobj*)&(*objCOPtr), objCOPtr.getRefCountPtr() );
    }

    template<typename Service_type, typename T_CO_Service>
      bool Base<Service_type, T_CO_Service>::getChunkIfInCache(ChunkID id, 
							       cht::shared_ptr<Chunk const> & objPtr) {

      // NOTE, OBS OBS OBS! : The mutex should be locked when calling
      // this function and may or may not be locked at return. If the
      // mutex is unlocked, the function returns false, otherwise it returns
      // true.
      std::string cid_class_id_str = T_CO_Service::instance().getChunkTypeIDStr(id.chunkTypeID);
      // Add to chunk type set
      setOfChunkTypes.insert( cid_class_id_str );
      typename AccessMap::iterator it = map_used.find(id);
      if ( it != map_used.end() ) {
	// Chunk found in map_used
	number_of_chunks_found_in_used++;
	map_stat_used[cid_class_id_str]++;
	CachedChunk* theCachedChunk = it->second;
	theCachedChunk->LockMutex();
	theCachedChunk->latestAccess = ++accessCounter;    
	UnlockMutex();
	// Note that theCachedChunk->mutex is locked here
	objPtr = theCachedChunk->objPtr;
	theCachedChunk->UnlockMutex();    
	checkMapSizes(6);
	return true;
      }
      it = map_dangling.find(id);
      if ( it != map_dangling.end() ) {
	// Chunk found in map_dangling
	number_of_chunks_found_in_dangling++;
	map_stat_dangling[cid_class_id_str]++;
	CachedChunk* theCachedChunk = it->second;
	theCachedChunk->LockMutex();
	if(map_dangling_latestAccess.erase(theCachedChunk->latestAccess) != 1)
	  throw std::runtime_error("Error: failed to erase from map_dangling_latestAccess.");
	map_dangling.erase(it);
	theCachedChunk->latestAccess = ++accessCounter;
	if(map_used.find(id) != map_used.end())
	  throw std::runtime_error("Error: id already found in map_used (2).");
	if(map_dangling.find(id) != map_dangling.end())
	  throw std::runtime_error("Error: id to be added to map_used found in map_dangling (2).");
	map_used[id] = theCachedChunk;
	if(theCachedChunk->objPtr == 0)
	  throw std::runtime_error("Error: (theCachedChunk->objPtr == 0)");
	totalMemoryUsageDangling -= theCachedChunk->objPtr->memoryUsage();
	UnlockMutex();
	// Note that theCachedChunk->mutex is locked here
	objPtr = theCachedChunk->objPtr;
	theCachedChunk->UnlockMutex();
	checkMapSizes(7);
	return true;    
      }
      // Object is not in cache.
      // NOT unlocking mutex!!!
      return false;
    }    

    template<typename Service_type, typename T_CO_Service>
      bool Base<Service_type, T_CO_Service>::getChunkIfLocal(ChunkID cid, 
							     cht::shared_ptr<Chunk const> & objPtr) {
      if (mode == extras::Cache::Disabled)
	return T_CO_Service::instance().getChunkIfLocal( cid, objPtr );
      if ( T_CO_Service::instance().getChunkIfLocal( cid, objPtr ) )
	return true;
      // Check if chunk is in cache, then we can return it directly.
      LockMutex();
      if ( getChunkIfInCache( cid, objPtr) ) 
	return true;
      UnlockMutex();
      return false;
    }

    template<typename Service_type, typename T_CO_Service>
      template<typename Tobj>
      bool Base<Service_type, T_CO_Service>::getChunkIfLocal(ChunkID id, 
							     cht::shared_ptr<Tobj const> & objPtr) {
      cht::shared_ptr<Chunk const> objCOPtr;
      if(getChunkIfLocal(id, objCOPtr) == false)
	return false;
      // Ugly cast needed here!!
      objPtr = cht::shared_ptr<Tobj const>( (Tobj*)&(*objCOPtr), objCOPtr.getRefCountPtr() );
      return true;
    }

    // Mutex must not be locked when calling this function
    template<typename Service_type, typename T_CO_Service>
      void Base<Service_type, T_CO_Service>::onRefCountChangeCallback( void* contextPtr, 
										size_t ptrId, 
										int ref_count, 
										int previous_ref_count) {
      if ( ref_count == 0 && previous_ref_count != 1)
	throw std::runtime_error("ref_count == 0 && previous_ref_count != 1 in onRefCountChangeCallback");
      if ( ref_count != 1 )
	return;
      if ( previous_ref_count != 2)
	throw std::runtime_error("previous_ref_count != 2 when ref_count == 1 in onRefCountChangeCallback");
      // ref_count == 1
      Base<Service_type, T_CO_Service>* COC_ptr = (Base<Service_type, T_CO_Service>*)contextPtr;
      Base<Service_type, T_CO_Service> & COC = *COC_ptr;
      //Base<T_CO_Service> & COC = Base<T_CO_Service>::instance();
      // Check if object is cache resident
      typename PtrMap::iterator it;
      COC.LockMutex();
      it = COC.map_all_cached.find(ptrId);
      if ( it == COC.map_all_cached.end() ) {
	COC.UnlockMutex();
	return;
      }
      // OK, object is in cache
      CachedChunk* theCachedChunk = it->second;
      theCachedChunk->LockMutex();
      size_t previousAccess = theCachedChunk->latestAccess;
      theCachedChunk->latestAccess = ++COC.accessCounter;
      typename AccessMap::iterator it_used = COC.map_used.find(theCachedChunk->id);
      if ( it_used == COC.map_used.end() ) {
	// We accept that the chunk may be in map_dangling.
	if(COC.map_dangling.find(theCachedChunk->id) == COC.map_dangling.end())
	  throw std::runtime_error("onRefCountCallback: CachedChunk not in map_used or map_dangling!");
	// Now we know the chunk is in map_dangling. We must remove and
	// re-add it to map_dangling_latestAccess since latestAccess value has
	// changed.
	if(COC.map_dangling_latestAccess.erase(previousAccess) != 1)
	  throw std::runtime_error("onRefCountCallback: failed to erase from map_dangling_latestAccess.");
	COC.map_dangling_latestAccess[theCachedChunk->latestAccess] = theCachedChunk;
	theCachedChunk->UnlockMutex();
	COC.UnlockMutex();
	return;
      }
      COC.map_used.erase(it_used);
      if(COC.map_dangling.find(theCachedChunk->id) != COC.map_dangling.end())
	throw std::runtime_error("Error: theCachedChunk already in map_dangling (1).");
      COC.map_dangling[theCachedChunk->id] = theCachedChunk;
      COC.map_dangling_latestAccess[theCachedChunk->latestAccess] = theCachedChunk;
      if(theCachedChunk->objPtr == 0)
	throw std::runtime_error("Error: (theCachedChunk->objPtr == 0)");
      COC.totalMemoryUsageDangling += theCachedChunk->objPtr->memoryUsage();
      COC.memoryUsageDanglingMaxUsed = 
	COC.memoryUsageDanglingMaxUsed > COC.totalMemoryUsageDangling ? COC.memoryUsageDanglingMaxUsed : COC.totalMemoryUsageDangling;
      theCachedChunk->UnlockMutex();
      COC.UnlockMutex();
      COC.freeMemoryIfNeeded();
    }

  }; // end namespace
}; // end namespace

#endif
