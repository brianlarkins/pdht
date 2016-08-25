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
#ifndef CHUNK_OBJ_SERVICE_HEADER
#define CHUNK_OBJ_SERVICE_HEADER

#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include "utilities/cht_utils.h"
#include "utilities/cht_threads.h"
#include "utilities/Singleton.h"
#include "services/services_utils.h"
#include "services/MPIWrapperInclude.h"
#include "services/chunk_obj_service/ChunkObject.h"
#include "services/output_service/OutputService.h"

namespace cht {
  /** A service for management and distribution of Chunk objects. */
  namespace ChunkObjService {    

    template<typename Service_type>
      class Base : public Service_type {
    public:
      static std::string object_type_id() {
	return "ChunkObjService";
      }
      std::string getChunkTypeIDStr(int const ChunkTypeIDInt) const;
      int getChunkTypeIDInt(std::string const ChunkTypeIDStr) const;
    private:
      typedef std::map<std::string, int> strToIntMap;
      typedef std::map<int, std::string> intToStrMap;
      strToIntMap chunkTypeIDToIntIDMap;
      intToStrMap intIDToChunkTypeIDMap;
    protected:
      void populateChunkIDMaps(std::list<std::string> const & strList) {
	cht::Service::populateMapsForGivenStrList(strList,
						  chunkTypeIDToIntIDMap,
						  intIDToChunkTypeIDMap);
      }
      void checkStrListAgainstIDMap(std::list<std::string> const & localStrList) {
	cht::Service::checkStrListAgainstMap(localStrList, chunkTypeIDToIntIDMap);
      }
      /* This constant determines the size of "small" writeChunkObj
	 messages; if the chunk data fits within a message of this size
	 then the data is included directly, otherwise two messages are
	 used.  */
      static const int MAX_WRITECHUNKOBJ_SMALL_MESSAGE_SIZE = 1000;
      // FIXME: this constant should probably be defined and checked elsewhere.
      static const int max_class_id_str_sz = 888;
      void writeChunkObj_communicate(cht::vector<char> & buf,
				     int bufSzExceptTag,
				     int destRank,
				     MPI_Comm* comm);
      void initTags();
      class TagManager {
      protected:
	Threads::Mutex mutex; 
	std::map<std::string, int> permanentTags;
	std::list<int> temporaryTags;
      public:
	~TagManager() {
	  if ( !temporaryTags.empty() )
	    throw std::runtime_error("Destruction of TagManager with "
				     "unreleased temporary tags");
	}
	int getPermanentTag(std::string tag_name) {
	  mutex.lock();
	  if (permanentTags.find(tag_name) == permanentTags.end())
	    throw std::runtime_error("getPermanentTag: requested tag not "
				     "registered, tag_name = '" + tag_name + "'.");
	  int tag = permanentTags[tag_name];
	  mutex.unlock();
	  return tag;
	}
	void registerPermanentTag(std::string tag_name, int tag) {
	  mutex.lock();
	  permanentTags[tag_name] = tag;
	  mutex.unlock();
	}
	int getTemporaryTag(bool odd) {
	  int tag = 0;
	  mutex.lock();
	  std::map<std::string, int>::iterator it = permanentTags.begin();
	  while ( it != permanentTags.end() ) {
	    tag = tag > (*it).second ? tag : (*it).second;
	    ++it;
	  }
	  tag++;
	  // Make sure tag is even or odd as requested.
	  if(odd) {
	    if(tag%2 == 0)
	      tag++; // Now tag is odd.
	  } else {
	    if(tag%2 == 1)
	      tag++; // Now tag is even.
	  }
	  while ( std::find( temporaryTags.begin(), temporaryTags.end(), tag ) != 
		  temporaryTags.end() )
	    tag += 2;
	  if(odd)
	    assert(tag%2 == 1);
	  else
	    assert(tag%2 == 0);
	  temporaryTags.push_back(tag);
	  mutex.unlock();
	  return tag;
	}
	/* We need to have two distinct sets of temporary tags, for incoming and outgoing messages. */
	int getTemporaryTag_incoming() { return getTemporaryTag(true); }
	int getTemporaryTag_outgoing() { return getTemporaryTag(false); }
	void releaseTemporaryTag(int tag) {
	  mutex.lock();
	  std::list<int>::iterator it = std::find( temporaryTags.begin(), 
						   temporaryTags.end(), tag );
	  if (it == temporaryTags.end())
	    throw std::runtime_error("Attempt to release tag that is not "
				     "managed by TagManager");
	  temporaryTags.erase(it);
	  mutex.unlock();
	}
      };
      TagManager tagManager;
    }; // end class Base

    template<typename Service_type>
      std::string Base<Service_type>::
      getChunkTypeIDStr(int const ChunkTypeIDInt) const {
      if(intIDToChunkTypeIDMap.find(ChunkTypeIDInt) == intIDToChunkTypeIDMap.end()) {
	std::stringstream ss;
	ss << "In getChunkTypeIDStr(): int " << ChunkTypeIDInt << " not found in map";
	throw std::runtime_error(ss.str());
      }
      return intIDToChunkTypeIDMap.find(ChunkTypeIDInt)->second;
    }
    template<typename Service_type>
      int Base<Service_type>::
      getChunkTypeIDInt(std::string const ChunkTypeIDStr) const {
      if(chunkTypeIDToIntIDMap.find(ChunkTypeIDStr) == chunkTypeIDToIntIDMap.end()) {
	std::stringstream ss;
	ss << "In getChunkTypeIDInt(): string '" << ChunkTypeIDStr << "' not found in map";
	throw std::runtime_error(ss.str());
      }
      return chunkTypeIDToIntIDMap.find(ChunkTypeIDStr)->second;
    }


    template<typename Service_type>
      void Base<Service_type>::initTags() {
      int tag = 1;
      tagManager.registerPermanentTag("TAG_Create_small_chunk",        tag++);
      tagManager.registerPermanentTag("TAG_Create_large_chunk_prep",   tag++);
      tagManager.registerPermanentTag("TAG_Create_large_chunk",        tag++);
      tagManager.registerPermanentTag("TAG_Create_chunk_ack",          tag++);
      tagManager.registerPermanentTag("TAG_Delete_chunk",              tag++);
      tagManager.registerPermanentTag("TAG_Copy_chunk",                tag++);
      tagManager.registerPermanentTag("TAG_Copy_chunk_ack",            tag++);
      tagManager.registerPermanentTag("TAG_Get_chunk_data",            tag++);
      tagManager.registerPermanentTag("TAG_Time_to_stop",              tag++);
      tagManager.registerPermanentTag("TAG_Ready_to_stop",             tag++);
      tagManager.registerPermanentTag("TAG_ChunkTypeID_map",           tag++);
      tagManager.registerPermanentTag("TAG_chunk_service_params",      tag++);
    }


    template<typename Service_type>
      void Base<Service_type>::writeChunkObj_communicate(cht::vector<char> & buf,
								    int bufSzExceptTag,
								    int destRank,
								    MPI_Comm* comm) {
      MPI_Wrapper& MW = MPI_Wrapper::instance();
      /* Now there are two different ways to send the data to the chunk
       * owner worker process: either as a single message (if small enough) or
       * as two messages where the first one gives the size. */
      int bufSz = buf.size();
      if(bufSz > MAX_WRITECHUNKOBJ_SMALL_MESSAGE_SIZE) {
	// Too big: use two messages.  First send an empty "prep" message
	// to let the receiver know that a large message is to be
	// expected.
	MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, destRank, 
		    tagManager.getPermanentTag("TAG_Write_data_to_chunk_prep"), *comm);
	// Now send the "real" message.
	int TMP_TAG_Write_data_to_chunk_ack = tagManager.getTemporaryTag();
	memcpy(&buf[bufSzExceptTag], 
	       &TMP_TAG_Write_data_to_chunk_ack, sizeof(int));
	MW._MPI_Send(&buf[0], bufSz, MPI_UNSIGNED_CHAR, destRank, 
		    tagManager.getPermanentTag("TAG_Write_data_to_chunk_large"), *comm);
	MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, destRank, 
		    TMP_TAG_Write_data_to_chunk_ack, *comm, MPI_STATUS_IGNORE);
	tagManager.releaseTemporaryTag(TMP_TAG_Write_data_to_chunk_ack);
      }
      else {
	// Small enough: use a single message.
	int TMP_TAG_Write_data_to_chunk_ack = tagManager.getTemporaryTag();
	memcpy(&buf[bufSzExceptTag], 
	       &TMP_TAG_Write_data_to_chunk_ack, sizeof(int));
	MW._MPI_Send(&buf[0], bufSz, MPI_UNSIGNED_CHAR, destRank, 
		    tagManager.getPermanentTag("TAG_Write_data_to_chunk_small"), *comm);
	MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, destRank, 
		    TMP_TAG_Write_data_to_chunk_ack, *comm, MPI_STATUS_IGNORE);
	tagManager.releaseTemporaryTag(TMP_TAG_Write_data_to_chunk_ack);
      }
    }

  }; // end namespace 
}; // end namespace 
#endif
