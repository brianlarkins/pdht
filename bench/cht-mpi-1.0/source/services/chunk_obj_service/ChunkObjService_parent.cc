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
#include "services/chunk_obj_service/ChunkObjService_parent.h"
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include "utilities/cht_utils.h"
#include "services/services_utils.h"

namespace cht {
  namespace ChunkObjService {

    void Parent::start_derived() {
      AccessKey key(this);
      std::list<std::string> strList;
      cht::obj_factory<Chunk>::instance().getObjectTypeIDStrings(strList);
      populateChunkIDMaps(strList);
      cht::Service::sendStrBufToWorkers(strList, 
					tagManager.getPermanentTag("TAG_ChunkTypeID_map"),
					key.comm_to_workers(),
					MW);

      // Send parameters to workers.
      for(int i = 0; i < key.n_workers(); i++)
	MW._MPI_Send(&probabilityToRegisterChunksLocally, sizeof(double), 
		    MPI_UNSIGNED_CHAR, i, 
		    tagManager.getPermanentTag("TAG_chunk_service_params"), *key.comm_to_workers());
    }
    
    void Parent::stop_derived() {
      AccessKey key(this);
      int n_workers = key.n_workers(); 
      MPI_Comm* comm_to_workers = key.comm_to_workers();
      // Send message to workers saying it is time to quit.
      for(int i = 0; i < n_workers; i++)
	MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, i, tagManager.getPermanentTag("TAG_Time_to_stop"), *comm_to_workers);
      // Receive signals from workers
      for(int i = 0; i < n_workers; i++)
	MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, i, 
		    tagManager.getPermanentTag("TAG_Ready_to_stop"), *comm_to_workers, MPI_STATUS_IGNORE);
      // send signal to workers: all workers are ready
      for(int i = 0; i < n_workers; i++)
	MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, i, tagManager.getPermanentTag("TAG_Ready_to_stop"), *comm_to_workers);  
    }

    Parent::Parent() 
      : OS( OutputService::Parent::instance() ),
	probabilityToRegisterChunksLocally(1.0)
    {
      initTags();
    }

    void Parent::deleteChunk(ChunkID cid) {
      AccessKey key(this);
      if (cid == CHUNK_ID_NULL)
	return;
      MW._MPI_Send(&cid, sizeof(ChunkID), MPI_UNSIGNED_CHAR, cid.ownerRank, 
		  tagManager.getPermanentTag("TAG_Delete_chunk"), *key.comm_to_workers());
      // NOTE: earlier we waited for an acknowledgment message here,
      // but now that has been removed. This means that this routine
      // can (and probably will) return before the chunk is actually
      // deleted.
    }

    bool Parent::getChunkIfExists(ChunkID cid, 
				  cht::shared_ptr<Chunk const> & objPtr) {
      AccessKey key(this);
      std::string cid_class_id_str = getChunkTypeIDStr(cid.chunkTypeID);
      int bufSz = sizeof(ChunkID) + 2*sizeof(int);
      cht::vector<char> buf(bufSz);
      memcpy(&buf[0], &cid, sizeof(ChunkID));
      memcpy(&buf[sizeof(ChunkID)], &cid.chunkTypeID, sizeof(int));
      int TMP_TAG_Chunk_data = tagManager.getTemporaryTag_incoming();
      memcpy(&buf[sizeof(ChunkID) + sizeof(int)], &TMP_TAG_Chunk_data, sizeof(int));
      MW._MPI_Send(&buf[0], bufSz, MPI_UNSIGNED_CHAR, cid.ownerRank, 
		  tagManager.getPermanentTag("TAG_Get_chunk_data"), 
		  *key.comm_to_workers());
      // Wait for reply from worker. We do not know what size to expect.
      MPI_Status status;
      MW._MPI_Probe(cid.ownerRank, TMP_TAG_Chunk_data, *key.comm_to_workers(), &status);
      int msgSizeTot;
      MW._MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &msgSizeTot);
      cht::vector<char> recvBuf(msgSizeTot);
      MW._MPI_Recv(&recvBuf[0], msgSizeTot, MPI_UNSIGNED_CHAR, cid.ownerRank, 
		  TMP_TAG_Chunk_data, *key.comm_to_workers(), MPI_STATUS_IGNORE);
      tagManager.releaseTemporaryTag(TMP_TAG_Chunk_data);
      // Now we have received the object data.
      // The received buffer contains first one bool saying if the chunk exists or not, and then the chunk data (if hte chunk exists).
      bool theBool = false;
      memcpy(&theBool, &recvBuf[0], sizeof(bool));
      if(theBool == false) {
	// This means the chunk does not exist. Probably it has been deleted.
	return false;
      }
      size_t chunkDataSize = msgSizeTot - sizeof(bool);
      const char* chunkDataPtr = &recvBuf[sizeof(bool)];
      Chunk* newObj = cht::obj_factory<Chunk>::instance().createObject(cid_class_id_str);
      newObj->assignFromBuffer (chunkDataPtr, chunkDataSize);
      objPtr = newObj; 		/* Transfer responsibility of new object to shared_ptr */  
      return true;
    }

    void Parent::getChunk(ChunkID cid, 
			  cht::shared_ptr<Chunk const> & objPtr) {
      bool success = getChunkIfExists(cid, objPtr);
      if(success == false)
	throw std::runtime_error("Error in Parent::getChunk: chunk does not seem to exist.");
    }

    ChunkID Parent::registerChunkDirectly(Chunk const * objPtr,
					std::string class_id_str) {
      cht::shared_ptr<Chunk const> objPtrDummy;
      return registerAndGetChunk(objPtr, objPtrDummy, class_id_str);
    }

    ChunkID Parent::registerAndGetChunk(Chunk const * objPtrInput,
					cht::shared_ptr<Chunk const> & objPtrOutput,
					std::string class_id_str) {
      AccessKey key(this);
      int chunkTypeID = getChunkTypeIDInt(class_id_str);
      // Create a buffer where the chunkTypeID is first, then a ChunkID which we set to NULL, then the data, then a tag to use for ack.
      size_t objDataSz = objPtrInput->getSize();
      int bufSz = sizeof(int) + sizeof(ChunkID) + objDataSz + sizeof(int);
      cht::vector<char> buf(bufSz);
      memcpy(&buf[0], &chunkTypeID, sizeof(int));
      ChunkID cid_null;
      memcpy(&buf[sizeof(int)], &cid_null, sizeof(ChunkID));
      objPtrInput->writeToBuffer(&buf[sizeof(int)+sizeof(ChunkID)], objDataSz);
      int TMP_TAG_Chunk_data = tagManager.getTemporaryTag_outgoing();
      int TMP_TAG_Write_data_to_chunk_ack = tagManager.getTemporaryTag_incoming();
      int bufSzExceptTag = bufSz - sizeof(int);
      memcpy(&buf[bufSzExceptTag], 
	     &TMP_TAG_Write_data_to_chunk_ack, sizeof(int));
      // Send message to any worker, randomly selected.
      int destRank = (int)(key.n_workers() * ( (double)rand() / RAND_MAX ));
      if(destRank < 0 || destRank >= key.n_workers())
	throw std::runtime_error("Error rand() !!!!.");      
      ChunkID cid;
      if(bufSz > MAX_WRITECHUNKOBJ_SMALL_MESSAGE_SIZE) {
	// Too big: use two messages. 
	// First send a TAG_Create_large_chunk_prep message containing 2 int values: the tag and the size of the large message that will come.
	int intsToSend[2];
	intsToSend[0] = TMP_TAG_Chunk_data;
	intsToSend[1] = bufSz;
	MW._MPI_Send(&intsToSend[0], 2*sizeof(int), MPI_UNSIGNED_CHAR, destRank, 
		    tagManager.getPermanentTag("TAG_Create_large_chunk_prep"), *key.comm_to_workers());
	// Now send the "real" message.
	MW._MPI_Send(&buf[0], bufSz, MPI_UNSIGNED_CHAR, destRank, 
		    TMP_TAG_Chunk_data, *key.comm_to_workers());
      }
      else {
	// Small enough: use a single message.
	MW._MPI_Send(&buf[0], bufSz, MPI_UNSIGNED_CHAR, destRank, 
		    tagManager.getPermanentTag("TAG_Create_small_chunk"), *key.comm_to_workers());
      }
      // Now we expect to receive an ack message containing a ChunkID and an int (tmp tag unused in this case).
      int recvBufSize = sizeof(ChunkID) + sizeof(int);
      char recvBuf[recvBufSize];
      MW._MPI_Recv(recvBuf, recvBufSize, MPI_UNSIGNED_CHAR, destRank, 
		  TMP_TAG_Write_data_to_chunk_ack, *key.comm_to_workers(), MPI_STATUS_IGNORE);
      memcpy(&cid, &recvBuf[0], sizeof(ChunkID));
      tagManager.releaseTemporaryTag(TMP_TAG_Chunk_data);
      tagManager.releaseTemporaryTag(TMP_TAG_Write_data_to_chunk_ack);
      objPtrOutput = objPtrInput;  /* Note that this assignment means transferring 
				      the responsibility for the object! */
      return cid;
    }
    
    bool Parent::chunkResidesLocally(ChunkID cid) {
      return false;
    }

    void Parent::setDebugParams(double probabilityToRegisterChunksLocally_) {
      // This routine may only be called when the service is stopped,
      // since the parameters are transferred to workers in start(). We
      // check that service is stopped by demanding that comm_to_workers
      // pointer is NULL.
      if ( serviceIsRunning() )
	throw std::runtime_error("Error! cht::TaskSchedulerService::Parent::setParams() called while service running.");
      // Be aware that service in principle could start to run while
      // setting the parameters below...
      probabilityToRegisterChunksLocally = probabilityToRegisterChunksLocally_;
    }


  }; // end namespace 
}; // end namespace 

