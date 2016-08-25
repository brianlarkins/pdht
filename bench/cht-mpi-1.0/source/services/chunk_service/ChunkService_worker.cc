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
#include <iomanip>
#include <vector>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include "utilities/cht_utils.h"
#include "services/output_service/OutputService_worker.h"
#include "services/chunk_service/ChunkService_worker.h"

namespace cht {
  namespace ChunkService {

    void* global_thread_func(void* arg)
    {
      Worker* p = (Worker*) arg;
      try {
	p->worker_thread_func();
      } catch ( std::runtime_error e ) {
	std::cerr << "Error! Exception caught in ChunkService global_thread_func()." 
		  << " what(): " << e.what() << std::endl;
      } catch ( std::exception e ) {
	std::cerr << "Error! Exception caught in ChunkService global_thread_func()." 
		  << " what(): " << e.what() << std::endl;
      }
      return NULL;
    }

    // When the service starts, we want to create a thread that will be
    // responsible for listening for MPI messages from other processes
    // (workers or parent) who may request data from this process.
    void Worker::start_derived() {
      time_to_stop = false;
      if (threadRequestHandler)
	throw std::runtime_error("threadRequestHandler != NULL in "
				 "cht::ChunkService::Worker::start_derived()");
      threadRequestHandler = new Threads::Thread("CS-request", global_thread_func, this);
    }

    void Worker::stop_derived() {
      // Output some statistics.
      std::stringstream s;
      LockMutex();
      s << "ChunkService::stop() (worker), totDataSizeBytes = " << totDataSizeBytes 
	<< ", totDataSizeBytesMax = " << totDataSizeBytesMax;
      UnlockMutex();
      OutputService::Worker::instance().output(Output::Info, s.str());
      // Signal exit to thread
      LockMutex();
      time_to_stop = true;
      UnlockMutex();
      // Join thread here!!
      delete threadRequestHandler;
      threadRequestHandler = NULL;
    }

    Worker::Worker() 
      : totDataSizeBytes(0),
	totDataSizeBytesMax(0),
	threadRequestHandler(NULL)
    {
    }

    void Worker::LockMutex()
    {
      mutex.lock();
    }

    void Worker::UnlockMutex()
    {
      mutex.unlock();
    }

    Worker::Chunk* Worker::FindChunkInList(ChunkID id)
    {
      LockMutex();
      ChunkListMap::iterator it = chunkList.find(id);
      if(it == chunkList.end())
	throw std::runtime_error("Error, chunk not found");
      Chunk* theChunk = (*it).second;
      UnlockMutex();
      return theChunk;
    }

    Worker::Chunk* Worker::GetChunkPtrFromIdInBuffer(const char* buf) {
      ChunkID id;
      memcpy(&id, buf, sizeof(ChunkID));
      return FindChunkInList(id);
    }

    ChunkID Worker::GetChunkIdFromBuffer(const char* buf) {
      ChunkID id;
      memcpy(&id, buf, sizeof(ChunkID));
      return id;
    }

    void Worker::worker_thread_func() {
      AccessKey key(this);
  
      // Prepare list of tags that are relevant for this thread to probe for.
      const int nTags = 5;
      int tagList[nTags];
      tagList[0] = TAG_Create_chunk;
      tagList[1] = TAG_Write_data_to_chunk;
      tagList[2] = TAG_Get_chunk_size;
      tagList[3] = TAG_Get_chunk_data;
      tagList[4] = TAG_Delete_chunk;
      while(1) {
	usleep(10000);

	// Check for message from parent.
	int flag;
	MPI_Status status;
	// Note: instead of probing for MPI_ANY_TAG here, we use a list of
	// specific tags. Otherwise the program can deadlock, if we "by
	// mistake" catch one of the tags intended for the other thread.
	MW.NonMPI_Iprobe_multiple(0, tagList, nTags, *key.comm_to_parent(), &flag, &status);
	if (flag == 1) {
	  int msgSize;
	  MW._MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &msgSize);
	  cht::vector<char> buf(msgSize);
	  MW._MPI_Recv(&buf[0], msgSize, MPI_UNSIGNED_CHAR, 0, 
		      status.MPI_TAG, *key.comm_to_parent(), &status);
	  // OK, now we have a message from parent stored in buf.
	  if(status.MPI_TAG == TAG_Create_chunk) {
	    ChunkID id = registerChunk(); // This creates a chunk locally at this process.
	    MW._MPI_Send(&id, sizeof(ChunkID), MPI_UNSIGNED_CHAR, 0, TAG_Chunk_id, *key.comm_to_parent());
	  }
	  else if(status.MPI_TAG == TAG_Write_data_to_chunk) {
	    writeDataToChunk(GetChunkIdFromBuffer(&buf[0]), &buf[sizeof(ChunkID)], msgSize - sizeof(ChunkID));
	    MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, TAG_Write_data_to_chunk_ack, *key.comm_to_parent());
	  }
	  else if(status.MPI_TAG == TAG_Get_chunk_size) {
	    int chunkSize = getChunkSize(GetChunkIdFromBuffer(&buf[0]));
	    MW._MPI_Send(&chunkSize, sizeof(int), MPI_UNSIGNED_CHAR, 0, TAG_Chunk_size, *key.comm_to_parent());
	  }
	  else if(status.MPI_TAG == TAG_Get_chunk_data) {
	    Chunk* theChunk = GetChunkPtrFromIdInBuffer(&buf[0]);
	    if(!theChunk->readOnly)
	      throw std::runtime_error("Error: not readonly.");
	    MW._MPI_Send(&theChunk->dataPtr[0], theChunk->dataPtr.size(), MPI_UNSIGNED_CHAR, 0, TAG_Chunk_data, *key.comm_to_parent());
	  }
	  else if(status.MPI_TAG == TAG_Delete_chunk) {
	    deleteChunk(GetChunkIdFromBuffer(&buf[0]));
	    MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, TAG_Delete_chunk_ack, *key.comm_to_parent());
	  }
	  else {
	    throw std::runtime_error("Unexpected tag received from parent.");
	  }
	}
	// Check for message from workers.
	// Note: instead of probing for MPI_ANY_TAG here, we use a list of
	// specific tags. Otherwise the program can deadlock, if we "by
	// mistake" catch one of the tags intended for the other thread.
	MW.NonMPI_Iprobe_multiple(MPI_ANY_SOURCE, tagList, nTags, *key.comm_worker_world(), &flag, &status);
	if (flag == 1) {
	  int msgSize;
	  MW._MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &msgSize);
	  cht::vector<char> buf(msgSize);
	  MW._MPI_Recv(&buf[0], msgSize, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, 
		      status.MPI_TAG, *key.comm_worker_world(), &status);
	  // OK, now we have a message from a worker stored in buf.
	  if(status.MPI_TAG == TAG_Get_chunk_size) {
	    int chunkSize = getChunkSize(GetChunkIdFromBuffer(&buf[0]));
	    MW._MPI_Send(&chunkSize, sizeof(int), MPI_UNSIGNED_CHAR, status.MPI_SOURCE, TAG_Chunk_size, *key.comm_worker_world());
	  }
	  else if(status.MPI_TAG == TAG_Get_chunk_data) {
	    Chunk* theChunk = GetChunkPtrFromIdInBuffer(&buf[0]);
	    MW._MPI_Send(&theChunk->dataPtr[0], theChunk->dataPtr.size(), MPI_UNSIGNED_CHAR, 
			status.MPI_SOURCE, TAG_Chunk_data, *key.comm_worker_world());
	  }
	  else if(status.MPI_TAG == TAG_Write_data_to_chunk) {
	    writeDataToChunk(GetChunkIdFromBuffer(&buf[0]), &buf[sizeof(ChunkID)], msgSize - sizeof(ChunkID));
	    MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, TAG_Write_data_to_chunk_ack, *key.comm_worker_world());
	  }
	  else if(status.MPI_TAG == TAG_Delete_chunk) {
	    deleteChunk(GetChunkIdFromBuffer(&buf[0]));
	    MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, TAG_Delete_chunk_ack, *key.comm_worker_world());
	  }
	  else {
	    throw std::runtime_error("Unexpected tag received from worker.");
	  }
	}
	LockMutex();
	if (time_to_stop) {
	  UnlockMutex();
	  break;
	}
	UnlockMutex();
      } // end while
    }

    ChunkID Worker::registerChunk() {
      AccessKey key(this);
      // This means a new chunk is to be created here, at this worker
      // process.
      Chunk* chunk = new Chunk;
      chunk->id.ownerRank = key.my_rank();
      chunk->id.localID = ++idCounter;
      LockMutex();
      chunkList.insert(ChunkListMap::value_type(chunk->id, chunk));
      UnlockMutex();
      return chunk->id;
    }

    void Worker::writeDataToChunk(ChunkID id, const char* dataPtr, int dataSize) {
      AccessKey key(this);
      if(id.ownerRank == key.my_rank()) {
	// This chunk is stored locally by this process.
	Chunk* theChunk = FindChunkInList(id);
	if(theChunk->readOnly)
	  throw std::runtime_error("Error: readonly.");
	theChunk->dataPtr.setsize(dataSize);
	memcpy(&theChunk->dataPtr[0], dataPtr, dataSize);
	theChunk->readOnly = true;
	LockMutex();
	totDataSizeBytes += dataSize;
	totDataSizeBytesMax = (totDataSizeBytes > totDataSizeBytesMax) ?  totDataSizeBytes : totDataSizeBytesMax;
	UnlockMutex();
	return;
      }
      // Need to communicate!
      int bufSz = sizeof(ChunkID) + dataSize;
      cht::vector<char> buf(bufSz);
      memcpy(&buf[0], &id, sizeof(ChunkID));
      memcpy(&buf[sizeof(ChunkID)], dataPtr, dataSize);
      MW._MPI_Send(&buf[0], bufSz, MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Write_data_to_chunk, *key.comm_worker_world());
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Write_data_to_chunk_ack, *key.comm_worker_world(), MPI_STATUS_IGNORE);
    }

    void Worker::deleteChunk(ChunkID id) {
      AccessKey key(this);
      if(id.ownerRank == key.my_rank()) {
	// This chunk is stored locally by this process.
	Chunk* theChunk = FindChunkInList(id);
	LockMutex();
	totDataSizeBytes -= theChunk->dataPtr.size();
	delete theChunk;
	chunkList.erase(id);
	UnlockMutex();
	return;
      }
      // Need to communicate!
      MW._MPI_Send(&id, sizeof(ChunkID), MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Delete_chunk, *key.comm_worker_world());
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Delete_chunk_ack, *key.comm_worker_world(), MPI_STATUS_IGNORE);
    }

    int Worker::getChunkSize(ChunkID id) {
      AccessKey key(this);
      if(id.ownerRank == key.my_rank()) {
	// This chunk is stored locally by this process.
	Chunk* theChunk = FindChunkInList(id);
	return theChunk->dataPtr.size();
      }
      // Need to communicate!
      MW._MPI_Send(&id, sizeof(ChunkID), MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Get_chunk_size, *key.comm_worker_world());
      int chunkSize;
      MW._MPI_Recv(&chunkSize, sizeof(int), MPI_UNSIGNED_CHAR, id.ownerRank, 
		  TAG_Chunk_size, *key.comm_worker_world(), MPI_STATUS_IGNORE);
      return chunkSize;
    }

    void Worker::getChunkData(ChunkID id, void* dataPtr, int dataSize) {
      AccessKey key(this);
      if(id.ownerRank == key.my_rank()) {
	clock_t startClock = clock();
	// This chunk is stored locally by this process.
	Chunk* theChunk = FindChunkInList(id);
	if(!theChunk->readOnly)
	  throw std::runtime_error("Error in ChunkService::getChunkData: not readonly.");
	if(theChunk->dataPtr.size() != dataSize)
	  throw std::runtime_error("Error in ChunkService::getChunkData: wrong size.");
	memcpy(dataPtr, &theChunk->dataPtr[0], dataSize);
	double secondsTaken = ((double)(clock() - startClock)) / CLOCKS_PER_SEC;
	std::stringstream s;
	s << "ChunkService::getChunkData (locally) took " << secondsTaken << " cpu s for dataSize = " << dataSize;
	OutputService::Worker::instance().outputDebug(s.str());
	return;
      }

      // Need to communicate!
      MW._MPI_Send(&id, sizeof(ChunkID), MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Get_chunk_data, *key.comm_worker_world());
      MW._MPI_Recv(dataPtr, dataSize, MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Chunk_data, *key.comm_worker_world(), MPI_STATUS_IGNORE);
    }


  }; // end namespace 
}; // end namespace 
