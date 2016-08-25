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
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include "utilities/cht_utils.h"
#include "services/output_service/OutputService_parent.h"
#include "services/chunk_service/ChunkService_parent.h"



namespace cht {
  namespace ChunkService {

    ChunkID Parent::registerChunk() {
      AccessKey key(this);
      OutputService::Parent::instance().outputDebug("entering ChunkService::Parent::registerChunk() (parent)");
      // Send message to any worker, randomly selected.
      int rank = (int)(key.n_workers() * ( (double)rand() / RAND_MAX ));
      if(rank < 0 || rank >= key.n_workers())
	throw std::runtime_error("Error rand() !!!!.");    
      MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, rank, TAG_Create_chunk, *key.comm_to_workers());
      OutputService::Parent::instance().outputDebug("ChunkService::Parent::registerChunk() after send.");
      // Wait for reply from worker.
      MPI_Status status;
      ChunkID chunkID;
      MW._MPI_Recv(&chunkID, sizeof(ChunkID), MPI_UNSIGNED_CHAR, rank, TAG_Chunk_id, *key.comm_to_workers(), &status);
      return chunkID;
    }

    void Parent::writeDataToChunk(ChunkID id, const char* dataPtr, int dataSize) {
      AccessKey key(this);
      OutputService::Parent::instance().outputDebug("entering ChunkService::writeDataToChunk() (parent)");
      // Create a buffer where the ChunkID is first, then the data.
      int bufSz = sizeof(ChunkID) + dataSize;
      cht::vector<char> buf(bufSz);
      memcpy(&buf[0], &id, sizeof(ChunkID));
      memcpy(&buf[sizeof(ChunkID)], dataPtr, dataSize);
      MW._MPI_Send(&buf[0], bufSz, MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Write_data_to_chunk, *key.comm_to_workers());
      // Wait for reply from worker.
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Write_data_to_chunk_ack, *key.comm_to_workers(), MPI_STATUS_IGNORE);
    }

    void Parent::deleteChunk(ChunkID id) {
      AccessKey key(this);
      MW._MPI_Send(&id, sizeof(ChunkID), MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Delete_chunk, *key.comm_to_workers());
      // Wait for reply from worker.
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Delete_chunk_ack, *key.comm_to_workers(), MPI_STATUS_IGNORE);
    }

    int Parent::getChunkSize(ChunkID id) {
      AccessKey key(this);
      MW._MPI_Send(&id, sizeof(ChunkID), MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Get_chunk_size, *key.comm_to_workers());
      // Wait for reply from worker.
      int chunkSize;
      MW._MPI_Recv(&chunkSize, sizeof(int), MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Chunk_size, *key.comm_to_workers(), MPI_STATUS_IGNORE);
      return chunkSize;
    }

    void Parent::getChunkData(ChunkID id, void* dataPtr, int dataSize) {
      AccessKey key(this);
      MW._MPI_Send(&id, sizeof(ChunkID), MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Get_chunk_data, *key.comm_to_workers());
      // Wait for reply from worker.
      MW._MPI_Recv(dataPtr, dataSize, MPI_UNSIGNED_CHAR, id.ownerRank, TAG_Chunk_data, *key.comm_to_workers(), MPI_STATUS_IGNORE);
    }


  }; // end namespace 
}; // end namespace 
