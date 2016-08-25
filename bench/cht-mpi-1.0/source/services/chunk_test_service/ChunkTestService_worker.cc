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
#include "services/output_service/OutputService_worker.h"
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <unistd.h>
#include "utilities/cht_utils.h"
#include "services/chunk_service/ChunkService_worker.h"
#include "services/chunk_test_service/ChunkTestService_worker.h"

namespace cht {
  namespace ChunkTestService {

    void Worker::start_derived() {
      AccessKey key(this);
      OutputService::Worker::instance().outputInfo("This is chunk test service in start() (worker)");

      // Recieve a chunk id from parent.
      int tag = 77;
      ChunkService::ChunkID cidFromParent;
      MW._MPI_Recv(&cidFromParent, sizeof(ChunkService::ChunkID), 
		  MPI_UNSIGNED_CHAR, 0, tag, *key.comm_to_parent(), MPI_STATUS_IGNORE);

      // Play around with chunk id from parent only if this is one particular rank.
      if(key.my_rank() == 1)
	{
	  // Write some data to the chunk
	  int int1234567 = 1234567;
	  ChunkService::Worker::instance().writeDataToChunk(cidFromParent, (char*)&int1234567, sizeof(int));

	  // Get chunk size
	  int foreignChunkSize = ChunkService::Worker::instance().getChunkSize(cidFromParent);
	  char sss[888];
	  sprintf(sss, "ChunkTestService::start(), got foreignChunkSize = %d", foreignChunkSize);
	  OutputService::Worker::instance().outputInfo(sss);
	  if(foreignChunkSize != sizeof(int))
	    throw std::runtime_error("Error, (foreignChunkSize != sizeof(int))");
  
	  // Get chunk data
	  int valueFromForeignChunk = 0;
	  ChunkService::Worker::instance().getChunkData(cidFromParent, &valueFromForeignChunk, sizeof(int));
	  if(valueFromForeignChunk != 1234567)
	    throw std::runtime_error("Error, (valueFromForeignChunk != 1234567)");
	}
  
      // Send message to parent to indicate we are done with the chunk.
      MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, tag, *key.comm_to_parent());
  

      // Create a chunk
      OutputService::Worker::instance().outputInfo("ChunkTestService::start: Registering a chunk...");
      ChunkService::ChunkID cid = ChunkService::Worker::instance().registerChunk();

      // Write some data to the newly created chunk
      OutputService::Worker::instance().outputInfo("ChunkTestService::start: Writing some data to chunk...");
      int dummyInt1 = 6666;
      ChunkService::Worker::instance().writeDataToChunk(cid, (char*)&dummyInt1, 
							    sizeof(int));

      // Get chunk size
      OutputService::Worker::instance().outputInfo("ChunkTestService::start: Getting chunk size...");
      int chunkSize = ChunkService::Worker::instance().getChunkSize(cid);
      if(chunkSize != sizeof(int))
	throw std::runtime_error("Error, (chunkSize != sizeof(int))");

      // Get chunk data
      OutputService::Worker::instance().outputInfo("ChunkTestService::start: Getting chunk data...");
      int dummyInt2 = 0;
      ChunkService::Worker::instance().getChunkData(cid, &dummyInt2, sizeof(int));
      if(dummyInt2 != dummyInt1)
	throw std::runtime_error("Error, (dummyInt2 != dummyInt1)");

      // Delete chunk
      OutputService::Worker::instance().outputInfo("ChunkTestService::start: Deleting chunk...");
      ChunkService::Worker::instance().deleteChunk(cid);

      // Test communication speed for some test messages between rank 0
      // and rank 1.
      if(key.n_workers() < 2)
	OutputService::Worker::instance().outputInfo("Skipping comm speed test since (n_workers < 2).");
      else {
	OutputService::Worker::instance().outputInfo("Doing comm speed test for some test "
							 "messages between rank 0 and rank 1.");
	int nTests = 10;
	int largeMsgSize = 10000000;
	if(key.my_rank() == 0) {
	  sleep(1);
	  // Send some small messages.
	  for(int i = 0; i < nTests; i++) {
	    double startTime = cht::get_wall_seconds();
	    // Send a message.
	    int msgInt = 7;
	    int otherRank = 1;
	    MW._MPI_Send(&msgInt, sizeof(int), MPI_UNSIGNED_CHAR, otherRank, 0, *key.comm_worker_world());
	    // Recive response.
	    int receivedInt;
	    MW._MPI_Recv(&receivedInt, sizeof(int), MPI_UNSIGNED_CHAR, 
			otherRank, 0, *key.comm_worker_world(), MPI_STATUS_IGNORE);
	    double secondsTaken = cht::get_wall_seconds() - startTime;
	    std::stringstream ss;
	    ss << "Received response, " << receivedInt << ", took " 
	       << std::fixed << std::setw(10) << std::setprecision(3)
	       << secondsTaken*1000000 << " microseconds.";
	    OutputService::Worker::instance().outputInfo(ss.str());
	    usleep(1000);
	  }
	  // Now some larger messages to test transfer rate in GB/second.
	  cht::vector<char> buf(largeMsgSize);
	  for(int k = 0; k < largeMsgSize; k++)
	    buf[k] = 'a';
	  for(int i = 0; i < nTests; i++) {
	    double startTime = cht::get_wall_seconds();
	    // Send a message.
	    int otherRank = 1;
	    MW._MPI_Send(&buf[0], largeMsgSize, MPI_UNSIGNED_CHAR, otherRank, 0, *key.comm_worker_world());
	    // Recive response.
	    char receivedChar;
	    MW._MPI_Recv(&receivedChar, sizeof(char), MPI_UNSIGNED_CHAR, otherRank, 
			0, *key.comm_worker_world(), MPI_STATUS_IGNORE);
	    double secondsTaken = cht::get_wall_seconds() - startTime;
	    double bytesPerSecond = largeMsgSize / secondsTaken;
	    std::stringstream ss;
	    ss << "Received response after large message, " << receivedChar << ", took " 
	       << std::fixed << std::setw(10) << std::setprecision(6)
	       << secondsTaken << " seconds  <--> " 
	       << std::fixed << std::setw(10) << std::setprecision(6)
	       << bytesPerSecond / 1000000000 << " GB/second.";
	    OutputService::Worker::instance().outputInfo(ss.str());
	    usleep(2000);
	  }
	  // Send some small messages again.
	  for(int i = 0; i < nTests; i++) {
	    double startTime = cht::get_wall_seconds();
	    // Send a message.
	    int msgInt = 7;
	    int otherRank = 1;
	    MW._MPI_Send(&msgInt, sizeof(int), MPI_UNSIGNED_CHAR, otherRank, 0, *key.comm_worker_world());
	    // Recive response.
	    int receivedInt;
	    MW._MPI_Recv(&receivedInt, sizeof(int), MPI_UNSIGNED_CHAR,
			otherRank, 0, *key.comm_worker_world(), MPI_STATUS_IGNORE);
	    double secondsTaken = cht::get_wall_seconds() - startTime;
	    std::stringstream ss;
	    ss << "Received response, " << receivedInt << ", took "
	       << std::fixed << std::setw(10) << std::setprecision(3)
	       << secondsTaken*1000000 << " microseconds.";
	    OutputService::Worker::instance().outputInfo(ss.str());
	    usleep(1000);
	  }
	} // end if rank 0
	if(key.my_rank() == 1) {
	  for(int i = 0; i < nTests; i++) {
	    int receivedInt;
	    MW._MPI_Recv(&receivedInt, sizeof(int), MPI_UNSIGNED_CHAR, 0, 0, *key.comm_worker_world(), MPI_STATUS_IGNORE);    
	    int intToSend = receivedInt * 2;
	    MW._MPI_Send(&intToSend, sizeof(int), MPI_UNSIGNED_CHAR, 0, 0, *key.comm_worker_world());
	  }
	  cht::vector<char> buf(largeMsgSize);
	  for(int i = 0; i < nTests; i++) {
	    MW._MPI_Recv(&buf[0], largeMsgSize, MPI_UNSIGNED_CHAR, 0, 0, *key.comm_worker_world(), MPI_STATUS_IGNORE);
	    char charToSend = buf[largeMsgSize-1];
	    MW._MPI_Send(&charToSend, sizeof(char), MPI_UNSIGNED_CHAR, 0, 0, *key.comm_worker_world());
	  }
	  for(int i = 0; i < nTests; i++) {
	    int receivedInt;
	    MW._MPI_Recv(&receivedInt, sizeof(int), MPI_UNSIGNED_CHAR, 0, 0, *key.comm_worker_world(), MPI_STATUS_IGNORE);
	    int intToSend = receivedInt * 2;
	    MW._MPI_Send(&intToSend, sizeof(int), MPI_UNSIGNED_CHAR, 0, 0, *key.comm_worker_world());
	  }
	} // end else rank 1
      }
  
      OutputService::Worker::instance().outputInfo("Chunk tests finished OK!");  
    }


    void Worker::stop_derived() {
      OutputService::Worker::instance().outputInfo("This is chunk test service in stop() (worker)");
    }

  }; // end namespace
}; // end namespace
