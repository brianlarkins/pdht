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

// The following files are normally not included in user code. Needed
// here because this test tests some internal stuff.
#include "services/chunk_service/ChunkService_parent.h"
#include "services/service_manager/ServiceManager.h"

#include <stdexcept>
#include <unistd.h>
#include <iostream>
#include "chunks_and_tasks.h"
#include "test_service_chunk_obj_implementations.h"
#include "test_tasksched_task_implementations_ABC.h"

int main() {

  int n_workers = 5;

  cht::extras::setNWorkers(n_workers);
  cht::start();

  std::cout<<"test_service: Creating A_chunk chunk"<<std::endl;
  A_chunk* tmpPtr = new A_chunk;
  tmpPtr->xx = 5;
  cht::ChunkID cid_ac_1 = cht::registerChunk<A_chunk>(tmpPtr);

  std::cout<<"test_service: Getting Chunk Obj"<<std::endl;
  cht::shared_ptr<A_chunk const> aPtr_2;
  cht::getChunk(cid_ac_1, aPtr_2);
  std::cout<<"test_service: After Chunk Obj"<<std::endl;
  if (aPtr_2 == 0)
    throw std::runtime_error("aPtr_2 == 0!!");
  if (aPtr_2->xx != 5)
    throw std::runtime_error("This doesn't work - lets go home ;(");
  std::cout << "SUCCESS!!! aPtr_2->xx == " << aPtr_2->xx <<std::endl;

  std::cout<<"test_service: deleting Chunk Obj"<<std::endl;
  cht::deleteChunk(cid_ac_1);


#if 1
  /// Cache stuff
  size_t oneGB = 1000000000;
  size_t cacheMemoryUsageLimit = oneGB;
  //  ChunkObjCacheService::Parent<ChunkObjService::Parent>::instance().setCacheParams(cht::Cache::Enabled, cacheMemoryUsageLimit);
  std::cout<<"test_service: Creating a chunk via cache"<<std::endl;

  A_chunk* tmpPtr_cache = new A_chunk;
  tmpPtr_cache->xx = 4;
  cht::ChunkID cid_1_cache = cht::registerChunk<A_chunk>(tmpPtr_cache);

  cht::shared_ptr<A_chunk const> sPtr_cache;
  cht::getChunk(cid_1_cache, sPtr_cache);
  if (sPtr_cache == 0)
    throw std::runtime_error("sPtr_cache == 0!!");
  if (sPtr_cache->xx != 4)
    throw std::runtime_error("This doesn't work - lets go home ;(");
  std::cout << "SUCCESS!!! sPtr_cache->xx == " << sPtr_cache->xx <<std::endl;
  std::cout<<"test_service: deleting Chunk Obj via cache"<<std::endl;
  cht::deleteChunk(cid_1_cache);
  /// End cache stuff
#endif

  std::cout<<"test_service: stopping  CHUNK OBJ service"<<std::endl;
  cht::ServiceManager::instance().stopService("ChunkObjService");
  std::cout<<"test_service: starting CHUNK OBJ service"<<std::endl;
  cht::ServiceManager::instance().startService("ChunkObjService");


  std::cout<<"test_service: starting CHUNK service"<<std::endl;
  cht::ServiceManager::instance().startService("ChunkService");

  // Create a chunk
  cht::ChunkService::ChunkID cid = cht::ChunkService::Parent::instance().registerChunk();
  // Write some data to the newly created chunk
  int dummyInt1 = 6666;
  cht::ChunkService::Parent::instance().writeDataToChunk(cid, (char*)&dummyInt1, 
							 sizeof(int));
  
  // Get chunk size
  int chunkSize = cht::ChunkService::Parent::instance().getChunkSize(cid);
  if(chunkSize != sizeof(int))
    throw std::runtime_error("Error, (chunkSize != sizeof(int))");
  
  // Get chunk data
  int dummyInt2 = 0;
  cht::ChunkService::Parent::instance().getChunkData(cid, &dummyInt2, sizeof(int));
  if(dummyInt2 != dummyInt1)
    throw std::runtime_error("Error, (dummyInt2 != dummyInt1)");
  
  // Delete chunk
  cht::ChunkService::Parent::instance().deleteChunk(cid);

  std::cout<<"test_service: starting CHUNK TEST service"<<std::endl;
  cht::ServiceManager::instance().startService("ChunkTestService");

  std::cout<<"test_service: stopping CHUNK TEST service"<<std::endl;
  cht::ServiceManager::instance().stopService("ChunkTestService");

  std::cout<<"test_service: stopping CHUNK service"<<std::endl;
  cht::ServiceManager::instance().stopService("ChunkService");
  std::cout<<"test_service: after stopping CHUNK service"<<std::endl;

  A_chunk* aaa = new A_chunk;
  aaa->xx = 4;
  cht::ChunkID cid_a = cht::registerChunk<A_chunk>(aaa);

  std::vector<cht::ChunkID> inputChunks(1);
  inputChunks[0] = cid_a;
  cht::executeMotherTask<TaskTypeA>(inputChunks);
  cht::executeMotherTask<TaskTypeA>(inputChunks);
  cht::executeMotherTask<TaskTypeA>(inputChunks);
  cht::deleteChunk(cid_a);

  std::cout<<"test_service: stopping TASK SCHEDULER service"<<std::endl;
  cht::ServiceManager::instance().stopService("TaskSchedulerService");

  std::cout<<"test_service: starting TASK SCHEDULER service"<<std::endl;
  cht::ServiceManager::instance().startService("TaskSchedulerService");

  // Stop cht services
  std::cout<<"test_service: Before cht::stop()"<<std::endl;
  cht::stop();
  std::cout<<"test_service: After cht::stop()"<<std::endl;

  std::cout<<"test_service: DONE."<<std::endl;
  
}
