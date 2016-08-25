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
#include "services/chunk_obj_service/ChunkObjService_worker.h"
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <sstream>
#include <iostream>
#include <cassert>
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include "utilities/cht_utils.h"
#include "services/services_utils.h"

namespace cht {
  namespace ChunkObjService {

    const int smallMsgSizeLimitBytes = 1000;

    void* global_thread_func(void* arg)
    {
      Worker* p = (Worker*) arg;
      try {
	p->worker_thread_func();
      } 
      catch ( std::exception & e ) {
	std::cerr << "Error  Exception caught in cht::ChunkObjService::Worker global_thread_func(), what(): " << e.what() <<std::endl;
	p->call_mpi_abort();
      }
      catch ( ... ) {
	std::cerr << "Error! Exception caught in cht::ChunkObjService::Worker global_thread_func()." << std::endl;
	p->call_mpi_abort();
      }
      return NULL;
    }

    void Worker::call_mpi_abort() {
      AccessKey key(this);
      MW._MPI_Abort(*key.comm_to_parent(), -1);
    }

    void Worker::start_derived() {
      AccessKey key(this);
      std::list<std::string> strList;
      cht::Service::receiveStrBufFromParent(strList, 
					    tagManager.getPermanentTag("TAG_ChunkTypeID_map"),
					    key.comm_to_parent(),  
					    MW);
      populateChunkIDMaps(strList);
      // Check that chunktype ID strings match the ones in local
      // Chunk object factory
      std::list<std::string> localStrList;
      cht::obj_factory<Chunk>::instance().getObjectTypeIDStrings(localStrList);
      checkStrListAgainstIDMap(localStrList);

      // Receive parameters from parent.
      MW._MPI_Recv(&probabilityToRegisterChunksLocally, sizeof(double), MPI_UNSIGNED_CHAR, 0, 
		  tagManager.getPermanentTag("TAG_chunk_service_params"), *key.comm_to_parent(), MPI_STATUS_IGNORE);
      if(probabilityToRegisterChunksLocally < 0 || probabilityToRegisterChunksLocally > 1)
	throw std::runtime_error("Error: bad value for probabilityToRegisterChunksLocally.");

      // When the service starts, we want to create a thread that will be
      // responsible for listening for MPI messages from other processes
      // (workers or parent) who may request data from this process.
      if (threadRequestHandler)
	throw std::runtime_error("threadRequestHandler != NULL in "
				 "cht::ChunkObjService::Worker::start_derived()");
      
      threadRequestHandler = new Threads::Thread("COS-request", global_thread_func, this);
    }

    void Worker::stop_derived() {
      AccessKey key(this);

      // It may happen that other workers are about to send
      // "deletechunk" messages to this worker to delete chunks that
      // this worker owns, but that those messages have not yet
      // arrived. Therefore, we wait here until all chunks have been
      // removed.
      int waitCounter = 0;
      int waitTimeTot_us = 0;
      while(1) {
	// Check how many chunks are left.
	LockMutex();
	int noOfChunks = chunkList.size();
	UnlockMutex();
	if(noOfChunks == 0)
	  break;
	const size_t waitTime_us = 10000;
	usleep(waitTime_us);
	waitTimeTot_us += waitTime_us;
	if(waitTimeTot_us > 10000000)
	  throw std::runtime_error("ChunkObjService::Worker::stop_derived() error: "
				   "waited too long (more than 10 seconds) for chunks to be removed. Chunk leaks? Extremely slow communication?");
      }

      // Output some statistics.
      std::stringstream s;
      LockMutex();
      s << "cht::ChunkObjService::Worker::stop() (worker), totDataSizeBytes = " << totDataSizeBytes 
	<< ", totDataSizeBytesMax = " << totDataSizeBytesMax;
      OS.outputInfo(s.str());
      {
	// Output memoryBufferPool statistics also
	s.str("");
	size_t noOfBuffers, noOfUsedBuffers, totAllocatedSize;
	memoryBufferPool.getStatistics(noOfBuffers, noOfUsedBuffers, totAllocatedSize);      
	s << "cht::ChunkObjService::Worker::stop() (worker), memoryBufferPool statistics: noOfBuffers = " << noOfBuffers << " noOfUsedBuffers = " << noOfUsedBuffers << " totAllocatedSize = " << totAllocatedSize;
	OS.outputInfo(s.str());
      }

      if(totDataSizeBytes > 0) {
	int leakCount = 0;
	std::stringstream s;
	s << "Leaked objects: " << std::endl;
	ChunkListMap::iterator it = chunkList.begin();
	while(it != chunkList.end()) {
	  int theIntId = it->first.chunkTypeID;
	  s << " '" << getChunkTypeIDStr(theIntId) << "'" << std::endl;
	  it++;
	  leakCount++;
	  if(leakCount > 100) {
	    s << "Breaking, cannot show so many leaks.  :-/" << std::endl;
	    break;
	  }
	}
	OS.outputInfo(s.str());
      }

      UnlockMutex();
      // Parent will send exit message that will cause the thread to exit,
      // so we just delete the thread here (the destructor will invoke join).  
      //
      // EMANUEL COMMENT: Perhaps it would be nicer to signal the thread
      // here? As it is now two "stop-messages" are sent from parent, one
      // for the base class and one for the derived class - could the stop 
      // be handled by the base class via this stop_derived function?
      delete threadRequestHandler;
      threadRequestHandler = NULL;
      // signal to parent: I am ready
      MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, 0, tagManager.getPermanentTag("TAG_Ready_to_stop"), *key.comm_to_parent());  
      // receive signal from parent
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, 0, 
		  tagManager.getPermanentTag("TAG_Ready_to_stop"), *key.comm_to_parent(), MPI_STATUS_IGNORE);
    }

    Worker::Worker() 
      : totDataSizeBytes(0),
	totDataSizeBytesMax(0),
	getChunk_total_count(0),
	getChunk_found_locally_count(0),
	getChunkIfLocal_total_count(0),
	getChunkIfLocal_found_locally_count(0),
	getChunk_small_totNoOfBytesReceived(0),
	getChunk_large_totNoOfBytesReceived(0),
	OS( OutputService::Worker::instance() ),
	threadRequestHandler(NULL)
    {
      initTags();
    }

    void Worker::LockMutex()
    {
      mutex.lock();
    }

    void Worker::UnlockMutex()
    {
      mutex.unlock();
    }

    bool Worker::TryToGetChunkFromList(ChunkID cid, 
				       cht::shared_ptr<Chunk const> & objPtr) {
      LockMutex();
      ChunkListMap::iterator it = chunkList.find(cid);
      if(it == chunkList.end()) {
	UnlockMutex();
	return false;
      }
      ChunkStruct* theChunk = (*it).second;
      if(theChunk == 0)
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::TryToGetChunkFromList: (theChunk == 0).");
      if(theChunk->objPtr == 0)
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::TryToGetChunkFromList: (theChunk->objPtr == 0).");
      objPtr = theChunk->objPtr;
      UnlockMutex();
      return true;
    }

    Worker::ChunkStruct* Worker::TryToFindChunkInList(ChunkID cid, 
						      bool objShouldExist) {
      LockMutex();
      ChunkListMap::iterator it = chunkList.find(cid);
      if(it == chunkList.end())
	return NULL;
      ChunkStruct* theChunk = (*it).second;
      if(theChunk == 0)
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::TryToFindChunkInList: (theChunk == 0).");
      if(objShouldExist && theChunk->objPtr == 0)
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::TryToFindChunkInList: (objShouldExist && objPtr == 0).");
      if(!objShouldExist && !(theChunk->objPtr == 0))
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::TryToFindChunkInList: (!objShouldExist && objPtr != 0).");
      UnlockMutex();
      return theChunk;
    }

    Worker::ChunkStruct* Worker::FindChunkInList(ChunkID cid, 
						 bool objShouldExist) {
      ChunkStruct* theChunk = TryToFindChunkInList(cid, objShouldExist);
      if(theChunk == 0)
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::FindChunkInList: (theChunk == 0).");
      return theChunk;
    }

    ChunkID Worker::GetChunkIdFromBuffer(const char* buf) {
      ChunkID cid;
      memcpy(&cid, buf, sizeof(ChunkID));
      return cid;
    }
    
    int Worker::getIntFromEndOfBuffer(const char* buf, int bufSz) {
      int i;
      memcpy(&i, &buf[bufSz-sizeof(int)], sizeof(int));
      return i;
    }




    void Worker::PostAsyncSendMessage(AsyncMsgStruct* msgStruct) {
      MW._MPI_Isend(msgStruct->bufPtr, msgStruct->bufSz, 
		   MPI_UNSIGNED_CHAR, msgStruct->rank, msgStruct->tag, *(msgStruct->comm), &msgStruct->request);
      LockMutex();
      sendMsgStructList.push_back(msgStruct);
      UnlockMutex();
    }

    void Worker::PostAsyncRecvMessage(AsyncMsgStruct* msgStruct) {
      MW._MPI_Irecv(msgStruct->bufPtr, msgStruct->bufSz, 
		   MPI_UNSIGNED_CHAR, msgStruct->rank, msgStruct->tag, *(msgStruct->comm), &msgStruct->request);
      LockMutex();
      recvMsgStructList.push_back(msgStruct);
      UnlockMutex();
    }

    bool Worker::sendMsgStructListIsEmpty() {
      LockMutex();
      bool listIsEmpty = sendMsgStructList.empty();
      UnlockMutex();
      return listIsEmpty;
    }

    void Worker::FreeMemoryForCompletedSendOperations() {
      cht::timer timer;
      LockMutex();
      std::list<AsyncMsgStruct*>::iterator it = sendMsgStructList.begin();
      while(it != sendMsgStructList.end()) {
	AsyncMsgStruct* msgStructPtr = (*it);
	int flag;
	MPI_Status st_dummy; // Bug in Open MPI when passing MPI_STATUS_IGNORE ?
	MW._MPI_Test(&msgStructPtr->request, &flag, &st_dummy);
	if(flag) {
	  it = sendMsgStructList.erase(it);
	  delete msgStructPtr;
	}
	else
	  it++;
      }
      FreeMemoryForCompletedSendOperations_stats.add_timings(timer);
      UnlockMutex();
    }


    /*
     * handle_request(). This routine is used to handle incoming requests
     * from parent or from other workers, depending on the comm argument
     * that may be set to either comm_to_parent or comm_worker_world.
     */
    void Worker::handle_request(const char* buf,
				int msgSize,
				int tag,
				int rank,
				MPI_Comm* comm,
				bool & timeToStopFlag) {
      AccessKey key(this);
      int my_rank = key.my_rank(); 
      if(tag == tagManager.getPermanentTag("TAG_Get_chunk_data")) {
	ChunkID cid;     // Used to find the chunk
	int chunkTypeID; // Only used for type check
	memcpy(&cid, buf, sizeof(ChunkID));
	memcpy(&chunkTypeID, &buf[sizeof(ChunkID)], sizeof(int));
	assert(cid.ownerRank == my_rank);
	cht::shared_ptr<Chunk const> objPtr;
	bool chunkExists = TryToGetChunkFromList(cid, objPtr);
	if(chunkExists == false) {
	  // This means the chunk does not exist; probably it has been deleted earlier.
	  // We need to send a message containing just one bool with the value false.
	  bool theFalseBool = false;
	  int TMP_TAG_Chunk_data = getIntFromEndOfBuffer(&buf[0], msgSize);
	  AsyncMsgStruct* sendMsgStruct = new AsyncMsgStruct(sizeof(bool), rank, TMP_TAG_Chunk_data, comm, memoryBufferPool);
	  memcpy(sendMsgStruct->bufPtr, &theFalseBool, sizeof(bool));
	  PostAsyncSendMessage(sendMsgStruct);
	}
	else {
	  // OK, now we know the chunk exists. We need to send a message containing first one bool with the value true, and then the chunk data.
	  // FIXME check that chunkTypeID is correct here? (old check was removed)
	  size_t chunkDataSize = objPtr->getSize();
	  int TMP_TAG_Chunk_data = getIntFromEndOfBuffer(&buf[0], msgSize);
	  size_t msgSizeToSend = sizeof(bool) + chunkDataSize;
	  AsyncMsgStruct* sendMsgStruct = new AsyncMsgStruct(msgSizeToSend, rank, TMP_TAG_Chunk_data, comm, memoryBufferPool);
	  bool theTrueBool = true;
	  memcpy(sendMsgStruct->bufPtr, &theTrueBool, sizeof(bool));
	  cht::timer timer_writeToBuffer;
	  objPtr->writeToBuffer( &sendMsgStruct->bufPtr[sizeof(bool)] , chunkDataSize );
	  LockMutex();
	  if(msgSizeToSend < smallMsgSizeLimitBytes)
	    getChunkReq_small_writeToBuffer_stats.add_timings(timer_writeToBuffer);
	  else 
	    getChunkReq_large_writeToBuffer_stats.add_timings(timer_writeToBuffer);
	  UnlockMutex();
	  PostAsyncSendMessage(sendMsgStruct);
	}
      }
      else if(tag == tagManager.getPermanentTag("TAG_Create_small_chunk")) {
	int chunkTypeID;
	memcpy(&chunkTypeID, buf, sizeof(int));
	std::string class_id_str = getChunkTypeIDStr(chunkTypeID);
	assert(msgSize > 2*sizeof(int)+sizeof(ChunkID));
	ChunkID cid_received;
	memcpy(&cid_received, &buf[sizeof(int)], sizeof(ChunkID));
	// We pass cid_received to createChunkObjectFromBuffer, so if cid_received is set it is used, otherwise a new ChunkID will be created.
	ChunkID cid = createChunkObjectFromBuffer(&buf[sizeof(int)+sizeof(ChunkID)], 
						  msgSize - 2*sizeof(int) - sizeof(ChunkID),
						  class_id_str,
						  cid_received);
	int TMP_TAG_Write_data_to_chunk_ack = getIntFromEndOfBuffer(&buf[0], msgSize);
	// We check the TMP_TAG_Write_data_to_chunk_ack in the message and use it if it was specified, otherwise use general ack tag.
	int tagToUseForAck = TMP_TAG_Write_data_to_chunk_ack;
	if(tagToUseForAck == -1)
	  tagToUseForAck = tagManager.getPermanentTag("TAG_Create_chunk_ack");
	// The ack message should contain a ChunkID and an int (tmp tag unused in this case, we set it to -1).
	int unusedInt = -1;
	AsyncMsgStruct* sendMsgStruct = new AsyncMsgStruct(sizeof(ChunkID)+sizeof(int), rank, tagToUseForAck, comm, memoryBufferPool);
	memcpy(&sendMsgStruct->bufPtr[              0], &cid      , sizeof(ChunkID));
	memcpy(&sendMsgStruct->bufPtr[sizeof(ChunkID)], &unusedInt, sizeof(int));
	PostAsyncSendMessage(sendMsgStruct);
      }
      else if(tag == tagManager.getPermanentTag("TAG_Create_large_chunk_prep")) {
	// The TAG_Create_large_chunk_prep message contains 2 int values: the tag and the size of the large message that will come.
	assert(msgSize == 2*sizeof(int));
	int tagForLargeChunkMsg, sizeOfLargeChunkMsg;
	memcpy(&tagForLargeChunkMsg, &buf[          0], sizeof(int));
	memcpy(&sizeOfLargeChunkMsg, &buf[sizeof(int)], sizeof(int));
	assert(sizeOfLargeChunkMsg > 0);
	// OK, now we know the tag and the size of the large message. We post a nonblocking receive for that and take care of the result later.
	AsyncMsgStruct* recvMsgStruct = new AsyncMsgStruct(sizeOfLargeChunkMsg, rank, tagForLargeChunkMsg, comm, memoryBufferPool);
	PostAsyncRecvMessage(recvMsgStruct);
      }
      else if(tag == tagManager.getPermanentTag("TAG_Delete_chunk")) {
	ChunkID cid = GetChunkIdFromBuffer(&buf[0]);
	if(cid.ownerRank != key.my_rank())
	  throw std::runtime_error("Worker::handle_request, TAG_Delete_chunk but (cid.ownerRank != key.my_rank()).");
	deleteChunk(cid);
      }
      else if(tag == tagManager.getPermanentTag("TAG_Copy_chunk")) {
	assert(msgSize == 2*sizeof(ChunkID));
	const char* p = &buf[0];
	ChunkID cid_org, cid_new;
	memcpy(&cid_org, p, sizeof(ChunkID)); p += sizeof(ChunkID);
	memcpy(&cid_new, p, sizeof(ChunkID)); p += sizeof(ChunkID);
	if(cid_org.ownerRank != key.my_rank())
	  throw std::runtime_error("Worker::handle_request, TAG_Copy_chunk but (cid_org.ownerRank != key.my_rank()).");
	std::string class_id_str = getChunkTypeIDStr(cid_org.chunkTypeID);
	createChunkStruct(class_id_str, cid_new.size, cid_new);
	CopyChunkLocally(cid_org, cid_new);
	AsyncMsgStruct* sendMsgStruct = new AsyncMsgStruct(sizeof(ChunkID), rank, tagManager.getPermanentTag("TAG_Copy_chunk_ack"), comm, memoryBufferPool);
	memcpy(&sendMsgStruct->bufPtr[0], &cid_new, sizeof(ChunkID));
	PostAsyncSendMessage(sendMsgStruct);
      }
      else if(tag == tagManager.getPermanentTag("TAG_Create_chunk_ack")) {
	// This message is expected to contain a ChunkID (for the created chunk) and an int (a temporary tag that should be released).
	assert(msgSize == sizeof(ChunkID)+sizeof(int));
	const char* p = &buf[0];
	ChunkID cid;
	memcpy(&cid, p, sizeof(ChunkID)); p += sizeof(ChunkID);
	int tmpTag;
	memcpy(&tmpTag, p, sizeof(int)); p += sizeof(int);
	LockMutex();
	// Check that cid exists in incompleteCreateOperations list
	if(incompleteRegisterChunkOperations.find(cid) == incompleteRegisterChunkOperations.end()) {
	  std::cerr << "Error in handle_request, TAG_Create_chunk_ack received but cid not found in incompleteRegisterChunkOperations list." << std::endl;
	  std::cerr << "cid = " << cid.str() << std::endl;
	  throw std::runtime_error("Error in handle_request, TAG_Create_chunk_ack received but cid not found in incompleteRegisterChunkOperations list.");
	}
	incompleteRegisterChunkOperations.erase(cid);
	UnlockMutex();
	if(tmpTag >= 0)
	  tagManager.releaseTemporaryTag(tmpTag);
      }
      else if(tag == tagManager.getPermanentTag("TAG_Copy_chunk_ack")) {
	assert(msgSize == sizeof(ChunkID));
	ChunkID cid = GetChunkIdFromBuffer(&buf[0]);
	LockMutex();
	// Check that cid exists in incompleteCopyChunkOperations list
	if(incompleteCopyChunkOperations.find(cid) == incompleteCopyChunkOperations.end())
	  throw std::runtime_error("Error in handle_request, TAG_Copy_chunk_ack received but cid not found in incompleteCopyChunkOperations list.");
	incompleteCopyChunkOperations.erase(cid);
	UnlockMutex();
      }

      else if(tag == tagManager.getPermanentTag("TAG_Time_to_stop")) {
	if(comm != key.comm_to_parent())
	  throw std::runtime_error("Error: TAG_Time_to_stop message received, but not from parent.");
	timeToStopFlag = true;
      }
      else {
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::handle_request: unexpected tag received.");
      }
    }




    void Worker::PostReceiveOperations(RecvMsgBufStruct & bufs,
				       MPI_Comm* comm,
				       MPI_Request* requests,
				       int rank) {
      // Number of MPI_Irecv calls here is supposed to be nReceiveTags. Use assert here to make sure.
      assert(nReceiveTags == 8);
      MW._MPI_Irecv(&bufs.buf_Delete_chunk[0], 
		   bufs.buf_Delete_chunk.size(), MPI_UNSIGNED_CHAR, rank,
		   tagManager.getPermanentTag("TAG_Delete_chunk"),
		   *comm, &requests[0]);
      MW._MPI_Irecv(&bufs.buf_Create_small_chunk[0], 
		   bufs.buf_Create_small_chunk.size(), MPI_UNSIGNED_CHAR, rank,
		   tagManager.getPermanentTag("TAG_Create_small_chunk"),
		   *comm, &requests[1]);
      MW._MPI_Irecv(&bufs.buf_Create_large_chunk_prep[0], 
		   bufs.buf_Create_large_chunk_prep.size(), MPI_UNSIGNED_CHAR, rank,
		   tagManager.getPermanentTag("TAG_Create_large_chunk_prep"),
		   *comm, &requests[2]);
      MW._MPI_Irecv(&bufs.buf_Get_chunk_data[0], 
		   bufs.buf_Get_chunk_data.size(), MPI_UNSIGNED_CHAR, rank,
		   tagManager.getPermanentTag("TAG_Get_chunk_data"),
		   *comm, &requests[3]);
      MW._MPI_Irecv(NULL, 
		   0, MPI_UNSIGNED_CHAR, rank,
		   tagManager.getPermanentTag("TAG_Time_to_stop"),
		   *comm, &requests[4]);
      MW._MPI_Irecv(&bufs.buf_Copy_chunk[0], 
		   bufs.buf_Copy_chunk.size(), MPI_UNSIGNED_CHAR, rank,
		   tagManager.getPermanentTag("TAG_Copy_chunk"),
		   *comm, &requests[5]);
      MW._MPI_Irecv(&bufs.buf_Copy_chunk_ack[0], 
		   bufs.buf_Copy_chunk_ack.size(), MPI_UNSIGNED_CHAR, rank,
		   tagManager.getPermanentTag("TAG_Copy_chunk_ack"),
		   *comm, &requests[6]);
      MW._MPI_Irecv(&bufs.buf_Create_chunk_ack[0], 
		   bufs.buf_Create_chunk_ack.size(), MPI_UNSIGNED_CHAR, rank,
		   tagManager.getPermanentTag("TAG_Create_chunk_ack"),
		   *comm, &requests[7]);
    }

    char* Worker::GetMessageBuf(RecvMsgBufStruct & bufs, 
				int i, 
				int tag_in,
				int & bufSz) {
      cht::vector<char>* buf;
      int tag = -1;
      switch(i) {
      case 0:
	tag = tagManager.getPermanentTag("TAG_Delete_chunk");
	buf = &bufs.buf_Delete_chunk;
	break;
      case 1:
	tag = tagManager.getPermanentTag("TAG_Create_small_chunk");
	buf = &bufs.buf_Create_small_chunk;
	break;
      case 2:
	tag = tagManager.getPermanentTag("TAG_Create_large_chunk_prep");
	buf = &bufs.buf_Create_large_chunk_prep;
	break;
      case 3:
	tag = tagManager.getPermanentTag("TAG_Get_chunk_data");
	buf = &bufs.buf_Get_chunk_data;
	break;
      case 4:
	tag = tagManager.getPermanentTag("TAG_Time_to_stop");
	buf = NULL;
	break;
      case 5:
	tag = tagManager.getPermanentTag("TAG_Copy_chunk");
	buf = &bufs.buf_Copy_chunk;
	break;
      case 6:
	tag = tagManager.getPermanentTag("TAG_Copy_chunk_ack");
	buf = &bufs.buf_Copy_chunk_ack;
	break;
      case 7:
	tag = tagManager.getPermanentTag("TAG_Create_chunk_ack");
	buf = &bufs.buf_Create_chunk_ack;
	break;
      default:
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::GetMessageInfo.");
      }
      if(tag != tag_in)
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::GetMessageInfo(): wrong tag.");
      if(buf) {
	bufSz = buf->size();
	return &((*buf)[0]);
      }
      bufSz = 0;
      return NULL;
    }

    void Worker::HandleCompletedLargeReceiveOperation(AsyncMsgStruct* recvMsgStructPtr) {
      int chunkTypeID;
      memcpy(&chunkTypeID, recvMsgStructPtr->bufPtr, sizeof(int));
      std::string class_id_str;
      try {
	class_id_str = getChunkTypeIDStr(chunkTypeID);
      }
      catch(std::runtime_error e) {
	std::cerr << "HandleCompletedLargeReceiveOperation: error in getChunkTypeIDStr" << std::endl;
	throw e;
      }
      ChunkID cid_received;
      memcpy(&cid_received, &recvMsgStructPtr->bufPtr[sizeof(int)], sizeof(ChunkID));
      // We pass cid_received to createChunkObjectFromBuffer, so if cid_received is set it is used, otherwise a new ChunkID will be created.
      ChunkID cid = createChunkObjectFromBuffer(&recvMsgStructPtr->bufPtr[sizeof(int)+sizeof(ChunkID)], 
						recvMsgStructPtr->bufSz - 2*sizeof(int) - sizeof(ChunkID),
						class_id_str,
						cid_received);
      int TMP_TAG_Write_data_to_chunk_ack = getIntFromEndOfBuffer(recvMsgStructPtr->bufPtr, recvMsgStructPtr->bufSz);
      // We check the TMP_TAG_Write_data_to_chunk_ack in the message and use it if it was specified, otherwise use general ack tag.
      int tagToUseForAck = TMP_TAG_Write_data_to_chunk_ack;
      if(tagToUseForAck == -1)
	tagToUseForAck = tagManager.getPermanentTag("TAG_Create_chunk_ack");
      AsyncMsgStruct* sendMsgStruct = new AsyncMsgStruct(sizeof(ChunkID) + sizeof(int), recvMsgStructPtr->rank, tagToUseForAck, recvMsgStructPtr->comm, memoryBufferPool);
      memcpy(&sendMsgStruct->bufPtr[0], &cid, sizeof(ChunkID));
      memcpy(&sendMsgStruct->bufPtr[sizeof(ChunkID)], &recvMsgStructPtr->tag, sizeof(int));
      PostAsyncSendMessage(sendMsgStruct);
      delete recvMsgStructPtr;
    }

    void Worker::HandleAnyCompletedLargeReceiveOperations() {
      // Go through recvMsgStructList to see if any of those operations have finished.
      std::list<AsyncMsgStruct*> finishedOperations;
      LockMutex();
      std::list<AsyncMsgStruct*>::iterator it = recvMsgStructList.begin();
      while(it != recvMsgStructList.end()) {
	AsyncMsgStruct* msgStructPtr = (*it);
	int flag;
	MPI_Status status;
	MW._MPI_Test(&msgStructPtr->request, &flag, &status);
	if(flag) {
	  // OK, this receive operation is done. We should handle this and remove it from the list.
	  int sourceRank = status.MPI_SOURCE;
	  int tag = status.MPI_TAG;
	  assert(sourceRank == msgStructPtr->rank && tag == msgStructPtr->tag);
	  finishedOperations.push_back(*it);
	  it = recvMsgStructList.erase(it);
	}
	else
	  it++;
      }
      UnlockMutex();      
      for(it = finishedOperations.begin(); it != finishedOperations.end(); it++) {
	AsyncMsgStruct* msgStructPtr = (*it);
	HandleCompletedLargeReceiveOperation(msgStructPtr);
      }
    }

    void Worker::worker_thread_func() {
      AccessKey key(this);
      // Here, we want to wait for the following events:
      // - A message from parent (several tags but not any tag)
      // - A message from another worker (several tags but not any tag)

      // To do this, we post a number of non-blocking receive operations
      // and then call MPI_Wait_some(). The use of MPI_Wait_some() is
      // supposed to be good because it helps avoid "starvation".

      // The nRecvsForWorkers constant determines how many simultaneous
      // Irecv() calls are made. A large number here should help avoid
      // "starvation".
      const int nRecvsForWorkers = 10;
      const int nRequestsTot = nReceiveTags * (1 + nRecvsForWorkers);
      MPI_Request requests[nRequestsTot];
  
      // First take care of non-blocking receive operations for messages
      // from parent.
      RecvMsgBufStruct parentBufs;  
      PostReceiveOperations(parentBufs, key.comm_to_parent(), &requests[0], 0);

      // Now take care of non-blocking receive operations for messages
      // from other workers.  
      RecvMsgBufStruct workerBufs[nRecvsForWorkers];
      for(int i = 0; i < nRecvsForWorkers; i++)
	PostReceiveOperations(workerBufs[i], 
			      key.comm_worker_world(), 
			      &requests[nReceiveTags + i*nReceiveTags],
			      MPI_ANY_SOURCE);
  

      cht::timer timer_loop; // Timer to measure loop time
      bool first_loop = true;
      bool timeToStopFlag = false;
      while(1) {
	// Wait for some kind of message(s) from either parent or other
	// workers.
	MPI_Status  statuses[nRequestsTot];
	int         indices [nRequestsTot];
	int ndone = -1;
	if (first_loop) {
	  first_loop = false;
	}
	else {
	  LockMutex();
	  worker_thread_func_loop_stats.add_timings(timer_loop);
	  UnlockMutex();
	}

	usleep(2000); // FIXME: how long should this thread sleep here?
	
	/* Now we want to wait until some kind of message has
	   arrived. Earlier this was done using MPI_Waitsome but that
	   is not possible now since we also have a list of
	   nonblocking receive operations that we also need to
	   check. So instead we do MPI_Testsome first and then check
	   the nonblocking receive operations one by one. */
	MW._MPI_Testsome( nRequestsTot, requests, &ndone, indices, statuses );

	timer_loop.reset();
	for (int i = 0; i < ndone; i++) {
	  int msgSize;
	  MW._MPI_Get_count(&statuses[i], MPI_UNSIGNED_CHAR, &msgSize);
	  int tag = statuses[i].MPI_TAG;
	  int sourceRank = statuses[i].MPI_SOURCE;
	  int j = indices[i];
	  if(j < nReceiveTags) {
	    // parent
	    int bufSz = -1;
	    char* bufPtr = GetMessageBuf(parentBufs, j, tag, bufSz);
	    cht::timer timer_handle_request;
	    handle_request(bufPtr, msgSize, tag, sourceRank, key.comm_to_parent(), timeToStopFlag);
	    LockMutex();
	    worker_thread_func_handle_request_parent_stats.add_timings(timer_handle_request);
	    UnlockMutex();
	    cht::timer timer_irecv;
	    MW._MPI_Irecv(bufPtr, bufSz, MPI_UNSIGNED_CHAR, 0, tag, *key.comm_to_parent(), &requests[j] );
	    LockMutex();
	    worker_thread_func_irecv_stats.add_timings(timer_irecv);
	    UnlockMutex();
	  }
	  else {
	    // worker
	    if(sourceRank == key.my_rank())
	      throw std::runtime_error("Error: (sourceRank == my_rank) for worker!!");
	    int idxInList = (j - nReceiveTags) / nReceiveTags;
	    int jj = j % nReceiveTags;
	    int bufSz = -1;
	    char* bufPtr = GetMessageBuf(workerBufs[idxInList], jj, tag, bufSz);
	    cht::timer timer_handle_request;
	    handle_request(bufPtr, msgSize, tag, sourceRank, key.comm_worker_world(), timeToStopFlag);
	    LockMutex();
	    worker_thread_func_handle_request_worker_stats.add_timings(timer_handle_request);
	    UnlockMutex();
	    cht::timer timer_irecv;
	    MW._MPI_Irecv(bufPtr, bufSz, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, 
			 tag, *key.comm_worker_world(), &requests[j] );
	    LockMutex();
	    worker_thread_func_irecv_stats.add_timings(timer_irecv);
	    UnlockMutex();
	  } // end else worker
	} // end for i
	// Check if any large receive operations have finished.
	HandleAnyCompletedLargeReceiveOperations();
	// Check if any send operations have finished so that we can free
	// memory for their buffers.
	FreeMemoryForCompletedSendOperations();
	if(timeToStopFlag) {
          LockMutex();
          int noOfChunks = chunkList.size();
          UnlockMutex();
          if(noOfChunks == 0)
	    break;
        }
      } // end while
      // Cleanup.
      while(!sendMsgStructListIsEmpty()) {
	usleep(10000);
	FreeMemoryForCompletedSendOperations();
      }
      for(int i = 0; i < nRequestsTot; i++) {
	MW._MPI_Cancel(&requests[i]);
	MW._MPI_Request_free(&requests[i]);
      }
    } // end Worker::worker_thread_func




    ChunkID Worker::createChunkStruct(std::string class_id_str, size_t size, ChunkID known_cid) {
      AccessKey key(this);
      ChunkStruct* chunk = new ChunkStruct;
      int chunkTypeID = getChunkTypeIDInt(class_id_str);
      chunk->id = known_cid;
      LockMutex(); // Must lock mutex before accessing idCounter
      if(chunk->id == CHUNK_ID_NULL) {
	// Generate new chunk id
	chunk->id.ownerRank = key.my_rank();
	chunk->id.creatorRank = key.my_rank();
	chunk->id.size = size;
	chunk->id.chunkTypeID = chunkTypeID;
	chunk->id.creatorLocalID = ++idCounter;
      }
      else {
	// Verify that given ChunkID seems reasonable.
	assert(chunk->id.ownerRank == key.my_rank());
	assert(chunk->id.chunkTypeID == chunkTypeID);
	assert(chunk->id.size == size);
      }
      chunk->class_id_str = class_id_str; // This somehow duplicates
					  // info in chunk ID (where
					  // type ID is stored as an
					  // integer). Perhaps
					  // class_id_str should be
					  // removed from the Chunk
					  // class?
      chunkList.insert(ChunkListMap::value_type(chunk->id, chunk));
      UnlockMutex();
      return chunk->id;
    }

    void Worker::CopyChunkLocally(ChunkID cid_org, ChunkID cid_new) {
      ChunkStruct* theOrgChunk = FindChunkInList(cid_org, true);
      ChunkStruct* theNewChunk = FindChunkInList(cid_new, false);
      theNewChunk->objPtr = theOrgChunk->objPtr;
      LockMutex();
      chunkCounters[&(*theNewChunk->objPtr)]++;
      UnlockMutex();
    }

    ChunkID Worker::getIdForCopyChunk(ChunkID cid) {
      AccessKey key(this);
      ChunkID cid_new;
      cid_new.ownerRank = cid.ownerRank;
      cid_new.creatorRank = key.my_rank();
      cid_new.size = cid.size;
      cid_new.chunkTypeID = cid.chunkTypeID;
      LockMutex(); // Must lock mutex before accessing idCounter
      cid_new.creatorLocalID = ++idCounter;
      incompleteCopyChunkOperations.insert(cid_new);
      UnlockMutex();
      return cid_new;
    }

    void Worker::startCopyChunk(ChunkID cid_old, ChunkID cid_new) {
      AccessKey key(this);
      assert(cid_old.ownerRank == cid_new.ownerRank);
      assert(cid_old.chunkTypeID == cid_new.chunkTypeID);
      assert(cid_old.size == cid_new.size);
      if(cid_old.ownerRank != key.my_rank()) {
	// The original chunk is owned by another worker.
	// Send a message to the owner worker, saying this chunk is supposed to be copied.
	AsyncMsgStruct* sendMsgStruct = new AsyncMsgStruct(2*sizeof(ChunkID), cid_new.ownerRank, tagManager.getPermanentTag("TAG_Copy_chunk"), key.comm_worker_world(), memoryBufferPool);
	char* p = &sendMsgStruct->bufPtr[0];
	memcpy(p, &cid_old, sizeof(ChunkID)); p += sizeof(ChunkID);
	memcpy(p, &cid_new, sizeof(ChunkID)); p += sizeof(ChunkID);
	// Note: we need to insert cid_new into
	// incompleteCopyChunkOperations list before sending the message,
	// otherwise the response ack message could come before we
	// inserted it, then the ack could not be handled properly.
	PostAsyncSendMessage(sendMsgStruct);
	return;
      }
      // Now we know the ownerRank is the same as my_rank. This means we can do the copy locally, no communication needed.
      std::string class_id_str = getChunkTypeIDStr(cid_old.chunkTypeID);
      createChunkStruct(class_id_str, cid_new.size, cid_new);
      CopyChunkLocally(cid_old, cid_new);
      LockMutex();
      incompleteCopyChunkOperations.erase(cid_new);
      UnlockMutex();
      return;
    }

    bool Worker::isCopyChunkStillInProgress(ChunkID cid) {
      LockMutex();
      bool found = false;
      if(incompleteCopyChunkOperations.find(cid) != incompleteCopyChunkOperations.end())
	found = true;
      UnlockMutex();
      return found;
    }

    ChunkID Worker::createChunkObjectFromBuffer(char const * dataBuffer, size_t const bufferSize, 
						std::string const class_id_str, ChunkID known_cid) {
      AccessKey key(this);
      // This chunk is stored locally by this process.
      ChunkID cid = createChunkStruct(class_id_str, bufferSize, known_cid);
      ChunkStruct* theChunk = FindChunkInList(cid, false);
      Chunk* newObjPtr = cht::obj_factory<Chunk>::instance().createObject(class_id_str);
      newObjPtr->assignFromBuffer(dataBuffer, bufferSize);
      theChunk->objPtr = newObjPtr;
      LockMutex();
      chunkCounters[newObjPtr]++;
      totDataSizeBytes += theChunk->objPtr->memoryUsage();
      totDataSizeBytesMax = (totDataSizeBytes > totDataSizeBytesMax) ?  totDataSizeBytes : totDataSizeBytesMax;
      UnlockMutex();
      return cid;
    }



    void Worker::deleteChunk(ChunkID cid) {
      AccessKey key(this);
      // Check if cid is "null"
      if(cid == CHUNK_ID_NULL)
	return;
      if(cid.ownerRank == key.my_rank()) {
	// This chunk is stored locally by this process.
	cht::timer timer_localwork_1; // do timing in two parts since we may call this routine recursively.
	ChunkStruct* theChunk = TryToFindChunkInList(cid, true);
	if(theChunk == NULL)
	  throw std::runtime_error("Error in ChunkObjService::Worker::deleteChunk: TryToFindChunkInList returned NULL. Trying to delete a chunk that was already deleted?");
	// Now we have the chunk. Before deleting it, call
	// getChildChunks function to make sure any "child chunks" are
	// also deleted. Note that this needs to be done before the
	// mutex lock since deleteChunk() (this routine) is called
	// recursively.
	// FIXME: think about risk of problems with using theChunk here without mutex being locked. Is that a problem?
	// Delete child chunks only if this is the last ChunkID referring to that object.
	LockMutex();	
	chunkCounters[&(*theChunk->objPtr)]--;
	bool thisWasFinalDelete = false;
	if(chunkCounters[&(*theChunk->objPtr)] == 0)
	  thisWasFinalDelete = true;
	deleteChunk_localwork_stats.add_timings(timer_localwork_1);
	UnlockMutex();
	if(thisWasFinalDelete) {
	  // Get list of child chunks and call deleteChunk for each of them.
	  std::list<ChunkID> childChunksIDs;
	  theChunk->objPtr->getChildChunks(childChunksIDs);
	  std::list<ChunkID>::iterator it;
	  for(it = childChunksIDs.begin(); it != childChunksIDs.end(); it++)
	    deleteChunk(*it);
	}
	cht::timer timer_localwork_2;
	LockMutex();
	size_t chunkMemUsage = theChunk->objPtr->memoryUsage();
	if(chunkMemUsage > totDataSizeBytes) {
	  std::cerr << "Error in cht::ChunkObjService::Worker::deleteChunk: (chunkMemUsage > totDataSizeBytes)." << std::endl
		    << "chunkMemUsage    = " << chunkMemUsage    << std::endl
		    << "totDataSizeBytes = " << totDataSizeBytes << std::endl;
	  throw std::runtime_error("Error in cht::ChunkObjService::Worker::deleteChunk: (chunkMemUsage > totDataSizeBytes).");
	}
	if(thisWasFinalDelete)
	  totDataSizeBytes -= chunkMemUsage;
	delete theChunk;
	chunkList.erase(cid);
	deleteChunk_localwork_stats.add_timings(timer_localwork_2);
	UnlockMutex();
	return;
      }
      /* 
	 ELIAS NOTE 2011-09-10 (1): There was a problem with "thread not
	 in any thread group" earlier, that happened when deleteChunk
	 was called for a chunk that had been created on another
	 worker. To avoid that problem, the code was changed so that
	 startingBlockingOperation() and blockingOperationFinished()
	 calls are only done if the thread belongs to any thread
	 group.
      */
      /* 
	 ELIAS NOTE 2011-09-10 (2): Now it turned out that the
	 blocking communication here anyway had to be removed since it
	 caused deadlock when deleteChunk() was called for a chunk
	 that has child-chunks on other workers and those child-chunks
	 have child-chunks on the first worker, or something like
	 that. So now the note (1) above is already obsolete, and the
	 threadBelongsToThreadGroup() function is no longer used.
      */
      // if(Threads::threadBelongsToThreadGroup())
      //	Threads::startingBlockingOperation();
      AsyncMsgStruct* sendMsgStruct = new AsyncMsgStruct(sizeof(ChunkID), cid.ownerRank, tagManager.getPermanentTag("TAG_Delete_chunk"), key.comm_worker_world(), memoryBufferPool);
      memcpy(&sendMsgStruct->bufPtr[0], &cid, sizeof(ChunkID));
      PostAsyncSendMessage(sendMsgStruct);
      // if(Threads::threadBelongsToThreadGroup())
      //	Threads::blockingOperationFinished();
    }


    void Worker::resetStatistics() {
      LockMutex();
      getChunk_total_count = 0;
      getChunk_found_locally_count = 0;
      getChunkIfLocal_total_count = 0;
      getChunkIfLocal_found_locally_count = 0;
      getChunk_small_communication_stats.clear();
      getChunk_small_assignFromBuffer_stats.clear();
      getChunk_large_communication_stats.clear();
      getChunk_large_assignFromBuffer_stats.clear();
      getChunkReq_small_writeToBuffer_stats.clear();
      getChunkReq_large_writeToBuffer_stats.clear();
      writeChunkObj_communication_stats.clear();
      writeChunkObj_writeToBuffer_stats.clear();
      deleteChunk_communication_stats.clear();
      deleteChunk_localwork_stats.clear();
      getChunk_small_totNoOfBytesReceived = 0;
      getChunk_large_totNoOfBytesReceived = 0;
      worker_thread_func_loop_stats.clear(); 
      worker_thread_func_handle_request_parent_stats.clear(); 
      worker_thread_func_handle_request_worker_stats.clear(); 
      FreeMemoryForCompletedSendOperations_stats.clear(); 
      worker_thread_func_irecv_stats.clear(); 
      UnlockMutex();
    }



    static void output_stats(std::stringstream & ss, const cht::work_statistics & stats, const char* name) {
      const int w = 14;
      const int p = 6;
      ss << "'" << name << "', count = " << stats.count << "\n";
      ss << "   wall       : "
	 << std::fixed << std::setw(w) << std::setprecision(p) << stats.wall_seconds_tot << "\n";
      ss << "   min        : "
	 << std::fixed << std::setw(w) << std::setprecision(p) << stats.wall_seconds_min << "\n";
      ss << "   max        : "
	 << std::fixed << std::setw(w) << std::setprecision(p) << stats.wall_seconds_max << "\n";
      if(stats.count > 0) 
	ss << "   average    : "
	   << std::fixed << std::setw(w) << std::setprecision(p) << stats.wall_seconds_tot / stats.count << "\n";
      else
	ss << "   average    : " << "\n";
    }

    void Worker::reportStatistics(std::string messageHeader) {
      std::stringstream ss;
      ss << messageHeader << "\n";
      LockMutex();
      ss << "getChunk (total)                = "
	 << std::setw(15) << getChunk_total_count <<"\n";
      ss << "getChunk (found locally)        = "
	 << std::setw(15) << getChunk_found_locally_count <<"\n";    
      ss << "getChunkIfLocal (total)         = "
	 << std::setw(15) << getChunkIfLocal_total_count <<"\n";    
      ss << "getChunkIfLocal (found locally) = "
	 << std::setw(15) << getChunkIfLocal_found_locally_count <<"\n";    
      ss << "getChunk_small_totNoOfBytesReceived = " 
	 << std::setw(15) << getChunk_small_totNoOfBytesReceived 
	 << " = " << std::fixed << std::setw(10) << std::setprecision(5)
	 << (double)getChunk_small_totNoOfBytesReceived/1000000000 << " GB\n";
      ss << "getChunk_large_totNoOfBytesReceived = " 
	 << std::setw(15) << getChunk_large_totNoOfBytesReceived 
	 << " = " << std::fixed << std::setw(10) << std::setprecision(5)
	 << (double)getChunk_large_totNoOfBytesReceived/1000000000 << " GB\n";
      output_stats(ss, getChunk_small_communication_stats, "getChunk_small_communication_stats");
      output_stats(ss, getChunk_small_assignFromBuffer_stats, "getChunk_small_assignFromBuffer_stats");
      output_stats(ss, getChunk_large_communication_stats, "getChunk_large_communication_stats");
      output_stats(ss, getChunk_large_assignFromBuffer_stats, "getChunk_large_assignFromBuffer_stats");
      output_stats(ss, getChunkReq_small_writeToBuffer_stats, "getChunkReq_small_writeToBuffer_stats");
      output_stats(ss, getChunkReq_large_writeToBuffer_stats, "getChunkReq_large_writeToBuffer_stats");
      output_stats(ss, writeChunkObj_communication_stats, "writeChunkObj_communication_stats");
      output_stats(ss, writeChunkObj_writeToBuffer_stats, "writeChunkObj_writeToBuffer_stats");
      output_stats(ss, deleteChunk_communication_stats, "deleteChunk_communication_stats");
      output_stats(ss, deleteChunk_localwork_stats, "deleteChunk_localwork_stats");
      output_stats(ss, worker_thread_func_loop_stats, "worker_thread_func_loop_stats");
      output_stats(ss, worker_thread_func_handle_request_parent_stats, "worker_thread_func_handle_request_parent_stats");
      output_stats(ss, worker_thread_func_handle_request_worker_stats, "worker_thread_func_handle_request_worker_stats");
      output_stats(ss, FreeMemoryForCompletedSendOperations_stats, "FreeMemoryForCompletedSendOperations_stats");
      output_stats(ss, worker_thread_func_irecv_stats, "worker_thread_func_irecv_stats");
      UnlockMutex();
      ss << "=========================================================";
      OS.outputInfo(ss.str());
    }

    ChunkID Worker::getIdForRegisterChunk(Chunk const * objPtr,
					std::string class_id_str) {
      AccessKey key(this);
      // Check if chunk should be created locally.
      double rand_number_0_to_1 = ((double)rand()) / RAND_MAX;
      bool createLocally = false;
      if(rand_number_0_to_1 < probabilityToRegisterChunksLocally || key.n_workers() == 1)
	createLocally = true;
      // Generate new chunk id
      ChunkID cid_new;
      if(createLocally)
	cid_new.ownerRank = key.my_rank();
      else {
	// Select random other worker rank to set as ownerRank.
	int nOtherWorkers = key.n_workers() - 1;
	int randomRank = (int)((((double)rand())*nOtherWorkers) / RAND_MAX);
	assert(randomRank >= 0 && randomRank < nOtherWorkers);
	if(randomRank == key.my_rank())
	  randomRank = key.n_workers() - 1;
	cid_new.ownerRank = randomRank;
      }
      cid_new.creatorRank = key.my_rank();
      cid_new.size = objPtr->getSize();
      cid_new.chunkTypeID = getChunkTypeIDInt(class_id_str);
      // Look at child chunk ids to determine childDepth, totalDeepSize, etc.
      std::list<ChunkID> childChunkIDs;
      objPtr->getChildChunks(childChunkIDs);
      int nChildChunks = childChunkIDs.size();
      int maxChildDepthForChildren = 0;
      size_t totalChunkCountForChildren = 0;
      size_t totalDeepSizeForChildren = 0;
      ChunkID::MajorOwnerInfo majorOwners_all[nChildChunks*ChunkID::MAX_NO_OF_OWNERS_IN_LIST + 1];
      int majorOwners_count = 0;
      for(std::list<ChunkID>::const_iterator it = childChunkIDs.begin(); it != childChunkIDs.end(); it++) {
	if(it->childDepth > maxChildDepthForChildren)
	  maxChildDepthForChildren = it->childDepth;
	totalChunkCountForChildren += it->totalChunkCount;
	totalDeepSizeForChildren += it->totalDeepSize;
	for(int i = 0; i < ChunkID::MAX_NO_OF_OWNERS_IN_LIST; i++) {
	  if(it->majorOwners[i].rank != -1) {
	    majorOwners_all[majorOwners_count] = it->majorOwners[i];
	    majorOwners_count++;
	  }
	}
      }
      majorOwners_all[majorOwners_count].rank = cid_new.ownerRank;
      majorOwners_all[majorOwners_count].size = cid_new.size;
      majorOwners_count++;
      if(nChildChunks == 0)
	cid_new.childDepth = 0;
      else
	cid_new.childDepth = 1 + maxChildDepthForChildren;
      cid_new.totalChunkCount = 1 + totalChunkCountForChildren;
      cid_new.totalDeepSize = objPtr->getSize() + totalDeepSizeForChildren;
      // Get list of unique "major owners"
      ChunkID::MajorOwnerInfo majorOwners_unique[majorOwners_count];
      int majorOwners_unique_count = 0;
      for(int i = 0; i < majorOwners_count; i++) {
	int foundIndex = -1;
	for(int j = 0; j < majorOwners_unique_count; j++) {
	  if(majorOwners_unique[j].rank == majorOwners_all[i].rank)
	    foundIndex = j;
	}
	if(foundIndex == -1) {
	  majorOwners_unique[majorOwners_unique_count] = majorOwners_all[i];
	  majorOwners_unique_count++;
	}
	else
	  majorOwners_unique[foundIndex].size += majorOwners_all[i].size;
      }
      // Sort list of unique "major owners" so that the largest sizes come first (bubble sort)
      for(int i = majorOwners_unique_count-1; i >= 0; i--) {
	for(int j = 0; j < i; j++) {
	  if(majorOwners_unique[j].size < majorOwners_unique[j+1].size) {
	    ChunkID::MajorOwnerInfo tmp = majorOwners_unique[j];
	    majorOwners_unique[j] = majorOwners_unique[j+1];
	    majorOwners_unique[j+1] = tmp;
	  }
	}
      }
      // Verify that list is sorted
      for(int i = 0; i < majorOwners_unique_count-1; i++)
	assert(majorOwners_unique[i].size >= majorOwners_unique[i+1].size);
      // Now set cid_new.majorOwners
      for(int i = 0; i < ChunkID::MAX_NO_OF_OWNERS_IN_LIST; i++) {
	if(i < majorOwners_unique_count)
	  cid_new.majorOwners[i] = majorOwners_unique[i];
      }
      // Set cid_new.creatorLocalID using idCounter
      LockMutex(); // Must lock mutex before accessing idCounter
      cid_new.creatorLocalID = ++idCounter;
      incompleteRegisterChunkOperations.insert(cid_new);
      UnlockMutex();
      return cid_new;
    }

    void Worker::startRegisterChunk(Chunk const * objPtr, ChunkID cid) {
      AccessKey key(this);
      if(cid.ownerRank == key.my_rank()) {
	// Create chunk locally.
	std::string class_id_str = getChunkTypeIDStr(cid.chunkTypeID);
	assert(class_id_str.length() > 0);
	createChunkStruct(class_id_str, objPtr->getSize(), cid);
	ChunkStruct* theChunk = FindChunkInList(cid, false);
	LockMutex();
	assert(theChunk->objPtr == 0);
	theChunk->objPtr = objPtr;
	chunkCounters[&(*theChunk->objPtr)]++;
	totDataSizeBytes += theChunk->objPtr->memoryUsage();
	totDataSizeBytesMax = (totDataSizeBytes > totDataSizeBytesMax) ?  totDataSizeBytes : totDataSizeBytesMax;
	incompleteRegisterChunkOperations.erase(cid);
	UnlockMutex();
      }
      else {
	// Create chunk on other worker.
	// Send a message to the owner worker, saying this chunk is supposed to be created.
	int chunkTypeID = cid.chunkTypeID;
	// Create a buffer where the chunkTypeID is first, then the ChunkID, then the data, then an unused int tag which we set to -1.
	size_t objDataSz = objPtr->getSize();
	int bufSz = sizeof(int) + sizeof(ChunkID) + objDataSz + sizeof(int);
	int TMP_TAG_Chunk_data = tagManager.getTemporaryTag_outgoing();
	AsyncMsgStruct* sendMsgStruct = new AsyncMsgStruct(bufSz, cid.ownerRank, TMP_TAG_Chunk_data, key.comm_worker_world(), memoryBufferPool);
	memcpy(&sendMsgStruct->bufPtr[0], &chunkTypeID, sizeof(int));
	memcpy(&sendMsgStruct->bufPtr[sizeof(int)], &cid, sizeof(ChunkID));
	objPtr->writeToBuffer(&sendMsgStruct->bufPtr[sizeof(int)+sizeof(ChunkID)], objDataSz);
	// The receiver expects an int with -1 when a tmp ack tag is not used.
	int tagNotUsed = -1;
	int bufSzExceptTag = bufSz - sizeof(int);
	memcpy(&sendMsgStruct->bufPtr[bufSzExceptTag], &tagNotUsed, sizeof(int));
	// Note: we need to insert cid into
	// incompleteRegisterChunkOperations list before sending the message,
	// otherwise the response ack message could come before we
	// inserted it, then the ack could not be handled properly.
	if(bufSz > MAX_WRITECHUNKOBJ_SMALL_MESSAGE_SIZE) {
	  // Too big: use two messages.  
	  // First send a TAG_Create_large_chunk_prep message containing 2 int values: the tag and the size of the large message that will come.
	  int intsToSend[2];
	  intsToSend[0] = TMP_TAG_Chunk_data;
	  intsToSend[1] = bufSz;
	  MW._MPI_Send(&intsToSend[0], 2*sizeof(int), MPI_UNSIGNED_CHAR, cid.ownerRank, 
		      tagManager.getPermanentTag("TAG_Create_large_chunk_prep"), *key.comm_worker_world());
	  // Now send the "real" message.
	  PostAsyncSendMessage(sendMsgStruct);
	  // FIXME: make sure the temporary tag is released when the ack message arrives (not here, somewhere else in the code).
	}
	else {
	  // Small enough: use a single message.
	  sendMsgStruct->tag = tagManager.getPermanentTag("TAG_Create_small_chunk");
	  PostAsyncSendMessage(sendMsgStruct);
	  // The tmp tag was not needed in this case, so we release it.
	  tagManager.releaseTemporaryTag(TMP_TAG_Chunk_data);
	}
      }
    }

    bool Worker::isRegisterChunkStillInProgress(ChunkID cid) {
      LockMutex();
      bool found = false;
      if(incompleteRegisterChunkOperations.find(cid) != incompleteRegisterChunkOperations.end())
	found = true;
      UnlockMutex();
      return found;
    }

    // Only used from parent; on worker, the CreateChunk operation is done in several steps instead, so that the first call can return directly.
    ChunkID Worker::registerChunkDirectly(Chunk const * objPtr,
					std::string class_id_str) {
      ChunkID cid = getIdForRegisterChunk(objPtr, class_id_str);
      startRegisterChunk(objPtr, cid);
      if(isRegisterChunkStillInProgress(cid))
	throw std::runtime_error("Error in registerChunkDirectly; operation still in progress.");
      return cid;
    }

    ChunkID Worker::registerAndGetChunk(Chunk const * objPtrInput,
				      cht::shared_ptr<Chunk const> & objPtrOutput,
				      std::string class_id_str) {
      AccessKey key(this);
      if(class_id_str.length() <= 0)
	throw std::runtime_error("cht::ChunkObjService::Worker::registerAndGetChunk error: (class_id_str.length() <= 0)");
      ChunkID cid = createChunkStruct(class_id_str, objPtrInput->getSize());
      ChunkStruct* theChunk = FindChunkInList(cid, false);
      LockMutex();
      if(!(theChunk->objPtr == 0))
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::registerAndGetChunk(...): objPtr already set.");
      objPtrOutput = objPtrInput;	/* Note that this assignment means transferring 
					   the responsibility for the object! */
      theChunk->objPtr = objPtrOutput;
      chunkCounters[&(*theChunk->objPtr)]++;
      totDataSizeBytes += theChunk->objPtr->memoryUsage();
      totDataSizeBytesMax = (totDataSizeBytes > totDataSizeBytesMax) ?  totDataSizeBytes : totDataSizeBytesMax;
      UnlockMutex();
      return cid;
    }

    bool Worker::chunkResidesLocally(ChunkID cid) {
      AccessKey key(this);
      if(cid.ownerRank < 0 || cid.ownerRank >= key.n_workers())
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::chunkResidesLocally: "
				 "(cid.ownerRank < 0 || cid.ownerRank >= n_workers)");
      return (cid.ownerRank == key.my_rank());
    }

    void Worker::getChunk(ChunkID cid, 
			     cht::shared_ptr<Chunk const> & objPtr) {
      bool success = getChunkIfExists(cid, objPtr);
      if(success == false)
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::getChunk: chunk does not seem to exist.");
    }

    bool Worker::getChunkIfExists(ChunkID cid, 
				     cht::shared_ptr<Chunk const> & objPtr) {
      assert(isRegisterChunkStillInProgress(cid) == false);
      AccessKey key(this);
      if(cid.ownerRank < 0 || cid.ownerRank >= key.n_workers())
	throw std::runtime_error("Error in cht::ChunkObjService::Worker::getChunkIfExists: "
				 "(cid.ownerRank < 0 || cid.ownerRank >= n_workers)");
      LockMutex();
      ++getChunk_total_count;
      UnlockMutex();
      std::string cid_class_id_str = getChunkTypeIDStr(cid.chunkTypeID);
      if(cid.ownerRank == key.my_rank()) {
	// This chunk is stored locally by this process.
	LockMutex();
	++getChunk_found_locally_count;
	UnlockMutex();
	bool chunkExists = TryToGetChunkFromList(cid, objPtr);
	if(!chunkExists) {
	  // This means the chunk does not exist; probably it has been deleted earlier.
	  return false;
	}
	// FIXME: check that cid_class_id_str is correct here? (old check was removed)
	return true;
      }
      // Need to communicate!
      size_t sendBufSz = sizeof(ChunkID) + 2*sizeof(int);
      cht::vector<char> sendBuf(sendBufSz);
      memcpy(&sendBuf[0], &cid, sizeof(ChunkID));
      memcpy(&sendBuf[sizeof(ChunkID)], &cid.chunkTypeID, sizeof(int));
      int TMP_TAG_Chunk_data = tagManager.getTemporaryTag_incoming();
      memcpy(&sendBuf[sizeof(ChunkID) + sizeof(int)], &TMP_TAG_Chunk_data, sizeof(int));
      cht::timer timer_MPI;
      MW._MPI_Send(&sendBuf[0], sendBufSz, MPI_UNSIGNED_CHAR, cid.ownerRank, 
		  tagManager.getPermanentTag("TAG_Get_chunk_data"), 
		  *key.comm_worker_world());  
      // Something should come back from owner worker, but before
      // receiving it we need to find out the size. 
      // FIXME: we actually already know the size since that info is included in the ChunkID. Could that info be used to make this more efficient?
      MPI_Status status;
      MW._MPI_Probe(cid.ownerRank, TMP_TAG_Chunk_data, *key.comm_worker_world(), &status);
      int msgSizeTot;
      MW._MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &msgSizeTot);
      char* recvBuffer = memoryBufferPool.getBufPtr(msgSizeTot);
      MW._MPI_Recv(&recvBuffer[0], msgSizeTot, MPI_UNSIGNED_CHAR, cid.ownerRank, 
		  TMP_TAG_Chunk_data, *key.comm_worker_world(), MPI_STATUS_IGNORE);
      tagManager.releaseTemporaryTag(TMP_TAG_Chunk_data);
      LockMutex();
      // Now the message is supposed to be first one bool saying if the chunk exists or not, and then the chunk data (if the chunk exists).
      bool chunkExists = false;
      memcpy(&chunkExists, &recvBuffer[0], sizeof(bool));
      if(chunkExists == false) {
	// This means the chunk does not exist; probably it has been deleted earlier.
	UnlockMutex();
	memoryBufferPool.releaseBufPtr(recvBuffer);
	return false;
      }
      if(msgSizeTot < smallMsgSizeLimitBytes) {
	getChunk_small_communication_stats.add_timings(timer_MPI);
	getChunk_small_totNoOfBytesReceived += msgSizeTot;
      }
      else {
	getChunk_large_communication_stats.add_timings(timer_MPI);
	getChunk_large_totNoOfBytesReceived += msgSizeTot;
      }
      UnlockMutex();
      // Now we have received the object data.
      int chunkDataSize = msgSizeTot - sizeof(bool);
      const char* chunkDataPtr = &recvBuffer[sizeof(bool)];
      Chunk* newObj = cht::obj_factory<Chunk>::instance().createObject(cid_class_id_str);
      cht::timer timer_assignFromBuffer;
      newObj->assignFromBuffer ( chunkDataPtr, chunkDataSize);
      LockMutex();
      if(msgSizeTot < smallMsgSizeLimitBytes)
	getChunk_small_assignFromBuffer_stats.add_timings(timer_assignFromBuffer);
      else 
	getChunk_large_assignFromBuffer_stats.add_timings(timer_assignFromBuffer);
      UnlockMutex();
      objPtr = newObj; 		/* Transfer responsibility of new object to shared_ptr */
      memoryBufferPool.releaseBufPtr(recvBuffer);
      return true;
    }

    bool Worker::getChunkIfLocal(ChunkID cid, 
				    cht::shared_ptr<Chunk const> & objPtr) {
      AccessKey key(this);
      LockMutex();
      ++getChunkIfLocal_total_count;
      if(cid.ownerRank == key.my_rank()) {
	// This chunk is stored locally by this process.
	++getChunkIfLocal_found_locally_count;
	UnlockMutex();
	ChunkStruct* theChunk = FindChunkInList(cid, true);
	if(theChunk->objPtr == 0)
	  throw std::runtime_error("Error in cht::ChunkObjService::Worker::getChunk(): objPtr is null.");
	objPtr = theChunk->objPtr;
	return true;
      }
      UnlockMutex();
      return false;
    }


  }; // end namespace 
}; // end namespace 

