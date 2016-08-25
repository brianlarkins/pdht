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
#ifndef CHUNK_OBJ_SERVICE_WORKER_HEADER
#define CHUNK_OBJ_SERVICE_WORKER_HEADER

#include "services/Service_worker.h"
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <cstdlib>
#include "utilities/cht_utils.h"
#include "utilities/cht_threads.h"
#include "services/output_service/OutputService_worker.h"
#include "services/chunk_obj_service/ChunkObjService.h"
#include "services/chunk_obj_service/MemoryBufferPool.h"

namespace cht {
  namespace ChunkObjService {

    static void* global_thread_func(void*);

    class Worker : 
    public Base<Service::Worker>,
      public Singleton<Worker> {
	friend class Singleton<Worker>;
    private:  
	struct AsyncMsgStruct {
	  char* bufPtr;
	  size_t bufSz;
	  int rank;
	  int tag;
	  MPI_Comm* comm;
	  MPI_Request request;
	  MemoryBufferPool & memBufPool;
	AsyncMsgStruct(size_t size, int rank_, int tag_, MPI_Comm* comm_, MemoryBufferPool & memBufPool_) : bufSz(size), rank(rank_), tag(tag_), comm(comm_), memBufPool(memBufPool_) { bufPtr = memBufPool.getBufPtr(size); }
	  ~AsyncMsgStruct() { memBufPool.releaseBufPtr(bufPtr); }
	};
	std::list<AsyncMsgStruct*> sendMsgStructList;
	std::list<AsyncMsgStruct*> recvMsgStructList;
	bool sendMsgStructListIsEmpty();
	void PostAsyncSendMessage(AsyncMsgStruct* msgStruct);
	void PostAsyncRecvMessage(AsyncMsgStruct* msgStruct);
	void FreeMemoryForCompletedSendOperations();
	void HandleAnyCompletedLargeReceiveOperations();
	void HandleCompletedLargeReceiveOperation(AsyncMsgStruct* recvMsgStructPtr);

	MemoryBufferPool memoryBufferPool;

	void writeMessageToFile(int my_rank, const char* message);
  
	struct RecvMsgBufStruct {
	  cht::vector<char> buf_Delete_chunk;
	  cht::vector<char> buf_Create_small_chunk;
	  cht::vector<char> buf_Create_large_chunk_prep;
	  cht::vector<char> buf_Get_chunk_data;
	  cht::vector<char> buf_Copy_chunk;
	  cht::vector<char> buf_Copy_chunk_ack;
	  cht::vector<char> buf_Create_chunk_ack;
	RecvMsgBufStruct() 
	: buf_Delete_chunk(sizeof(ChunkID) + sizeof(int)),
	    buf_Create_small_chunk(MAX_WRITECHUNKOBJ_SMALL_MESSAGE_SIZE),
	    buf_Create_large_chunk_prep(2*sizeof(int)),
	    buf_Get_chunk_data(max_class_id_str_sz+sizeof(ChunkID)+sizeof(int)),
            buf_Copy_chunk(2*sizeof(ChunkID)), 
	    buf_Copy_chunk_ack(sizeof(ChunkID)),
	    buf_Create_chunk_ack(sizeof(ChunkID)+sizeof(int))
	  { }
	};
	char* GetMessageBuf(RecvMsgBufStruct & bufs, int i, int tag_in, int & bufSz);

    private:
	static int const nReceiveTags = 8;
	void PostReceiveOperations(RecvMsgBufStruct & bufs,
				   MPI_Comm* comm,
				   MPI_Request* requests,
				   int rank);
	int getChunk_total_count;
	int getChunk_found_locally_count;
	int getChunkIfLocal_total_count;
	int getChunkIfLocal_found_locally_count;
	cht::work_statistics getChunk_small_communication_stats;
	cht::work_statistics getChunk_small_assignFromBuffer_stats;
	cht::work_statistics getChunk_large_communication_stats;
	cht::work_statistics getChunk_large_assignFromBuffer_stats;
	cht::work_statistics getChunkReq_small_writeToBuffer_stats;
	cht::work_statistics getChunkReq_large_writeToBuffer_stats;
	cht::work_statistics writeChunkObj_communication_stats;
	cht::work_statistics writeChunkObj_writeToBuffer_stats;
	cht::work_statistics deleteChunk_communication_stats;
	cht::work_statistics deleteChunk_localwork_stats;
	size_t getChunk_small_totNoOfBytesReceived;
	size_t getChunk_large_totNoOfBytesReceived;
	cht::work_statistics worker_thread_func_loop_stats;
	cht::work_statistics worker_thread_func_handle_request_parent_stats;
	cht::work_statistics worker_thread_func_handle_request_worker_stats;
	cht::work_statistics FreeMemoryForCompletedSendOperations_stats;
	cht::work_statistics worker_thread_func_irecv_stats;

	double probabilityToRegisterChunksLocally;
  
	struct ChunkStruct {
	  ChunkID id;
	  cht::shared_ptr<Chunk const> objPtr;
	  std::string class_id_str;
          ChunkStruct() : objPtr(0) { }
	};
	friend void* global_thread_func(void*);
	virtual void start_derived();
  	virtual void stop_derived();

    public:
	void deleteChunk(ChunkID cid);
	void getChunk(ChunkID cid, cht::shared_ptr<Chunk const> & objPtr);
	bool getChunkIfLocal(ChunkID cid, 
			     cht::shared_ptr<Chunk const> & objPtr);
	bool getChunkIfExists(ChunkID cid, 
			      cht::shared_ptr<Chunk const> & objPtr);
	template<typename Tobj>
	  void getChunks(ChunkID cid1, cht::shared_ptr<Tobj const> & objPtr1,
			 ChunkID cid2, cht::shared_ptr<Tobj const> & objPtr2);
	void resetStatistics();
	void reportStatistics(std::string messageHeader);

	ChunkID getIdForRegisterChunk(Chunk const * objPtr,
				      std::string class_id_str);
	void startRegisterChunk(Chunk const * objPtr, ChunkID cid);
	bool isRegisterChunkStillInProgress(ChunkID cid);
	ChunkID registerChunkDirectly(Chunk const * objPtr,
				      std::string class_id_str);

	ChunkID getIdForCopyChunk(ChunkID cid_old);
	void startCopyChunk(ChunkID cid_old, ChunkID cid_new);
	bool isCopyChunkStillInProgress(ChunkID cid);

	ChunkID registerAndGetChunk(Chunk const * objPtrInput,
				    cht::shared_ptr<Chunk const> & objPtrOutput,
				    std::string class_id);
	bool chunkResidesLocally(ChunkID cid);
	void call_mpi_abort();
    protected:
	ChunkID createChunkStruct(std::string class_id_str, size_t size, ChunkID known_cid = CHUNK_ID_NULL);
	ChunkID createChunkObjectFromBuffer(char const * dataBuffer, size_t const bufferSize, 
					    std::string const class_id_str, ChunkID known_cid = CHUNK_ID_NULL);
    private:
	
	int getIntFromEndOfBuffer(const char* buf, int bufSz);

	int idCounter;
	// A more optimized map could possibly speed up lookups if needed
	typedef std::map<ChunkID, ChunkStruct*> ChunkListMap;
	ChunkListMap chunkList; // Must be protected by mutex!
	// We use a separate map to keep track of how many ChunkIDs are referring to the same object.
	typedef std::map<Chunk const *, int> ChunkCounterMap;
	ChunkCounterMap chunkCounters; // Must be protected by mutex!
	// Keep track of ongoing registerChunk and copyChunk and operations using sets of ChunkIDs
	std::set<ChunkID> incompleteRegisterChunkOperations;
	std::set<ChunkID> incompleteCopyChunkOperations;
	size_t totDataSizeBytes;
	size_t totDataSizeBytesMax;
	Threads::Thread* threadRequestHandler; /**< Request handler thread pointer. @warning Not protected by mutex. */
	Threads::Mutex mutex;
	void LockMutex();
	void UnlockMutex();
	bool TryToGetChunkFromList(ChunkID cid, cht::shared_ptr<Chunk const> & objPtr); // Returns true if chunk existed, false otherwise.
	ChunkStruct* TryToFindChunkInList(ChunkID cid, bool objShouldExist); // Returns NULL if chunk does not exist in list.
	ChunkStruct* FindChunkInList(ChunkID cid, bool objShouldExist); // Should only be used if it is known that the chunk exists in list.
	ChunkID GetChunkIdFromBuffer(const char* buf);
	void CopyChunkLocally(ChunkID cid_org, ChunkID cid_new);

	void handle_request(const char* buf,
			    int msgSize,
			    int tag,
			    int rank,
			    MPI_Comm* comm,
			    bool & timeToStopFlag);

	void worker_thread_func();

    private:
	OutputService::Worker& OS;
	Worker();
      };


    // FIXME: remove this if not used/needed.
    template<typename Tobj>
      void Worker::getChunks(ChunkID cid1, cht::shared_ptr<Tobj const> & objPtr1,
			     ChunkID cid2, cht::shared_ptr<Tobj const> & objPtr2) {
      getChunk(cid1, objPtr1);
      getChunk(cid2, objPtr2);
    }
    
  }; // end namespace 
}; // end namespace 

#endif
