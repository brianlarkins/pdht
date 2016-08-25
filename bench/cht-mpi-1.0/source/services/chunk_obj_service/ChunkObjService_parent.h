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
#ifndef CHUNK_OBJ_SERVICE_PARENT_HEADER
#define CHUNK_OBJ_SERVICE_PARENT_HEADER

#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <cstdlib>
#include "services/Service_parent.h"
#include "services/chunk_obj_service/ChunkObjService.h"
#include "services/output_service/OutputService_parent.h"

namespace cht {
  namespace ChunkObjService {

    static void* global_thread_func(void*);

    class Parent : 
    public Base<Service::Parent>,
      public Singleton<Parent> {
	friend class Singleton<Parent>;
    private:    
	virtual void start_derived();
	virtual void stop_derived();
    public:
	void deleteChunk(ChunkID cid);

	bool getChunkIfExists(ChunkID cid, cht::shared_ptr<Chunk const> & objPtr);
	void getChunk(ChunkID cid, cht::shared_ptr<Chunk const> & objPtr);
	bool getChunkIfLocal(ChunkID cid, cht::shared_ptr<Chunk const> & objPtr) { return false; }

	ChunkID registerChunkDirectly(Chunk const * objPtr,
				      std::string class_id_str);

	ChunkID registerAndGetChunk(Chunk const * objPtrInput,
				    cht::shared_ptr<Chunk const> & objPtrOutput,
				    std::string class_id);
		
	bool chunkResidesLocally(ChunkID cid);

	void setDebugParams(double probabilityToRegisterChunksLocally_);

    protected:
	ChunkID registerChunk(std::string class_id_str);
	void createChunkObjectFromBuffer(ChunkID cid, 
					 char const * dataBuffer, size_t const bufferSize, 
					 std::string const class_id_str);  
    private:
	OutputService::Parent& OS;

	/* The following member variables exist for both parent and worker,
	 * but for the parent the vaiables are only used in order to send
	 * the values to the workers. */
	double probabilityToRegisterChunksLocally; /* Mainly for debugging. */

    private:
	Parent();
      };

  }; // end namespace 
}; // end namespace 


#endif
