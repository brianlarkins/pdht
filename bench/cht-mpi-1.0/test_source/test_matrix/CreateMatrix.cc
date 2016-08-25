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
#include "CreateMatrix.h"
#include "CreateMatrixFromIds.h"
#include "MatrixElementValues.h"
#include <cmath>

CHT_TASK_TYPE_IMPLEMENTATION((CreateMatrix));
cht::ID CreateMatrix::execute(CInt const & matSize,
			      CInt const & baseIdx1,
			      CInt const & baseIdx2,
			      CInt const & matType) {
  int n = matSize;
  if(n <= CMatrix::BLOCK_SIZE) {
    // Lowest level
    CMatrix* A = new CMatrix();
    A->n = n;
    A->elements.resize(n*n);
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++) {
	int idx1 = baseIdx1 + i;
	int idx2 = baseIdx2 + j;
	A->elements[i*matSize+j] = matElementFunc(matType, idx1, idx2);
      }
    return registerChunk(A, cht::persistent);
  }
  else {
    // Not lowest level
    if(matSize % 2 != 0)
      throw std::runtime_error("Error in CreateMatrix::execute: matSize not divisible by 2.");
    int nHalf = matSize / 2;
    cht::ChunkID cid_nHalf = registerChunk( new CInt(nHalf) );
    cht::ID childTaskIDs[4];
    for(int i1 = 0; i1 < 2; i1++) {
      cht::ChunkID cid_baseIdx_i1 = registerChunk( new CInt(baseIdx1+i1*nHalf) );
      for(int i2 = 0; i2 < 2; i2++) {
	cht::ChunkID cid_baseIdx_i2 = registerChunk( new CInt(baseIdx2+i2*nHalf) );
	childTaskIDs[i1*2+i2] = registerTask<CreateMatrix>(cid_nHalf, cid_baseIdx_i1, cid_baseIdx_i2, getInputChunkID(matType));
      }
    }
    return registerTask<CreateMatrixFromIds>(getInputChunkID(matSize), childTaskIDs[0], childTaskIDs[1], childTaskIDs[2], childTaskIDs[3], cht::persistent);
  }
} // end execute
