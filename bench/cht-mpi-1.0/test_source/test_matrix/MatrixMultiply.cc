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
#include "MatrixMultiply.h"
#include "MatrixAdd.h"
#include "CreateMatrixFromIds.h"
#include <cstring>

CHT_TASK_TYPE_IMPLEMENTATION((MatrixMultiply));
cht::ID MatrixMultiply::execute(CMatrix const & A, CMatrix const & B) {
  int nA = A.n;
  int nB = B.n;
  if(nA != nB)
    throw std::runtime_error("Error in MatrixMultiply::execute: (nA != nB).");
  int n = nA;
  if(n <= CMatrix::BLOCK_SIZE) {
    // Lowest level
    CMatrix* C = new CMatrix();
    C->n = n;
    C->elements.resize(n*n);
    memset(&C->elements[0], 0, n*n*sizeof(double));
    for(int i = 0; i < n; i++)
      for(int k = 0; k < n; k++)
	for(int j = 0; j < n; j++) {
	  double Aik = A.elements[i*n+k];
	  double Bkj = B.elements[k*n+j];
	  C->elements[i*n+j] += Aik * Bkj;
	}
    return registerChunk(C, cht::persistent);
  }
  else {
    for(int i = 0; i < 4; i++) {
      if(A.children[i] == cht::CHUNK_ID_NULL)
	throw std::runtime_error("Error in MatrixMultiply::execute: CHUNK_ID_NULL found for A.");
      if(B.children[i] == cht::CHUNK_ID_NULL)
	throw std::runtime_error("Error in MatrixMultiply::execute: CHUNK_ID_NULL found for B.");
    }
    // Not lowest level
    cht::ID childTaskIDs[4];
    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 2; j++) {
	cht::ID childTaskIDsForSum[2];
	for(int k = 0; k < 2; k++)
	  childTaskIDsForSum[k] = registerTask<MatrixMultiply>(A.children[i*2+k], B.children[k*2+j]);
	childTaskIDs[i*2+j] = registerTask<MatrixAdd>(childTaskIDsForSum[0], childTaskIDsForSum[1]);
      }
    cht::ChunkID cid_n = registerChunk( new CInt(n) );
    return registerTask<CreateMatrixFromIds>(cid_n, childTaskIDs[0], childTaskIDs[1], childTaskIDs[2], childTaskIDs[3], cht::persistent);
  }
} // end execute
