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
#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include "chunks_and_tasks.h"
#include "CInt.h"
#include "CDouble.h"
#include "CMatrix.h"
#include "CreateMatrix.h"
#include "GetMatrixElement.h"
#include "MatrixMultiply.h"
#include "MatrixAdd.h"
#include "MatrixElementValues.h"

static double compute_product_matrix_element(int N, int i, int j) {
  double sum = 0;
  for(int k = 0; k < N; k++) {
    double Aik = matElementFunc(MATRIX_TYPE_A, i, k);
    double Bkj = matElementFunc(MATRIX_TYPE_B, k, j);
    sum += Aik * Bkj;
  }
  return sum;
}

int main(int argc, char* const  argv[])
{
  try {
    if(argc != 3) {
      std::cout << "Please give 2 arguments: N nWorkerProcs" << std::endl;
      return -1;
    }
    long int N = atoi(argv[1]);
    int nWorkerProcs = atoi(argv[2]);
    std::cout << "CMatrix::BLOCK_SIZE = " << CMatrix::BLOCK_SIZE << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "nWorkerProcs = " << nWorkerProcs << std::endl;
    size_t size_of_matrix_in_bytes = N*N*sizeof(double);
    double size_of_matrix_in_GB = (double)size_of_matrix_in_bytes / 1000000000;
    std::cout << "size_of_matrix_in_GB = " << size_of_matrix_in_GB << std::endl;
    int nThreads = 1;
    cht::extras::setNWorkers(nWorkerProcs);
    cht::setOutputMode(cht::Output::AllInTheEnd);
    cht::extras::setNoOfWorkerThreads(nThreads);

    cht::extras::Cache::Mode cache_mode = cht::extras::Cache::Enabled;
    cht::extras::setCacheMode(cache_mode);
    size_t cacheMemoryUsageLimit = 100000000;
    cht::extras::setCacheSize(cacheMemoryUsageLimit);

    cht::start();

    cht::ChunkID cid_baseIdx1 = cht::registerChunk<CInt>(new CInt(0));
    cht::ChunkID cid_baseIdx2 = cht::registerChunk<CInt>(new CInt(0));
    cht::ChunkID cid_n = cht::registerChunk<CInt>(new CInt(N));
    cht::ChunkID cid_matType_A = cht::registerChunk<CInt>(new CInt(1));
    cht::ChunkID cid_matType_B = cht::registerChunk<CInt>(new CInt(2));

    std::cout << "Calling executeMotherTask() for CreateMatrix for A..." << std::endl;
    cht::ChunkID cid_matrix_A = cht::executeMotherTask<CreateMatrix>(cid_n, cid_baseIdx1, cid_baseIdx2, cid_matType_A);

    std::cout << "Calling executeMotherTask() for CreateMatrix for B..." << std::endl;
    cht::ChunkID cid_matrix_B = cht::executeMotherTask<CreateMatrix>(cid_n, cid_baseIdx1, cid_baseIdx2, cid_matType_B);

    std::cout << "Calling executeMotherTask() for MatrixMultiply to compute C = A * B ..." << std::endl;
    cht::ChunkID cid_matrix_C = cht::executeMotherTask<MatrixMultiply>(cid_matrix_A, cid_matrix_B);

    int nElementsToVerify = 20;
    std::cout << "Verifying result by checking " << nElementsToVerify << " C matrix elements..." << std::endl;
    double max_abs_diff = 0;
    for(int i = 0; i < nElementsToVerify; i++) {
      int idx1 = rand() % N;
      int idx2 = rand() % N;
      cht::ChunkID cid_idx1 = cht::registerChunk<CInt>(new CInt(idx1));
      cht::ChunkID cid_idx2 = cht::registerChunk<CInt>(new CInt(idx2));      
      cht::ChunkID cid_value = cht::executeMotherTask<GetMatrixElement>(cid_matrix_C, cid_idx1, cid_idx2);
      cht::shared_ptr<CDouble const> valuePtr;
      cht::getChunk(cid_value, valuePtr);
      double value = *valuePtr;
      // Compute expected value for this C matrix element
      double value_expected = compute_product_matrix_element(N, idx1, idx2);
      double absdiff = std::fabs(value - value_expected);
      //      std::cout << "Checking C matrix element ( " << idx1 << " , " << idx2 << " ) : value = " << value << " , value_expected = " << value_expected << " , absdiff = " << absdiff << std::endl;
      if(absdiff > max_abs_diff)
	max_abs_diff = absdiff;
      cht::deleteChunk(cid_idx1);
      cht::deleteChunk(cid_idx2);
      cht::deleteChunk(cid_value);
    }
    if(max_abs_diff > 1e-8) {
      std::cout << "Error: absdiff too large, result seems wrong, max_abs_diff = " << max_abs_diff << "." << std::endl;
      return -1;
    }
    std::cout << "OK, result seems correct, max_abs_diff = " << max_abs_diff << "." << std::endl;
    
    std::cout << "Cleaning up..." << std::endl;
    cht::deleteChunk(cid_baseIdx1);
    cht::deleteChunk(cid_baseIdx2);
    cht::deleteChunk(cid_n);
    cht::deleteChunk(cid_matrix_A);
    cht::deleteChunk(cid_matrix_B);
    cht::deleteChunk(cid_matrix_C);
    cht::deleteChunk(cid_matType_A);
    cht::deleteChunk(cid_matType_B);

    // Stop cht services
    cht::stop();

    std::cout << "Done, test_matrix finished OK." << std::endl;

  } catch ( std::exception & e) {
    std::cerr << "Exception caught in test main! What: "<< e.what() << std::endl;
    return 1;
  }
  
  return 0;
}

