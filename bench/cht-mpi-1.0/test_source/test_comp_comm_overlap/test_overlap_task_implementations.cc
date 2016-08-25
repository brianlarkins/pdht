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
#include <sstream>
#include "ChunkObjDummyData.h"
#include "ChunkObjInputParams.h"
#include "test_overlap_task_implementations.h"
// Normally not included:
#include "utilities/cht_time.h"


class TaskTypeOverlapTestChild;

CHT_TASK_TYPE_IMPLEMENTATION((TaskTypeOverlapTestMother));


cht::ID TaskTypeOverlapTestMother::execute(ChunkObjInputParams const & inp) {
  // Get input params (n_child_tasks and data_size) 
  n_child_tasks = inp.x.n_child_tasks;
  data_size = inp.x.data_size;

  // CREATE CHUNK
  cid = registerChunk<ChunkObjDummyData>(new ChunkObjDummyData(data_size), cht::persistent);

  // REGISTER TASKS
  childTaskInputParamChunks.resize(n_child_tasks);
  for ( int ind = 0; ind < n_child_tasks; ++ind ) {
    ChunkObjInputParams child_params = inp;
    cht::ChunkID cid_param_child = 
      registerChunk<ChunkObjInputParams>(new ChunkObjInputParams(child_params));
    std::vector<cht::ID> childTaskInputAndDepend(2);
    childTaskInputAndDepend[0] = cid_param_child;
    childTaskInputAndDepend[1] = cid;
    cht::TaskID tid = registerTask<TaskTypeOverlapTestChild>(childTaskInputAndDepend);
    childTaskInputParamChunks[ind] = cid_param_child;
  }
  return cid;
}

/////// Child task type

class TaskTypeOverlapTestChild : public cht::Task {
  size_t data_size;
public:
  cht::ID execute(ChunkObjInputParams const &,
		  ChunkObjDummyData const &);
  CHT_TASK_INPUT((ChunkObjInputParams,ChunkObjDummyData));
  CHT_TASK_OUTPUT((ChunkObjDummyData));
  CHT_TASK_TYPE_DECLARATION;
};

CHT_TASK_TYPE_IMPLEMENTATION((TaskTypeOverlapTestChild));

cht::ID TaskTypeOverlapTestChild::execute(ChunkObjInputParams const & inp,
					  ChunkObjDummyData const & dummy) {
  double work_seconds  = inp.x.work_seconds;
  data_size = inp.x.data_size;

  // Get input params 
  std::stringstream s;
  s << " TID: "<< pthread_self() <<" Child task before getElement<1>";
  output(s.str());
  cht::timer computePretendTimer;
  std::stringstream s2;
  s2 << " TID: "<< pthread_self() <<" Child task before compute";
  output(s2.str());
  while ( computePretendTimer.get_elapsed_wall_seconds() < work_seconds ) {}
  std::stringstream s3;
  s3 << " TID: "<< pthread_self() <<" Child task after compute";
  output(s3.str());
  
  // CREATE CHUNK
  cht::ChunkID cid = registerChunk<ChunkObjDummyData>(new ChunkObjDummyData(data_size), cht::persistent);
  return cid;
}


