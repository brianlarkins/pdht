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
#include <cstdio>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <unistd.h>
#include "chunks_and_tasks.h"
#include "ChunkObjDummyData.h"
#include "ChunkObjInputParams.h"
#include "test_overlap_task_implementations.h"

int main(int argc, char* const  argv[])
{
  try {
    int nWorkerProcs   = 2;
    int nThreads       = 2;
    int n_child_tasks = 4;
    size_t data_size = 1000000;
    double work_seconds = 1.0;
    cht::extras::Cache::Mode cache_mode = cht::extras::Cache::Disabled;
    bool verbose = false;
    // INPUT PARSING
    int c;
    opterr = 0;
    while ((c = getopt (argc, argv, "w:t:c:s:x:ev")) != -1)
      switch (c)
	{
	case 'w': // 
	  nWorkerProcs = atoi(optarg);
	  break;
	case 't':
	  nThreads = atoi(optarg);
	  break;
	case 'c':
	  n_child_tasks = atoi(optarg);
	  break;
	case 's':
	  data_size = (size_t)1000000*(double)atof(optarg);
	  break;
	case 'x':
	  work_seconds = atof(optarg);
	  break;
	case 'e':
	  cache_mode = cht::extras::Cache::Enabled;
	  break;
	case 'v':
	  verbose = true;
	  break;
	case '?':
	  if (optopt == 'w' || optopt == 't' || optopt == 'c' || optopt == 's')
	    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	  else if (isprint (optopt))
	    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	  else
	    fprintf (stderr,
		     "Unknown option character `\\x%x'.\n",
		     optopt);
	  return 1;
	default:
	  abort ();
	}
    if (argc - optind != 0) {
      std::cout<<"Usage: "<<argv[0]<<" [options] \n";
      std::cout<<"Available options:\n";
      std::cout<<"  -w workers          : Number of MPI worker processes\n";
      std::cout<<"  -t threads          : Number of threads for each worker process\n";
      std::cout<<"  -c child_tasks      : Total number of child tasks\n";
      std::cout<<"  -s data_size        : Size of data to be communicated for each task (MB)\n";
      std::cout<<"  -x work_seconds     : Time that each task simulates 'work' (wall seconds)\n";
      std::cout<<"  -e                  : Enable cache for chunk objects. \n";
      std::cout<<"  -v                  : Run program with verbose output. \n";
      std::cout<<"Examples:\n  "
	       << argv[0]<<" -w 2 -t 2 -a 1"<<std::endl;
      std::exit(1);
    }
    
    printf("nWorkerProcs   = %5d\n", nWorkerProcs);
    printf("nThreads       = %5d\n", nThreads);
    printf("child tasks    = %5d\n", n_child_tasks);
    printf("data size      = %5d MB\n", (int)(data_size/1000000));
    printf("work sim time  = %5.1f seconds\n", work_seconds);
    if ( cache_mode == cht::extras::Cache::Enabled )
      printf("Cache enabled\n");
    else
      printf("Cache disabled\n");
    
    size_t cacheMemoryUsageLimit = 1000000000;
    cht::extras::setNWorkers(nWorkerProcs);
    cht::setOutputMode(cht::Output::AllInTheEnd);
    cht::setOutputLevel(cht::Output::Debug);
    cht::extras::setCacheMode(cache_mode);
    cht::extras::setCacheSize(cacheMemoryUsageLimit);
    cht::extras::setNoOfWorkerThreads(nThreads);
    cht::start();

    cht::output("Testing output message from parent!");

    std::cout << "Calling executeMotherTask() for TaskTypeOverlapTestMother..." << std::endl;
    cht::resetStatistics();
    // Create chunk with input data to TaskTypeOverlapTestMother.    
    ChunkObjInputParams inputParams;
    inputParams.x.n_child_tasks = n_child_tasks;
    inputParams.x.data_size = data_size;
    inputParams.x.work_seconds = work_seconds;
    cht::ChunkID cid_input_params = cht::registerChunk<ChunkObjInputParams>(new ChunkObjInputParams(inputParams));
    std::vector<cht::ChunkID> inputChunks(1);
    inputChunks[0] = cid_input_params;
    cht::ChunkID cid_result = cht::executeMotherTask<TaskTypeOverlapTestMother>(inputChunks);
    cht::reportStatistics();

    cht::deleteChunk(cid_input_params);
    cht::deleteChunk(cid_result);

    // Stop cht services
    cht::stop();

    std::cout << "test_comp_comm_overlap OK, after cht::stop()." << std::endl;

#if 0
    // Check result
    std::cout << "Checking result..." << std::endl;    
    size_t serialResult = Fib(n);
    std::cout << "Chunks&Tasks computed result :  Fib(" << n << ") =  " << result->x << std::endl;
    std::cout << "Serially computed result     :  Fib(" << n << ") =  " << serialResult << std::endl;
    if (result->x == serialResult)
      std::cout << "Result OK." << std::endl;
    else      
      std::cout << "Test failed: wrong result" << std::endl;
#endif
  } catch ( std::exception & e) {
    std::cerr << "Exception caught in test main! What: "<< e.what() << std::endl;
    return 1;
  }
  
  return 0;
}

