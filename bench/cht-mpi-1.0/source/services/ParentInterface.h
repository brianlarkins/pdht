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
#ifndef PARENTINTERFACE_HEADER
#define PARENTINTERFACE_HEADER
#include "services/MPIWrapperInclude.h"
#include <string>
#include "utilities/Singleton.h"
#include "utilities/cht_threads.h"
#include "services/chunks_and_tasks_params.h"

namespace cht {
  class ParentInterface : public Singleton<ParentInterface> {
    friend class Singleton<ParentInterface>;
  protected:
    ParentInterface();
    ~ParentInterface();
    std::string service_worker_program_name;
    unsigned int n_workers;
    MPI_Comm intraComm; // Only needed for single program case
    MPI_Comm comm_to_workers;
    Threads::Thread* mainThread;
    bool access_granted;
    void verifyInstanceAccess();    
  public:
    void i_am_worker();
    void start();
    void stop();
    void setNWorkers(unsigned int n_workers);    
    void setOutputMode(Output::Mode mode);
    void setOutputLevel(Output::Priority prio);
    void setCacheMode(extras::Cache::Mode mode);
    void setCacheSize(size_t cache_memory_usage_limit);
    void setNoOfWorkerThreads(unsigned int n_worker_threads);
    void setNoOfPrefetchTasks(unsigned int n_prefetch_tasks);
    void setChildChunkFetchLevel(int child_chunk_fetch_level);
    void setDebugParams(bool, bool, bool, int, int, bool, double);
  };
} // end namespace cht
#endif