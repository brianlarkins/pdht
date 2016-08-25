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
#include "MPIWrapperSerialized.h"
#include <stdexcept>
#include <pthread.h>
#include <unistd.h>



MPIWrapperSerialized::MPIWrapperSerialized()
{
  pthread_mutex_init(&mutex, NULL);
}

MPIWrapperSerialized& MPIWrapperSerialized::instance()
{
  static MPIWrapperSerialized theInstance; // local static variable
  return theInstance;
}  


void MPIWrapperSerialized::SleepShortTime() 
{
  // return; // return here to skip all usleep() calls.
  const int sleep_time_useconds = 1000; // set sleep time here
  usleep(sleep_time_useconds);  
}

void MPIWrapperSerialized::LockMutex()
{
  pthread_mutex_lock(&mutex);
}

void MPIWrapperSerialized::UnlockMutex()
{
  pthread_mutex_unlock(&mutex);
}


int MPIWrapperSerialized::MPI_Comm_rank ( MPI_Comm comm, int *rank )
{
  LockMutex();
  int rc = ::MPI_Comm_rank(comm, rank);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Recv(void *buf, int count, MPI_Datatype dtype, 
				   int src, int tag, MPI_Comm comm, MPI_Status *stat)
{
  int rc;
  while(1) {
    MPI_Status status; int flag;
    LockMutex();
    ::MPI_Iprobe(src, tag, comm, &flag, &status);
    if(flag == 1) {
      rc = ::MPI_Recv(buf, count, dtype, src, tag, comm, stat);
      UnlockMutex();
      break;
    }  
    UnlockMutex();
    SleepShortTime();
  }
  return rc;
}


int MPIWrapperSerialized::MPI_Send(void* buf, int count, MPI_Datatype datatype, 
				   int dest, int tag, MPI_Comm comm) 
{
  int rc;
  MPI_Request request;
  LockMutex();
  ::MPI_Isend(buf, count, datatype, dest, tag, comm, &request);
  UnlockMutex();
  // Wait for send to finish.
  while(1) {
    MPI_Status status; int flag;
    LockMutex();
    rc = ::MPI_Test(&request, &flag, &status);
    UnlockMutex();
    if(flag == 1) break;
    SleepShortTime();
  }
  return rc;
}

int MPIWrapperSerialized::MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status) {
  int rc;
  // Wait for message to arrive.
  while(1) {
    int flag;
    LockMutex();
    rc = ::MPI_Iprobe(source, tag, comm, &flag, status);
    UnlockMutex();
    if(flag == 1) break;
    SleepShortTime();
  }
  return rc; 
}


int MPIWrapperSerialized::MPI_Iprobe(int src, int tag, MPI_Comm comm, 
				     int *flag, MPI_Status *stat)
{
  LockMutex();
  int rc = ::MPI_Iprobe(src, tag, comm, flag, stat);
  UnlockMutex();
  return rc;
}


int MPIWrapperSerialized::MPI_Comm_size(MPI_Comm comm, int *psize)
{
  LockMutex();
  int rc = ::MPI_Comm_size(comm, psize);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Comm_remote_size(MPI_Comm comm, int *psize)
{
  LockMutex();
  int rc = ::MPI_Comm_remote_size(comm, psize);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm) 
{
  LockMutex();
  int rc = ::MPI_Comm_dup(comm, newcomm);
  UnlockMutex();
  return rc;
}


int MPIWrapperSerialized::MPI_Comm_free(MPI_Comm *comm)
{
  LockMutex();
  int rc = ::MPI_Comm_free(comm);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Isend(void *buf, int count, MPI_Datatype dtype, 
				    int dest, int tag, MPI_Comm comm, 
				    MPI_Request *req)
{
  LockMutex();
  int rc = ::MPI_Isend(buf, count, dtype, dest, tag, comm, req);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Irecv(void *buf, int count, MPI_Datatype dtype, 
				    int src, int tag, MPI_Comm comm, 
				    MPI_Request *req)
{
  LockMutex();
  int rc = ::MPI_Irecv(buf, count, dtype, src, tag, comm, req);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Waitall(int count, MPI_Request *array_of_requests, 
				      MPI_Status *array_of_statuses) {
  int rc;
  while(1) {
    LockMutex();
    int flag;
    rc = ::MPI_Testall(count, array_of_requests, &flag, 
		       array_of_statuses);
    UnlockMutex();
    if(flag) break;
    SleepShortTime();
  }
  return rc; 
}
int MPIWrapperSerialized::MPI_Waitany(int count, MPI_Request *array_of_requests, 
				      int *index, MPI_Status *status) {
  int rc;
  while(1) {
    LockMutex();
    int flag;
    rc = ::MPI_Testany(count, array_of_requests, index, 
		       &flag, status);
    UnlockMutex();
    if(flag) break;
    SleepShortTime();
  }
  return rc; 
}


int MPIWrapperSerialized::MPI_Waitsome(int count, MPI_Request* reqs, int *ndone, 
				       int *indices, MPI_Status *stats)
{
  int rc;
  while(1) {
    LockMutex();
    rc = ::MPI_Testsome(count, reqs, ndone, indices, stats);
    UnlockMutex();
    if(*ndone) break;
    SleepShortTime();
  }
  return rc; 
}

int MPIWrapperSerialized::MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count)
{
  LockMutex();
  int rc = ::MPI_Get_count(status, datatype, count);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Testsome(int incount, MPI_Request *array_of_requests, int *outcount, 
				       int *array_of_indices, MPI_Status *array_of_statuses)
{
  LockMutex();
  int rc = ::MPI_Testsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Get_processor_name(char *name, int *resultlen) {
  LockMutex();
  int rc = ::MPI_Get_processor_name(name, resultlen);
  UnlockMutex();
  return rc;  
}

int MPIWrapperSerialized::MPI_Cancel(MPI_Request *request) {
  LockMutex();
  int rc = ::MPI_Cancel(request);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Request_free(MPI_Request *request) {
  LockMutex();
  int rc = ::MPI_Request_free(request);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Test(MPI_Request *request, int *flag, MPI_Status *status) {
  LockMutex();
  int rc = ::MPI_Test(request, flag, status);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) {
  LockMutex();
  int rc = ::MPI_Comm_split(comm, color, key, newcomm);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, 
					       MPI_Comm bridge_comm, int remote_leader, 
					       int tag, MPI_Comm *newintercomm) {
  LockMutex();
  int rc = ::MPI_Intercomm_create(local_comm, local_leader, 
				  bridge_comm, remote_leader, 
				  tag, newintercomm);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Init_thread(int *argc, char ***argv, int required, int *provided) {
  LockMutex();
  int rc = ::MPI_Init_thread(argc, argv, required, provided);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Init(int *argc, char ***argv) {
  LockMutex();
  int rc = ::MPI_Init(argc, argv);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Finalize() {
  LockMutex();
  int rc = ::MPI_Finalize();
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Abort(MPI_Comm comm, int errorcode) {
  LockMutex();
  int rc = ::MPI_Abort(comm, errorcode);
  UnlockMutex();
  return rc;
}


#if MPI_VERSION_MPI_WRAPPER_SERIALIZED>1

int MPIWrapperSerialized::MPI_Comm_get_parent(MPI_Comm *parent)
{
  LockMutex();
  int rc = ::MPI_Comm_get_parent(parent);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Comm_spawn(char* command, char** argv, int maxprocs, MPI_Info info, 
					 int root, MPI_Comm comm, MPI_Comm *intercomm, 
					 int *errcodes)
{
  LockMutex();
  int rc = ::MPI_Comm_spawn(command, argv, maxprocs, info, root, comm, intercomm, errcodes);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status)
{
  LockMutex();
  int rc = ::MPI_Request_get_status(request, flag, status);
  UnlockMutex();
  return rc;
}

int MPIWrapperSerialized::MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag)
{
  LockMutex();
  int rc = ::MPI_Comm_get_attr(comm, comm_keyval, attribute_val, flag);
  UnlockMutex();
  return rc;
}

#endif



// Extra non-MPI routine: probe for any of the tags in given list. 
int MPIWrapperSerialized::NonMPI_Iprobe_multiple(int src, const int tagList[], int nTags, MPI_Comm comm, 
						 int *flag, MPI_Status *stat)
{
  int rc;
  for(int i = 0; i < nTags; i++)
    {
      LockMutex();
      rc = ::MPI_Iprobe(src, tagList[i], comm, flag, stat);
      UnlockMutex();
      if(*flag) return rc;
    }
  return rc;
}

