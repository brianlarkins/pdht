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
#include "determine_process_roles.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cstring>

namespace cht {

  void determine_process_roles(int nWorkerProcs, 
			       ProcessRole & resultProcessRole,
			       MPI_Comm & resultIntraComm, 
			       MPI_Comm & resultInterComm) {

    // Find out what role this process should take: parent, worker, or
    // idle.

    MPI_Wrapper & MW = MPI_Wrapper::instance();

    int nProcsTot;
    MW._MPI_Comm_size(MPI_COMM_WORLD, &nProcsTot);
    int myRank;
    MW._MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // Check that there are enough processes.
    // There must be at least nWorkerProcs+1 processes (workers+parent).
    if(nProcsTot < nWorkerProcs+1)
      throw std::runtime_error("Error in determine_process_roles: (nProcsTot < nWorkerProcs+1). "
			       "There must be at least nWorkerProcs+1 processes (workers+parent).");

    char nameList[nProcsTot][MPI_MAX_PROCESSOR_NAME];
    int resultlen;
    MW._MPI_Get_processor_name(nameList[myRank], &resultlen);
    // Send this processor_name to all other processes.
    MPI_Request sendRequests[nProcsTot-1];
    int count = 0;
    for(int i = 0; i < nProcsTot; i++) {
      if(i == myRank)
	continue;
      int tag = 0;
      MW._MPI_Isend(nameList[myRank], MPI_MAX_PROCESSOR_NAME, 
		   MPI_UNSIGNED_CHAR, i, tag, MPI_COMM_WORLD, &sendRequests[count]);
      count++;
    }
    // Receive processor name from all other processes.
    MPI_Request recvRequests[nProcsTot-1];
    count = 0;
    for(int i = 0; i < nProcsTot; i++) {
      if(i == myRank)
	continue;
      int tag = 0;
      MW._MPI_Irecv(nameList[i], MPI_MAX_PROCESSOR_NAME, 
		   MPI_UNSIGNED_CHAR, i, tag, MPI_COMM_WORLD, &recvRequests[count]);
      count++;
    }
    // Wait for all sends and recvs to complete.
    MPI_Status statuses_send[nProcsTot-1];
    MW._MPI_Waitall(nProcsTot-1, sendRequests, statuses_send);
    MPI_Status statuses_recv[nProcsTot-1];
    MW._MPI_Waitall(nProcsTot-1, recvRequests, statuses_recv);
    // At this point all processes should have all the names.
#if 0
    const int ROLE_UNDEFINED = 1;
    const int ROLE_PARENT    = 2;
    const int ROLE_WORKER    = 3;
    const int ROLE_IDLE      = 4;
#endif
    // Find out how many nodes there are.
    char uniqueNameList[nProcsTot][MPI_MAX_PROCESSOR_NAME];
    int noOfUniqueNames = 0;
    for(int i = 0; i < nProcsTot; i++) {
      // Check if this name is already in uniqueNameList.
      bool found = false;
      for(int j = 0; j < noOfUniqueNames; j++) {
	if(strcmp(nameList[i], uniqueNameList[j]) == 0) {
	  found = true;
	  break;
	}
      }
      if(!found) {
	strcpy(uniqueNameList[noOfUniqueNames], nameList[i]);
	noOfUniqueNames++;
      }
    }
    int noOfNodes = noOfUniqueNames;
    // Determine max number of allowed workers per node.
    int maxNoOfWorkersPerNode = nProcsTot / noOfNodes;
    if(nProcsTot%noOfNodes)
      maxNoOfWorkersPerNode++;
    ProcessRole roleList[nProcsTot];
    for(int i = 0; i < nProcsTot; i++)
      roleList[i] = Undefined;
    bool parentFound = false;
    int workerCount = 0;
    for(int i = 0; i < nProcsTot; i++) {
      // Determine role of process i.
      if(parentFound == false) {
	// Check if this process can be the parent.
	int noOfMatches = 0;
	for(int j = 0; j < nProcsTot; j++) {
	  if(j == i)
	    continue;
	  if(strcmp(nameList[i], nameList[j]) == 0)
	    noOfMatches++;
	}
        // We let this process be the parent if there is at least one other process on the same node.
        // If nProcsTot > nWorkerProcs we also let this node be parent; then the parent process will be alone on that node.
	if(noOfMatches > 0 || nProcsTot > nWorkerProcs) {
	  // OK, this can be the parent!
	  roleList[i] = Parent;
	  parentFound= true;
	}
	else {
	  // Only one process on this node, must be worker.
	  roleList[i] = Worker;
	  workerCount++;
	}
      }
      else {
	// Parent already found. This process must be either worker or
	// idle. 
	if(workerCount == nWorkerProcs) {
	  roleList[i] = Idle;
	}
	else {
	  // Check how many workers have already been assigned to this
	  // node.
	  int workerCountLocal = 0;
	  for(int j = 0; j < nProcsTot; j++) {
	    if(j == i)
	      continue;
	    if(strcmp(nameList[i], nameList[j]) == 0) {
	      if(roleList[j] == Worker)
		workerCountLocal++;
	    }
	  }
	  if(workerCountLocal < maxNoOfWorkersPerNode) {
	    // OK, this can be a worker!
	    roleList[i] = Worker;
	    workerCount++;
	  }
	  else {
	    roleList[i] = Idle;
	  }
	}
      }
    }
    // Now we have assigned roles to all processes.
    if(myRank == 0) {
      for(int i = 0; i < nProcsTot; i++) {
	std::cout << "MPI process with original rank " << std::setw(5) << i << " role: ";
	switch(roleList[i]) {
	case Parent: std::cout << "PARENT"; break;
	case Worker: std::cout << "WORKER"; break;
	case Idle  : std::cout << "IDLE  "  ; break;
	default: throw std::runtime_error("Error: undefined role for process.");
	}
	std::cout << " (node " << nameList[i] << ")" << std::endl;
      } // end for
    } // end if rank 0
    if(parentFound == false)
      throw std::runtime_error("Error: (parentFound == false)");
    workerCount = 0;
    for(int i = 0; i < nProcsTot; i++) {
      if(roleList[i] == Worker)
	workerCount++;
    }
    if(workerCount != nWorkerProcs)
      throw std::runtime_error("Error: (workerCount != nWorkerProcs)");
    ProcessRole myRole = roleList[myRank];
    // Create new intra-communicators.
    int leaderRank_parent = 0;
    int leaderRank_worker = -1;
    int leaderRank_idle = -1;
    for(int i = 0; i < nProcsTot; i++) {
      if(roleList[i] == Worker && leaderRank_worker == -1)
	leaderRank_worker = i;
      if(roleList[i] == Idle && leaderRank_idle == -1)
	leaderRank_idle = i;
    }
    // Note: by passing myRank as "key" argument to MPI_Comm_split()
    // we make sure that the ranks in resultIntraComm will be in the same
    // order as they were in MPI_COMM_WORLD, so we know the local
    // "leader" process has rank 0 locally.
    MW._MPI_Comm_split(MPI_COMM_WORLD, myRole, myRank, &resultIntraComm);

    int tag = 0;
    if(myRole == Parent) {
      MW._MPI_Intercomm_create( resultIntraComm, 0, MPI_COMM_WORLD, 
			       leaderRank_worker, tag, &resultInterComm);
    }
    else if(myRole == Worker) {
      if(resultIntraComm == MPI_COMM_NULL)
	throw std::runtime_error("Error: (resultIntraComm == MPI_COMM_NULL).");
      MW._MPI_Intercomm_create( resultIntraComm, 0, MPI_COMM_WORLD, 
			       leaderRank_parent, tag, &resultInterComm);
    }
    else if(myRole == Idle) {
      // Do nothing here.
    }
    else {
      throw std::runtime_error("Error: bad 'role' value.");
    }
    // Now we have the two communicators ready!

    resultProcessRole = myRole;

  }

}; // end namespace
