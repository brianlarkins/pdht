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
#include "services/output_service/OutputService_parent.h"
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <cstring>
#include "utilities/cht_utils.h"
#include "utilities/cht_threads.h"

namespace cht {
  namespace OutputService {

    void Parent::LockMutex()
    {
      mutex.lock();
    }

    void Parent::UnlockMutex()
    {
      mutex.unlock();
    }

    void Parent::outputInfo(std::string outText) {
      output(Output::Info, outText);
    }

    void Parent::outputDebug(std::string outText) {
      output(Output::Debug, outText);
    }



    void* global_thread_func(void* arg)
    {
      Parent* p = (Parent*) arg;
      try {
	p->parent_thread_func();
      } catch ( std::runtime_error e ) {
	std::cerr << "========= ERROR ERROR ERROR =========================" << std::endl;
	std::cerr << "Error! Exception std::runtime_error caught in "
	  "cht::OutputService::Parent global_thread_func()." << std::endl;
	std::cerr << "what(): " << e.what() << std::endl;
	std::cerr << "=====================================================" << std::endl;
      } catch ( std::exception e ) {
	std::cerr << "========= ERROR ERROR ERROR =========================" << std::endl;
	std::cerr << "Error! Exception std::exception caught in "
	  "cht::OutputService::Parent global_thread_func()." << std::endl;
	std::cerr << "what(): " << e.what() << std::endl;
	std::cerr << "=====================================================" << std::endl;
      } catch ( ... ) {
	std::cerr << "========= ERROR ERROR ERROR =========================" << std::endl;
	std::cerr << "Error! Exception caught in "
	  "cht::OutputService::Parent global_thread_func()." << std::endl;
	std::cerr << "=====================================================" << std::endl;
      }
      return NULL;
    }

    void Parent::start_derived() {
      AccessKey key(this);
      outFileForText.open("output.txt");
      if ( !outFileForText.good() )
	throw std::runtime_error("Output file stream to file output.txt is not good.");
      time_to_stop = false;
      if (threadPtr)
	throw std::runtime_error("threadPtr != NULL in "
				 "cht::OutputService::Parent::start_derived()");
      threadPtr = new Threads::Thread("OS-parent", global_thread_func, this);
      // Send mode and priority to workers.
      for(int i = 0; i < key.n_workers(); i++) {
	MW._MPI_Send(&mode, sizeof(Output::Mode), 
		    MPI_UNSIGNED_CHAR, i, 
		    Tag_mode, *key.comm_to_workers());
	MW._MPI_Send(&priority, sizeof(Output::Priority), 
		    MPI_UNSIGNED_CHAR, i, 
		    Tag_priority, *key.comm_to_workers());
      }
    }

    void Parent::stop_derived() {
      // Do not close file directly since the other thread may still need
      // to write some messages.
      // Signal exit to executor thread
      LockMutex();
      time_to_stop = true;
      UnlockMutex();
      // Join thread here!!
      delete threadPtr;
      threadPtr = NULL;
      // Now we can close file.
      outFileForText.close();
    }

    void Parent::write_output_to_file(const char* p) {
      // Check if it is a multi-line string.
      if(strchr(p, '\n')) {
	// Multi-line output message.
	if ( !outFileForText.good() )
	  throw std::runtime_error("Multi-line output message: Output file stream to file output.txt is not good.");
	outFileForText << " --- MULTILINE OUTPUT ---" << std::endl;
	while(1) {
	  // Now p points to beginning of new output line. The line
	  // could end with a newline character or with a
	  // end-of-string character.
	  std::string lineString;
	  const char* pp = strchr(p, '\n');
	  if(pp) {
	    lineString = std::string(p, pp-p);
	    p = pp + 1;
	  } 
	  else {
	    lineString = std::string(p);
	    p += lineString.length();
	  }
	  outFileForText << "         " << lineString << std::endl;
	  if(*p == '\0')
	    break;
	}
      }
      else {
	// Single-line output message.
	outFileForText << " " << p << std::endl;
      }
    }

    void Parent::output(Output::Priority prio, std::string outText) {
      if(mode == Output::Disabled)
	return;
      if(prio < priority)
	return;
      LockMutex();
      if(time_to_stop)
	throw std::runtime_error("OutputService::output() (parent) error: "
				 "output() called while service is stopping.");
      std::string timeStr = GetTimeStrWithDecimals();
      if ( !outFileForText.good() )
	throw std::runtime_error("cht::OutputService::Parent::output(...): Output file stream to file output.txt is not good.");

      UnlockMutex();
      size_t threadID = Threads::get_ID_number();
      std::string threadString = Threads::get_ID_string();
      LockMutex();
      outFileForText << "PARENT  " << " "
		     << std::right << std::setw(14) << threadString << "-" 
		     << std::left  << std::setw(7) << threadID
		     << " " << timeStr;
      write_output_to_file(outText.c_str());
      UnlockMutex();
    }

    void Parent::parent_thread_func() {
      AccessKey key(this);   
      int n_workers = key.n_workers(); 
      MPI_Comm* comm_to_workers = key.comm_to_workers();
  
      std::vector<int> messageCounts(n_workers);
      std::vector<int> messageCountsRequired(n_workers);
      for(int i = 0; i < n_workers; i++) {
	messageCounts[i] = 0;
	messageCountsRequired[i] = -1;
      }
      bool finalPhase = false;
      while(1) {
	// Check if time to quit
	if(finalPhase) {
	  // Compare messageCounts to messageCountsRequired. At the same
	  // time, check that the number of received messages does not
	  // exceed the number of required.
	  bool allEqual = true;
	  for(int k = 0; k < n_workers; k++) {
	    if(messageCounts[k] > messageCountsRequired[k])
	      throw std::runtime_error("OutputService::parent_thread_func() "
				       "error: (messageCounts[k] > messageCountsRequired[k]).");
	    if(messageCounts[k] < messageCountsRequired[k])
	      allEqual = false;
	  }
	  if(allEqual) {
	    // Send acknowledgements back to workers.
	    for(int k = 0; k < n_workers; k++) 
	      MW._MPI_Send(NULL, 0, MPI_UNSIGNED_CHAR, k, 
			  Tag_parent_final_confirmation, *comm_to_workers);
	    // OK, we are done!
	    outFileForText << "OutputService::parent_thread_func() message counts match: exiting now!" << std::endl;
	    break;
	  }
	}
	int flag;
	MPI_Status status;
	MW._MPI_Iprobe( MPI_ANY_SOURCE, Tag_flush_output, *comm_to_workers, &flag, 
		       &status);
	if(flag == 0) {
	  // No messages arrived. Check if time to quit.
	  LockMutex();
	  if (time_to_stop && !finalPhase) {
	    // Before closing we need to get message count from each
	    // worker, then continue processing output until the given
	    // number of output messages have been processed for each
	    // worker. However, we cannot use blocking receive to get the
	    // message counts since output messages may still arrive and
	    // then we must still process those.
	    for(int i = 0; i < n_workers; i++) {
	      int flag2;
	      MW._MPI_Iprobe( i, Tag_message_count, *comm_to_workers, &flag2, MPI_STATUS_IGNORE);
	      if(flag2)
		{
		  MW._MPI_Recv(&messageCountsRequired[i], sizeof(int), MPI_UNSIGNED_CHAR, i, 
			      Tag_message_count, *comm_to_workers, MPI_STATUS_IGNORE);
		}
	    }
	    // Check if all messageCounts have arrived.
	    bool allArrived = true;
	    for(int i = 0; i < n_workers; i++) {
	      if(messageCountsRequired[i] == -1)
		allArrived = false;
	    }
	    if(allArrived)
	      finalPhase = true;
	  }
	  UnlockMutex();
	  // Not time to quit yet, just wait a while before checking for
	  // new messages again.
	  usleep(10000);
	}
	else {
	  // Something has arrived from some worker.
	  int size;
	  MW._MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &size);
	  std::vector<char> buf(size);
	  MW._MPI_Recv(&buf[0], size, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, 
		      Tag_flush_output, *comm_to_workers, &status);
	  // Now the message is in buf. It may contain one or more output
	  // messages. Double-check that last char is 0.
	  if(buf[size-1] != '\0')
	    throw std::runtime_error("Error: (buf[size-1] != '\0')");
	  const char* p = &buf[0];
	  LockMutex();
	  while(p < &buf[size]) {
	    // Now p points to beginning of new output message string.
	    if ( !outFileForText.good() )
	      throw std::runtime_error("Output message: Output file stream to file output.txt is not good.");
	    outFileForText << "WRK " << std::setw(4) << status.MPI_SOURCE;
	    write_output_to_file(p);
	    p += strlen(p) + 1;
	    messageCounts[status.MPI_SOURCE]++;
	  }
	  UnlockMutex();
	}
      } // end while
    }

    // Note: the mode and priority variables are not mutex protected,
    // instead we ensure that setMode() is called before the service is
    // started.
    void Parent::setMode(Output::Mode mode_) {
      if ( serviceIsRunning() )
	throw std::runtime_error("Call to cht::OutputService::Parent::setMode while service is running.");
      mode = mode_;
    }
    void Parent::setLevel(Output::Priority priority_) {
      if ( serviceIsRunning() )
	throw std::runtime_error("Call to cht::OutputService::Parent::setLevel while service is running.");
      priority = priority_;
    }

    Parent::Parent() 
      : time_to_stop(false), 
	mode(Output::Standard),
	priority(Output::Info),
	threadPtr(NULL)
    {
    }


  }; // end namespace
}; // end namespace

