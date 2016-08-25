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
#include "services/output_service/OutputService_worker.h"
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <cstring>
#include "utilities/cht_utils.h"
#include "services/output_service/OutputService.h"

namespace cht {
  namespace OutputService {

    void Worker::LockMutex()
    {
      mutex.lock();
    }

    void Worker::UnlockMutex()
    {
      mutex.unlock();
    }

    void Worker::outputInfo(std::string outText) {
      output(Output::Info, outText);
    }

    void Worker::outputDebug(std::string outText) {
      output(Output::Debug, outText);
    }


    void Worker::start_derived() {
      AccessKey key(this);
      messageCounter = 0;
      worker_stopping_flag = false;
      // Receive mode and priority from parent.
      MW._MPI_Recv(&mode, sizeof(Output::Mode), MPI_UNSIGNED_CHAR, 0, 
		  Tag_mode, *key.comm_to_parent(), MPI_STATUS_IGNORE);
      MW._MPI_Recv(&priority, sizeof(Output::Priority), MPI_UNSIGNED_CHAR, 0, 
		  Tag_priority, *key.comm_to_parent(), MPI_STATUS_IGNORE);
      outputInfo("cht::OutputService::Worker::start(): this is probably the first output message from this worker.");
      std::stringstream s;
      s << "Mode: ";
      switch (mode) {
      case Output::Standard:
	s << "Standard";
	break;
      case Output::AllInTheEnd:
	s << "AllInTheEnd";
	break;
      case Output::Disabled:
	s << "Disabled";
	break;
      default:
	throw std::runtime_error("Unknown mode");
      }
      s << "   Priority: ";
      switch (priority) {
      case Output::Debug:
	s << "Debug";
	break;
      case Output::Info:
	s << "Info";
	break;
      case Output::High:
	s << "High";
	break;
      default:
	throw std::runtime_error("Unknown mode");
      }
      outputInfo(s.str());
      char procName[MPI_MAX_PROCESSOR_NAME];
      int resultlen;
      MW._MPI_Get_processor_name(procName, &resultlen);
      std::stringstream ss;
      ss << "cht::OutputService::Worker::start(): MPI_Get_processor_name() returned '" << procName << "'";
      outputInfo(ss.str());
    }

    void Worker::stop_derived() {
      AccessKey key(this);
      // Before exit, make sure that all messages have been output by
      // parent. Worker sends a message to parent indicating how many
      // output messages have been sent. Parent is expected to send a
      // confirmation message when that number of output messages have
      // really been processed.
      std::stringstream s;
      LockMutex();
      s << "Message count (including this final message): " << messageCounter+1;
      UnlockMutex();
      outputInfo(s.str());
      LockMutex();
      // Set flag to indicate to other threads that stop() has been called
      // so they know that any call to output() is now illegal.
      worker_stopping_flag = true;
      int noOfMessagesInList = messageList.size();
      UnlockMutex();
      if(noOfMessagesInList > 0)
	flushOutputToParent();
      LockMutex();
      MW._MPI_Send(&messageCounter, sizeof(int), 
		  MPI_UNSIGNED_CHAR, 0, 
		  Tag_message_count, *key.comm_to_parent());
      UnlockMutex();
      MW._MPI_Recv(NULL, 0, MPI_UNSIGNED_CHAR, 0, 
		  Tag_parent_final_confirmation, *key.comm_to_parent(), MPI_STATUS_IGNORE);
    }

    // Note: it is important that the whole flushOutputToParent()
    // implementation is inside the mutex lock, otherwise it can happen
    // that messageList is flushed again by another thread before it has
    // been cleared, leading to mismatch between the message count and the
    // actual number of messages sent to the parent.
    void Worker::flushOutputToParent() {
      AccessKey key(this);
      LockMutex();
      // Compute total buffer size needed.
      std::list<std::string>::iterator it1 = messageList.begin();
      int bufSz = 0;
      while(it1 != messageList.end()) {
	bufSz += (*it1).length() + 1;
	it1++;
      }
      if(bufSz == 0) {
	// This can happen if someone else has flushed the output
	// already. In that case we return here. This avoids problem with
	// trying to access msgBuf[bufSz-1] (negative index) which
	// sometimes caused program crash before this check was added.
	UnlockMutex();
	return;
      }
      // OK, now we have bufSz. Create buffer.
      std::vector<char> msgBuf(bufSz);
      // Copy strings into buffer.
      char* p = &msgBuf[0];
      std::list<std::string>::iterator it2 = messageList.begin();
      while(it2 != messageList.end()) {
	strcpy(p, (*it2).c_str());
	p += (*it2).length() + 1;
	it2++;
      }
      // Double check end of buffer.
      if(msgBuf[bufSz-1] != '\0')
	throw std::runtime_error("Error: (msgBuf[bufSz-1] != '\0')");
      MW._MPI_Send(&msgBuf[0], bufSz, 
		  MPI_UNSIGNED_CHAR, 0, 
		  Tag_flush_output, *key.comm_to_parent());
      messageList.clear();
      UnlockMutex();
    }

    // public version of output function
    void Worker::output(Output::Priority prio, std::string outText_) {
      AccessKey key(this);
      if(mode == Output::Disabled)
	return;
      if(prio < priority)
	return;    
      LockMutex();
      if(worker_stopping_flag)
	throw std::runtime_error("cht::OutputService::Worker::output() called while stopping. Message: " + outText_);
      UnlockMutex();
      std::string timeStr = GetTimeStrWithDecimals();

      size_t threadID = Threads::get_ID_number();
      std::string threadString = Threads::get_ID_string();
      std::stringstream ss;
      ss  << std::right << std::setw(14) << threadString << "-" 
	  << std::left  << std::setw(7)  << threadID
	  << " " << timeStr << " " << outText_;
      std::string msgString = ss.str();
      // Always store in list, regadless of mode.
      LockMutex();
      messageList.push_back(msgString);
      UnlockMutex();
      // Now send message directly or keep in list, depending on mode.
      if(mode == Output::Standard) {
	flushOutputToParent();
      }
      else {
	// Do nothing, it will hopefully be flushed later on.
      }
      LockMutex();
      messageCounter++;
      UnlockMutex();
    }

    Worker::Worker() 
      : mode(Output::Standard),
	priority(Output::Info)
    {
    }

  }; // end namespace
}; // end namespace

