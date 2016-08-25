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
#include "services_utils.h"
#include "utilities/cht_utils.h"
#include <cstring>

namespace cht {
  namespace Service {
    void populateMapsForGivenStrList(std::list<std::string> const & strList,
				     strToIntMap & strToInt,
				     intToStrMap & intToStr) {
      strToInt.clear();
      intToStr.clear();
      int intValue = 0;
      std::list<std::string>::const_iterator it = strList.begin();
      while( it != strList.end() ) {
	std::string strValue = *it;
	strToInt.insert(strToIntMap::value_type(strValue, intValue));
	intToStr.insert(intToStrMap::value_type(intValue, strValue));
	intValue++;
	it++;
      }      
    }

    void sendStrBufToWorkers(std::list<std::string> const & strList, 
			     int const tag,
			     MPI_Comm* comm_to_workers,
			     MPI_Wrapper & MW) {
      // First compute bufSz
      size_t bufSz = 0;
      std::list<std::string>::const_iterator it = strList.begin();
      while( it != strList.end() ) {
	std::string strValue = *it;
	bufSz += strValue.length()+1;
	it++;
      }      
      it = strList.begin();
      cht::vector<char> buf(bufSz);
      int count = 0;
      for (;it != strList.end();it++){
	std::string strValue = *it;
	strcpy(&buf[count], strValue.c_str());
	count += strValue.length()+1;
      }
      // Send buffer to workers
      int n_workers; 
      MW._MPI_Comm_remote_size(*comm_to_workers, &n_workers);
      // FIXME: Change to BCast
      for(int i = 0; i < n_workers; i++)
	MW._MPI_Send(&buf[0], bufSz, MPI_UNSIGNED_CHAR, i, tag, *comm_to_workers);
    }
    
    void receiveStrBufFromParent(std::list<std::string> & strList, 
				 int const tag,
				 MPI_Comm* comm_to_parent,
				 MPI_Wrapper & MW) {
      // Probe to get message size
      MPI_Status status;
      MW._MPI_Probe(0, tag, *comm_to_parent, &status);
      int msgSize;
      MW._MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &msgSize);
      cht::vector<char> buf(msgSize);
      MW._MPI_Recv(&buf[0], msgSize, MPI_UNSIGNED_CHAR, 0, 
		  tag, *comm_to_parent, &status);      
      if (buf[msgSize-1] != 0)
	throw std::runtime_error("receiveStrBufFromParent(...) : Buffer with strings should end with 0 (the end of the last string).");
      strList.clear();
      int count = 0;
      while (count < msgSize) {
	std::string strValue = &buf[count];
	strList.push_back(strValue);
	count += strValue.length()+1;
      }
    }

    void checkStrListAgainstMap(std::list<std::string> const & strList,
				strToIntMap const & strToInt) {
      if ( strList.size() != strToInt.size() )
	throw std::runtime_error("checkStrListAgainstMap(...) : Number of strings in list does not match number of strings in map.");
      std::list<std::string>::const_iterator it = strList.begin();
      while( it != strList.end() ) {
	std::string strValue = *it;
	if (strToInt.find(strValue) == strToInt.end())
	  throw std::runtime_error("checkStrListAgainstMap(...) : Strings in list do not match strings in map.");
	it++;
      }
      strToIntMap ::const_iterator it2 = strToInt.begin();
      while ( it2 != strToInt.end() ) {
	std::string strValue = it2->first;
	it = strList.begin();
	while( *it != strValue ) {
	  it++;
	  if (it == strList.end())
	    throw std::runtime_error("checkStrListAgainstMap(...) : Strings in list do not match strings in map (second loop).");
	}
	it2++;
      }
    }
    
  }; // end namespace
}; // end namespace
