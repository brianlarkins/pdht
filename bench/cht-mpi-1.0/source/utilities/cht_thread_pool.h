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
#ifndef CHT_THREAD_POOL_HEADER
#define CHT_THREAD_POOL_HEADER

#include <map>
#include <queue>
#include <list>
#include "Singleton.h"
#include "cht_time.h"
#include "cht_vector.h"

namespace cht {

  class thread_pool : public Singleton<thread_pool> {
  public:
    struct State_stats {
      int counter;
      double time_acc;
      void clear() {
	counter = 0;
	time_acc = 0;
      }
      State_stats() : counter(0), time_acc(0) {}
    };

    static const int beforeJoin = 0;
    static const int afterLeave = 1;
    static const int doingBlockingOperation = 2;
    static const int readyAfterBlockingOperation = 3;
    static const int working = 4;
    static const int waiting = 5;
    static const int numberOfPossibleStates = 6;

    static std::string get_thread_state_name(int ts) {
      switch (ts) {
      case beforeJoin:
	return "before-join               ";
      case afterLeave:
	return "after-leave               ";
      case doingBlockingOperation:
	return "doing-blocking-op         ";
      case readyAfterBlockingOperation:
	return "ready-after-blocking-op   ";
      case working:
	return "working                   ";
      case waiting:
	return "waiting                   ";
      };
      throw std::runtime_error("Invalid thread state");
    }
  private:
    friend class Singleton<thread_pool>;
    
    class ThreadInfo {
    public:
    private:
      size_t pool_id;
      int myState; 
      double timer_checkpoint;
      State_stats stats[numberOfPossibleStates];
    public:
      //pthread_mutex_t mutex;
      pthread_cond_t cond;
      
      size_t getThreadID() {
	return pool_id;
      }
      void changeState(int state) {
	// thread pool mutex must be locked when calling this function!!
	updateStats(); // Note: important to update before changing state.
	stats[state].counter++;
	myState = state;
      }
      void updateStats() {
	// thread pool mutex must be locked when calling this function!!
	stats[myState].time_acc += (cht::get_wall_seconds() - timer_checkpoint);
	timer_checkpoint = cht::get_wall_seconds();
      }
      void getStats(cht::vector<State_stats> & threadStats) {
	// thread pool mutex must be locked when calling this function!!
	updateStats();
	threadStats.setsize(numberOfPossibleStates);
	for (int i = 0; i < numberOfPossibleStates; ++i) 
	  threadStats[i] = stats[i];
      }
      int getState() const { return myState; }
    ThreadInfo(size_t pool_id_) 
      : pool_id(pool_id_), 
	myState(beforeJoin) {
	  //pthread_mutex_init(&mutex, NULL);
	  pthread_cond_init(&cond, NULL);
	  for (int i = 0; i < numberOfPossibleStates; ++i) {
	    stats[i].clear();
	  }
	  stats[beforeJoin].counter++;
	  timer_checkpoint = cht::get_wall_seconds();
	}
    }; // end class ThreadInfo
    
  public:
    typedef std::map<pthread_t, ThreadInfo*> ThreadInfoMap;
    typedef std::queue<pthread_t> ThreadQueue;

    void set_n_working_max(int const n_max); // Call only when infoMap is empty
    void join();
    void leave();
    void yield();
    void startingBlockingOperation();
    void blockingOperationFinished();
    void getThreadStats(pthread_t id, cht::vector<State_stats> & threadStats) {
      pthread_mutex_lock(&mutex_pool); 
      infoMap[id]->getStats(threadStats);
      pthread_mutex_unlock(&mutex_pool); 
    }
    void getAllStats(std::list<cht::vector<State_stats> > & statsList) {
      statsList.clear();
      pthread_mutex_lock(&mutex_pool); 
      for(ThreadInfoMap::iterator it = infoMap.begin(); it != infoMap.end(); it++) {
	cht::vector<State_stats> threadStats;
	(*it).second->getStats(threadStats);
	statsList.push_back(threadStats);
      }
      pthread_mutex_unlock(&mutex_pool); 
    }

    size_t getThreadID() {
      pthread_t myID = pthread_self();
      size_t tmp;      
      pthread_mutex_lock(&mutex_pool); 
      if (infoMap.find(myID) == infoMap.end())
	tmp = 0;
      else
	tmp = infoMap[myID]->getThreadID(); 
      pthread_mutex_unlock(&mutex_pool); 
      return tmp;
    }
  protected:
    size_t threadCounter; // Total number of threads that has joined 
                          // thread pool ever
    int n_working(); // thread pool mutex must be locked when 
                     // calling this function!!
    int n_working_max;
    pthread_mutex_t mutex_pool;
    ThreadInfoMap infoMap;
    ThreadQueue readyAfterBlockingQueue;
    void tryToActivateOneWorker(); // thread pool mutex must be locked when 
                                   // calling this function!!
    thread_pool();
  };

}

#endif
