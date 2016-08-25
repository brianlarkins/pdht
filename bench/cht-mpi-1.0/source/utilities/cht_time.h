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
#ifndef CHT_TIME_HEADER
#define CHT_TIME_HEADER

#include <sys/time.h>
#include <sys/resource.h>
#include <stdexcept>
#include <time.h>

namespace cht {
  
  static double get_wall_seconds() {
    struct timeval tv;
    if(gettimeofday(&tv, NULL) != 0)
      throw std::runtime_error("Error in cht::get_wall_seconds(), in gettimeofday().");
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
  }

  static double get_cpu_seconds_clock() {
    return (double)clock() / CLOCKS_PER_SEC;
  }

  static double get_cpu_seconds_user() {
    struct rusage ru;
    if(getrusage(RUSAGE_SELF, &ru) != 0)
      throw std::runtime_error("Error in cht::get_cpu_seconds_user(), in getrusage().");
    double seconds = ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
    return seconds;
  }

  static double get_cpu_seconds_system() {
    struct rusage ru;
    if(getrusage(RUSAGE_SELF, &ru) != 0)
      throw std::runtime_error("Error in cht::get_cpu_seconds_system(), in getrusage().");
    double seconds = ru.ru_stime.tv_sec + (double)ru.ru_stime.tv_usec / 1000000;
    return seconds;
  }

  struct timer {
  private:
    double wall_seconds;
    double cpu_seconds_clock;
    double cpu_seconds_user;
    double cpu_seconds_system;
  public:
    void reset() {
      wall_seconds = get_wall_seconds();
      cpu_seconds_clock = get_cpu_seconds_clock();
      cpu_seconds_user = get_cpu_seconds_user();
      cpu_seconds_system = get_cpu_seconds_system();
    }
    timer() { 
      reset();
    }
    double get_elapsed_wall_seconds() const { 
      return get_wall_seconds() - wall_seconds; 
    }
    double get_elapsed_cpu_seconds_clock() const { 
      return get_cpu_seconds_clock() - cpu_seconds_clock;
    }
    double get_elapsed_cpu_seconds_user() const { 
      return get_cpu_seconds_user() - cpu_seconds_user;
    }
    double get_elapsed_cpu_seconds_system() const {
      return get_cpu_seconds_system() - cpu_seconds_system;
    }
  };

  struct simple_statistics {
    int count;
    double time_min;
    double time_max;
    double time_tot;
    void clear() {
      count = 0;
      time_min = -1; // does not matter, set first time.
      time_max = -1; // does not matter, set first time.
      time_tot = 0;
    }
    void add_time(double t) {
      time_tot += t;
      if(count == 0) {
	time_min = t;
	time_max = t;
      }
      else {
	time_min = (t < time_min) ? t : time_min;
	time_max = (t > time_max) ? t : time_max;
      }
      count++;
    }
    simple_statistics() {
      clear();
    }
  };

  struct work_statistics {
    int count;
    double wall_seconds_tot;
    double cpu_clock_seconds_tot;
    double cpu_user_seconds_tot;
    double cpu_system_seconds_tot;
    double wall_seconds_min;
    double wall_seconds_max;
    work_statistics() { clear(); }
    void clear() {
      count = 0; 
      wall_seconds_tot = 0;
      cpu_clock_seconds_tot = 0;
      cpu_user_seconds_tot = 0;
      cpu_system_seconds_tot = 0;
      wall_seconds_min = 0;
      wall_seconds_max = 0;
    }
    void add_timings(const cht::timer & timer) {
      count++;
      double wallsecs = timer.get_elapsed_wall_seconds();
      wall_seconds_tot += wallsecs;
      if(count == 1) {
	// First time
	wall_seconds_min = wallsecs;
	wall_seconds_max = wallsecs;
      }
      else {
	wall_seconds_min = (wallsecs < wall_seconds_min) ? wallsecs : wall_seconds_min;
	wall_seconds_max = (wallsecs > wall_seconds_max) ? wallsecs : wall_seconds_max;
      }
      cpu_clock_seconds_tot += timer.get_elapsed_cpu_seconds_clock();
      cpu_user_seconds_tot += timer.get_elapsed_cpu_seconds_user();
      cpu_system_seconds_tot += timer.get_elapsed_cpu_seconds_system();
    }
  };


}; // end namespace cht

#endif
