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
#ifndef CHT_TASK_TYPE_MACROS 
#define CHT_TASK_TYPE_MACROS 
#include "utilities/cht_static_check.h"

namespace cht {
  namespace internal {
    bool dummy_register_task_types_fun(bool xx);
    template <typename T> struct get_type;
    template <typename T>
      struct get_type<void (*)(T)> {
      typedef T type;
    };
  }
} 

#define CHT_TASK_INPUT(input_types)		\
  typedef CHT_TYPELIST(input_types) inputTypes; 
#define CHT_TASK_OUTPUT(output_type)					\
  typedef typename cht::internal::get_type< void (*)output_type >::type outputType; 


#define CHT_TASK_TYPE_DECLARATION					\
  public:								\
    static std::string get_class_id();					\
  private:								\
  using Task::execute;							\
    static bool const cht_internal_registered;				\
    cht::ID internal_execute() {					\
      return Exec<inputTypes::length>::exec(this);			\
    }	
	


#define CHT_TASK_TYPE_IMPLEMENTATION(tasktype)				\
  bool const cht::internal::get_type< void (*)tasktype >::type::cht_internal_registered = \
		 cht::registerTaskType<cht::internal::get_type< void (*)tasktype >::type>(); \
  std::string cht::internal::get_type< void (*)tasktype >::type::get_class_id() { \
    cht::internal::dummy_register_task_types_fun(cht_internal_registered); \
    return #tasktype;} 

#define CHT_TASK_TYPE_SPECIALIZATION_IMPLEMENTATION(tasktype)		\
  template<> bool const cht::internal::get_type< void (*)tasktype >::type::cht_internal_registered = \
    cht::registerTaskType< cht::internal::get_type< void (*)tasktype >::type >(); \
  template<> std::string cht::internal::get_type< void (*)tasktype >::type::get_class_id() { \
    cht::internal::dummy_register_task_types_fun(cht_internal_registered); \
    return #tasktype;} 

////// Alternative way of solving the comma problem:
//#define CHT_TASK_TYPE_SPECIALIZATION_IMPLEMENTATION(...)	   \
//  template<> bool const __VA_ARGS__::cht_internal_registered =	\
//    cht::registerTaskType< __VA_ARGS__ >();			\
//  template<> std::string __VA_ARGS__::get_class_id() {			\
//    cht::internal::dummy_register_task_types_fun(cht_internal_registered); \
//    return #__VA_ARGS__;} 

#endif
