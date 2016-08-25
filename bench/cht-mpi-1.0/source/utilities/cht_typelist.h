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
#ifndef CHT_TYPELIST_HEADER
#define CHT_TYPELIST_HEADER

#include <list>
#include <vector>
#include <string>
#include <cassert>
#include <typeinfo>
#include "cht_shared_ptr.h"
namespace cht {

  // TYPELIST
  struct null_type {};
  struct null_typelist {
    static int const length = 0;
    static bool type_check(std::list<std::string> x) {
      return x.empty();
    }
    static void get_type_strings(std::list<std::string> & x) {
      return; 
    }
    template<int i>
    struct get_type;

    template<typename Type>
    static bool type_check(int ind) {
      throw std::runtime_error("Call to bool typelist::typecheck(int) with too large integer input.");
    }
  };

  template<int i>
    struct null_typelist::get_type {
    typedef null_type type;
  };


  template<typename H, typename T>
    struct typelist {
      static int const length = T::length+1;;
      typedef H head;
      typedef T tail;
      static bool type_check(std::list<std::string> & x) {
	bool head_ok = head::get_class_id() == x.front();
	x.pop_front();
	return head_ok && tail::type_check(x);
      }
      static void get_type_strings(std::list<std::string> & x) {
	x.push_back(head::get_class_id());
	return tail::get_type_strings(x);
      }
      template<int i, int dummy = 0>
	struct get_type;

      template<typename Type>
      static bool type_check(int ind) {
	if (ind == 0)
	  return head::get_class_id() == Type::get_class_id();
	else
	  return tail::template type_check<Type>(ind-1);
      }

      template<int i, typename base_type>
	static typename get_type<i>::type const & 
	getElement(std::vector<base_type const * > const & objList) {
	// Cannot do assert on typeid since base class may not be virtual 
	// assert(typeid(*objList.get_element(i)) == typeid(typename get_type<i>::type));
	return *((typename get_type<i>::type*)&(*(objList[i])));
      }
    }; // end struct typelist

  template<typename H, typename T>
    template<int i, int dummy>
    struct typelist<H,T>::get_type {
    typedef typename T::template get_type<i-1>::type type;
  };
  template<typename H, typename T>
    template<int dummy>
    struct typelist<H,T>::get_type<0,dummy> {
    typedef H type;
  };

  template <typename Fun> struct cons;

  template <typename T1> 
    struct cons<void (*)(T1)> {
    typedef typelist<T1, null_typelist> type;
  };

  template <typename T1, typename T2> 
    struct cons<void (*)(T1,T2)> {
    typedef typelist<T1, typelist<T2,null_typelist> > type;
  };

  template <typename T1, typename T2, typename T3> 
    struct cons<void (*)(T1,T2,T3)> {
    typedef typelist<T1, typelist<T2,typelist<T3,null_typelist> > > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4> 
    struct cons<void (*)(T1,T2,T3,T4)> {
    typedef typelist<T1, typelist<T2,typelist<T3,typelist<T4,null_typelist> > > > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5> 
    struct cons<void (*)(T1,T2,T3,T4,T5)> {
    typedef typelist<T1, typelist<T2,typelist<T3,typelist<T4,typelist<T5,null_typelist> > > > > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6> 
    struct cons<void (*)(T1,T2,T3,T4,T5,T6)> {
    typedef typelist<T1, typelist<T2,typelist<T3,typelist<T4,typelist<T5,typelist<T6,null_typelist> > > > > > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7> 
    struct cons<void (*)(T1,T2,T3,T4,T5,T6,T7)> {
    typedef typelist<T1, typelist<T2,typelist<T3,typelist<T4,typelist<T5,typelist<T6,typelist<T7,null_typelist> > > > > > > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8> 
    struct cons<void (*)(T1,T2,T3,T4,T5,T6,T7,T8)> {
    typedef typelist<T1, typelist<T2,typelist<T3,typelist<T4,typelist<T5,typelist<T6,typelist<T7,typelist<T8,null_typelist> > > > > > > > type;
  };
    
  
} // end namespace cht

#define CHT_TYPELIST(a) typename cht::cons< void (*)a >::type 

#endif
