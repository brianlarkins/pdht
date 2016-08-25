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
#ifndef CHT_SHARED_POINTER_HEADER
#define CHT_SHARED_POINTER_HEADER

#include <memory>
#include "cht_threads.h"

namespace cht {
  typedef void (*OnRefCountChangeCallback)( void*, size_t, int, int );
  
  class Ref_count {
  public:
  Ref_count() : count(0), callbackFun(0) {
    }
    void count_increment() {
      mutex.lock();
      count++;
      int tmpCount = count;
      void* tmpContextPtr = contextPtrForCallback;
      size_t tmpId = idForCallback;
      bool callbackFunExists = callbackFun ? true : false;
      mutex.unlock();
      if (callbackFunExists)
        callbackFun(tmpContextPtr, tmpId, tmpCount, tmpCount - 1);
    }
    void count_decrement() {
      mutex.lock();
      count--;
      int tmpCount = count;
      void* tmpContextPtr = contextPtrForCallback;
      size_t tmpId = idForCallback;
      bool callbackFunExists = callbackFun ? true : false;
      mutex.unlock();
      if (callbackFunExists)
        callbackFun(tmpContextPtr, tmpId, tmpCount, tmpCount + 1);
    }
    int get_count() {
      mutex.lock();
      int tmpCount = count;
      mutex.unlock();
      return tmpCount;
    }
    void setRefCountChangeCallback(OnRefCountChangeCallback fun, 
				   void* contextPtr,
				   size_t id) {
      mutex.lock();
      if(callbackFun)
        throw std::runtime_error("Error in cht Ref_count: "
				 "callbackFun already set."); 
      callbackFun = fun;
      idForCallback = id;
      contextPtrForCallback = contextPtr;
      mutex.unlock();
    }
    size_t get_callback_id() {
      mutex.lock();
      size_t rc = idForCallback;
      mutex.unlock();
      return rc;
    }
  private:
    Ref_count(Ref_count const &);
    Ref_count & operator= (Ref_count const & other);
  
    Threads::Mutex mutex;
    int count;
    OnRefCountChangeCallback callbackFun;
    size_t idForCallback;
    void* contextPtrForCallback;
  };



  template <typename T>
    class shared_ptr {
  public:
    
    template<typename U>
      operator shared_ptr<U>() {
      return shared_ptr<U>(pointer, ref);
    }

    Ref_count* getRefCountPtr() const {
      return ref;
    }
    
    explicit shared_ptr (T * p, Ref_count * ref_)
      : pointer(p), ref(ref_) {
      //pointer = p;
      if (pointer != 0) { ref->count_increment(); }
    }
  
    explicit shared_ptr (T * p = 0)
      : pointer(p), ref(0) {
      if (pointer != 0) { 
	ref = new Ref_count;
	ref->count_increment(); 
      }
    }

    ~shared_ptr () {
      if (pointer != 0) {
	if (ref->get_count() <= 1) {
	  /* Last reference, we need to delete the object. */
	  delete pointer;
	  delete ref;
	}
	else {
	  /* There are still other references, we only decrement the reference
	   * count. */
	  ref->count_decrement();
	}
      }
    }

    T & operator*  () const { return *pointer; }
    T * operator-> () const { return  pointer; }

    /** Copy constructor. */
  shared_ptr(shared_ptr const & other)
    : pointer(other.pointer), ref(other.ref) {
      /* Increment reference counter. */
      if ( pointer )
	ref->count_increment();
    }
  
    /** Assignment operator, assigning shared_ptr to shared_ptr.
     */
    shared_ptr & operator= 
      (shared_ptr const & other) {
      
      /* Begin with incrementing "other"'s reference counter.
	 (important in case of self-assignment) */
      if ( other.pointer )
	other.ref->count_increment();
    
      /* Decrement reference counter. */
      if ( pointer ) {
	if (ref->get_count() <= 1) {
	  /* No more references, delete object. */
	  delete pointer;
	  delete ref;
	}
	else { 
	  ref->count_decrement();
	}
      }
      
      /* Point to other. */
      pointer = other.pointer;
      ref = other.ref;
      return *this;
    }

    /** Assignment operator, assigning ptr to shared_ptr.
     */
    shared_ptr & operator= 
      (T * const ptr) {
      
      /* Decrement reference counter. */
      if ( pointer ) {
	if (ref->get_count() <= 1) {
	  /* No more references, delete object. */
	  delete pointer;
	  delete ref;
	}
	else { 
	  ref->count_decrement();
	  ref = 0;
	}
      }
      pointer = ptr;
      if (pointer != 0) { 
	ref = new Ref_count;
	ref->count_increment();
      }
      return *this;
    }

    bool operator==(T const * const ptr) const {
      return (pointer == ptr);
    }
    
    void setRefCountChangeCallback(OnRefCountChangeCallback fun, 
				   void* contextPtr,
				   size_t id) {
      if (ref == 0)
	throw std::runtime_error("Attempt to set RefCountChangeCallback "
				 "when ref == 0 ");
      ref->setRefCountChangeCallback(fun, contextPtr, id);
    }

    T* get() const { return pointer; } 
    size_t get_callback_id() const { return ref->get_callback_id(); }
  protected:
    Ref_count* ref;
    T* pointer;
  };


}; // end namespace

#endif
