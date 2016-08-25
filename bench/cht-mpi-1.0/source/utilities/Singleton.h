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
#ifndef SINGLETON_HEADER
#define SINGLETON_HEADER

#include <pthread.h>

/*
  Example usage:

  class A : public Singleton<A> {
    friend class Singleton<A>;
   private:
    A();

    // Implement your A here
   public:
    void foo();
  };

  Now, only a single A obejct can be accessed via
  A.instance() which returns a reference to the singleton. 
  
  Call to foo():
    A.instance().foo();

*/

template <typename T>
class Singleton  {
  // Hidden stuff
 private:
  static void create();
  static T* ptr_to_instance;
  // All values allowed for char - 
  // may be important for double checked locking pattern in instance()
  static volatile char ptr_to_instance_is_valid; 
  static pthread_mutex_t singleInstantiationMutex; 
  Singleton(Singleton const &);
 protected:
  Singleton() {} // No instances of Singleton ever!
  
  // Function that can be overridden in the implementation of T,
  // e.g. to check if access is allowed.
  void verifyInstanceAccess() {}
 public:
  static T& instance();
  template<typename Base_class_to_T>
    static Base_class_to_T* instance_ptr();
};



template <typename T>
void Singleton<T>::create() {
  static T theInstance;
  ptr_to_instance = &theInstance;
}

template <typename T>
T& Singleton<T>::instance() {
  if (!ptr_to_instance_is_valid) {
    pthread_mutex_lock(&singleInstantiationMutex);  
    if (!ptr_to_instance_is_valid) {
      create();
      ptr_to_instance_is_valid = 1;
    }
    pthread_mutex_unlock(&singleInstantiationMutex);  
  }
  ptr_to_instance->verifyInstanceAccess();
  return *ptr_to_instance;
}

template <typename T>
template<typename Base_class_to_T>
Base_class_to_T* Singleton<T>::instance_ptr() {
  return &instance();
}


// Initialization of static members

template <typename T>
T* Singleton<T>::ptr_to_instance = 0;

template <typename T>
volatile char Singleton<T>::ptr_to_instance_is_valid = 0;

template <typename T>
pthread_mutex_t Singleton<T>::singleInstantiationMutex = 
  PTHREAD_MUTEX_INITIALIZER;




#endif
