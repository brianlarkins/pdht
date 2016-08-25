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
#ifndef GENERIC_FACTORY_HEADER
#define GENERIC_FACTORY_HEADER

#include <map>
#include <list>
#include <string>
#include <sstream>
#include <stdexcept>

namespace cht {

  template<typename TObject>
    class obj_factory {
    template<typename ObjectType> 
      static TObject* createObject() { return new ObjectType; }
  public:
    // createObjectCallback: pointer to a function returning TObject* 
    typedef TObject* (*createObjectCallback)(); 
  private:
    // A more optimized map could possibly speed up lookups if needed
    typedef std::map<std::string, createObjectCallback> CallbackMap;
  public:
    // New TObject types must be registered in startup code  
    // The TObject type ID should be a string with the name of the TObject class
    // Returns true if the registration was succesful
    bool registerObjectType(std::string ObjectTypeID, 
			    createObjectCallback createFunction);
    template<typename TObjectDerived>
      bool registerObjectType();

    template<typename TObjectDerived>
      bool registerSingletonObjectType();

    // Returns true if the Object type was registered before
    bool unregisterObjectType(std::string ObjectTypeID);
    TObject* createObject(std::string ObjectTypeID);
    void getObjectTypeIDStrings(std::list<std::string> & strList) const;
  private:
    CallbackMap callbacks;
    // Stuff needed to make this a singleton.
  private:
    // Prevent creation of undesired object factories
    obj_factory();
    obj_factory(obj_factory const &);
    // Assignment doesn't make sense for a singleton
    obj_factory & operator=(obj_factory const &);
  public:
    // Better return a reference 
    static obj_factory& instance();
  };

  template<typename TObject>
    template<typename TObjectDerived>
    bool obj_factory<TObject>::registerObjectType() {
    std::string objectTypeID = TObjectDerived::get_class_id();
    if ( callbacks.find(objectTypeID) != callbacks.end() )
      throw std::runtime_error("Attempt to register already registered object type '" + objectTypeID + "' to obj_factory<" + TObject::object_base_type_id() + ">");
    return callbacks.insert
      ( typename CallbackMap::value_type(objectTypeID, createObject<TObjectDerived> ) ).second;
  }

  template<typename TObject>
    template<typename TObjectDerived>
    bool obj_factory<TObject>::registerSingletonObjectType() {
    if ( callbacks.find(TObjectDerived::object_type_id()) != callbacks.end() )
      throw std::runtime_error("Attempt to register already registered object type '" + TObjectDerived::object_type_id() + "' to obj_factory<" + TObject::object_base_type_id() + ">");
    return callbacks.insert
      ( typename CallbackMap::value_type(TObjectDerived::object_type_id(), 
					 TObjectDerived::template instance_ptr<TObject>) ).second;
  }

  template<typename TObject>
    bool obj_factory<TObject>::registerObjectType(std::string objectTypeID, 
						  createObjectCallback createFunction) 
    {
      if ( callbacks.find(objectTypeID) != callbacks.end() )
	throw std::runtime_error("Attempt to register already registered object type '" + objectTypeID + "' to obj_factory<" + TObject::object_base_type_id() + ">");
      //typedef typename CallbackMap::value_type myVT;
      return callbacks.insert
	( typename CallbackMap::value_type(objectTypeID, createFunction) ).second;
    }
  
  template<typename TObject>
    bool obj_factory<TObject>::unregisterObjectType(std::string objectTypeID)
    {
      return callbacks.erase( objectTypeID ) == 1;
    }

  // This function will typically be called very often 
  template<typename TObject>
    TObject* obj_factory<TObject>::createObject(std::string objectTypeID) 
    {
      typename CallbackMap::const_iterator i = callbacks.find( objectTypeID );
      if ( i == callbacks.end() )
	{
	  // Object type not found
	  // std::cout << "Error in obj_factory::createObject: Unknown Object type '" << objectTypeID << "'\n";
	  std::stringstream ss;
	  ss << "Error in obj_factory::createObject: Unknown Object type '" << objectTypeID << "'";
	  throw std::runtime_error(ss.str());
	}
      // Invoke the callback creation function
      return (i->second)(); 
    }

  template<typename TObject>
    void obj_factory<TObject>::getObjectTypeIDStrings(std::list<std::string> & strList) const {
    strList.clear();
    typename CallbackMap::const_iterator it;
    for (it = callbacks.begin();it != callbacks.end();it++)
      strList.push_back(it->first);      
  }


  template<typename TObject>
    obj_factory<TObject>::obj_factory() 
    {  
  
    }

  template<typename TObject>
    obj_factory<TObject>& obj_factory<TObject>::instance()
    {
      static obj_factory<TObject> theInstance; // local static variable
      return theInstance;
    }

} // end namespace cht

#endif
