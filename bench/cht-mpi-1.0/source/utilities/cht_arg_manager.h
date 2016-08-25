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
#ifndef ARG_MANAGER_HEADER
#define ARG_MANAGER_HEADER

#include <iostream>
#include <map>
#include <list>
#include <string>
#include <sstream>
#include <stdexcept>

namespace cht {

  template<typename TObject>
    class arg_manager {
  public:
  private:
    typedef std::map<std::string, std::list<std::string> > TypeArgMap;
    typedef std::map<std::string, std::string > OutputTypeMap;
  public:
    bool registerArgumentTypes(std::string const & type_str, std::list<std::string> const & args, std::string const & output_type_str);

    // Check if argument number ind of type type_str has type arg_str
    bool check(std::string const & type_str, std::string const & arg_str, int const ind);

    // Check if output type of type type_str has type output_type_str
    bool checkOutputType(std::string const & type_str, std::string const & output_type_str);

  private:
    TypeArgMap arg_map;
    OutputTypeMap output_type_map;
    // Stuff needed to make this a singleton.
  private:
    // Prevent creation of undesired arg_manager objects.
    arg_manager();
    arg_manager(arg_manager const &);
    // Assignment doesn't make sense for a singleton
    arg_manager & operator=(arg_manager const &);
  public:
    // Better return a reference 
    static arg_manager& instance();
  };

  template<typename TObject>
    bool arg_manager<TObject>::registerArgumentTypes(std::string const & type_str, 
						     std::list<std::string> const & args,
						     std::string const & output_type_str) {
    if (arg_map.find(type_str) != arg_map.end())
      throw std::runtime_error("Attempt to register argument types for type already"
			       " registered in arg_manager");
    return arg_map.insert( typename TypeArgMap::value_type(type_str, args)).second
      && output_type_map.insert( typename OutputTypeMap::value_type(type_str, output_type_str)).second;
  }

  template<typename TObject>
    bool arg_manager<TObject>::check(std::string const & type_str, 
				     std::string const & arg_str, 
				     int const ind) {
    if (arg_map.find(type_str) == arg_map.end())
      throw std::runtime_error("Attempt to check type not registered"
			       " in arg_manager");
    std::list<std::string> const & arg_strings = arg_map[type_str];
    if (ind >= arg_strings.size() || ind < 0)
      throw std::runtime_error("Bad index in arg_manager<TObject>::check with type_str = " + type_str + " and arg_str = " + arg_str);
    std::list<std::string>::const_iterator it = arg_strings.begin();
    for (int i = 0; i<ind;i++)
      it++;    
    return (*it == arg_str);
  }

  template<typename TObject>
    bool arg_manager<TObject>::checkOutputType(std::string const & type_str, 
					       std::string const & output_type_str) {
    if (output_type_map.find(type_str) == output_type_map.end())
      throw std::runtime_error("checkOutputType error: Attempt to check type not registered in arg_manager");
    return (output_type_map[type_str] == output_type_str);
  }

  template<typename TObject>
    arg_manager<TObject>::arg_manager() 
    {    
    }

  template<typename TObject>
    arg_manager<TObject>& arg_manager<TObject>::instance()
    {
      static arg_manager<TObject> theInstance; // local static variable
      return theInstance;
    }

} // end namespace cht

#endif
