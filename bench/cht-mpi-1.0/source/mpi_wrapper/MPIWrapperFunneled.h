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
#ifndef MPI_WRAPPER_FUNNELED_HEADER
#define MPI_WRAPPER_FUNNELED_HEADER

#include <list>
#include <mpi.h>
#include <pthread.h>

#define MPI_VERSION_MPI_WRAPPER_FUNNELED MPI_VERSION
// To ensure use of only mpi version 1 functionality
//#define MPI_VERSION_MPI_WRAPPER_FUNNELED 1

struct Proxy_MPI;
struct Proxy_MPI_Wait_calls;
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
struct Proxy_MPI_Non_blocking_send_calls;
#endif

static void* global_thread_func(void*);

/** MPIWrapperFunneled 
 *  Currently this wrapper funnels all MPI calls to a single pthread.
 *  
 *  REQUIREMENT: MPI implementation supporting at least 
 *               thread level MPI_THREAD_FUNNELED 
 *
 */
class MPIWrapperFunneled
{
  friend void* global_thread_func(void*);

 private:
  int mpi_was_already_initialized;
  pthread_t thread;

  static int const mutex_id_init                     = 0;
  static int const mutex_id_finalize                 = 1;
  static int const mutex_id_non_blocking_list        = 2;
  static int const mutex_id_waitcalls_list           = 3;
  static int const mutex_id_non_blocking_send_list   = 4;
  static int const mutex_id_executor_wakeup_flag     = 5;
  static int const number_of_mutexes                 = 6;
  pthread_mutex_t mutexes[number_of_mutexes];
  pthread_cond_t executor_wakeup_cond;

  void LockMutex(int const mutex_id);
  void UnlockMutex(int const mutex_id);
  
  void executor_thread_func();
  
  // Stuff needed to make this a singleton.
 private:
  MPIWrapperFunneled();
  ~MPIWrapperFunneled();
  MPIWrapperFunneled(MPIWrapperFunneled const &);
  // Assignment doesn't make sense for a singleton
  MPIWrapperFunneled & operator=(MPIWrapperFunneled const &);
 public:
  static MPIWrapperFunneled& instance();

 private:
  Proxy_MPI* init_proxy;
  Proxy_MPI* finalize_proxy;
  std::list<Proxy_MPI*> pendingNonBlocking;
  std::list<Proxy_MPI_Wait_calls*> pendingWaitCalls; 
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  std::list<Proxy_MPI_Non_blocking_send_calls*> pendingNonBlockingSendCalls[2];
#endif
  bool executor_wakeup_flag;
  //  MPI_Request req; // For now don't use request to wake up executor
  inline void wake_up_executor() {
    LockMutex(mutex_id_executor_wakeup_flag);
    executor_wakeup_flag = true;
    pthread_cond_signal(&executor_wakeup_cond);
    UnlockMutex(mutex_id_executor_wakeup_flag);
    //  MPI_Grequest_complete( req );
  }
  

 public:
  
  // Wrappers for all needed MPI routines

  // Blocking calls handled specifically
  
  int _MPI_Bsend(void *buf, int count, MPI_Datatype datatype, 
		int dest, int tag, MPI_Comm comm);
  int _MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, 
	       int tag, MPI_Comm comm, MPI_Status *status);
  int _MPI_Rsend(void *ibuf, int count, MPI_Datatype datatype, int dest, 
		int tag, MPI_Comm comm);
  int _MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, 
	       int tag, MPI_Comm comm);
  int _MPI_Ssend(void *buf, int count, MPI_Datatype datatype, int dest, 
		int tag, MPI_Comm comm);

  // Blocking and non blocking calls not handled in a special way
  
  int _MPI_Abort(MPI_Comm comm, int errorcode);

  // Scali MPI mpi.h definition:
  // int MPI_Accumulate(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,  MPI_Op, MPI_Win);


#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Accumulate(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
		     int target_rank, MPI_Aint target_disp, int target_count,
		     MPI_Datatype target_datatype, MPI_Op op, MPI_Win win); 
  int _MPI_Add_error_class(int *errorclass);
  int _MPI_Add_error_code(int errorclass, int *errorcode);
  int _MPI_Add_error_string(int errorcode, char *string);
#endif
  int _MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
		    void *recvbuf, int recvcount, 
		    MPI_Datatype recvtype, MPI_Comm comm);
  int _MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
		     void *recvbuf, int *recvcounts, 
		     int *displs, MPI_Datatype recvtype, MPI_Comm comm);
  int _MPI_Alloc_mem(MPI_Aint size, MPI_Info info, 
		    void *baseptr);
  int _MPI_Allreduce(void *sendbuf, void *recvbuf, int count, 
		    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm); 
  int _MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
		   void *recvbuf, int recvcount, 
		   MPI_Datatype recvtype, MPI_Comm comm);
  int _MPI_Alltoallv(void *sendbuf, int *sendcounts, int *sdispls, 
		    MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
		    int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Alltoallw(void *sendbuf, int *sendcounts, int *sdispls, MPI_Datatype *sendtypes, 
		    void *recvbuf, int *recvcounts, int *rdispls, MPI_Datatype *recvtypes,
		    MPI_Comm comm);
#endif
  int _MPI_Barrier(MPI_Comm comm);
  int _MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, 
		int root, MPI_Comm comm);
  int _MPI_Bsend_init(void *buf, int count, MPI_Datatype datatype, 
		     int dest, int tag, MPI_Comm comm, MPI_Request *request); 
  int _MPI_Buffer_attach(void *buffer, int size);
  int _MPI_Buffer_detach(void *buffer, int *size);
  int _MPI_Cancel(MPI_Request *request);
  int _MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int *coords);
  int _MPI_Cart_create(MPI_Comm old_comm, int ndims, int *dims, 
		      int *periods, int reorder, MPI_Comm *comm_cart);
  int _MPI_Cart_get(MPI_Comm comm, int maxdims, int *dims, 
		   int *periods, int *coords);
  int _MPI_Cart_map(MPI_Comm comm, int ndims, int *dims, 
		   int *periods, int *newrank);
  int _MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank);
  int _MPI_Cart_shift(MPI_Comm comm, int direction, int disp, 
		     int *rank_source, int *rank_dest);
  int _MPI_Cart_sub(MPI_Comm comm, int *remain_dims, MPI_Comm *new_comm);
  int _MPI_Cartdim_get(MPI_Comm comm, int *ndims);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Close_port(char *port_name);
  int _MPI_Comm_accept(char *port_name, MPI_Info info, int root, 
		      MPI_Comm comm, MPI_Comm *newcomm);
  MPI_Fint _MPI_Comm_c2f(MPI_Comm comm);
  int _MPI_Comm_call_errhandler(MPI_Comm comm, int errorcode);
#endif
  int _MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Comm_connect(char *port_name, MPI_Info info, int root, 
		       MPI_Comm comm, MPI_Comm *newcomm);
#endif
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn, 
			     MPI_Comm_delete_attr_function *comm_delete_attr_fn, 
			     int *comm_keyval, void *extra_state);
#endif
  int _MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval);
  int _MPI_Comm_disconnect(MPI_Comm *comm);
#endif
  int _MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  MPI_Comm _MPI_Comm_f2c(MPI_Fint comm);
  int _MPI_Comm_free_keyval(int *comm_keyval);
#endif
  int _MPI_Comm_free(MPI_Comm *comm);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, 
			void *attribute_val, int *flag);
#endif
  int _MPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *erhandler);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen);
  int _MPI_Comm_get_parent(MPI_Comm *parent);
#endif
  int _MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Comm_join(int fd, MPI_Comm *intercomm);
#endif
  int _MPI_Comm_rank(MPI_Comm comm, int *rank);
  int _MPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group);
  int _MPI_Comm_remote_size(MPI_Comm comm, int *size);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val);
#endif
  int _MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Comm_set_name(MPI_Comm comm, char *comm_name);
#endif
  int _MPI_Comm_size(MPI_Comm comm, int *size);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Comm_spawn(char *command, char **argv, int maxprocs, MPI_Info info, 
		     int root, MPI_Comm comm, MPI_Comm *intercomm, 
		     int *array_of_errcodes);
  int _MPI_Comm_spawn_multiple(int count, char **array_of_commands, char ***array_of_argv, 
			      int *array_of_maxprocs, MPI_Info *array_of_info, 
			      int root, MPI_Comm comm, MPI_Comm *intercomm, 
			      int *array_of_errcodes);
#endif
  int _MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
  int _MPI_Comm_test_inter(MPI_Comm comm, int *flag);
  int _MPI_Dims_create(int nnodes, int ndims, int *dims);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  MPI_Fint _MPI_Errhandler_c2f(MPI_Errhandler errhandler);
#endif
  int _MPI_Errhandler_create(MPI_Handler_function *function, 
			    MPI_Errhandler *errhandler);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  MPI_Errhandler _MPI_Errhandler_f2c(MPI_Fint errhandler);
#endif
  int _MPI_Errhandler_free(MPI_Errhandler *errhandler);
  int _MPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler);
  int _MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler);
  int _MPI_Error_class(int errorcode, int *errorclass);
  int _MPI_Error_string(int errorcode, char *string, int *resultlen);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Exscan(void *sendbuf, void *recvbuf, int count, 
		 MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
  MPI_Fint _MPI_File_c2f(MPI_File file);
  MPI_File _MPI_File_f2c(MPI_Fint file);
  int _MPI_File_call_errhandler(MPI_File fh, int errorcode);
  int _MPI_File_set_errhandler( MPI_File file, MPI_Errhandler errhandler);
  int _MPI_File_get_errhandler( MPI_File file, MPI_Errhandler *errhandler);
  int _MPI_File_open(MPI_Comm comm, char *filename, int amode,
		    MPI_Info info, MPI_File *fh);
  int _MPI_File_close(MPI_File *fh);
  int _MPI_File_delete(char *filename, MPI_Info info);
  int _MPI_File_set_size(MPI_File fh, MPI_Offset size);
  int _MPI_File_preallocate(MPI_File fh, MPI_Offset size);
  int _MPI_File_get_size(MPI_File fh, MPI_Offset *size);
  int _MPI_File_get_group(MPI_File fh, MPI_Group *group);
  int _MPI_File_get_amode(MPI_File fh, int *amode);
  int _MPI_File_set_info(MPI_File fh, MPI_Info info);
  int _MPI_File_get_info(MPI_File fh, MPI_Info *info_used);
  int _MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype,
			MPI_Datatype filetype, char *datarep, MPI_Info info);
  int _MPI_File_get_view(MPI_File fh, MPI_Offset *disp,
			MPI_Datatype *etype, 
			MPI_Datatype *filetype, char *datarep);
  int _MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf,
		       int count, MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
			   int count, MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_write_at(MPI_File fh, MPI_Offset offset, void *buf,
			int count, MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, void *buf,
			    int count, MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_iread_at(MPI_File fh, MPI_Offset offset, void *buf,
			int count, MPI_Datatype datatype, MPI_Request *request);
  int _MPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, void *buf,
			 int count, MPI_Datatype datatype, MPI_Request *request);
  int _MPI_File_read(MPI_File fh, void *buf, int count,
		    MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_read_all(MPI_File fh, void *buf, int count,
			MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_write(MPI_File fh, void *buf, int count,
		     MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_write_all(MPI_File fh, void *buf, int count,
			 MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_iread(MPI_File fh, void *buf, int count,
		     MPI_Datatype datatype, MPI_Request *request);
  int _MPI_File_iwrite(MPI_File fh, void *buf, int count,
		      MPI_Datatype datatype, MPI_Request *request);
  int _MPI_File_seek(MPI_File fh, MPI_Offset offset, int whence);
  int _MPI_File_get_position(MPI_File fh, MPI_Offset *offset);
  int _MPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset,
			       MPI_Offset *disp);
  int _MPI_File_read_shared(MPI_File fh, void *buf, int count,
			   MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_write_shared(MPI_File fh, void *buf, int count,
			    MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_iread_shared(MPI_File fh, void *buf, int count,
			    MPI_Datatype datatype, MPI_Request *request);
  int _MPI_File_iwrite_shared(MPI_File fh, void *buf, int count,
			     MPI_Datatype datatype, MPI_Request *request);
  int _MPI_File_read_ordered(MPI_File fh, void *buf, int count,
			    MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_write_ordered(MPI_File fh, void *buf, int count,
			     MPI_Datatype datatype, MPI_Status *status);
  int _MPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence);
  int _MPI_File_get_position_shared(MPI_File fh, MPI_Offset *offset);
  int _MPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf,
				 int count, MPI_Datatype datatype);
  int _MPI_File_read_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
  int _MPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf,
				  int count, MPI_Datatype datatype);
  int _MPI_File_write_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
  int _MPI_File_read_all_begin(MPI_File fh, void *buf, int count,
			      MPI_Datatype datatype);
  int _MPI_File_read_all_end(MPI_File fh, void *buf, MPI_Status *status);
  int _MPI_File_write_all_begin(MPI_File fh, void *buf, int count,
			       MPI_Datatype datatype);
  int _MPI_File_write_all_end(MPI_File fh, void *buf, MPI_Status *status);
  int _MPI_File_read_ordered_begin(MPI_File fh, void *buf, int count,
				  MPI_Datatype datatype);
  int _MPI_File_read_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
  int _MPI_File_write_ordered_begin(MPI_File fh, void *buf, int count,
				   MPI_Datatype datatype);
  int _MPI_File_write_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
  int _MPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype,
			       MPI_Aint *extent);
  int _MPI_File_set_atomicity(MPI_File fh, int flag);
  int _MPI_File_get_atomicity(MPI_File fh, int *flag);
  int _MPI_File_sync(MPI_File fh);
  /*
   * file functions end
   */
#endif
  int _MPI_Finalize(void);
  int _MPI_Finalized(int *flag);
  int _MPI_Free_mem(void *base);
  int _MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
		 void *recvbuf, int recvcount, MPI_Datatype recvtype, 
		 int root, MPI_Comm comm);
  int _MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
		  void *recvbuf, int *recvcounts, int *displs, 
		  MPI_Datatype recvtype, int root, MPI_Comm comm);
  int _MPI_Get_address(void *location, MPI_Aint *address);
  int _MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);
  int _MPI_Get_elements(MPI_Status *status, MPI_Datatype datatype, int *count);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Get(void *origin_addr, int origin_count, 
	      MPI_Datatype origin_datatype, int target_rank, 
	      MPI_Aint target_disp, int target_count, 
	      MPI_Datatype target_datatype, MPI_Win win);
#endif
  int _MPI_Get_processor_name(char *name, int *resultlen);
  int _MPI_Get_version(int *version, int *subversion);
  int _MPI_Graph_create(MPI_Comm comm_old, int nnodes, int *index, 
		       int *edges, int reorder, MPI_Comm *comm_graph);
  int _MPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, 
		    int *index, int *edges);
  int _MPI_Graph_map(MPI_Comm comm, int nnodes, int *index, int *edges, 
		    int *newrank);
  int _MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors);
  int _MPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, 
			  int *neighbors);
  int _MPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Grequest_complete(MPI_Request request);
  int _MPI_Grequest_start(MPI_Grequest_query_function *query_fn,
			 MPI_Grequest_free_function *free_fn,
			 MPI_Grequest_cancel_function *cancel_fn,
			 void *extra_state, MPI_Request *request);
  MPI_Fint _MPI_Group_c2f(MPI_Group group);
#endif
  int _MPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result);
  int _MPI_Group_difference(MPI_Group group1, MPI_Group group2, 
			   MPI_Group *newgroup);
  int _MPI_Group_excl(MPI_Group group, int n, int *ranks, 
		     MPI_Group *newgroup);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  MPI_Group _MPI_Group_f2c(MPI_Fint group);
#endif
  int _MPI_Group_free(MPI_Group *group);
  int _MPI_Group_incl(MPI_Group group, int n, int *ranks, 
		     MPI_Group *newgroup);
  int _MPI_Group_intersection(MPI_Group group1, MPI_Group group2, 
			     MPI_Group *newgroup);
  int _MPI_Group_rank(MPI_Group group, int *rank);
  int _MPI_Group_size(MPI_Group group, int *size);
  int _MPI_Group_translate_ranks(MPI_Group group1, int n, int *ranks1, 
				MPI_Group group2, int *ranks2);
  int _MPI_Group_union(MPI_Group group1, MPI_Group group2, 
		      MPI_Group *newgroup);
  int _MPI_Ibsend(void *buf, int count, MPI_Datatype datatype, int dest, 
		 int tag, MPI_Comm comm, MPI_Request *request);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  MPI_Fint _MPI_Info_c2f(MPI_Info info);
#endif
  int _MPI_Info_create(MPI_Info *info);
  int _MPI_Info_delete(MPI_Info info, char *key);
  int _MPI_Info_dup(MPI_Info info, MPI_Info *newinfo);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  MPI_Info _MPI_Info_f2c(MPI_Fint info);
#endif
  int _MPI_Info_free(MPI_Info *info);
  int _MPI_Info_get(MPI_Info info, char *key, int valuelen, 
		   char *value, int *flag);
  int _MPI_Info_get_nkeys(MPI_Info info, int *nkeys);
  int _MPI_Info_get_nthkey(MPI_Info info, int n, char *key);
  int _MPI_Info_get_valuelen(MPI_Info info, char *key, int *valuelen, 
			    int *flag);
  int _MPI_Info_set(MPI_Info info, char *key, char *value);
  int _MPI_Init(int *argc, char ***argv);
  int _MPI_Initialized(int *flag);
  int _MPI_Init_thread(int *argc, char ***argv, int required, 
		      int *provided);
  int _MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, 
			   MPI_Comm bridge_comm, int remote_leader, 
			   int tag, MPI_Comm *newintercomm);
  int _MPI_Intercomm_merge(MPI_Comm intercomm, int high, 
			  MPI_Comm *newintercomm);
  int _MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, 
		 MPI_Status *status);
  int _MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, 
		int tag, MPI_Comm comm, MPI_Request *request);
  int _MPI_Irsend(void *buf, int count, MPI_Datatype datatype, int dest, 
		 int tag, MPI_Comm comm, MPI_Request *request);
  int _MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, 
		int tag, MPI_Comm comm, MPI_Request *request);
  int _MPI_Issend(void *buf, int count, MPI_Datatype datatype, int dest, 
		 int tag, MPI_Comm comm, MPI_Request *request);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Is_thread_main(int *flag);
#endif
  int _MPI_Keyval_create(MPI_Copy_function *copy_fn, 
			MPI_Delete_function *delete_fn, 
			int *keyval, void *extra_state);
  int _MPI_Keyval_free(int *keyval);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Lookup_name(char *service_name, MPI_Info info, char *port_name);
  MPI_Fint _MPI_Op_c2f(MPI_Op op); 
#endif
  int _MPI_Op_create(MPI_User_function *function, int commute, MPI_Op *op);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Open_port(MPI_Info info, char *port_name);
  MPI_Op _MPI_Op_f2c(MPI_Fint op);
#endif
  int _MPI_Op_free(MPI_Op *op);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Pack_external(char *datarep, void *inbuf, int incount,
			MPI_Datatype datatype, void *outbuf,
			MPI_Aint outsize, MPI_Aint *position);
  int _MPI_Pack_external_size(char *datarep, int incount, 
			     MPI_Datatype datatype, MPI_Aint *size);
#endif
  int _MPI_Pack(void *inbuf, int incount, MPI_Datatype datatype, 
	       void *outbuf, int outsize, int *position, MPI_Comm comm);
  int _MPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, 
		    int *size);
  int _MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Publish_name(char *service_name, MPI_Info info, 
		       char *port_name);
  int _MPI_Put(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, 
	      int target_rank, MPI_Aint target_disp, int target_count, 
	      MPI_Datatype target_datatype, MPI_Win win);
  int _MPI_Query_thread(int *provided);
#endif
  int _MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source,
		    int tag, MPI_Comm comm, MPI_Request *request);
  int _MPI_Reduce(void *sendbuf, void *recvbuf, int count, 
		 MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
  int _MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts, 
			 MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Register_datarep(char *datarep, 
			   MPI_Datarep_conversion_function *read_conversion_fn,
			   MPI_Datarep_conversion_function *write_conversion_fn,
			   MPI_Datarep_extent_function *dtype_file_extent_fn,
			   void *extra_state);
  MPI_Fint _MPI_Request_c2f(MPI_Request request);
  MPI_Request _MPI_Request_f2c(MPI_Fint request);
#endif
  int _MPI_Request_free(MPI_Request *request);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Request_get_status(MPI_Request request, int *flag, 
			     MPI_Status *status);
#endif
  int _MPI_Rsend_init(void *buf, int count, MPI_Datatype datatype, 
		     int dest, int tag, MPI_Comm comm, 
		     MPI_Request *request);
  int _MPI_Scan(void *sendbuf, void *recvbuf, int count, 
	       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
  int _MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
		  void *recvbuf, int recvcount, MPI_Datatype recvtype, 
		  int root, MPI_Comm comm);
  int _MPI_Scatterv(void *sendbuf, int *sendcounts, int *displs, 
		   MPI_Datatype sendtype, void *recvbuf, int recvcount, 
		   MPI_Datatype recvtype, int root, MPI_Comm comm);
  int _MPI_Send_init(void *buf, int count, MPI_Datatype datatype, 
		    int dest, int tag, MPI_Comm comm, 
		    MPI_Request *request);
  int _MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
		   int dest, int sendtag, void *recvbuf, int recvcount,
		   MPI_Datatype recvtype, int source, int recvtag, 
		   MPI_Comm comm,  MPI_Status *status);
  int _MPI_Sendrecv_replace(void * buf, int count, MPI_Datatype datatype, 
			   int dest, int sendtag, int source, int recvtag,
			   MPI_Comm comm, MPI_Status *status);
  int _MPI_Ssend_init(void *buf, int count, MPI_Datatype datatype, 
		     int dest, int tag, MPI_Comm comm, 
		     MPI_Request *request);
  int _MPI_Start(MPI_Request *request);
  int _MPI_Startall(int count, MPI_Request *array_of_requests);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Status_c2f(MPI_Status *c_status, MPI_Fint *f_status);
  int _MPI_Status_f2c(MPI_Fint *f_status, MPI_Status *c_status);
  int _MPI_Status_set_cancelled(MPI_Status *status, int flag);
  int _MPI_Status_set_elements(MPI_Status *status, MPI_Datatype datatype,
			      int count);
#endif
  int _MPI_Testall(int count, MPI_Request array_of_requests[], int *flag, 
		  MPI_Status array_of_statuses[]);
  int _MPI_Testany(int count, MPI_Request array_of_requests[], int *index, 
		  int *flag, MPI_Status *status);
  int _MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
  int _MPI_Test_cancelled(MPI_Status *status, int *flag);
  int _MPI_Testsome(int incount, MPI_Request array_of_requests[], 
		   int *outcount, int array_of_indices[], 
		   MPI_Status array_of_statuses[]);
  int _MPI_Topo_test(MPI_Comm comm, int *status);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  MPI_Fint _MPI_Type_c2f(MPI_Datatype datatype);
#endif
  int _MPI_Type_commit(MPI_Datatype *type);
  int _MPI_Type_contiguous(int count, MPI_Datatype oldtype, 
			  MPI_Datatype *newtype);
  int _MPI_Type_create_darray(int size, int rank, int ndims, 
			     int gsize_array[], int distrib_array[], 
			     int darg_array[], int psize_array[],
			     int order, MPI_Datatype oldtype, 
			     MPI_Datatype *newtype);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Type_create_f90_complex(int p, int r, MPI_Datatype *newtype);
  int _MPI_Type_create_f90_integer(int r, MPI_Datatype *newtype);
  int _MPI_Type_create_f90_real(int p, int r, MPI_Datatype *newtype);
#endif
  int _MPI_Type_create_hindexed(int count, int array_of_blocklengths[], 
			       MPI_Aint array_of_displacements[], 
			       MPI_Datatype oldtype, 
			       MPI_Datatype *newtype);
  int _MPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, 
			      MPI_Datatype oldtype, 
			      MPI_Datatype *newtype);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn, 
			     MPI_Type_delete_attr_function *type_delete_attr_fn, 
			     int *type_keyval, void *extra_state);
  int _MPI_Type_create_indexed_block(int count, int blocklength,
				    int array_of_displacements[],
				    MPI_Datatype oldtype,
				    MPI_Datatype *newtype);
#endif
  int _MPI_Type_create_struct(int count, int array_of_block_lengths[], 
			     MPI_Aint array_of_displacements[], 
			     MPI_Datatype array_of_types[], 
			     MPI_Datatype *newtype);
  int _MPI_Type_create_subarray(int ndims, int size_array[], int subsize_array[], 
			       int start_array[], int order, 
			       MPI_Datatype oldtype, MPI_Datatype *newtype);
  int _MPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb, 
			      MPI_Aint extent, MPI_Datatype *newtype); 
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Type_delete_attr(MPI_Datatype type, int type_keyval);
#endif
  int _MPI_Type_dup(MPI_Datatype type, MPI_Datatype *newtype);
  int _MPI_Type_extent(MPI_Datatype type, MPI_Aint *extent);
  int _MPI_Type_free(MPI_Datatype *type);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Type_free_keyval(int *type_keyval);
  MPI_Datatype _MPI_Type_f2c(MPI_Fint datatype);
  int _MPI_Type_get_attr(MPI_Datatype type, int type_keyval, 
			void *attribute_val, int *flag);
#endif
  int _MPI_Type_get_contents(MPI_Datatype mtype, int max_integers, 
			    int max_addresses, int max_datatypes, 
			    int array_of_integers[], 
			    MPI_Aint array_of_addresses[], 
			    MPI_Datatype array_of_datatypes[]);
  int _MPI_Type_get_envelope(MPI_Datatype type, int *num_integers, 
			    int *num_addresses, int *num_datatypes, 
			    int *combiner);
  int _MPI_Type_get_extent(MPI_Datatype type, MPI_Aint *lb, 
			  MPI_Aint *extent);
  int _MPI_Type_get_name(MPI_Datatype type, char *type_name, 
			int *resultlen);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb, 
			       MPI_Aint *true_extent);
#endif
  int _MPI_Type_hindexed(int count, int array_of_blocklengths[], 
			MPI_Aint array_of_displacements[], 
			MPI_Datatype oldtype, MPI_Datatype *newtype);
  int _MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, 
		       MPI_Datatype oldtype, MPI_Datatype *newtype);
  int _MPI_Type_indexed(int count, int array_of_blocklengths[], 
		       int array_of_displacements[], 
		       MPI_Datatype oldtype, MPI_Datatype *newtype);
  int _MPI_Type_lb(MPI_Datatype type, MPI_Aint *lb);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Type_match_size(int typeclass, int size, MPI_Datatype *type);
  int _MPI_Type_set_attr(MPI_Datatype type, int type_keyval, 
			void *attr_val);
#endif
  int _MPI_Type_set_name(MPI_Datatype type, char *type_name);
  int _MPI_Type_size(MPI_Datatype type, int *size);
  int _MPI_Type_struct(int count, int array_of_blocklengths[], 
		      MPI_Aint array_of_displacements[], 
		      MPI_Datatype array_of_types[], 
		      MPI_Datatype *newtype);
  int _MPI_Type_ub(MPI_Datatype mtype, MPI_Aint *ub);
  int _MPI_Type_vector(int count, int blocklength, int stride, 
		      MPI_Datatype oldtype, MPI_Datatype *newtype);
  int _MPI_Unpack(void *inbuf, int insize, int *position, 
		 void *outbuf, int outcount, MPI_Datatype datatype, 
		 MPI_Comm comm);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  int _MPI_Unpublish_name(char *service_name, MPI_Info info, char *port_name);
  int _MPI_Unpack_external (char *datarep, void *inbuf, MPI_Aint insize,
			   MPI_Aint *position, void *outbuf, int outcount,
			   MPI_Datatype datatype);
#endif
  int _MPI_Waitall(int count, MPI_Request *array_of_requests, 
		  MPI_Status *array_of_statuses);
  int _MPI_Waitany(int count, MPI_Request *array_of_requests, 
		  int *index, MPI_Status *status);
  int _MPI_Wait(MPI_Request *request, MPI_Status *status);
  int _MPI_Waitsome(int incount, MPI_Request *array_of_requests, 
		   int *outcount, int *array_of_indices, 
		   MPI_Status *array_of_statuses);
#if MPI_VERSION_MPI_WRAPPER_FUNNELED>1
  MPI_Fint _MPI_Win_c2f(MPI_Win win);
  int _MPI_Win_call_errhandler(MPI_Win win, int errorcode);
  int _MPI_Win_complete(MPI_Win win);
  int _MPI_Win_create(void *base, MPI_Aint size, int disp_unit, 
		     MPI_Info info, MPI_Comm comm, MPI_Win *win);
  int _MPI_Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn, 
			    MPI_Win_delete_attr_function *win_delete_attr_fn, 
			    int *win_keyval, void *extra_state);
  int _MPI_Win_delete_attr(MPI_Win win, int win_keyval);
  MPI_Win _MPI_Win_f2c(MPI_Fint win);
  int _MPI_Win_fence(int assert, MPI_Win win);
  int _MPI_Win_free(MPI_Win *win);
  int _MPI_Win_free_keyval(int *win_keyval);
  int _MPI_Win_get_attr(MPI_Win win, int win_keyval, 
		       void *attribute_val, int *flag);
  int _MPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler);
  int _MPI_Win_get_group(MPI_Win win, MPI_Group *group);
  int _MPI_Win_get_name(MPI_Win win, char *win_name, int *resultlen);
  int _MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win);
  int _MPI_Win_post(MPI_Group group, int assert, MPI_Win win);
  int _MPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val);
  int _MPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler);
  int _MPI_Win_set_name(MPI_Win win, char *win_name);
  int _MPI_Win_start(MPI_Group group, int assert, MPI_Win win);
  int _MPI_Win_test(MPI_Win win, int *flag);
  int _MPI_Win_unlock(int rank, MPI_Win win);
  int _MPI_Win_wait(MPI_Win win);
#endif
  double _MPI_Wtick(void);
  double _MPI_Wtime(void);


  // Extra non-MPI routines
  int NonMPI_Iprobe_multiple(int src, const int tagList[], 
			     int nTags, MPI_Comm comm, 
			     int *flag, MPI_Status *stat);

  
};

// Location mpi.h: /usr/lib64/openmpi/include/mpi.h

#endif
