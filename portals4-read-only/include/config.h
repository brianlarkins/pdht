/* include/config.h.  Generated from config.h.in by configure.  */
/* include/config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if bitfields are in forward order */
/* #undef BITFIELD_ORDER_FORWARD */

/* Define if bitfields are in reverse order */
#define BITFIELD_ORDER_REVERSE 1

/* The cacheline width */
#define CACHELINE_WIDTH 64

/* Define if `calloc'ing more than 16 bytes always returns a 16-byte-aligned
   address (common practice on MacOS X). */
#define HAVE_16ALIGNED_CALLOC 1

/* Define if `malloc'ing more than 16 bytes always returns a 16-byte-aligned
   address (common practice on MacOS X). */
#define HAVE_16ALIGNED_MALLOC 1

/* Define if `calloc'ing more than 8 bytes always returns a 8-byte-aligned
   address (common practice on most libcs). */
#define HAVE_8ALIGNED_CALLOC 1

/* Define if `malloc'ing more than 8 bytes always returns a 8-byte-aligned
   address (common practice on most libcs). */
#define HAVE_8ALIGNED_MALLOC 1

/* Define to 1 if you have the <arpa/inet.h> header file. */
#define HAVE_ARPA_INET_H 1

/* if the compiler and cpu can both handle the 128-bit CMPXCHG16B instruction
   */
#define HAVE_CMPXCHG16B 1

/* Define to 1 if you have the declaration of `strerror_r', and to 0 if you
   don't. */
#define HAVE_DECL_STRERROR_R 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `dlsym' function. */
#define HAVE_DLSYM 1

/* Define to 1 if you have the <endian.h> header file. */
#define HAVE_ENDIAN_H 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the `fork' function. */
#define HAVE_FORK 1

/* Define to 1 if you have the `ftruncate' function. */
#define HAVE_FTRUNCATE 1

/* Define to 1 if you have the `getpagesize' function. */
#define HAVE_GETPAGESIZE 1

/* define if clock_gettime is available */
#define HAVE_GETTIME_TIMER 1

/* Define if hwloc is available */
/* #undef HAVE_HWLOC */

/* Define to 1 if you have the <ia32intrin.h> header file. */
/* #undef HAVE_IA32INTRIN_H */

/* Define to 1 if you have the <ia64intrin.h> header file. */
/* #undef HAVE_IA64INTRIN_H */

/* Define to 1 if you have the `inet_ntoa' function. */
#define HAVE_INET_NTOA 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Defined to 1 if kitten support requested */
/* #undef HAVE_KITTEN */

/* Define to 1 if you have the <knem_io.h> header file. */
/* #undef HAVE_KNEM_IO_H */

/* Define to 1 if you have the `bsd-compat' library (-lbsd-compat). */
/* #undef HAVE_LIBBSD_COMPAT */

/* Define to 1 if you have the `dl' library (-ldl). */
#define HAVE_LIBDL 1

/* Define to 1 if you have the `xpmem' library (-lxpmem). */
/* #undef HAVE_LIBXPMEM */

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the `linux/ioctl.h' function. */
/* #undef HAVE_LINUX_IOCTL_H */

/* Define to 1 if you have the <linux/mmtimer.h> header file. */
#define HAVE_LINUX_MMTIMER_H 1

/* define is mach_absolute_time is available */
/* #undef HAVE_MACH_TIMER */

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the `memalign' function. */
#define HAVE_MEMALIGN 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 if you have a working `mmap' system call. */
#define HAVE_MMAP 1

/* Define to 1 if you have the `munmap' function. */
#define HAVE_MUNMAP 1

/* Define to 1 if you have the <netinet/in.h> header file. */
#define HAVE_NETINET_IN_H 1

/* Define to 1 if you have the `posix_memalign' function. */
/* #undef HAVE_POSIX_MEMALIGN */

/* Define if PTHREAD_PROCESS_SHARED attribute on mutexes and cond variables
   works */
#define HAVE_PTHREAD_SHMEM_LOCKS 1

/* Define if pthread supports spinlocks */
#define HAVE_PTHREAD_SPIN_INIT 1

/* Define to 1 if you have the `select' function. */
#define HAVE_SELECT 1

/* Define to 1 if you have the `setenv' function. */
#define HAVE_SETENV 1

/* Define to 1 if you have the <sn/mmtimer.h> header file. */
/* #undef HAVE_SN_MMTIMER_H */

/* Define to 1 if you have the `socket' function. */
#define HAVE_SOCKET 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strerror' function. */
#define HAVE_STRERROR 1

/* Define to 1 if you have the `strerror_r' function. */
#define HAVE_STRERROR_R 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strtol' function. */
#define HAVE_STRTOL 1

/* Define to 1 if you have the `strtoul' function. */
#define HAVE_STRTOUL 1

/* Define to 1 if you have the `syscall' function. */
#define HAVE_SYSCALL 1

/* Define to 1 if you have the <syscall.h> header file. */
#define HAVE_SYSCALL_H 1

/* Define to 1 if you have the <sys/file.h> header file. */
#define HAVE_SYS_FILE_H 1

/* Define to 1 if you have the <sys/ioctl.h> header file. */
#define HAVE_SYS_IOCTL_H 1

/* Define to 1 if you have the <sys/mman.h> header file. */
#define HAVE_SYS_MMAN_H 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/posix_shm.h> header file. */
/* #undef HAVE_SYS_POSIX_SHM_H */

/* Define to 1 if you have the <sys/socket.h> header file. */
#define HAVE_SYS_SOCKET_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the `tdestroy' function. */
#define HAVE_TDESTROY 1

/* Define to 1 if the system has the type `uint_fast32_t'. */
#define HAVE_UINT_FAST32_T 1

/* Define to 1 if the system has the type `uint_fast64_t'. */
#define HAVE_UINT_FAST64_T 1

/* Define to 1 if the system has the type `uint_fast8_t'. */
#define HAVE_UINT_FAST8_T 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vfork' function. */
#define HAVE_VFORK 1

/* Define to 1 if you have the <vfork.h> header file. */
/* #undef HAVE_VFORK_H */

/* Define to 1 if `fork' works. */
#define HAVE_WORKING_FORK 1

/* Define if pthread_*attr_setpshared(attr, PTHREAD_PROCESS_SHARED) works. */
#define HAVE_WORKING_PTHREAD_PROCESS_SHARED 1

/* Define if `valloc'ed memory can be `free'd. */
#define HAVE_WORKING_VALLOC 1

/* Define to 1 if `vfork' works. */
#define HAVE_WORKING_VFORK 1

/* Define to 1 if you have the <xpmem.h> header file. */
/* #undef HAVE_XPMEM_H */

/* Define to 1 if you have the `__mmap' function. */
/* #undef HAVE___MMAP */

/* Define to 1 if you have the `__munmap' function. */
/* #undef HAVE___MUNMAP */

/* libev headers in a weird spot */
/* #undef LIBEV_INC_PREFIXED */

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#define LT_OBJDIR ".libs/"

/* Define to disable argument checking */
/* #undef NO_ARG_VALIDATION */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "portals4-devel@googlegroups.com"

/* Define to the full name of this package. */
#define PACKAGE_NAME "portals4"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "portals4 1.0a1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "portals4"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0a1"

/* Defined to 1 if PMI implementation is Portals4. */
/* #undef PMI_PORTALS4 */

/* Defined to 1 if PMI implementation is SLURM. */
/* #undef PMI_SLURM */

/* if the compiler supports __attribute__((noreturn)) */
#define Q_NORETURN 

/* most gcc compilers know a function __attribute__((unused)) */
#define Q_UNUSED __attribute__((unused))

/* Define that makes XFE memory registration happen at MDBind time, rather
   than at data movement time. */
/* #undef REGISTER_ON_BIND */

/* specifying data alignment is allowed */
#define SANDIA_ALIGNEDDATA_ALLOWED 1

/* if the compiler supports __sync_val_compare_and_swap */
#define SANDIA_ATOMIC_BUILTINS 1

/* if the compiler supports __sync_val_compare_and_swap */
#define SANDIA_BUILTIN_CAS 1

/* if the compiler supports __sync_val_compare_and_swap on 128-bit ints */
#define SANDIA_BUILTIN_CAS128 1

/* if the compiler supports __sync_val_compare_and_swap on 32-bit ints */
#define SANDIA_BUILTIN_CAS32 1

/* if the compiler supports __sync_val_compare_and_swap on 64-bit ints */
#define SANDIA_BUILTIN_CAS64 1

/* if the compiler supports __sync_val_compare_and_swap on pointers */
#define SANDIA_BUILTIN_CAS_PTR 1

/* if the compiler supports __sync_fetch_and_add */
#define SANDIA_BUILTIN_INCR 1

/* if this header is necessary for builtin atomics */
/* #undef SANDIA_NEEDS_INTEL_INTRIN */

/* Defined to 1 if should use the asm pause instruction */
#define SANDIA_WANT_ASM_PAUSE 1

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if strerror_r returns char *. */
#define STRERROR_R_CHAR_P 1

/* if the compiler supports __builtin_unreachable())) */
#define UNREACHABLE 

/* Define to use KNEM */
/* #undef USE_KNEM */

/* Enable extensions on AIX 3, Interix.  */
#ifndef _ALL_SOURCE
# define _ALL_SOURCE 1
#endif
/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif
/* Enable threading extensions on Solaris.  */
#ifndef _POSIX_PTHREAD_SEMANTICS
# define _POSIX_PTHREAD_SEMANTICS 1
#endif
/* Enable extensions on HP NonStop.  */
#ifndef _TANDEM_SOURCE
# define _TANDEM_SOURCE 1
#endif
/* Enable general extensions on Solaris.  */
#ifndef __EXTENSIONS__
# define __EXTENSIONS__ 1
#endif


/* Define to enable PPE */
/* #undef WITH_PPE */

/* Define to enable RUDP */
/* #undef WITH_RUDP */

/* Define to enable IB support */
#define WITH_TRANSPORT_IB 1

/* Define to enable shared memory support */
/* #undef WITH_TRANSPORT_SHMEM */

/* Define to enable UDP support */
/* #undef WITH_TRANSPORT_UDP */

/* Define to enable triggered match list entry operations */
#define WITH_TRIG_ME_OPS 1

/* Define to enable MOFED V2.2 or greater or Qlogic */
/* #undef WITH_ZERO_MRS */

/* Enable large inode numbers on Mac OS X 10.5.  */
#ifndef _DARWIN_USE_64_BIT_INODE
# define _DARWIN_USE_64_BIT_INODE 1
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define to 1 if on MINIX. */
/* #undef _MINIX */

/* Define to 2 if the system does not provide POSIX.1 features except with
   this defined. */
/* #undef _POSIX_1_SOURCE */

/* Define to 1 if you need to in order for `stat' and other things to work. */
/* #undef _POSIX_SOURCE */

/* Define for Solaris 2.5.1 so the uint32_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT32_T */

/* Define for Solaris 2.5.1 so the uint64_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT64_T */

/* Define for Solaris 2.5.1 so the uint8_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT8_T */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `long int' if <sys/types.h> does not define. */
/* #undef off_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#define restrict __restrict
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef ssize_t */

/* Define to the type of an unsigned integer type of width exactly 16 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint16_t */

/* Define to the type of an unsigned integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint32_t */

/* Define to the type of an unsigned integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint64_t */

/* Define to the type of an unsigned integer type of width exactly 8 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint8_t */

/* Define as `fork' if `vfork' does not work. */
/* #undef vfork */

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */
