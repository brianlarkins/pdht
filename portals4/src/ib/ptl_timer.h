#ifndef PTL_INTERNAL_TIMER_H
#define PTL_INTERNAL_TIMER_H

#ifdef HAVE_MACH_TIMER
# include <mach/mach_time.h>
# define TIMER_TYPE uint64_t
# define MARK_TIMER(x)          do { x = mach_absolute_time(); } while (0)
# define TIMER_INTS(x)          (x)
static inline uint64_t MILLI_TO_TIMER_INTS(uint64_t x)
{
    return x * 1000000;
}

#elif defined(HAVE_GETTIME_TIMER)
# ifndef _POSIX_C_SOURCE
#  define _POSIX_C_SOURCE 199309L
# endif
# include <time.h>
# define TIMER_TYPE struct timespec
# define MARK_TIMER(x)          clock_gettime(CLOCK_MONOTONIC, &x)
# define TIMER_INTS(x)          (x.tv_sec * 1000000000 + x.tv_nsec)
static inline uint64_t MILLI_TO_TIMER_INTS(uint64_t x)
{
    return x * 1000000;
}
#else /* ifdef HAVE_MACH_TIMER */
# error No timers are defined!
#endif /* ifdef HAVE_MACH_TIMER */

#endif /* ifndef PTL_INTERNAL_TIMER_H */
/* vim:set expandtab: */
