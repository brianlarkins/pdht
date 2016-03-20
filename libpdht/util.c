#include "cfg.h"

/**
 *  eprintf - error printing wrapper
 *    @return number of bytes written to stderr
 */
int eprintf(const char *format, ...) {
  va_list ap;
  int ret;

  if (c.rank == 0) {
    va_start(ap, format);
    ret = vfprintf(stdout, format, ap);
    va_end(ap);
    return ret;
  }
  else
    return 0;
}
