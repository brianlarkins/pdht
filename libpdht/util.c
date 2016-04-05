#include <pdht_impl.h>

/**
 *  eprintf - error printing wrapper
 *    @return number of bytes written to stderr
 */
int eprintf(const char *format, ...) {
  va_list ap;
  int ret;

  if (c->rank == 0) {
    va_start(ap, format);
    ret = vfprintf(stdout, format, ap);
    va_end(ap);
    return ret;
  }
  else
    return 0;
}



/**
 *  pdht_dbg_printf - optionally compiled debug printer
 *    @return number of bytes written to stderr
 */
int pdht_dbg_printf(const char *format, ...) {
  va_list ap;
  int ret;

  va_start(ap, format);
  fprintf(stdout,"%d: ", c->rank);
  ret = vfprintf(stdout, format, ap);
  va_end(ap);
  fflush(stdout);
  return ret;
}


/**
 * pdht_ptl_error() - returns string for Portals error codes 
 * @param error_code number
 * @returns string for error code
 */
char *pdht_ptl_error(int error_code)
{                                      
    switch (error_code) {
        case PTL_OK:
            return "PTL_OK";

        case PTL_ARG_INVALID:
            return "PTL_ARG_INVALID";

        case PTL_CT_NONE_REACHED:
            return "PTL_CT_NONE_REACHED";

        case PTL_EQ_DROPPED:
            return "PTL_EQ_DROPPED";

        case PTL_EQ_EMPTY:
            return "PTL_EQ_EMPTY";

        case PTL_FAIL:
            return "PTL_FAIL";

        case PTL_IN_USE:
            return "PTL_IN_USE";

        case PTL_IGNORED:
            return "PTL_IGNORED";

        case PTL_INTERRUPTED:
            return "PTL_INTERRUPTED";

        case PTL_LIST_TOO_LONG:
            return "PTL_LIST_TOO_LONG";

        case PTL_NO_INIT:
            return "PTL_NO_INIT";

        case PTL_NO_SPACE:
            return "PTL_NO_SPACE";

        case PTL_PID_IN_USE:
            return "PTL_PID_IN_USE";

        case PTL_PT_FULL:
            return "PTL_PT_FULL";

        case PTL_PT_EQ_NEEDED:
            return "PTL_PT_EQ_NEEDED";

        case PTL_PT_IN_USE:
            return "PTL_PT_IN_USE";

        default:
            return "Unknown Portals return code";
    }
}
