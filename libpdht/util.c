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
    fflush(stdout);
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



/**
 * pdht_event_to_string - returns string for a Portals event type
 * @param evtype Portals4 event type
 * @returns character string for event
 */
char *pdht_event_to_string(ptl_event_kind_t evtype) {
  char *ret = "(unmatched)";
  switch (evtype) {
    case PTL_EVENT_GET:
      ret = "PTL_EVENT_GET";
      break;
    case PTL_EVENT_GET_OVERFLOW:
      ret = "PTL_EVENT_GET_OVERFLOW";
      break;
    case PTL_EVENT_PUT:
      ret = "PTL_EVENT_PUT";
      break;
    case PTL_EVENT_PUT_OVERFLOW:
      ret = "PTL_EVENT_PUT_OVERFLOW";
      break;
    case PTL_EVENT_ATOMIC:
      ret = "PTL_EVENT_ATOMIC";
      break;
    case PTL_EVENT_ATOMIC_OVERFLOW:
      ret = "PTL_EVENT_ATOMIC_OVERFLOW";
      break;
    case PTL_EVENT_FETCH_ATOMIC:
      ret = "PTL_EVENT_FETCH_ATOMIC";
      break;
    case PTL_EVENT_FETCH_ATOMIC_OVERFLOW:
      ret = "PTL_EVENT_FETCH_ATOMIC_OVERFLOW";
      break;
    case PTL_EVENT_REPLY:
      ret = "PTL_EVENT_REPLY";
      break;
    case PTL_EVENT_SEND:
      ret = "PTL_EVENT_SEND";
      break;
    case PTL_EVENT_ACK:
      ret = "PTL_EVENT_ACK";
      break;
    case PTL_EVENT_PT_DISABLED:
      ret = "PTL_EVENT_PT_DISABLED";
      break;
    case PTL_EVENT_LINK:
      ret = "PTL_EVENT_LINK";
      break;
    case PTL_EVENT_AUTO_UNLINK:
      ret = "PTL_EVENT_AUTO_UNLINK";
      break;
    case PTL_EVENT_AUTO_FREE:
      ret = "PTL_EVENT_AUTO_FREE";
      break;
    case PTL_EVENT_SEARCH:
      ret = "PTL_EVENT_SEARCH";
      break;
  }
  return ret;
}


/**
 * pdht_dump_event - prints out entire event
 *  @param ev Portals event
 */
void pdht_dump_event(ptl_event_t *ev) {
   char *fail = "(unset)";
   pdht_dprintf("event type: %s\n", pdht_event_to_string(ev->type));
   pdht_dprintf("\tstart: %p\n", ev->start);
   pdht_dprintf("\tuser_ptr: %p\n", ev->user_ptr);
   pdht_dprintf("\tmatch: %lu\n", ev->match_bits);
   pdht_dprintf("\trlength: %lu\n", ev->rlength);
   pdht_dprintf("\tmlength: %lu\n", ev->mlength);
   pdht_dprintf("\tremote: %lu\n", ev->remote_offset);
   pdht_dprintf("\tuid: %lu\n", ev->uid);
   pdht_dprintf("\tinitator: %lu\n", ev->initiator.rank);
   pdht_dprintf("\tpt_index: %lu\n", ev->pt_index);
   switch (ev->ni_fail_type) {
   case PTL_NI_OK:
     fail = "PTL_NI_OK";
     break;
   case PTL_NI_UNDELIVERABLE:
     fail = "PTL_NI_UNDELIVERABLE";
     break;
   case PTL_NI_PT_DISABLED:
     fail = "PTL_NI_PT_DISABLED";
     break;
   case PTL_NI_DROPPED:
     fail = "PTL_NI_DROPPED";
     break;
   case PTL_NI_PERM_VIOLATION:
     fail = "PTL_NI_PERM_VIOLATION";
     break;
   case PTL_NI_OP_VIOLATION:
     fail = "PTL_NI_OP_VIOLATION";
     break;
   case PTL_NI_SEGV:
     fail = "PTL_NI_SEGV";
     break;
   case PTL_NI_NO_MATCH:
     fail = "PTL_NI_NO_MATCH";
     break;
   }
   pdht_dprintf("\tni_fail: %s\n", fail);
   
}
