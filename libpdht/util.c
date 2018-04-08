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



#ifdef DEPRECATED
/**
 *  pdht_get_wtime - get a wall clock time for performance analysis
 */
double pdht_get_wtime() {
  double t;
  struct timeval tv;

  gettimeofday(&tv, NULL);

  t = (tv.tv_sec*1000000LL + tv.tv_usec)/1000000.0;

  return t;
}
#endif // moved to inline code, uses clock_gettime()


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
 *  pdht_lvl_dbg_printf - optionally compiled debug printer with priority level
 *    @param lvl priority level
 *    @param format format string
 *    @return number of bytes written to stderr
 */
int pdht_lvl_dbg_printf(int lvl, const char *format, ...) {
  va_list ap;
  int ret;

  if (lvl <= c->dbglvl) {
    va_start(ap, format);
    fprintf(stdout,"%d: ", c->rank);
    ret = vfprintf(stdout, format, ap);
    va_end(ap);
    fflush(stdout);
  }
  return ret;
}



/**
 *  pdht_lvl_dbg_printf - optionally compiled debug printer with priority level
 *    @param lvl priority level
 *    @param format format string
 *    @return number of bytes written to stderr
 */
int pdht_lvl_dbg_eprintf(int lvl, const char *format, ...) {
  va_list ap;
  int ret;

  if ((c->rank == 0) && (lvl <= c->dbglvl)) {
    va_start(ap, format);
    ret = vfprintf(stdout, format, ap);
    va_end(ap);
    fflush(stdout);
  }
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



/**
 * pdht_print_active
 */
void pdht_print_active(pdht_t *dht, void kprinter(void *key), void vprinter(void *val)) {
  char *iter;
  int pending = 0;
  _pdht_ht_trigentry_t *hte;
  long *key;

  iter = (char *)dht->ht;

  for (int i=0; i < dht->maxentries; i++) {
    hte = (_pdht_ht_trigentry_t *)iter;
    key = (long *)hte->key;
    if (hte->ame != PTL_INVALID_HANDLE) {
      pdht_dprintf("elem %d: mbits: %12"PRIx64" ptr: %p ", i, hte->me.match_bits, &hte->key);
      //pdht_dprintf(" pkey: %ld ", key[0]);
      //pdht_dprintf(" pkey: [%ld,%ld,%ld@%ld] ", key[0],key[1],key[2],key[3]); // MADNESS
      kprinter(hte->key);
      vprinter(&hte->data);
      printf("\n");
    } else {
      pending++;
    }
    iter += dht->entrysize;
  }
  pdht_dprintf("pending: %d\n", pending);
}



/**
 * pdht_print_stats - prints out runtime statistics
 */
void pdht_print_stats(pdht_t *dht) {
  u_int64_t ilocal[6];
  u_int64_t tilocal[6];
  u_int64_t isum[6];
  u_int64_t imin[6];
  u_int64_t imax[6];
  double    dlocal[8];
  double    dsum[8];
  double    dmin[8];
  double    dmax[8];
  double    tdlocal[8];
 

  ilocal[0] = dht->stats.puts;
  ilocal[1] = dht->stats.gets;
  ilocal[2] = dht->stats.collisions;
  ilocal[3] = dht->stats.notfound;
  ilocal[4] = dht->stats.updates;
  ilocal[5] = dht->stats.inserts;
 
  memcpy(tilocal,ilocal,sizeof(ilocal));

  dlocal[0] = PDHT_READ_TIMER(dht, ptimer);
  dlocal[1] = PDHT_READ_TIMER(dht, gtimer);
  dlocal[2] = PDHT_READ_TIMER(dht, t1);
  dlocal[3] = PDHT_READ_TIMER(dht, t2);
  dlocal[4] = PDHT_READ_TIMER(dht, t3);
  dlocal[5] = PDHT_READ_TIMER(dht, t4);
  dlocal[6] = PDHT_READ_TIMER(dht, t5);
  dlocal[7] = PDHT_READ_TIMER(dht, t6);
  
  memcpy(tdlocal,dlocal,sizeof(dlocal)); //have to set temp locals because they are manipulated for pdht_allreduce
  
  pdht_allreduce(tilocal, isum, PdhtReduceOpSum, LongType, 6);
  memcpy(tilocal,ilocal,sizeof(ilocal));
  pdht_allreduce(tilocal, imin, PdhtReduceOpMin, LongType, 6);
  memcpy(tilocal,ilocal,sizeof(ilocal));
  pdht_allreduce(tilocal, imax, PdhtReduceOpMax, LongType, 6);

  pdht_allreduce(tdlocal, dsum, PdhtReduceOpSum, DoubleType, 8);
  memcpy(tdlocal,dlocal,sizeof(dlocal));
  pdht_allreduce(tdlocal, dmin, PdhtReduceOpMin, DoubleType, 8);
  memcpy(tdlocal,dlocal,sizeof(dlocal));
  pdht_allreduce(tdlocal, dmax, PdhtReduceOpMax, DoubleType, 8);

  if (c->rank == 0) {
    printf("pdht global stats: \n");    

    printf("\tputs:       min: %12"PRIu64"\tmax: %12"PRIu64"\t total: %12"PRIu64"\n", imin[0], imax[0], isum[0]);
    printf("\tupdates:    min: %12"PRIu64"\tmax: %12"PRIu64"\t total: %12"PRIu64"\n", imin[4], imax[4], isum[4]);
    printf("\tinserts:    min: %12"PRIu64"\tmax: %12"PRIu64"\t total: %12"PRIu64"\n", imin[5], imax[5], isum[5]);
    printf("\tgets:       min: %12"PRIu64"\tmax: %12"PRIu64"\t total: %12"PRIu64"\n", imin[1], imax[1], isum[1]);
    printf("\tcollisions: min: %12"PRIu64"\tmax: %12"PRIu64"\t total: %12"PRIu64"\n", imin[2], imax[2], isum[2]);
    printf("\tnotfound:   min: %12"PRIu64"\tmax: %12"PRIu64"\t total: %12"PRIu64"\n", imin[3], imax[3], isum[3]);
    printf("\tputtime:    min: %10.4f sec\t max:%10.4f sec avg: %10.4f\n", 
                  dmin[0]/(double)1e9, dmax[0]/(double)1e9, dsum[0]/(double)(c->size * 1e9));
    printf("\tgettime:    min: %10.4f sec\t max:%10.4f sec avg: %10.4f\n", 
                  dmin[1]/(double)1e9, dmax[1]/(double)1e9, dsum[1]/(double)(c->size * 1e9));
    printf("\tt1:    min: %10.4f sec\t max:%10.4f sec avg: %10.4f\n", 
                  dmin[2]/(double)1e9, dmax[2]/(double)1e9, dsum[2]/(double)(c->size * 1e9));
    printf("\tt2:    min: %10.4f sec\t max:%10.4f sec avg: %10.4f\n", 
                  dmin[3]/(double)1e9, dmax[3]/(double)1e9, dsum[3]/(double)(c->size * 1e9));
    printf("\tt3:    min: %10.4f sec\t max:%10.4f sec avg: %10.4f\n", 
                  dmin[4]/(double)1e9, dmax[4]/(double)1e9, dsum[4]/(double)(c->size * 1e9));
    printf("\tt4:    min: %10.4f sec\t max:%10.4f sec avg: %10.4f\n", 
                  dmin[5]/(double)1e9, dmax[5]/(double)1e9, dsum[5]/(double)(c->size * 1e9));
    printf("\tt5:    min: %10.4f sec\t max:%10.4f sec avg: %10.4f\n", 
                  dmin[6]/(double)1e9, dmax[6]/(double)1e9, dsum[6]/(double)(c->size * 1e9));
    printf("\tt6:    min: %10.4f sec\t max:%10.4f sec avg: %10.4f\n", 
                  dmin[7]/(double)1e9, dmax[7]/(double)1e9, dsum[7]/(double)(c->size * 1e9));
  }
}



/**
 * pdht_print_distribution - prints the distribution of puts over the processor ranks
 * @param dht - DHT of the hash table of interest
 */
void pdht_print_distribution(pdht_t *dht) {
  u_int64_t iglobal[PDHT_MAX_RANKS];

  // hash distribution data
  pdht_allreduce(dht->stats.rankputs, iglobal, PdhtReduceOpSum, LongType, PDHT_MAX_RANKS);
  
  if (c->rank == 0) {
    printf("  put distribution: \n");
    for (int i=0; i < c->size; i++)
      printf("    rank[%d] : %12"PRIu64"\n", i, iglobal[i]);
  }
  pdht_barrier();
}
