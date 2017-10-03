/********************************************************/
/*                                                      */
/*    alloc.c                                           */
/*                                                      */
/*    author: d. brian larkins                          */
/*    id: $Id: alloc.c 1628 2009-03-05 19:45:08Z brian $m */
/*                                                      */
/********************************************************/


#include <stdio.h>
#include <stdlib.h>

int memalloc = 0;

/**
 *  talloc - custom memory allocator (zero-initialized)
 *    @param  size size of block to allocate
 *    @return pointer to new memory allocation
 *    @return NULL if error
 */
void *talloc(size_t size) {
  void *ret;
  memalloc += size;
  ret =  malloc(size);
  if (!ret) 
    printf("out of memory. good-bye.\n");
  return ret;
}



/**
 *  talloc - custom memory allocator (zero-initialized)
 *    @param  size size of block to allocate
 *    @return pointer to new memory allocation
 *    @return NULL if error
 */
void *tcalloc(size_t size) {
  memalloc += size;
  return calloc(1, size);
}



/**
 *  tfree - custom deallocator
 *    @return pointer of memory to deallocate
 */
void tfree(void *p) {
  free(p);
}
