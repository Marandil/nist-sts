#include <stdlib.h>
#include <stdint.h>
#include <sys/mman.h>

#include <errno.h>
#include <stdio.h>
#include <string.h>

typedef struct {
  uintptr_t ptr;
  size_t size;
} _map_entry_t;

#define ALLOC_BUFFER_SIZE 16384
_map_entry_t alloc_buffer[ALLOC_BUFFER_SIZE] = {0};
size_t alloc_size = 0;

void *mcalloc(size_t nmemb, size_t size) {
  if ( ALLOC_BUFFER_SIZE <= alloc_size ) return NULL;

  size_t total_size = nmemb * size;
  void* addr = mmap(NULL, total_size, PROT_READ | PROT_WRITE,
                    MAP_ANONYMOUS | MAP_PRIVATE, -1, 0);
  printf("MMAP %zd, returned %p, errno %s(%x).\n", total_size, addr, strerror(errno), errno);
  if ( MAP_FAILED == addr )
    return NULL;

  /* record length in the alloc_buffer */
  alloc_buffer[alloc_size++] = (_map_entry_t) { .ptr = (uintptr_t)addr, .size = size } ;

  return addr;
}

void mfree(void* addr) {
  int    retval;
  size_t size, idx;

  /* Find the right index in the table */
  for ( idx = alloc_size; idx > 0; idx--)
  {
    if (alloc_buffer[idx-1].ptr == (uintptr_t)addr) {
      break;
    }
  }
  /* if idx reached 0, no suitable buffer was found. */
  if ( 0 == idx ) {
    fprintf(stderr, "Potential double-mfree detected: Couldn't deallocate slot %p\n", addr);
    return;
  }
  idx--;
  size = alloc_buffer[idx].size;
  /* swap data at idx (deprecated) with the last used cell (alloc_size-1) */
  alloc_size--;
  alloc_buffer[idx] = alloc_buffer[alloc_size];
  alloc_buffer[alloc_size] = (_map_entry_t) {0, 0};

  retval = munmap(addr, size);
  printf("MUNMAP %p, returned %d, errno %s(%x).\n", addr, retval, strerror(errno), errno);
}
