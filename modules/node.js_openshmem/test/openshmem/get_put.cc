
#include<stdio.h>
#include<shmem.h>

int main() {
  //printf("size %d %d %d %d %d\n", sizeof(int), sizeof(long), sizeof(long long), sizeof(float), sizeof(double));
  shmem_init();
  int rank =  shmem_my_pe();
  printf("rank %d\n", rank);
  int size = sizeof(long);
  long *ptr = (long*)shmem_malloc(size);

  if(rank == 0) {
    shmem_long_p(ptr, 34, 1);
    shmem_barrier_all();
    long val = shmem_long_g(ptr, 0);
  }
  else if(rank == 1) {
    shmem_barrier_all();
    printf("value for on 1 from 0 %d\n", *ptr);
  }
  else {
    shmem_barrier_all();
  }
  shmem_finalize();
}
