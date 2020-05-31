
#ifndef MYSHMEM_H
#define MYSHMEM_H

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

#include <time.h>
#include <sys/time.h>

double wall_seconds() {
  struct timeval tp;
  int retVal = gettimeofday(&tp,NULL);
  if (retVal == -1) { perror("gettimeofday:"); fflush(stderr); }
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

#define T0_fprintf if(MYTHREAD==0) fprintf

#define myshmem_global_malloc(nblock_on_all_threads,blocksize) shmem_malloc( ( ( (size_t)(nblock_on_all_threads)+(shmem_n_pes())-1)/(shmem_n_pes()) ) *(blocksize) ) \

void myshmem_int64_p(int64_t *addr, size_t index, int64_t val) {
    int64_t pe = index % shmem_n_pes();
    int64_t local_index = index / shmem_n_pes();
    shmem_int64_p(addr+local_index, val, pe);
}

int64_t myshmem_int64_g(int64_t *addr, size_t index) {
    int64_t pe = index % shmem_n_pes();
    int64_t local_index = index / shmem_n_pes();
    return shmem_int64_g(addr+local_index, pe);
}

static void *setup_shmem_reduce_workdata(long **psync, size_t xsize) {
  void *work;
  int i;

  work=shmem_malloc(_SHMEM_REDUCE_MIN_WRKDATA_SIZE*xsize);
  *psync=(long *)shmem_malloc(_SHMEM_REDUCE_SYNC_SIZE*sizeof(long));
  for(i=0;i<_SHMEM_REDUCE_SYNC_SIZE;++i) {
    (*psync)[i]=_SHMEM_SYNC_VALUE;
  }
  shmem_barrier_all();
  return work;
}

double myshmem_reduce_add(double myval){
  static double *buff=NULL, *work;
  static long *sync;
  if (buff==NULL) {
    buff=(double*)shmem_malloc(2*sizeof(double));
    work=(double*)setup_shmem_reduce_workdata(&sync,sizeof(double));
  }
  buff[0]=myval;

  shmem_double_sum_to_all(&buff[1],buff,1,0,0,shmem_n_pes(),work,sync);

  shmem_barrier_all();
  return buff[1];
}
#endif //MYSHMEM_H

