
#include <stdio.h>
#include <stdlib.h>
#include <shmem.h>
#define NUM_ELEMS 4
//static long pSync[SHMEM_BCAST_SYNC_SIZE];
//static long source[NUM_ELEMS], dest[NUM_ELEMS];

int main(void) {
    int i, me, npes;
    shmem_init();
    me = shmem_my_pe();
    npes = shmem_n_pes();

    long *pReduce = (long*)shmem_malloc(SHMEM_REDUCE_SYNC_SIZE*sizeof(long));
    long *ipWrk = (long*)shmem_malloc(SHMEM_REDUCE_MIN_WRKDATA_SIZE*sizeof(long));
    long *source = (long*)shmem_malloc(NUM_ELEMS*sizeof(long));
    long *dest = (long*)shmem_malloc(NUM_ELEMS*sizeof(long));

    //if (me == 0)
        for (i = 0; i < NUM_ELEMS; i++)
            source[i] = i+22;
    for (i=0; i < SHMEM_REDUCE_SYNC_SIZE; i++) {
        pReduce[i] = SHMEM_SYNC_VALUE;
    }
    shmem_barrier_all(); /* Wait for all PEs to initialize pSync */
    shmem_long_sum_to_all(dest, source, NUM_ELEMS, 0, 0, npes, ipWrk, pReduce);
    printf("%d: %ld", me, dest[0]);
    for (i = 1; i < NUM_ELEMS; i++)
        printf(", %ld", dest[i]);
    printf("\n");
    return 0;
}

