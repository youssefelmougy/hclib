
/******************************************************************
//
//
//  Copyright(C) 2019, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//  This material may be reproduced by or for the US Government
//  pursuant to the copyright license under the clauses at DFARS
//  252.227-7013 and 252.227-7014.
//
//
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the copyright holder nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT HOLDER NOR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
//  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.
//
 *****************************************************************/

/*! \file histo_conveyor.upc
 * \brief A conveyor implementation of histogram.
 */
 extern "C" {
#include <convey.h>
}
#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <stdlib.h>
#include <shmem.h>
#include "mytime.h"

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()


/*!
 * \brief This routine implements histogram using conveyors
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param T the length of the pcindx array
 * \param *lcounts localized pointer to the count array.
 * \return average run time
 *
 */

double histo_conveyor(int64_t *pckindx, int64_t T,  int64_t *lcounts) {
  int ret;
  int64_t i;
  double tm;
  int64_t pe, col;
  int64_t pop_col;

  //convey_t* conveyor = convey_new(SIZE_MAX, 0, NULL, convey_opt_SCATTER);
  convey_t* conveyor = convey_new(SIZE_MAX, 0, NULL, 0);
  if(!conveyor){printf("ERROR: histo_conveyor: convey_new failed!\n"); return(-1.0);}

  ret = convey_begin(conveyor, sizeof(int64_t));
  if(ret < 0){printf("ERROR: histo_conveyor: begin failed!\n"); return(-1.0);}

  tm = wall_seconds();
  i = 0UL;
  while(convey_advance(conveyor, i == T)) {
    for(; i< T; i++){
      col = pckindx[i] >> 16;
      pe  = pckindx[i] & 0xffff;
      if( !convey_push(conveyor, &col, pe))
	    break;
    }
    while( convey_pull(conveyor, &pop_col, NULL) == convey_OK)
    //while( (pop_col = convey_pull(conveyor, NULL)) != NULL)
      lcounts[pop_col] += 1;
  }
  tm = wall_seconds() - tm;
  printf("time %f\n", tm);

  convey_free(conveyor);
  return 0;
}

int main(int argc, char * argv[]) {

  shmem_init();
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  printf("Hostname: %s rank: %d\n", hostname, MYTHREAD);

  int64_t buf_cnt = 1024;
  int64_t models_mask = 0; // run all the programing models
  int64_t l_num_ups  = 10; //1000000;     // per thread number of requests (updates)
  int64_t lnum_counts = 4; //1000;       // per thread size of the table
  int64_t cores_per_node = 0;       // Default to 0 so it won't give misleading bandwidth numbers

  int64_t i;

  int printhelp = 0;
  int opt;
  while( (opt = getopt(argc, argv, "hb:M:n:c:T:")) != -1 ) {
    switch(opt) {
    case 'h': printhelp = 1; break;
    case 'b': sscanf(optarg,"%ld" ,&buf_cnt);  break;
    case 'M': sscanf(optarg,"%ld" ,&models_mask);  break;
    case 'n': sscanf(optarg,"%ld" ,&l_num_ups);  break;
    case 'T': sscanf(optarg,"%ld" ,&lnum_counts);  break;
    case 'c': sscanf(optarg,"%ld" ,&cores_per_node); break;
    default:  break;
    }
  }

  fprintf(stderr,"Running histo on %d threads\n", THREADS);
  fprintf(stderr,"buf_cnt (number of buffer pkgs)      (-b)= %ld\n", buf_cnt);
  fprintf(stderr,"Number updates / thread              (-n)= %ld\n", l_num_ups);
  fprintf(stderr,"Table size / thread                  (-T)= %ld\n", lnum_counts);
  fprintf(stderr,"models_mask                          (-M)= %ld\n", models_mask);
  fflush(stderr);


  // Allocate and zero out the counts array
  int64_t num_counts = lnum_counts*THREADS;
  int64_t * counts = (int64_t *)shmem_malloc(lnum_counts* sizeof(int64_t)); assert(counts != NULL);
  int64_t *lcounts = counts;
  for(i = 0; i < lnum_counts; i++)
    lcounts[i] = 0L;

  // index is a local array of indices into the shared counts array.
  // This is used by the _agi version.
  // To avoid paying the UPC tax of computing index[i]/THREADS and index[i]%THREADS
  // when using the exstack and conveyor models
  // we also store a packed version that holds the pe (= index%THREADS) and lindx (=index/THREADS)
  int64_t *index   = (int64_t *) calloc(l_num_ups, sizeof(int64_t)); assert(index != NULL);
  int64_t *pckindx = (int64_t *) calloc(l_num_ups, sizeof(int64_t)); assert(pckindx != NULL);
  int64_t indx, lindx, pe;

  srand(MYTHREAD + 120348);
  for(i = 0; i < l_num_ups; i++) {
    //indx = i % num_counts;          //might want to do this for debugging
    indx = rand() % num_counts;
    index[i] = indx;
    lindx = indx / THREADS;
    pe  = indx % THREADS;
    pckindx[i]  =  (lindx << 16L) | (pe & 0xffff);
  }
  double volume_per_node = (8*l_num_ups*cores_per_node)*(1.0E-9);

  shmem_barrier_all();

  int64_t use_model;
  double laptime = 0.0;
  double injection_bw = 0.0;
  int64_t num_models = 0L;               // number of models that are executed


      fprintf(stderr,"Conveyors: ");
      laptime = histo_conveyor(pckindx, l_num_ups, lcounts);

    //injection_bw = volume_per_node / laptime;
    //fprintf(stderr,"  %8.3lf seconds  %8.3lf GB/s injection bandwidth\n", laptime, injection_bw);
 
  if(printhelp) {
    for(int i=0;i<lnum_counts; i++)
        printf("%ld ", lcounts[i]);
    printf("\n");
  }

  shmem_barrier_all();

  // Check the results
  // Assume that the atomic add version will correctly zero out the counts array
  //for(i = 0; i < l_num_ups; i++) {
  //  lgp_atomic_add(counts, index[i], -num_models);
  //}
  //shmem_barrier_all();

  //int64_t num_errors = 0, totalerrors = 0;
  //for(i = 0; i < lnum_counts; i++) {
  //  if(lcounts[i] != 0L) {
  //    num_errors++;
  //    if(num_errors < 5)  // print first five errors, report number of errors below
  //      fprintf(stderr,"ERROR: Thread %d error at %ld (= %ld)\n", MYTHREAD, i, lcounts[i]);
  //  }
  //}
  //totalerrors = shmem_reduce_add_l(num_errors);
  //if(totalerrors) {
  //   fprintf(stderr,"FAILED!!!! total errors = %ld\n", totalerrors);
  //}

  shmem_free(counts);
  free(index);
  free(pckindx);
  shmem_finalize();

  return 0;
}
