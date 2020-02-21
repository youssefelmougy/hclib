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

/*! \file toposort_conveyor.upc
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 */

#include <iostream>
#include <math.h>
#include <shmem.h>

extern "C" {

#include <exstack.h>
#include <spmat.h>
#include <libgetput.h>
#include <convey.h>

}

#include "mytime.h"

typedef struct pkg_tri_t {
    int64_t w;    
    int64_t vj;
} pkg_tri_t;

/*!
 * \brief pop routine to handle the pushes
 * \param *c a place to return the number of hits.
 * \param *conv the conveyor
 * \param *mat the input sparse matrix 
     NB. The nonzero within a row must be increasing
 * \param done the signal to convey_advance that this thread is done
 * \return the return value from convey_advance
 */
static int64_t tri_convey_push_process(int64_t* c, convey_t* conv, sparsemat_t* mat, int64_t done) {
    int64_t k, cnt = 0;
    struct pkg_tri_t pkg;
    
    while (convey_pull(conv, &pkg, NULL) == convey_OK) {
    //Is pkg.w on row pkg.vj
        for (k = mat->loffset[pkg.vj]; k < mat->loffset[pkg.vj + 1]; k++) {
            if (pkg.w == mat->lnonzero[k]) {
                cnt ++;
                break;
            }
            if ( pkg.w < mat->lnonzero[k]) // requires that nonzeros are increasing
            break;
        }
    }
    *c += cnt;

    return convey_advance(conv, done);
}

/*!
* \brief This routine implements the exstack variant of triangle counting. 
 *   NB. The column indices of the nonzeros must be in increasing order within each row.
 * \param *count  a place to return the counts from each thread
 * \param *sr a place to return the number of shared references
 * \param *L lower triangle of the input matrix
 * \param *U upper triangle of the input matrix
 * \param alg 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
 * \return average run time
 */
double triangle_convey_push(int64_t* count, int64_t* sr, sparsemat_t* L, sparsemat_t* U, int64_t alg) {
    convey_t * conv = convey_new(SIZE_MAX, 0, NULL, 0);
    if (conv == NULL) return(-1);
    if (convey_begin(conv, sizeof(pkg_tri_t)) != convey_OK) return(-1);

    int64_t cnt = 0;
    int64_t numpushed = 0;
    double t1 = wall_seconds();

    pkg_tri_t pkg;
    int64_t k,kk, pe;
    int64_t l_i, L_i, L_j;

    if (alg == 0) {
        // foreach nonzero (i,j) in L
        for (l_i=0; l_i < L->lnumrows; l_i++) { 
            for (k=L->loffset[l_i]; k< L->loffset[l_i + 1]; k++) {
                L_i = l_i * THREADS + MYTHREAD;
                L_j = L->lnonzero[k];
        
                pe = L_j % THREADS;
                pkg.vj = L_j / THREADS;
                for (kk = L->loffset[l_i]; kk < L->loffset[l_i + 1]; kk++) {
                    pkg.w = L->lnonzero[kk]; 
                    if( pkg.w > L_j) 
                        break;

                    numpushed++;
                    if (convey_push(conv, &pkg, pe) != convey_OK) {
                        tri_convey_push_process(&cnt, conv, L, 0); 
                        kk--;
                        numpushed--;
                    }
                }
            }
        }
        while (tri_convey_push_process(&cnt, conv, L, 1)) {}// keep popping til all threads are done
    } else {
        for (l_i=0; l_i < L->lnumrows; l_i++) { 
            for (k=L->loffset[l_i]; k< L->loffset[l_i + 1]; k++) {
                L_i = l_i * THREADS + MYTHREAD;
                L_j = L->lnonzero[k];
        
                pe = L_j % THREADS;
                pkg.vj = L_j / THREADS;
                
                for (kk = U->loffset[l_i]; kk < U->loffset[l_i + 1]; kk++) {
                    pkg.w = U->lnonzero[kk]; 
                    numpushed++;
          
                    if (convey_push(conv, &pkg, pe) != convey_OK) {
                        tri_convey_push_process(&cnt, conv, U, 0); 
                        kk--;
                        numpushed--;
                    }
                }
            }
        }
        while (tri_convey_push_process(&cnt, conv, U, 1)) {} // keep popping til all threads are done
    }

    lgp_barrier();
    *sr = numpushed;
    *count = cnt;
    minavgmaxD_t stat[1];
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d(stat, t1, THREADS);
  
    return(stat->avg);
}

/* the pull version of triangle counting */
double triangle_convey_pull(int64_t *count, int64_t *sr, sparsemat_t *mat) {

    int64_t cnt = 0, row, col, i, ii, j, k, rowcnt, pos, pos2;
    int64_t L_i, l_i, L_j, pe, frompe;

    // determine the largest nnz in a row
    int64_t max_row = 0;
    for (i = 0; i < mat->lnumrows; i++)
        max_row = ((max_row < mat->loffset[i+1] - mat->loffset[i]) ? mat->loffset[i+1] - mat->loffset[i]: max_row);
    
    max_row = lgp_reduce_max_l(max_row);

    struct pkg_tri_t pkg;
    int64_t* buf = (int64_t*)calloc(max_row + 2, sizeof(int64_t));
  
    convey_t* conv_req = convey_new(SIZE_MAX, 0, NULL, 0);
    convey_t* conv_resp = convey_new_elastic(sizeof(int64_t), (max_row + 2) * sizeof(int64_t), SIZE_MAX, 0, NULL, 0);
  
    if (!conv_req || !conv_resp) return(-1);
    if (convey_begin(conv_req, sizeof(pkg_tri_t)) != convey_OK) return -1;
    if (convey_begin(conv_resp, sizeof(int64_t)) != convey_OK) return -1;
  
    int64_t numpushed = 0;
    bool more, new_pull = 1;
    double t1 = wall_seconds();
    convey_item_t* convitem = (convey_item_t*)calloc(1, sizeof(convey_item_t));
    int64_t lastrow = -1, start;
  
    lgp_barrier();
  
    k = 0;
    l_i = 0;
    L_i = MYTHREAD;
  
    while (more = convey_advance(conv_req, k == mat->lnnz), more | convey_advance(conv_resp, !more)) {

        /* push requests: each request is a (i,j) tuple where i is the row 
        this PE will need to look at once the response is sent back. 
        j is the row this PE is asking for. The response will be the 
        nonzeros in local row j on the PE the request is sent to */
    
        while (k < mat->loffset[l_i + 1]) {
            L_j = mat->lnonzero[k]; // shared name for col j
            assert( L_i > L_j );
            pkg.w = l_i;
            pkg.vj = L_j/THREADS;
            pe = L_j % THREADS;
            if (convey_push(conv_req, &pkg, pe) != convey_OK)
                break;      
            numpushed++;
            k++;
        }

        /* we just hit the end of a row */
        if (k == mat->loffset[l_i + 1]) {
            l_i++;
            L_i += THREADS;
        }

        /* pop requests and send responses */
        while(1) {
            if (new_pull) {
                if (convey_pull(conv_req, &pkg, &frompe) != convey_OK)
                    break;
                pos = 0;
                row = pkg.vj;
                buf[pos++] = pkg.w;
                buf[pos++] = mat->loffset[row + 1] - mat->loffset[row];
                for (j = mat->loffset[row]; j < mat->loffset[row+1]; j++) {
                    buf[pos++] = mat->lnonzero[j];
                }
            }

            // this push can still fail, so we have to save some state, or unpop.
            // how much code messiness is the elastic conveyor really saving us?
            if (convey_epush(conv_resp, pos * sizeof(int64_t), buf, frompe) != convey_OK) {
                new_pull = 0;
                break;
            }
            new_pull = 1;
            numpushed+=pos;
        }

        /* pop responses and search for collisions */
        // this depends on the row nonzeros being sorted
        while (convey_epull(conv_resp, convitem) == convey_OK) {
            int64_t * data = (int64_t *)convitem->data;
            pos2 = 0;
            row = data[pos2++];
            rowcnt = data[pos2++];
            start = mat->loffset[row];
            for (i = 0; i < rowcnt; i++, pos2++) {
                for (ii = start; ii < mat->loffset[row+1]; ii++) {
                    if(mat->lnonzero[ii] == data[pos2]) {
                        cnt++;
                        start = ii + 1;
                        break;
                    } else if(mat->lnonzero[ii] > data[pos2]) {
                        start = ii; // since pulled row is sorted, we can change start
                        break;
                    }
                }
            }
        }
    }
  
  
    *sr = numpushed;
    *count = cnt;
    minavgmaxD_t stat[1];
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d(stat, t1, THREADS);

    free(buf);
    free(convitem);
    convey_free(conv_req);
    convey_free(conv_resp);
  
    return stat->avg;
}

sparsemat_t* generate_kronecker_graph(
    int64_t* B_spec,
    int64_t B_num,
    int64_t* C_spec,
    int64_t C_num,
    int mode)
{

    T0_fprintf(stderr, "Generating Mode %d Kronecker Product graph (A = B X C) with parameters:  ", mode);
        for(int i = 0; i < B_num; i++) T0_fprintf(stderr, "%ld ", B_spec[i]);
    T0_fprintf(stderr, "X ");
        for(int i = 0; i < C_num; i++) T0_fprintf(stderr, "%ld ", C_spec[i]);   
    T0_fprintf(stderr, "\n");

    sparsemat_t* B = gen_local_mat_from_stars(B_num, B_spec, mode);
    sparsemat_t* C = gen_local_mat_from_stars(C_num, C_spec, mode);   
    if(!B || !C) {
        T0_fprintf(stderr,"ERROR: triangles: error generating input!\n"); lgp_global_exit(1);
    }

    T0_fprintf(stderr, "B has %ld rows/cols and %ld nnz\n", B->numrows, B->lnnz);
    T0_fprintf(stderr, "C has %ld rows/cols and %ld nnz\n", C->numrows, C->lnnz);

    sparsemat_t* A = kron_prod_dist(B, C, 1);
  
    return A;
}

int main (int argc, char* argv[]) {

    lgp_init(argc, argv);

    int64_t buf_cnt = 1024;
    int64_t models_mask = ALL_Models;  // default is running all models
    int64_t l_numrows = 10000;         // number of a rows per thread
    int64_t nz_per_row = 35;           // target number of nonzeros per row (only for Erdos-Renyi)
    int64_t read_graph = 0L;           // read graph from a file
    char filename[64];
    int64_t cores_per_node = 0;
  
    double t1;
    int64_t i, j;
    int64_t alg = 0;
    int64_t gen_kron_graph = 0L;
    int kron_graph_mode = 0;
    char * kron_graph_string;
    double erdos_renyi_prob = 0.0;
  
    int printhelp = 0;
    int opt;
    while ((opt = getopt(argc, argv, "hb:c:M:n:f:a:e:K:")) != -1) {
        switch (opt) {
            case 'h': printhelp = 1; break;
            case 'b': sscanf(optarg,"%ld", &buf_cnt);  break;
            case 'c': sscanf(optarg,"%ld" ,&cores_per_node); break;
            case 'M': sscanf(optarg,"%ld", &models_mask);  break;
            case 'n': sscanf(optarg,"%ld", &l_numrows); break;
            case 'f': read_graph = 1; sscanf(optarg,"%s", filename); break;

            case 'a': sscanf(optarg,"%ld", &alg); break;
            case 'e': sscanf(optarg,"%lg", &erdos_renyi_prob); break;
            case 'K': gen_kron_graph = 1; kron_graph_string = optarg; break;
            default:  break;
        }
    }

    int64_t numrows = l_numrows * THREADS;
    if (erdos_renyi_prob == 0.0) { // use nz_per_row to get erdos_renyi_prob
        erdos_renyi_prob = (2.0 * (nz_per_row - 1)) / numrows;
        if (erdos_renyi_prob > 1.0) erdos_renyi_prob = 1.0;
    } else {                     // use erdos_renyi_prob to get nz_per_row
        nz_per_row = erdos_renyi_prob * numrows;
    }
  
    T0_fprintf(stderr,"Running triangle on %d threads\n", THREADS);
    if (!read_graph && !gen_kron_graph) {
        T0_fprintf(stderr,"Number of rows per thread   (-N)   %ld\n", l_numrows);
        T0_fprintf(stderr,"Erdos Renyi prob (-e)   %g\n", erdos_renyi_prob);
    }
    T0_fprintf(stderr,"Model mask (M) = %ld (should be 1,2,4,8,16 for agi, exstack, exstack2, conveyors, alternates\n", models_mask);  
    T0_fprintf(stderr,"algorithm (a) = %ld (0 for L & L*U, 1 for L & U*L)\n", alg);
  
    if (printhelp)
        lgp_global_exit(0);

    // alg = 0 only needs L
    // alg = 1 needs both U and L

    double correct_answer = -1;
  
    sparsemat_t *A, *L, *U;
    if (read_graph) {
        A = read_matrix_mm_to_dist(filename);
        if (!A) lgp_global_exit(1);

        T0_fprintf(stderr,"Reading file %s...\n", filename);
        T0_fprintf(stderr, "A has %ld rows/cols and %ld nonzeros.\n", A->numrows, A->nnz);

        // we should check that A is symmetric!
    
        if (!is_lower_triangular(A, 0)) { //if A is not lower triangular... make it so.      
            T0_fprintf(stderr, "Assuming symmetric matrix... using lower-triangular portion...\n");
            tril(A, -1);
            L = A;
        } else {
            L = A;
        }
    
        sort_nonzeros(L);

    }  else if (gen_kron_graph) {
        // string should be <mode> # # ... #
        // we will break the string of numbers (#s) into two groups and create
        // two local kronecker graphs out of them.
        int num;
        char* ptr = kron_graph_string;
        int64_t* kron_specs = (int64_t*)calloc(32, sizeof(int64_t*));
    
        // read the mode
        int ret = sscanf(ptr, "%d ", &kron_graph_mode);
        if (ret == 0) ret = sscanf(ptr, "\"%d ", &kron_graph_mode);
        if (ret == 0) { T0_fprintf(stderr, "ERROR reading kron graph string!\n"); return(1); }
        T0_fprintf(stderr,"kron string: %s return = %d\n", ptr, ret);
        T0_fprintf(stderr,"kron mode: %d\n", kron_graph_mode);
        ptr += 2;
        int mat, num_ints = 0;
        while (sscanf(ptr, "%d %n", &num, &mat) == 1) {
            T0_fprintf(stderr,"%s %d\n", ptr, mat);
            kron_specs[num_ints++] = num;
            ptr+=mat;
        }

        if (num_ints <= 1) {
            T0_fprintf(stderr,"ERROR: invalid kronecker product string (%s): must contain at least three integers\n", kron_graph_string); return(-1);
        }

        /* calculate the number of triangles */
        if (kron_graph_mode == 0) {
            correct_answer = 0.0;
        } else if (kron_graph_mode == 1) {
            correct_answer = 1;
            for (i = 0; i < num_ints; i++)
                correct_answer *= (3 * kron_specs[i] + 1);
            correct_answer *= 1.0 / 6.0;
            double x = 1;
            for (i = 0; i < num_ints; i++)
                x *= (kron_specs[i] + 1);
            correct_answer = correct_answer - 0.5 * x + 1.0 / 3.0;
        } else if (kron_graph_mode == 2) {
            correct_answer = (1.0 / 6.0) * pow(4, num_ints) - pow(2.0, (num_ints - 1)) + 1.0 / 3.0;
        }

        correct_answer = round(correct_answer);
        T0_fprintf(stderr, "Pre-calculated answer = %ld\n", (int64_t)correct_answer);

        int64_t half = num_ints / 2;

        L = generate_kronecker_graph(kron_specs, half, &kron_specs[half], num_ints - half, kron_graph_mode);
    } else {
        L = gen_erdos_renyi_graph_dist(numrows, erdos_renyi_prob, 0, 1, 12345);
    }

    lgp_barrier();

    if(alg == 1)
        U = transpose_matrix(L);

    lgp_barrier();

    T0_fprintf(stderr,"L has %ld rows/cols and %ld nonzeros.\n", L->numrows, L->nnz);

    if (!is_lower_triangular(L, 0)) {
        T0_fprintf(stderr,"ERROR: L is not lower triangular!\n");
        lgp_global_exit(1);
    }

    T0_fprintf(stderr,"Run triangle counting ...\n");
    int64_t tri_cnt;           // partial count of triangles on this thread
    int64_t total_tri_cnt;     // the total number of triangles on all threads
    int64_t sh_refs;         // number of shared reference or pushes
    int64_t total_sh_refs;

  
    SHARED int64_t * cc = (int64_t*)lgp_all_alloc(L->numrows, sizeof(int64_t));
    int64_t * l_cc = lgp_local_part(int64_t, cc);
    for (i = 0; i < L->lnumrows; i++)
        l_cc[i] = 0;
  
    lgp_barrier();
  
    /* calculate col sums */
    for (i = 0; i < L->lnnz; i++) {
        lgp_fetch_and_inc(cc, L->lnonzero[i]);
    }
  
    lgp_barrier();
  
    int64_t rtimesc_calc = 0;
    for (i = 0; i < L->lnumrows; i++) {
        int64_t deg = L->loffset[i + 1] - L->loffset[i];        
        rtimesc_calc += deg * l_cc[i];
    }

    /* calculate sum (r_i choose 2) */
    int64_t rchoose2_calc = 0;
    for (i = 0; i < L->lnumrows; i++) {
        int64_t deg = L->loffset[i + 1] - L->loffset[i];
        rchoose2_calc += deg * (deg - 1) / 2;
    }
  
    /* calculate sum (c_i choose 2) */
    int64_t cchoose2_calc = 0;
    for (i = 0; i < L->lnumrows; i++) {
        int64_t deg = l_cc[i];
        cchoose2_calc += deg * (deg - 1) / 2;
    }

    int64_t pulls_calc = 0;
    int64_t pushes_calc = 0;
    if (alg == 0) {
        pulls_calc = lgp_reduce_add_l(rtimesc_calc);
        pushes_calc = lgp_reduce_add_l(rchoose2_calc);
    } else {
        pushes_calc = lgp_reduce_add_l(rtimesc_calc);
        pulls_calc = lgp_reduce_add_l(cchoose2_calc);
    }

    lgp_all_free(cc);
  
    T0_fprintf(stderr,"Calculated: Pulls = %ld\n            Pushes = %ld\n\n", pulls_calc, pushes_calc);
  
    int64_t use_model;
    double laptime = 0.0;

    tri_cnt = 0;
    total_tri_cnt = 0;
    sh_refs = 0;
    total_sh_refs = 0;
    T0_fprintf(stderr, "      Conveyor: ");
    laptime = triangle_convey_push(&tri_cnt, &sh_refs, L, U, alg);
    
    lgp_barrier();
    total_tri_cnt = lgp_reduce_add_l(tri_cnt);
    total_sh_refs = lgp_reduce_add_l(sh_refs);
    T0_fprintf(stderr,"  %8.3lf seconds: %16ld triangles", laptime, total_tri_cnt);
    T0_fprintf(stderr,"%16ld shared refs\n", total_sh_refs);
    if ((correct_answer >= 0) && (total_tri_cnt != (int64_t)correct_answer)) {
        T0_fprintf(stderr, "ERROR: Wrong answer!\n");
    }
  
    if(correct_answer == -1) {
        correct_answer = total_tri_cnt;
    }
  
    lgp_barrier();
    lgp_finalize();

    return 0;
}