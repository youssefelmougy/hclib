#include <math.h>
#include <shmem.h>
#include <stdio.h>

extern "C" {

#include <spmat.h>

}

#include "hclib_bale_actor.h"
#include "mytime.h"
#include "selector.h"


#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

typedef struct TrianglePkt {
    int64_t w;
    int64_t vj;
} TrianglePkt;

enum MailBoxType {REQUEST};

class TriangleSelector: public hclib::Selector<1, TrianglePkt> {
public:
    TriangleSelector(int64_t* cnt, sparsemat_t* mat) : cnt_(cnt), mat_(mat) {
        mb[REQUEST].process = [this] (TrianglePkt pkt, int sender_rank) { 
            this->req_process(pkt, sender_rank);
        };
    }

private:
    //shared variables
    int64_t* cnt_;
    sparsemat_t* mat_;

    void req_process(TrianglePkt pkg, int sender_rank) {
        int64_t tempCount = 0;

        for (int64_t k = mat_->loffset[pkg.vj]; k < mat_->loffset[pkg.vj + 1]; k++) {
            if (pkg.w == mat_->lnonzero[k]) {
                tempCount++;
                break;
            }
            if (pkg.w < mat_->lnonzero[k]) { // requires that nonzeros are increasing
                break;
            }
        }

        // update the share variable
        *cnt_ += tempCount;
    }
};


double triangle_selector(int64_t* count, int64_t* sr, sparsemat_t* L, sparsemat_t* U, int64_t alg) {
    int64_t numpushed = 0;

    if (!L) {
        T0_printf("ERROR: triangle_selector: NULL L!\n");
    }

    // Start timing
    double t1 = wall_seconds();
  
    if (alg == 0) {
        TriangleSelector* triSelector = new TriangleSelector(count, L);

        hclib::finish([=, &numpushed]() {
            hclib::selector::finish(triSelector, [=, &numpushed]() {
                int64_t k,kk, pe;
                int64_t l_i, L_i, L_j;

                TrianglePkt pkg;
                //foreach nonzero (i, j) in L
                for (l_i = 0; l_i < L->lnumrows; l_i++) { 
                    for (k = L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
                        L_i = l_i * THREADS + MYTHREAD;
                        L_j = L->lnonzero[k];

                        pe = L_j % THREADS;
                        pkg.vj = L_j / THREADS;
                        for (kk = L->loffset[l_i]; kk < L->loffset[l_i + 1]; kk++) {
                            pkg.w = L->lnonzero[kk];

                            if (pkg.w > L_j) {
                                break;
                            }

                            numpushed++;
                            // if (convey_push(conv, &pkg, pe) != convey_OK) {
                            //     tri_convey_push_process(&cnt, conv, L, 0); 
                            //     kk--;
                            //     numpushed--;
                            // }

                            // block above is the original sending mechanism
                            // assume sending will not fail
                            triSelector->send(REQUEST, pkg, pe);
                        }
                    }    
                }

                // Indicate that we are done with sending messages to the REQUEST mailbox
                triSelector->done(REQUEST);
            });
        });
        
    } else {
        if (!U) {
            T0_printf("ERROR: triangle_selector: NULL U!\n");
            
            return -1;
        }

        TriangleSelector* triSelector = new TriangleSelector(count, U);

        hclib::finish([=, &numpushed] () {
            hclib::selector::finish(triSelector, [=, &numpushed] () {
                int64_t k,kk, pe;
                int64_t l_i, L_i, L_j;
                TrianglePkt pkg;

                //foreach nonzero (i, j) in L
                for (l_i = 0; l_i < L->lnumrows; l_i++) { 
                    for (k = L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
                        L_i = l_i * THREADS + MYTHREAD;
                        L_j = L->lnonzero[k];
                
                        pe = L_j % THREADS;
                        pkg.vj = L_j / THREADS;

                        for (kk = U->loffset[l_i]; kk < U->loffset[l_i + 1]; kk++) {
                            pkg.w = U->lnonzero[kk]; 
                            numpushed++;

                            // if (convey_push(conv, &pkg, pe) != convey_OK){
                            //     tri_convey_push_process(&cnt, conv, U, 0); 
                            //     kk--;
                            //     numpushed--;
                            // }

                            // block above is the original sending mechanism
                            // assume sending will not fail
                            triSelector->send(REQUEST, pkg, pe);
                        }
                    }
                }

                // Indicate that we are done with sending messages to the REQUEST mailbox
                triSelector->done(REQUEST);
            });
        });
    }

    shmem_barrier_all();
    *sr = numpushed;
    // Updating count is not necessary since that is taken cared by the mailbox logic
    // *count = cnt;
    t1 = wall_seconds() - t1;

    return t1;
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

static void* setup_shmem_reduce_workdata(long** psync, size_t xsize) {
    long* work;
    long i;

    work = (long*)shmem_malloc(_SHMEM_REDUCE_MIN_WRKDATA_SIZE * xsize);
    *psync = (long*)shmem_malloc(_SHMEM_REDUCE_SYNC_SIZE * sizeof(long));

    for (i = 0; i <_SHMEM_REDUCE_SYNC_SIZE; ++i) {
        (*psync)[i] = _SHMEM_SYNC_VALUE;
    }

    shmem_barrier_all();
    return work;
}

long reduce_add_l(long myval) {
    static long* buff = nullptr;
    static long* work;
    static long* sync;

    if (!buff) {
        buff = (long*)shmem_malloc(2 * sizeof(long));
        work = (long*)setup_shmem_reduce_workdata(&sync, sizeof(long));
    }

    buff[0] = myval;
    shmem_long_sum_to_all(&buff[1], buff, 1, 0, 0, THREADS, work, sync);
    shmem_barrier_all();
    return buff[1];
}

int main(int argc, char* argv[]) {
    const char *deps[] = { "system", "bale_actor" };
    hclib::launch(deps, 2, [=] {

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

        // if (printhelp) usage(); // Skipping print help
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

        double correct_answer = -1;

        sparsemat_t *A, *L, *U;
        
        if (read_graph) {
            A = read_matrix_mm_to_dist(filename);
            if (!A) return -1;
            
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

        } else if (gen_kron_graph) {
            // string should be <mode> # # ... #
            // we will break the string of numbers (#s) into two groups and create
            // two local kronecker graphs out of them.
            int num;
            char* ptr = kron_graph_string;
            int64_t* kron_specs = (int64_t*)calloc(32, sizeof(int64_t *));
    
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
                T0_fprintf(stderr, "ERROR: invalid kronecker product string (%s): must contain at least three integers\n", kron_graph_string); return(-1);
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
                for (i = 0; i < num_ints; i++) {
                    x *= (kron_specs[i] + 1);
                }

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

        shmem_barrier_all();
        if (alg == 1)
            U = transpose_matrix(L);   

        shmem_barrier_all();

        T0_fprintf(stderr, "L has %ld rows/cols and %ld nonzeros.\n", L->numrows, L->nnz);
  
        if (!is_lower_triangular(L, 0)) {
            T0_fprintf(stderr,"ERROR: L is not lower triangular!\n");
            return -1;
        }
    
        T0_fprintf(stderr, "Run triangle counting ...\n");
        int64_t tri_cnt;           // partial count of triangles on this thread
        int64_t total_tri_cnt;     // the total number of triangles on all threads
        int64_t sh_refs;         // number of shared reference or pushes
        int64_t total_sh_refs;

        int64_t* cc = (int64_t*)shmem_malloc(L->numrows * sizeof(int64_t));
        int64_t* l_cc = cc;
        for (i = 0; i < L->lnumrows; i++)
            l_cc[i] = 0;
            
        shmem_barrier_all();
  
        /* calculate col sums */
        for (i = 0; i < L->lnnz; i++) {
            long lindex = L->lnonzero[i] / THREADS;
            long pe = L->lnonzero[i] % THREADS;
            shmem_long_finc(&cc[lindex], pe);
        }
  
        shmem_barrier_all();
  
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
            pulls_calc = reduce_add_l(rtimesc_calc);
            pushes_calc = reduce_add_l(rchoose2_calc);
        } else {
            pushes_calc = reduce_add_l(rtimesc_calc);
            pulls_calc = reduce_add_l(cchoose2_calc);
        }

        shmem_free(cc);

        T0_fprintf(stderr,"Calculated: Pulls = %ld\n            Pushes = %ld\n\n", pulls_calc, pushes_calc);

        int64_t use_model;
        double laptime = 0.0;

        tri_cnt = 0;
        total_tri_cnt = 0;
        sh_refs = 0;
        total_sh_refs = 0;

        // only running selector model
        T0_fprintf(stderr, "Running Selector: ");
        laptime = triangle_selector(&tri_cnt, &sh_refs, L, U, alg);
        shmem_barrier_all();

        total_tri_cnt = reduce_add_l(tri_cnt);
        total_sh_refs = reduce_add_l(sh_refs);
        T0_fprintf(stderr, "  %8.3lf seconds: %16ld triangles", laptime, total_tri_cnt);
        T0_fprintf(stderr, "%16ld shared refs\n", total_sh_refs);
        
        if ((correct_answer >= 0) && (total_tri_cnt != (int64_t)correct_answer)) {
            T0_fprintf(stderr, "ERROR: Wrong answer!\n");
        }
  
        if(correct_answer == -1) {
            correct_answer = total_tri_cnt;
        }

        shmem_barrier_all();
        shmem_finalize();
        
        return 0;
    });
}