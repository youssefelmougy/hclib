
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <cmath>
#include <iostream>
#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include <chrono>
#include <thread>
#include "spmv_tile.h"

/*
  Code that mocks up A PART OF CONJUGATE GRADIENT ONLY (Not Whole).
  Reference: https://en.wikipedia.org/wiki/Conjugate_gradient_method 
  Launches iterations spanning only the two inner kernels:

    1) [alpha]_k = dot([r]_k, [r]_k) / dot([p]_k, [Ap]_k)
    2) [x]_k+1   = [x]_k + [alpha]_k * [p]_k

  For reference, the additional kernels part of CG, not spawned in
  this iteration are:

    3) [r]_k+1   = [r]_k - [alpha]_k * [Ap]_k
    4) [beta]_k  = dot([r]_k+1, [r]_k+1) / dot([r]_k, [r]_k)   
    5) [p]_k+1   = [r]_k+1 + [beta]_k * [p]_k
 */



/*
  Struct used to create the array of promises over ntiles*niters
 */
typedef struct {
#if defined(USE_REPLAY)
  hclib::resilience::replay::promise_t<DoubleVec*>* promise;
#elif defined(USE_REPLICATION)
  hclib::resilience::diamond::promise_t<DoubleVec*>* promise;
#else
  hclib::ref_count::promise_t<DoubleVec*>* promise;
#endif
#if defined(DEBUG) || defined(DUMP_TIMINGS)
  unsigned long long task_idx;
#endif
} Tile_t;

typedef hclib::ref_count::promise_t<double*>* double_promise_t;


void initialize_tile(Tile_t &rvec_tile, Tile_t &pvec_tile, Tile_t &xvec_tile, int tile_idx, const spmv_tile * row_tile, const std::vector<double> * eig_vec,
		  const double & eig_val )
{

  hclib::async([=]{

      //compute initial conditions over tile
      DoubleVec *workspace1 = new DoubleVec;
      assert(workspace1);
      workspace1->resize( row_tile->getEndRow() - row_tile->getStartRow() + 1);

      DoubleVec *workspace2 = new DoubleVec;
      assert(workspace2);
      workspace2->resize( row_tile->getEndRow() - row_tile->getStartRow() + 1);

      DoubleVec *workspace3 = new DoubleVec;
      assert(workspace3);
      workspace3->resize( row_tile->getEndRow() - row_tile->getStartRow() + 1);      

      //We want to set b_vector (in the solution of Ax=b) to eigen_value * eigen_vector
      //This makes the algorithm seek the solution of x = eigen vector.

      // So, b = eig_val * eig_vector
      //==>  x_init = 0
      //==>  r_init = b - A*x_init = b = eig_val*eig_vector
      //==>  p_init = r_init = eig_val*eig_vector
      
      for(int k = row_tile->getStartRow(); k<= row_tile->getEndRow(); k++){
	(*workspace1)[k-row_tile->getStartRow()] = eig_val * (*eig_vec)[k];// r_vec
	(*workspace2)[k-row_tile->getStartRow()] = eig_val * (*eig_vec)[k];// p_vec
	(*workspace3)[k-row_tile->getStartRow()] = 0; // x_vec	
      }

      //printf("On tile %d allocated vector of size %d\n", tile_idx, workspace->size() ) ;
      // put initial conditions
      rvec_tile.promise->put(workspace1);
      pvec_tile.promise->put(workspace2);
      xvec_tile.promise->put(workspace3);      

   });
}



void launch_iters(Tile_t **rvec_tile_array, Tile_t **pvec_tile_array, Tile_t **xvec_tile_array, Tile_t **Apvec_tile_array,
		  double_promise_t **pAp_result_tile_array, double_promise_t **rr_result_tile_array,
		  double_promise_t **rr_next_result_tile_array, double_promise_t * alpha_iter, double_promise_t * beta_iter,
		  int n_tiles, int n_iters,
		  const std::vector<spmv_tile> * Amat_row_tiles_ptr)
{

  //One async task per tile per iteration
  for (int tt=0; tt < n_iters; ++tt){

    //----------------------START SPMV PHASE i.e. A*p --------------------------------
    for (int j=0; j<n_tiles; ++j){
      //assemble the list of futures needed for this tile
      //based on the neighbouring tiles depended on
      auto futures = new std::vector<hclib_future_t*>();
      //for (std::set<int>::iterator it = Amat_row_tiles[j].get_dependent_tiles_set().begin(); it!=Amat_row_tiles[j].get_dependent_tiles_set().end(); it++) {
      //	futures->push_back( pvec_tile_array[tt][*it].promise->get_future() );
      //}
      for (auto i : (*Amat_row_tiles_ptr)[j].get_dependent_tiles_set() ) {
	futures->push_back( pvec_tile_array[tt][i].promise->get_future() );
      }

      //advance_tile( futures, pvec_tile_array[tt+1][j], &Amat_row_tiles[j] );
      //Launch async await
      hclib::ref_count::async_await([=]{

#ifdef DUMP_TIMINGS
	  struct timespec begin;
	  clock_gettime(CLOCK_MONOTONIC, &begin);//gettimeofday(&begin,0);
#endif	  
	  std::map<int, DoubleVec*> vec_tiles_I_seek;

	  //for(std::set<int>::iterator it = Amat_row_tiles[j].get_dependent_tiles_set().begin(); it!=Amat_row_tiles[j].get_dependent_tiles_set().end(); it++) {
	  //  vec_tiles_I_seek[*it] = pvec_tile_array[tt][*it].promise->get_future()->get();	  
	  //}
	  for (auto i : (*Amat_row_tiles_ptr)[j].get_dependent_tiles_set() ) {
	    vec_tiles_I_seek[i] = pvec_tile_array[tt][i].promise->get_future()->get();	 
	  }


	  DoubleVec *workspace = new DoubleVec;
	  assert(workspace);
	  (*Amat_row_tiles_ptr)[j].do_tile_matvec(vec_tiles_I_seek, workspace);

          // //Normalise by the eigen value
          // for(std::vector<double>::iterator it = (*workspace).begin(); it != (*workspace).end(); ++it) {
          //   *it /= normalising_val;
          // }  
	  

	  Apvec_tile_array[tt][j].promise->put(workspace);
      	  //pvec_tile_array[tt+1][j].promise->put(workspace);

	  //std::cout<< "On iter " << tt << " tile " << j << std::endl;
	  
#ifdef DUMP_TIMINGS
	  struct timespec end;
	  clock_gettime(CLOCK_MONOTONIC, &end);//gettimeofday(&end,0);
	  double dur = ((end.tv_sec - begin.tv_sec)*1000000000
			+(end.tv_nsec - begin.tv_nsec))*1.0/1000000000;
	  printf("TASKDUR %llu %g\n", j, dur);
#endif	  

	}, futures );     
	
    }

    //----------------------END SPMV PHASE i.e. A*p --------------------------------

      
    //----------------------START alpha DOT PRODUCT PHASE i.e. dot(p,Ap), dot(r,r)-----------------

    //First loop over tiles and let each tile compute and put its local result

    for (int j=0; j<n_tiles; ++j){

      hclib::ref_count::async_await([=]{

	  DoubleVec* Apvec_ptr = Apvec_tile_array[tt][j].promise->get_future()->get();

	  DoubleVec* pvec_ptr = pvec_tile_array[tt][j].promise->get_future()->get();
	  
	  assert(Apvec_ptr->size() == pvec_ptr->size());

	  //do tile-local-dot-product and put in promise
	  double *local_dot_result = new double; *local_dot_result  = 0.0;
	  for(int i = 0; i < pvec_ptr->size(); ++i) {
	    (*local_dot_result) += (*Apvec_ptr)[i] * (*pvec_ptr)[i] ;
	  }

	  pAp_result_tile_array[tt][j]->put(local_dot_result);

	  double *rr_local_dot_result = new double; *rr_local_dot_result = 0.0;
	  DoubleVec* rvec_ptr = rvec_tile_array[tt][j].promise->get_future()->get();
	  for(std::vector<double>::iterator it = rvec_ptr->begin();
	      it != rvec_ptr->end(); ++it) {
	    (*rr_local_dot_result) += (*it) * (*it);
	  }
	  rr_result_tile_array[tt][j]->put(rr_local_dot_result);	 

	  
	}, pvec_tile_array[tt][j].promise->get_future(), Apvec_tile_array[tt][j].promise->get_future(),
	   rvec_tile_array[tt][j].promise->get_future() );


    }
    
    //Now that local results are put, one task to aggregate them all
    //First we need to accumulate all the futures
    auto alpha_res_futures = new std::vector<hclib_future_t*>(2*n_tiles);
    for (int j=0; j<n_tiles; ++j){
      (*alpha_res_futures)[j] = rr_result_tile_array[tt][j]->get_future() ; //first ntiles futures are rr_local_results
      (*alpha_res_futures)[n_tiles+j] = pAp_result_tile_array[tt][j]->get_future() ;//next ntiles futures are pAp_local_results
    }
    
    hclib::ref_count::async_await([=]{
	
	double *alpha_this_iter = new double; *alpha_this_iter = 0.0;
	double r_dot_r = 0.0; double p_dot_Ap = 0.0;
	for (int j=0; j<n_tiles; ++j){
	  r_dot_r += *( rr_result_tile_array[tt][j]->get_future()->get() ) ;
	}
	for (int j=0; j<n_tiles; ++j){
	  p_dot_Ap += *( pAp_result_tile_array[tt][j]->get_future()->get() ) ;
	}

	*alpha_this_iter = r_dot_r / p_dot_Ap ;
	
	alpha_iter[tt]->put(alpha_this_iter);
	
    }, alpha_res_futures ); 
    
    //----------------------END alpha DOT PRODUCT PHASE i.e. dot(p,Ap), dot(r,r) -----------------

    //----------------------START PHASE OF UPDATE OF x_vec, r_vec---------------------------------
    for (int j=0; j<n_tiles; ++j){

      auto update_futures = new std::vector<hclib_future_t*>(5);
      (*update_futures)[0] = xvec_tile_array[tt][j].promise->get_future();
      (*update_futures)[1] = pvec_tile_array[tt][j].promise->get_future();
      (*update_futures)[2] = rvec_tile_array[tt][j].promise->get_future();
      (*update_futures)[3] = Apvec_tile_array[tt][j].promise->get_future();
      (*update_futures)[4] = alpha_iter[tt]->get_future() ;
      
      hclib::ref_count::async_await([=]{

	  DoubleVec* x_next = new DoubleVec; DoubleVec* r_next = new DoubleVec;

	  double alpha = *( alpha_iter[tt]->get_future()->get() );

	  // Update x i.e. x_k+1 = x_k + alpha_k * p_k
	  DoubleVec* xvec_ptr = xvec_tile_array[tt][j].promise->get_future()->get();
	  DoubleVec* pvec_ptr = pvec_tile_array[tt][j].promise->get_future()->get();

	  for(int i = 0; i < xvec_ptr->size(); ++i) {
	    double daxpy_this_elem = (*xvec_ptr)[i] + alpha * (*pvec_ptr)[i] ;
	    x_next->push_back(daxpy_this_elem) ;
	  }	  

	  xvec_tile_array[tt+1][j].promise->put( x_next );

	  // Update r i.e. r_k+1 = r_k - alpha_k * Ap_k
	  DoubleVec* rvec_ptr = rvec_tile_array[tt][j].promise->get_future()->get();
	  DoubleVec* Apvec_ptr = Apvec_tile_array[tt][j].promise->get_future()->get();

	  for(int i = 0; i < rvec_ptr->size(); ++i) {
	    double daxpy_this_elem = (*rvec_ptr)[i] - alpha * (*Apvec_ptr)[i] ;
	    r_next->push_back(daxpy_this_elem) ;
	  }	  

	  rvec_tile_array[tt+1][j].promise->put( r_next );
	  
      }, update_futures);

    }
    
    //----------------------END PHASE OF UPDATE OF x_vec, r_vec---------------------------------

    //----------------------START beta DOT PRODUCT PHASE i.e. dot(r_k+1,r_k+1)-----------------
    //First loop over tiles and let each tile compute and put its local result

    for (int j=0; j<n_tiles; ++j){

      hclib::ref_count::async_await([=]{

	  double *rr_local_dot_result = new double; *rr_local_dot_result = 0.0;
	  DoubleVec* rvec_ptr = rvec_tile_array[tt+1][j].promise->get_future()->get(); //Note that this is tt+1
	  for(std::vector<double>::iterator it = rvec_ptr->begin();
	      it != rvec_ptr->end(); ++it) {
	    (*rr_local_dot_result) += (*it) * (*it);
	  }

	  rr_next_result_tile_array[tt][j]->put(rr_local_dot_result);

	}, rvec_tile_array[tt+1][j].promise->get_future() );
    }

    //Now that local results are put, one task to aggregate them all
    //First we need to accumulate all the futures
    auto beta_res_futures = new std::vector<hclib_future_t*>(2*n_tiles);
    for (int j=0; j<n_tiles; ++j){
      (*beta_res_futures)[j] = rr_result_tile_array[tt][j]->get_future() ; //first ntiles futures are rr_local_results
      (*beta_res_futures)[n_tiles+j] = rr_next_result_tile_array[tt][j]->get_future() ;//next ntiles futures are rr_next_local_results
    }

    hclib::ref_count::async_await([=]{
	
	double *beta_this_iter = new double; *beta_this_iter = 0.0;
	double r_dot_r = 0.0; double r_dot_r_next = 0.0;
	for (int j=0; j<n_tiles; ++j){
	  r_dot_r += *( rr_result_tile_array[tt][j]->get_future()->get() ) ;
	}
	for (int j=0; j<n_tiles; ++j){
	  r_dot_r_next += *( rr_next_result_tile_array[tt][j]->get_future()->get() ) ;
	}

	*beta_this_iter = r_dot_r_next / r_dot_r ;
	
	beta_iter[tt]->put(beta_this_iter);
	
    }, beta_res_futures ); 
    
    //----------------------END beta DOT PRODUCT PHASE i.e. dot(r_k+1,r_k+1)-----------------

    
    //--------------------Final task that update p vector and puts p_k+1---------------------
    for (int j=0; j<n_tiles; ++j){
      hclib::ref_count::async_await([=]{

	  DoubleVec* p_next = new DoubleVec;	  

	  double beta = *( beta_iter[tt]->get_future()->get() );

	  // Update p i.e. p_k+1 = r_k+1 + beta_k * p_k
	  DoubleVec* rvec_ptr = rvec_tile_array[tt+1][j].promise->get_future()->get();
	  DoubleVec* pvec_ptr = pvec_tile_array[tt][j].promise->get_future()->get();

	  for(int i = 0; i < pvec_ptr->size(); ++i) {
	    double daxpy_this_elem = (*rvec_ptr)[i] + beta * (*pvec_ptr)[i] ;
	    p_next->push_back(daxpy_this_elem) ;
	  }	  
	  
	  pvec_tile_array[tt+1][j].promise->put( p_next );

	}, rvec_tile_array[tt+1][j].promise->get_future(), pvec_tile_array[tt][j].promise->get_future(),
	   beta_iter[tt]->get_future() );

    }
    //-------------------End p vector update task-----------------------------------------------------
    
  }          
    
}
  
int main (int argc, char *argv[])
{

 printf("using no resilience\n");
 const char *deps[] = { "system" };
 hclib::launch(deps, 1, [&]() {

   std::string filename;
   filename = argv[1];
   lflr::CSRMatrix  matrix = readMatrixMarketFile( filename );

 assert( checkColIndices(matrix) );

 //Check that the matrix is square
 assert(matrix.getRowSize() == matrix.getColSize() );

 int n_tiles = 4; //default value
 int n_iters = 10; //default value

 n_tiles = atoi(argv[2]);
 n_iters = atoi(argv[3]);

 assert(n_tiles < matrix.getRowSize() );

 //Read the file that has the eigenvalue/vector for reference soln
 std::string eigfile, line;
 eigfile = argv[4];

 std::ifstream eFile;
 eFile.open(eigfile);
 std::getline(eFile, line);//read first line header

 double eig_value;
 std::getline(eFile, line);//read eigen value
 eig_value = std::stod(line);
 std::cout << "Eigen value is " << eig_value << std::endl;

 std::vector<double> eig_vector;
 while( std::getline(eFile, line) )//read eigen vector entries
   {
     eig_vector.push_back(std::stod(line));
   }
 //std::cout << eig_vector[0] << std::endl;
 eFile.close();
 
 //construct the tiles now
 std::vector<spmv_tile> Amat_row_tiles;
 for(int i = 0; i<n_tiles; i++) {
   Amat_row_tiles.push_back( spmv_tile(&matrix, n_tiles, i) );
 }

 //setup dependencies between tiles vectors
 for(int i = 0; i<n_tiles; i++) {
   Amat_row_tiles[i].setup_vec_tiles_dependency();
   //std::cout << "Tile " << i << " row extents: " << Amat_row_tiles[i].getStartRow() <<"-->" << Amat_row_tiles[i].getEndRow()
   //	     << "="<< Amat_row_tiles[i].getEndRow()-Amat_row_tiles[i].getStartRow()+1 << std::endl;
 }


 //Create the arrays of promises
 Tile_t** rvec_tile_array = (Tile_t **) malloc(sizeof(Tile_t*)*(n_iters+1));
 Tile_t** pvec_tile_array = (Tile_t **) malloc(sizeof(Tile_t*)*(n_iters+1));
 Tile_t** xvec_tile_array = (Tile_t **) malloc(sizeof(Tile_t*)*(n_iters+1));

 Tile_t** Apvec_tile_array = (Tile_t **) malloc(sizeof(Tile_t*)*(n_iters+1));

 assert(rvec_tile_array);assert(pvec_tile_array);
 assert(xvec_tile_array); assert(Apvec_tile_array);

 for (int tt=0; tt < n_iters+1; ++tt) {
   rvec_tile_array[tt] = (Tile_t *) malloc(sizeof(Tile_t)*(n_tiles));
   assert(rvec_tile_array[tt]);
   
   pvec_tile_array[tt] = (Tile_t *) malloc(sizeof(Tile_t)*(n_tiles));
   assert(pvec_tile_array[tt]);

   xvec_tile_array[tt] = (Tile_t *) malloc(sizeof(Tile_t)*(n_tiles));
   assert(xvec_tile_array[tt]);

   Apvec_tile_array[tt] = (Tile_t *) malloc(sizeof(Tile_t)*(n_tiles));
   assert(Apvec_tile_array[tt]);

   

   for (int j=0; j < n_tiles; ++j) {
     //Find out how many times the vector of tile j needs to accessed for reference counting
     //The symmetry of matrix, and dependency ensures that
     //  "the number of tiles that depend on me = the number I depend on"        
     int matvec_accesses = Amat_row_tiles[j].get_dependent_tiles_set().size();
     //Total #accesses is num for matvec + one for dot product pAp + one for x/r update task + one for p update
     int num_accesses = matvec_accesses + 3;
// #if defined(USE_REPLAY)
//         pvec_tile_array[tt][j].promise
// 	  = new hclib::resilience::replay::promise_t<DoubleVec*>(tt < n_iters ? num_accesses : 1);
// #elif defined(USE_REPLICATION)
//         pvec_tile_array[tt][j].promise
// 	  = new hclib::resilience::diamond::promise_t<DoubleVec*>(tt < n_iters ? num_accesses : 1);
// #else
       rvec_tile_array[tt][j].promise
	 = new hclib::ref_count::promise_t<DoubleVec*>(4); //Needed only for (1) dot(r,r), (2) r_next = r - alpha * Ap
                                                           // (3) p_next = r_next + beta * p
       pvec_tile_array[tt][j].promise
	 = new hclib::ref_count::promise_t<DoubleVec*>(tt < n_iters ? num_accesses : 1);       
       xvec_tile_array[tt][j].promise
	 = new hclib::ref_count::promise_t<DoubleVec*>(1); //Needed only for (1) x_next = x + alpha * Ap 
       Apvec_tile_array[tt][j].promise
	 = new hclib::ref_count::promise_t<DoubleVec*>(2); //Accessed only for (1) the dot product, (2) r_next = r - alpha * Ap	  
//#endif
       	 assert(rvec_tile_array[tt][j].promise);
       	 assert(pvec_tile_array[tt][j].promise);
	 assert(xvec_tile_array[tt][j].promise);
	 assert(Apvec_tile_array[tt][j].promise);
	
   }
 }

 double_promise_t** pAp_dot_tile_result_array = (double_promise_t **)malloc( sizeof(double_promise_t*)*(n_iters) );
 double_promise_t** rr_dot_tile_result_array = (double_promise_t **)malloc( sizeof(double_promise_t*)*(n_iters) );
 double_promise_t** rr_next_dot_tile_result_array = (double_promise_t **)malloc( sizeof(double_promise_t*)*(n_iters) );

 double_promise_t* alpha_iter = (double_promise_t *)malloc( sizeof(double_promise_t)*(n_iters) );
 double_promise_t* beta_iter = (double_promise_t *)malloc( sizeof(double_promise_t)*(n_iters) );

 for (int tt=0; tt < n_iters; ++tt) {

   pAp_dot_tile_result_array[tt] = (double_promise_t *)malloc(sizeof(double_promise_t)*n_tiles);
   for (int j=0; j < n_tiles; ++j) {
     pAp_dot_tile_result_array[tt][j]
       = ( new hclib::ref_count::promise_t<double*>(1) );
    }

   rr_dot_tile_result_array[tt] = (double_promise_t *)malloc(sizeof(double_promise_t)*n_tiles);
   for (int j=0; j < n_tiles; ++j) {
     rr_dot_tile_result_array[tt][j]
       = ( new hclib::ref_count::promise_t<double*>(2) );
    }

   rr_next_dot_tile_result_array[tt] = (double_promise_t *)malloc(sizeof(double_promise_t)*n_tiles);
   for (int j=0; j < n_tiles; ++j) {
     rr_next_dot_tile_result_array[tt][j]
       = ( new hclib::ref_count::promise_t<double*>(1) );
    }   

   alpha_iter[tt]
     = ( new hclib::ref_count::promise_t<double*>(n_tiles) );

   beta_iter[tt]
     = ( new hclib::ref_count::promise_t<double*>(n_tiles) );
 }


 
 printf("Allocated tile array\n");

 //having the lambda capture the Amat_row_tiles object was a major
 //performance bottleneck (thanks Sriraj) since each lambda makes a copy
 //use  pointer to the Amat_row_tiles vector instead
 auto Amat_row_tiles_ptr = &Amat_row_tiles;
 
 struct timeval begin,end;
 
 //Initialise the vector tiles
 hclib::finish([rvec_tile_array, pvec_tile_array, xvec_tile_array,
		n_tiles, Amat_row_tiles_ptr, &eig_vector, eig_value](){
   for (int j=0; j < n_tiles; ++j){
     initialize_tile(rvec_tile_array[0][j], pvec_tile_array[0][j], xvec_tile_array[0][j],
		     j, &((*Amat_row_tiles_ptr)[j]), &eig_vector, eig_value );
   }
 });

 gettimeofday(&begin,0);
 
 //Launch iterations
 hclib::finish([rvec_tile_array, pvec_tile_array, xvec_tile_array, Apvec_tile_array,
		pAp_dot_tile_result_array, rr_dot_tile_result_array, rr_next_dot_tile_result_array,
		alpha_iter, beta_iter, n_tiles, n_iters, Amat_row_tiles_ptr]() {

    launch_iters(rvec_tile_array, pvec_tile_array, xvec_tile_array, Apvec_tile_array,
		pAp_dot_tile_result_array, rr_dot_tile_result_array, rr_next_dot_tile_result_array,
		alpha_iter, beta_iter, n_tiles, n_iters, Amat_row_tiles_ptr);

 });  

 gettimeofday(&end,0);

 printf("The computation took %f seconds\n",
	((end.tv_sec - begin.tv_sec)*1000000
	 +(end.tv_usec - begin.tv_usec))*1.0/1000000);
 
#ifdef PRINT_SOLN 
 for (int j=0; j < n_tiles; ++j) {
   DoubleVec* soln = (DoubleVec*)pvec_tile_array[n_iters][j].promise->get_future()->get();
   //soln->atomic_print(j);
   for(std::vector<double>::iterator it = soln->begin(); it != soln->end(); it++) {
     std::cout << *it << std::endl;
   }
 }   
#endif
 
 for (int j=0; j < n_tiles; ++j) {
   pvec_tile_array[n_iters][j].promise->get_future()->release();
 }   
 for (int tt=0; tt < n_iters+1; ++tt) {
   free(pvec_tile_array[tt]);
 }   
 free(pvec_tile_array);

 });

 

 return 0;
 
}
