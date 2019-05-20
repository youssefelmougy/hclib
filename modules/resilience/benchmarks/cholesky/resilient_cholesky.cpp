
#include<iostream>
#include<cassert>
#include<numeric>
#include<cmath>
#include "mkl.h"
#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include "hclib_resilience_refstore.h"

//#define USE_IN_PLACE

enum TASK_STATE {NON_LEAF, LEAF};

//#define USE_REF_CT
//#define USE_REPLICATION
//#define USE_REPLAY
//#define USE_ABFT
//#define USE_CHECKPOINT

#define USE_SYRK
//#define PRINT_M

int size_global = -1;
int rate = 1;

//max fault rate allowed is 50%
void add_error(int task_index, int play_index, double *arr) {
    if(play_index == 0 && task_index % rate == 1)
      arr[2] += 3;
}

#ifdef USE_REPLICATION

  #define REF_CT 10000
  namespace diamond = hclib::resilience::diamond;
  namespace resilient {
  template<typename T>
  using promise_t = diamond::promise_t<T>;
  
  template<typename T>
  using future_t = diamond::future_t<T>;
  };

  #define START_TASK \
    hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();\
    diamond::async_await_check<LEAF>([=]() {

  #define END_TASK(...) \
    add_error(num_task, diamond::get_index(), out_arr);\
    }, prom_res, __VA_ARGS__);\
    check_success(prom_res);

#elif defined(USE_REPLAY)

  #define REF_CT 10000
  namespace replay = hclib::resilience::replay;
  namespace resilient {
  template<typename T>
  using promise_t = replay::promise_t<T>;
  
  template<typename T>
  using future_t = replay::future_t<T>;
  };

  #define START_TASK \
    hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();\
    check_params<double> *check_data = new check_params<double>();\
    replay::async_await_check<LEAF>([=]() {

  #define END_TASK(...) \
    add_error(num_task, replay::get_index(), out_arr);\
    }, prom_res, check_error<double>, check_data, __VA_ARGS__);\
    check_success(prom_res, check_data);

#elif defined(USE_ABFT)

  #define REF_CT 10000
  namespace abft = hclib::resilience::abft;
  namespace resilient {
  template<typename T>
  using promise_t = abft::promise_t<T>;

  template<typename T>
  using future_t = abft::future_t<T>;
  };

  #define START_TASK \
    hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();\
    check_params<double> *check_data = new check_params<double>();\
    abft::async_await_check([=]() {

  #define END_TASK(...) \
    add_error(num_task, abft::get_index(), out_arr);\
    }, prom_res, check_error<double>, check_data, \
    [=]() { recover_error<double>(check_data); }, \
    __VA_ARGS__);\
    check_success(prom_res, check_data);

#elif defined(USE_CHECKPOINT)

  #define REF_CT 10000
  namespace checkpoint = hclib::resilience::checkpoint;
  namespace resilient {
  template<typename T>
  using promise_t = checkpoint::promise_t<T>;
  
  template<typename T>
  using future_t = checkpoint::future_t<T>;
  };

  #define START_TASK \
    hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();\
    check_params<double> *check_data = new check_params<double>();\
    checkpoint::async_await_check<LEAF,1>([=]() {

  #define END_TASK(...) \
    add_error(num_task, checkpoint::get_index(), out_arr);\
    }, prom_res, check_error<double>, check_data, __VA_ARGS__);\
    check_success(prom_res, check_data);

#elif defined(USE_REF_CT)

  #define REF_CT 10000
  namespace ref_count = hclib::ref_count;
  namespace resilient {
  template<typename T>
  using promise_t = ref_count::promise_t<T>;

  template<typename T>
  using future_t = ref_count::future_t<T>;
  };

  #define START_TASK \
    hclib::async_await([=]() {

  #define END_TASK(...) \
    }, __VA_ARGS__);

#else

  #define REF_CT 
  namespace resilient {
  template<typename T>
  using promise_t = hclib::promise_t<T>;
  
  template<typename T>
  using future_t = hclib::future_t<T>;
  };

  #define START_TASK \
    hclib::async_await([=]() {

  #define END_TASK(...) \
    }, __VA_ARGS__);

#endif

double *unit_m, *weight_m;

bool double_equal(double i, double j) {
      return (i-j)*(i-j) < 0.000001;
}

#ifdef USE_CHECKPOINT
template<typename T>
struct objArr : public checkpoint::obj {
    T* arr;
    int size;

    objArr() { arr = nullptr; }
    objArr(T* ptr, int size_tmp = size_global) { arr = ptr; size=size_tmp; }
    objArr(checkpoint::archive_obj* ar_ptr) { 
        size = ar_ptr->size/sizeof(T);
        arr = new T[size];
        memcpy(arr, ar_ptr->data, ar_ptr->size);
    }
    ~objArr() { delete arr; arr=nullptr; }

    checkpoint::archive_obj* serialize() {
        auto ar_ptr = new checkpoint::archive_obj();
        ar_ptr->size = sizeof(T)*size;
        ar_ptr->data = malloc(ar_ptr->size);
        memcpy(ar_ptr->data, arr, ar_ptr->size);
        return ar_ptr;
    }
};

#else
template<typename T>
struct objArr : public hclib::resilience::obj {
    T* arr;
    int size;

    objArr() { arr = nullptr; }
    objArr(T* ptr, int size_tmp = size_global) { arr = ptr; size=size_tmp; }
    ~objArr() { delete arr; arr=nullptr; }

    bool equals(obj* obj2){
        bool ret = std::equal(arr, arr+size, ((objArr*)obj2)->arr, double_equal);
        return ret;
    }
};
#endif


template<typename T>
struct check_params {
    T* arr;
    int block_size;
    int block_size_cs;
    T* checksum;
    //char type;

    ~check_params() { delete checksum; }
};

template<typename T>
int check_error(void *ptr) {
    check_params<T> *data = (check_params<T>*)(ptr);
    
    data->checksum = new T[data->block_size];
    cblas_dgemv(CblasRowMajor, CblasTrans, data->block_size, data->block_size, 
            1, data->arr, data->block_size_cs, unit_m, 1, 0, data->checksum, 1);
    T* actual = &data->arr[data->block_size * data->block_size_cs];

    bool ret = std::equal(actual, actual+data->block_size, data->checksum, double_equal);
    //if(ret == false) {
    //    //printf("type %c ", data->type);
    //    for(int i=0;i<data->block_size; i++)
    //        printf("(%f %f) ", actual[i], data->checksum[i]);
    //    printf("\n");
    //}
    return ret;
}

template<typename T>
void recover_error(check_params<T> *data) {

    T* actual = &data->arr[data->block_size * data->block_size_cs];
    for(int i=0; i<data->block_size; i++)
        if(!double_equal(actual[i], data->checksum[i])) {
            T w_actual = actual[i + data->block_size_cs];
            T w_checksum = cblas_ddot(data->block_size, weight_m, 1, &(data->arr[i]), data->block_size_cs);
            if(!double_equal(w_actual, w_checksum)) {
                int err_index = round((w_checksum - w_actual)/(data->checksum[i] - actual[i]))-1;
                int err_index_lin = err_index*data->block_size_cs + i;
                data->arr[err_index_lin] += (actual[i] - data->checksum[i]);
            }
            else
                actual[i] = data->checksum[i];
        }
}

template<typename T>
struct Matrix {
    size_t row_blocks, col_blocks, rows, cols;
    int is_checksum = 0;
    objArr<T> ***ptr;

    Matrix(size_t total_rows, size_t total_cols, size_t row, size_t col) {
        row_blocks = total_rows/row;
        col_blocks = total_cols/col;
        rows = row;
        cols = col;

        assert(row_blocks*rows == total_rows);
        assert(col_blocks*cols == total_cols);

        ptr = new objArr<T>**[row_blocks];
        for(int i=0; i<row_blocks; i++) {
            ptr[i] = new objArr<T>*[col_blocks]{};
        }
    }

    ~Matrix() {
        for(int i=0; i<row_blocks; i++)
            delete ptr[i];
        delete ptr;
    }

    double* create_matrix() {
        double* arr_in = new double[row_blocks*rows*col_blocks*cols];
        double* arr_out = new double[row_blocks*rows*col_blocks*cols];
        int loc = 0;
        for(int i=0;i<row_blocks*rows;i++)  {
          for(int j=0;j<col_blocks*cols;j++) {
              arr_in[loc++] = (double)((loc%10)+2);
              if(j>i) arr_in[loc-1] = 0;
          }
        }
        auto c = col_blocks*cols;
        //char str[] = "MKL_NUM_THREADS=8";
        //putenv(str);
        //printf("%d ",unsetenv("MKL_NUM_THREADS"));
        //printf("%d ",setenv("MKL_NUM_THREADS", "16", 1));
        mkl_set_num_threads(16);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, c, c, c,
                      1, arr_in, c, arr_in, c, 0, arr_out, c);


        delete[] arr_in;
        return arr_out;
    }

    void read_matrix() {
        size_t rows_cs = rows + is_checksum;
        size_t cols_cs = cols + is_checksum;
        int loc =0;
        double * arr_out = create_matrix();

        for(size_t i=0; i<row_blocks; i++) {
          for(size_t j=0; j<rows; j++) {
            for(size_t k=0; k<col_blocks; k++) {
              if(ptr[i][k] == nullptr) {
                  ptr[i][k] = new objArr<T>(new T[rows_cs * cols_cs]);
              }
              for(size_t l=0; l<cols; l++){
                  //std::cin >> ptr[i][k]->arr[j*cols_cs + l];
                  ptr[i][k]->arr[j*cols_cs + l] = arr_out[loc++];
              }
            }
          }
        }
        delete[] arr_out;

        //checksum unit sum row
        if(is_checksum >= 1)
          for(size_t i=0; i<row_blocks; i++) {
            for(size_t k=0; k<col_blocks; k++) {
                cblas_dgemv(CblasRowMajor, CblasTrans, rows, cols, 1, ptr[i][k]->arr, cols_cs, unit_m, 1, 0, &(ptr[i][k]->arr[rows*(cols_cs)]), 1);

                //checksum unit sum column
                //cblas_dgemv(CblasRowMajor, CblasNoTrans, rows, cols, 1, ptr[i][k], cols+1, unit_m, 1, 0, &(ptr[i][k][cols]), cols+1);
                //checksum total i.e at location [n][n]
                //ptr[i][k][(rows+1) * (cols+1)-1] = cblas_ddot(cols, &(ptr[i][k][rows*(cols+1)]), 1, unit_m, 1);
            }
          }

        //checksum weighted sum row
        if(is_checksum >= 2)
          for(size_t i=0; i<row_blocks; i++) {
            for(size_t k=0; k<col_blocks; k++) {
                cblas_dgemv(CblasRowMajor, CblasTrans, rows, cols, 1, ptr[i][k]->arr, cols_cs, weight_m, 1, 0, &(ptr[i][k]->arr[(rows+1)*(cols_cs)]), 1);
            }
          }
    }

    void print_matrix() {
        size_t rows_cs = rows + is_checksum;
        size_t cols_cs = cols + is_checksum;
        printf("----\n");
        for(size_t i=0; i<row_blocks; i++) {
          for(size_t j=0; j<rows_cs; j++) {
            for(size_t k=0; k<col_blocks; k++) {
              if(ptr[i][k] == nullptr) break;
              for(size_t l=0; l<cols; l++){
              //for(size_t l=0; l<cols_cs; l++){
                  printf("%10.3f ", ptr[i][k]->arr[j*cols_cs + l]);
              }
            }
            printf("\n");
          }
        }
    }

    //void zero_upper() {
    //    size_t rows_cs = rows + is_checksum;
    //    size_t cols_cs = cols + is_checksum;
    //    for(size_t i=0; i<row_blocks; i++)
    //      for(size_t j=0; j<rows; j++)
    //        for(size_t l=j+1; l<cols; l++)
    //            ptr[i][i]->arr[j*cols_cs + l] = 0;
    //}

};

void print_matrix(size_t length, double *arr) {
    printf("----\n");
    int index = 0;
    for(int i=0;i<length; i++) {
        for(int j=0;j<length;j++) {
            printf("%11.6f ", arr[index++]);
        }
        std::cout<< std::endl;
    }
}

inline resilient::future_t<objArr<double>*>* create_future_helper(objArr<double> *data) {
    auto prom = new resilient::promise_t<objArr<double>*>(REF_CT);
    prom->put(data);
    return prom->get_future();
}

inline void check_success(hclib::promise_t<int>* prom_res, check_params<double>* check_data = nullptr) {
  hclib::async_await([=]{
    int res = prom_res->get_future()->get();
    if (res == 0){
      fprintf(stderr,"Fatal error: resilience failed to save you\n");
      exit(EXIT_FAILURE);
    }
    delete check_data;
  }, prom_res->get_future());
}

int main(int argc, char* argv[]) {
    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC, &start);
  assert(argc > 2);
  const char *deps[] = { "system" };
  hclib::launch(deps, 1, [=]() {
    int size = atoi(argv[1]);
    int block_size = atoi(argv[2]), res;
    if(argc > 3) {
        float rate_f = atof(argv[3]);
        rate = 1/rate_f;
        printf("fault rate = %f percent\n", rate_f*100);
    }

    assert(size%block_size == 0);

    //fill with value 1
    unit_m = new double[block_size];
    std::fill_n(unit_m, block_size, 1);
    //fill with natural numbers 1..n
    weight_m = new double[block_size];
    std::iota(weight_m, weight_m+block_size, 1);

    Matrix<double> *in_mat = new Matrix<double>(size, size, block_size, block_size);
#ifdef USE_REPLICATION
    in_mat->is_checksum = 0;
    printf("using replication\n");
#elif defined(USE_REPLAY)
    in_mat->is_checksum = 1;
    printf("using replay\n");
#elif defined(USE_ABFT)
    in_mat->is_checksum = 2;
    printf("using abft\n");
#elif defined(USE_CHECKPOINT)
    in_mat->is_checksum = 1;
    checkpoint::set_archive_store(new checkpoint::ref_store());
    printf("using checkpoint\n");
#elif defined(USE_REF_CT)
    in_mat->is_checksum = 0;
    printf("using only ref_count\n");
#else
    in_mat->is_checksum = 0;
    printf("using no resilience\n");
#endif

    size_global = (block_size+in_mat->is_checksum) * (block_size+in_mat->is_checksum);
    in_mat->read_matrix();
#ifdef PRINT_M
    in_mat->print_matrix();
#endif
    //in_mat->zero_upper();

    int num_blocks = size/block_size;
    int block_elems = in_mat->rows * in_mat->cols;
    resilient::promise_t<objArr<double>*> **potrf_prom  = new resilient::promise_t<objArr<double>*>*[num_blocks];
    resilient::promise_t<objArr<double>*> ***trsm_prom  = new resilient::promise_t<objArr<double>*>**[num_blocks];
    resilient::promise_t<objArr<double>*> ***syrk_prom  = new resilient::promise_t<objArr<double>*>**[num_blocks];
    resilient::promise_t<objArr<double>*> ****gemm_prom = new resilient::promise_t<objArr<double>*>***[num_blocks];
    for(int i=0;i<num_blocks;i++) {
        potrf_prom[i] = new resilient::promise_t<objArr<double>*>(REF_CT);
        trsm_prom[i]  = new resilient::promise_t<objArr<double>*>*[num_blocks];
        syrk_prom[i]  = new resilient::promise_t<objArr<double>*>*[num_blocks];
        gemm_prom[i]  = new resilient::promise_t<objArr<double>*>**[num_blocks];
        for(int j=0;j<num_blocks;j++) {
            trsm_prom[i][j] = new resilient::promise_t<objArr<double>*>(REF_CT);
            syrk_prom[i][j] = new resilient::promise_t<objArr<double>*>(REF_CT);
            gemm_prom[i][j] = new resilient::promise_t<objArr<double>*>*[num_blocks];

            for(int k=0;k<num_blocks;k++) {
                gemm_prom[i][j][k] = new resilient::promise_t<objArr<double>*>(REF_CT);
            }
        }
    }

  //block_size++;
  struct timespec begin, end;
  clock_gettime(CLOCK_MONOTONIC, &begin);
  int block_size_cs = block_size + in_mat->is_checksum;
  int arr_size_cs = block_size_cs * block_size_cs;
  hclib::finish([=]() {

        //printf("%d ",unsetenv("MKL_NUM_THREADS"));
        //printf("%d ",setenv("MKL_NUM_THREADS", "1", 1));

        mkl_set_num_threads(1);
        int num_task = 0;

    for(int k=0; k<num_blocks; k++) {

        resilient::future_t<objArr<double>*>* potr_dep = nullptr;
        if(k>0) potr_dep = syrk_prom[k][k-1]->get_future();
        else    potr_dep = create_future_helper(in_mat->ptr[k][k]);
        num_task++;
        START_TASK
            double* data_potr_dep = potr_dep->get()->arr;
#ifdef USE_IN_PLACE
            double* out_potr = data_potr_dep;
#else
            double* out_potr = new double[arr_size_cs];
            memcpy(out_potr, data_potr_dep, sizeof(double)* arr_size_cs);
#endif
            LAPACKE_dpotrf2(LAPACK_ROW_MAJOR, 'L', block_size_cs, out_potr, block_size_cs);
            potrf_prom[k]->put(new objArr<double>(out_potr));

#if defined(USE_REPLAY) || defined(USE_ABFT) || defined(USE_CHECKPOINT)
            //zeroing the upper triangle
            for(int i=0; i<block_size; i++)
                for(int j=i+1; j<block_size; j++)
                    out_potr[i*block_size_cs + j] = 0;
            check_data->arr = out_potr;
            check_data->block_size = block_size;
            check_data->block_size_cs = block_size_cs;
            //check_data->type = 'p';
#endif
            auto out_arr = out_potr;
        END_TASK(potr_dep)

        for(int m=k+1; m<num_blocks; m++) {

            resilient::future_t<objArr<double>*>* trsm_dep1 = potrf_prom[k]->get_future();
            resilient::future_t<objArr<double>*>* trsm_dep2 = nullptr;
            if(k>0) trsm_dep2 = gemm_prom[m][k][k-1]->get_future();
            else    trsm_dep2 = create_future_helper(in_mat->ptr[m][k]);
            num_task++;
            START_TASK
                double* data_trsm_dep1 = trsm_dep1->get()->arr;
                double* data_trsm_dep2 = trsm_dep2->get()->arr;
#ifdef USE_IN_PLACE
                double* out_trsm = data_trsm_dep2;
#else
                double* out_trsm = new double[arr_size_cs];
                memcpy(out_trsm, data_trsm_dep2, sizeof(double)* arr_size_cs);
#endif
                cblas_dtrsm(CblasRowMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, block_size_cs,
                  block_size, 1, data_trsm_dep1, block_size_cs, out_trsm, block_size_cs);
                trsm_prom[m][k]->put(new objArr<double>(out_trsm));

#if defined(USE_REPLAY) || defined(USE_ABFT) || defined(USE_CHECKPOINT)
                check_data->arr = out_trsm;
                check_data->block_size = block_size;
                check_data->block_size_cs = block_size_cs;
                //check_data->type = 't';
#endif
                auto out_arr = out_trsm;
            END_TASK(trsm_dep1, trsm_dep2)
        }

        for(int n=k+1; n<num_blocks; n++) {

            resilient::future_t<objArr<double>*>* syrk_dep1 = trsm_prom[n][k]->get_future();
            resilient::future_t<objArr<double>*>* syrk_dep2 = nullptr;
            if(k>0) syrk_dep2 = syrk_prom[n][k-1]->get_future();
            else    syrk_dep2 = create_future_helper(in_mat->ptr[n][n]);
            num_task++;
            START_TASK
                double* data_syrk_dep1 = syrk_dep1->get()->arr;
                double* data_syrk_dep2 = syrk_dep2->get()->arr;
#ifdef USE_IN_PLACE
                double* out_syrk = data_syrk_dep2;
#else
                double* out_syrk = new double[arr_size_cs];
                memcpy(out_syrk, data_syrk_dep2, sizeof(double)* arr_size_cs);
#endif

#ifdef USE_SYRK
                cblas_dsyrk(CblasRowMajor, CblasLower, CblasNoTrans, block_size_cs, block_size, -1,
                  data_syrk_dep1, block_size_cs, 1, out_syrk, block_size_cs);
#else
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, block_size_cs, block_size_cs, block_size,
                  -1, data_syrk_dep1, block_size_cs, data_syrk_dep1, block_size_cs, 1, out_syrk, block_size_cs);
#endif
                syrk_prom[n][k]->put(new objArr<double>(out_syrk));

#if defined(USE_REPLAY) || defined(USE_ABFT) || defined(USE_CHECKPOINT)
#ifdef USE_SYRK
                for(int i=0; i<block_size; i++)
                    for(int j=i+1; j<block_size; j++)
                      out_syrk[i*block_size_cs + j] = out_syrk[j*block_size_cs + i];
#endif
                check_data->arr = out_syrk;
                check_data->block_size = block_size;
                check_data->block_size_cs = block_size_cs;
                //check_data->type = 's';
#endif
                auto out_arr = out_syrk;
            END_TASK(syrk_dep1, syrk_dep2)

            for (int m=n+1; m<num_blocks; m++) {
                resilient::future_t<objArr<double>*>* gemm_dep1 = trsm_prom[m][k]->get_future();
                resilient::future_t<objArr<double>*>* gemm_dep2 = trsm_prom[n][k]->get_future();
                resilient::future_t<objArr<double>*>* gemm_dep3 = nullptr;
                if(k>0) gemm_dep3 = gemm_prom[m][n][k-1]->get_future();
                else    gemm_dep3 = create_future_helper(in_mat->ptr[m][n]);
                num_task++;
                START_TASK
                    double* data_gemm_dep1 = gemm_dep1->get()->arr;
                    double* data_gemm_dep2 = gemm_dep2->get()->arr;
                    double* data_gemm_dep3 = gemm_dep3->get()->arr;
#ifdef USE_IN_PLACE
                    double* out_gemm = data_gemm_dep3;
#else
                    double* out_gemm = new double[arr_size_cs];
                    memcpy(out_gemm, data_gemm_dep3, sizeof(double)* arr_size_cs);
#endif
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, block_size_cs, block_size_cs, block_size,
                      -1, data_gemm_dep1, block_size_cs, data_gemm_dep2, block_size_cs, 1, out_gemm, block_size_cs);
                    gemm_prom[m][n][k]->put(new objArr<double>(out_gemm));
                    //for(int i=0;i<block_size;i++)
                    //    out_gemm[(block_size-1)*block_size_cs + i] -= 1;
#if defined(USE_REPLAY) || defined(USE_ABFT) || defined(USE_CHECKPOINT)
                    check_data->arr = out_gemm;
                    check_data->block_size = block_size;
                    check_data->block_size_cs = block_size_cs;
                    //check_data->type = 'g';
#endif
                    auto out_arr = out_gemm;
                END_TASK(gemm_dep1, gemm_dep2, gemm_dep3, nullptr)
            }
        }
    }
    printf("num_task %d\n", num_task);
  });

  clock_gettime(CLOCK_MONOTONIC, &end);
  double dur = ((end.tv_sec - begin.tv_sec)*1000000000
               +(end.tv_nsec - begin.tv_nsec))*1.0/1000000000;
  printf("The computation took %f seconds\n", dur);
    double dur2 = ((end.tv_sec - start.tv_sec)*1000000000
                  +(end.tv_nsec - start.tv_nsec))*1.0/1000000000;
    //printf("The computation full took %f seconds\n", dur2);

#ifdef PRINT_M
    std::cout<<"Result"<<std::endl;
    //assemble output
    Matrix<double> *out_mat = new Matrix<double>(size, size, block_size, block_size);
    out_mat->is_checksum = in_mat->is_checksum;
    for(int i=0; i<num_blocks; i++) {
        out_mat->ptr[i][i] = potrf_prom[i]->get_future()->get();
        for(int j=0; j<i; j++)
            out_mat->ptr[i][j] = trsm_prom[i][j]->get_future()->get();
    }
    ////in_mat->print_matrix();
    out_mat->print_matrix();
#endif

  });
}

