/*
 * Copyright 2017 Rice University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// This version will compute multiple time steps per task if n_steps_task > 1.
// This is a variant of Jackson's 3-to-1 tasks approach where the task may
// not advance the maximum allowable number of steps given by the tile size.
// This allows us to find the best tradeoff between task granularity and
// redundant work. This code is used to analyze performance degradations
// introduced by the chosen resilience mechanism.

#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <cmath>
#include "hclib_cpp.h"
#include "hclib_resilience.h"
#include <chrono>
#include <thread>
#include <numeric>
#include <iostream>
#include <unordered_map>
#include <cstring>

//#define DEBUG

// this will dump raw time stamps, not elapsed time
//#define DUMP_TIMINGS
//#define DUMP_DURATION
//#define TIMING_INTERVAL 8

// instead of dumping at regular intervals, only dump timings for
// tasks that we classify as failed
//#define DUMP_FAIL_DURATION

// skip injecting the failure so we can see difference in timings
// between failing and not failing same task
//#define MOCK_FAIL

#ifdef DUMP_FAIL_DURATION
#ifdef TIMING_INTERVAL
#undef TIMING_INTERVAL
#endif
#ifndef DUMP_DURATION
#define DUMP_DURATION
#endif
#ifdef DUMP_TIMINGS
#undef DUMP_TIMINGS
#endif
#endif

#if defined(DUMP_TIMINGS) && defined(DUMP_DURATION)
#error Define at most one of DUMP_TIMINGS and DUMP_DURATION
#endif

#ifdef DEBUG
#include <sstream>
#endif

#if defined(USE_REPLAY) && defined(USE_REPLICATION)
#error Define at most one of USE_REPLAY and USE_REPLICATION
#endif

#define PI 3.14159265358979323846  /* pi */

enum TASK_STATE {NON_LEAF, LEAF};

typedef int int_t;

typedef struct {
  // default values defined below; can be overridden on command line
  int_t n_tiles = 64;
  int_t tilesize = 128;
  int_t n_steps_total = 128;
  int_t n_steps_task = 16;
  bool do_checksums = true;
  bool do_injection = false;
  std::string failfile;
} params_t;

class DoubleVec : public hclib::resilience::obj,
                 public std::vector<double> {
public:
  virtual ~DoubleVec() {}

  bool equals(obj* obj2){
    DoubleVec *o2 = (DoubleVec*)obj2;
    return *this == *o2;
  }

#ifdef DEBUG
  void atomic_print(unsigned long long task_idx){
    std::stringstream ss;
    ss << task_idx << ": ";
    for (auto it=cbegin(); it!=cend(); ++it){
      ss << *it << " ";
    }
    ss << std::endl;
    std::cout << ss.str();
  }
#endif

  void compute_chk(int_t n_steps_task){
    // need checksums over three potentially-distinct sets
    // when n_steps_task == tilesize, these will be identical
    if ((size_t)n_steps_task == size()){
      left_chk_ = right_chk_ = chk_ = std::accumulate(begin(), end(), 0.0);
    } else {
      chk_ = std::accumulate(begin(), end(), 0.0);
      left_chk_ = std::accumulate(begin(), begin()+n_steps_task, 0.0);
      right_chk_ = std::accumulate(end()-n_steps_task, end(), 0.0);
    }
  }

  double chk() const { return chk_; }
  double chk_left() const { return left_chk_; }
  double chk_right() const { return right_chk_; }

  bool chk_within_tol(double independent_chk) const {
    return (std::abs(chk_ - independent_chk) <= checksum_tol);
  }
 
private:
  double chk_ = 0.0;
  double left_chk_ = 0.0;
  double right_chk_ = 0.0;

  static constexpr double checksum_tol = 1e-5;
};

#if defined(USE_REPLICATION)
class deterministic_replication_failure_log {
 public:
  void push_failure(unsigned long long task_idx, int replica_idx){
    assert(replica_idx >= 0 && replica_idx < 3);
    int &bits = failures_[task_idx];
    bits |= 2<<replica_idx;
  }

  int should_fail(unsigned long long task_idx) const {
    auto it = failures_.find(task_idx);
    if (it == failures_.end()){
      return 0;
    } else {
      return it->second;
    }
  }

  static bool replica_should_fail(int fails, int replica_idx){
    int mask = 2<<replica_idx;
    return ((fails & mask) == mask);
  }

 private:
  std::unordered_map<unsigned long long, int/*replica_bits*/> failures_;
};

typedef deterministic_replication_failure_log failure_generator_t;
#else
class deterministic_replay_failure_log {
 public:
  void push_failure(unsigned long long task_idx){
    int &count = failures_[task_idx];
    ++count;
  }

  int should_fail(unsigned long long task_idx) const {
    auto it = failures_.find(task_idx);
    if (it == failures_.end()){
      return 0;
    } else {
      return it->second;
    }
  }

  static bool replay_should_fail(int fails, int replay_idx){
    return (replay_idx < fails);
  }

 private:
  std::unordered_map<unsigned long long, int/*failure_count*/> failures_;
};

typedef deterministic_replay_failure_log failure_generator_t;
#endif

typedef struct {
#if defined(USE_REPLAY)
  hclib::resilience::replay::promise_t<DoubleVec*>* promise;
#elif defined(USE_REPLICATION)
  hclib::resilience::diamond::promise_t<DoubleVec*>* promise;
#else
  hclib::ref_count::promise_t<DoubleVec*>* promise;
#endif
#if defined(DEBUG) || defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
  unsigned long long task_idx;
#endif
} Tile_t;

#if defined(USE_REPLAY)
int check(void *args)
{
  bool succeeded = *(bool*)args;
  return (succeeded ? 1 : 0);
}
#endif

double initial_conditions(double x)
{
  return std::sin(2.0 * PI * x);
}

// Lax-Wendroff stencil for
// \partial u / \partial t + \partial u / \partial x = 0
double stencil(double left, double center, double right, double cfl)
{
  const double cflsqr = cfl * cfl;
  return 0.5 * (cfl + cflsqr) * left
         + (1.0 - cflsqr) * center
         + 0.5 * (cflsqr - cfl) * right;
}

// flux near left edge of domain that decreases in size each iter
double left_flux(double left, double center, double cfl)
{
  const double cflsqr = cfl * cfl;
  return (1.0 - 0.5 * cfl - 0.5 * cflsqr) * left
         + 0.5 * (cflsqr - cfl) * center;
}

// flux near right edge of domain that decreases in size each iter
double right_flux(double center, double right, double cfl)
{
  const double cflsqr = cfl * cfl;
  return 0.5 * (cfl + cflsqr) * center
         + (1.0 + 0.5 * cfl - 0.5 * cflsqr) * right;
}

void initialize_tile(Tile_t &tile, int_t tile_idx, int_t tilesize,
                     int_t domainsize, int_t n_steps_task,
                     bool do_checksums)
{
#if defined(USE_REPLAY)
  bool *succeeded = new bool;
  hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();
  hclib::resilience::replay::async_await_check<LEAF,3>([=]{
#elif defined(USE_REPLICATION)
  hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();
  hclib::resilience::diamond::async_await_check<LEAF>([=]{
#else
  hclib::async([=]{
#endif
    double checksum = 0.0;
    // compute initial conditions over tile
    DoubleVec *workspace = new DoubleVec;
    assert(workspace);
    workspace->resize(tilesize);
    for (int_t k=0; k<tilesize; ++k){
      double m = tile_idx * tilesize + k;
      (*workspace)[k] = initial_conditions(m/domainsize);
      if (do_checksums)
        checksum += (*workspace)[k];
    }

    if (do_checksums){
      workspace->compute_chk(n_steps_task);
      bool pass = workspace->chk_within_tol(checksum);
#if defined(USE_REPLAY)
      *succeeded = pass;
#endif
      if (!pass){
#if !defined(DUMP_TIMINGS) && !defined(DUMP_DURATION)
        //fprintf(stderr,"checksums do not match: %g %g %g\n",
        //        checksum, workspace->chk(), checksum-workspace->chk());
#endif
#if !defined(USE_REPLAY) && !defined(USE_REPLICATION)
        exit(EXIT_FAILURE);
#endif
      }
    }

    // put initial conditions
    tile.promise->put(workspace);
#if defined(USE_REPLAY)
  }, prom_res, check, (void*)(succeeded), (hclib_future_t*)NULL);
#elif defined(USE_REPLICATION)
  }, prom_res, (hclib_future_t*)NULL);
#else
  });
#endif

#if defined(USE_REPLAY) || defined(USE_REPLICATION)
  hclib::async_await([=]{
    int res = prom_res->get_future()->get();
    if (res == 0){
      fprintf(stderr,"Fatal error: resilience failed to save you\n");
      exit(EXIT_FAILURE);
    }
#if defined(USE_REPLAY)
    delete succeeded;
#endif
  }, prom_res->get_future());
#endif
} 

void advance_tile(Tile_t &left_tile_prev, Tile_t &right_tile_prev,
                  Tile_t &center_tile_prev, Tile_t &center_tile_updated,
                  int_t n_steps_task, double cfl, bool do_checksums,
                  int fails
#if defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
                  ,bool do_timing
#endif
                  )
{
#if defined(USE_REPLAY)
  bool *succeeded = new bool;
  hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();
  hclib::resilience::replay::async_await_check<LEAF,3>([=]{
#elif defined(USE_REPLICATION)
  hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();
  hclib::resilience::diamond::async_await_check<LEAF>([=]{
#else
  hclib::ref_count::async_await([=]{
#endif

#if defined(DUMP_DURATION)
    struct timespec begin;
    if (do_timing)
      clock_gettime(CLOCK_MONOTONIC, &begin);
#endif

    DoubleVec* prev_ghost_to_left
        = (DoubleVec*)left_tile_prev.promise->get_future()->get();
    DoubleVec* prev_ghost_to_right
        = (DoubleVec*)right_tile_prev.promise->get_future()->get();
    DoubleVec* prev_soln
        = (DoubleVec*)center_tile_prev.promise->get_future()->get();

    // assuming same number of ghost zones to right as left
    int_t n_ghost = prev_ghost_to_left->size();
    int_t n_extra = n_ghost - n_steps_task;
    int_t tilesize = prev_soln->size();
    int_t extended_size = tilesize+2*n_steps_task;

    DoubleVec *workspace = new DoubleVec;
    assert(workspace);
    workspace->reserve(extended_size);

    // put left ghost, main tile, right ghost all contiguously in buffer
    workspace->insert(workspace->end(), prev_ghost_to_left->begin()+n_extra,
                                        prev_ghost_to_left->end());
    workspace->insert(workspace->end(), prev_soln->begin(),
                                        prev_soln->end());
    workspace->insert(workspace->end(), prev_ghost_to_right->begin(),
                                        prev_ghost_to_right->end()-n_extra);

    double checksum = 0.0;
    if (do_checksums)
      checksum = prev_ghost_to_left->chk_right() + prev_soln->chk()
               + prev_ghost_to_right->chk_left();

    // advance by as many iterations as ghost zones allow (left-justifying soln)
    for (int_t t=0; t < n_steps_task; ++t){
      if (do_checksums){
        checksum -= left_flux((*workspace)[0],
                              (*workspace)[1], cfl);
        checksum -= right_flux((*workspace)[extended_size-2-2*t],
                               (*workspace)[extended_size-1-2*t], cfl);
      }
      for (int_t k=0; k < extended_size-2-2*t; ++k){
        (*workspace)[k] = stencil((*workspace)[k],
                                  (*workspace)[k+1],
                                  (*workspace)[k+2], cfl);
      }
    }


    if (fails){
      int idx = 0;
#if defined(USE_REPLAY)
      idx = hclib::resilience::replay::get_index();
      if (failure_generator_t::replay_should_fail(fails, idx)){
#elif defined(USE_REPLICATION)
      idx = hclib::resilience::diamond::get_index();
      if (failure_generator_t::replica_should_fail(fails, idx)){
#endif
#ifndef MOCK_FAIL
        //printf("injecting failure\n");
        double *val = &((*workspace)[tilesize/2]);
        uint64_t tmp;
        std::memcpy(&tmp, val, sizeof(double));
        tmp ^= 2ULL<<(60-idx);
        std::memcpy(val, &tmp, sizeof(double));
#endif
#if defined(USE_REPLAY) || defined(USE_REPLICATION)
      }
#endif
    }

    // cut off the excess to keep only the fully advanced main tile values 
    workspace->resize(tilesize);

    if (do_checksums){
      workspace->compute_chk(n_steps_task);
      bool pass = workspace->chk_within_tol(checksum);
#if defined(USE_REPLAY)
      *succeeded = pass;
#endif
      if (!pass){
#if !defined(DUMP_TIMINGS) && !defined(DUMP_DURATION)
        //fprintf(stderr,"checksums do not match: %g %g %g\n",
        //        checksum, workspace->chk(), checksum-workspace->chk());
#endif
#if !defined(USE_REPLAY) && !defined(USE_REPLICATION)
        exit(EXIT_FAILURE);
#endif
      }
    }

    // put the advanced values
    center_tile_updated.promise->put(workspace);

#if defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
    if (do_timing){
      struct timespec end;
      clock_gettime(CLOCK_MONOTONIC, &end);
#if defined(DUMP_DURATION)
      double dur = ((end.tv_sec - begin.tv_sec)*1000000000
                   +(end.tv_nsec - begin.tv_nsec))*1.0/1000000000;
      printf("TASKDUR %llu %g\n", center_tile_prev.task_idx, dur);
#else
      printf("TASKN %llu %ld %ld\n", center_tile_prev.task_idx, end.tv_sec, end.tv_nsec);
#endif
    }
#endif
  },
#if defined(USE_REPLAY)
    prom_res, check, (void*)(succeeded),
#elif defined(USE_REPLICATION)
    prom_res,
#endif
    left_tile_prev.promise->get_future(),
    right_tile_prev.promise->get_future(),
    center_tile_prev.promise->get_future(), NULL);

#if defined(USE_REPLAY) || defined(USE_REPLICATION)
  hclib::async_await([=]{
    int res = prom_res->get_future()->get();
    if (res == 0){
      fprintf(stderr,"Fatal error: resilience failed to save you\n");
      exit(EXIT_FAILURE);
    }
#if defined(USE_REPLAY)
    delete succeeded;
#endif
  }, prom_res->get_future());
#endif
} 

void launch_iters(Tile_t **tile_array, int_t n_tiles, int_t n_rounds,
                  int_t n_steps_task, double cfl, bool do_checksums,
                  const failure_generator_t &failgen)
{
  // one task per tile per round
  for (int_t tt=0; tt < n_rounds; ++tt){
    for (int_t j=0; j<n_tiles; ++j){
      int fails = failgen.should_fail(1+j + tt*n_tiles);
      int_t left_nbr = (n_tiles+j-1)%n_tiles;
      int_t right_nbr = (j+1)%n_tiles;
      advance_tile(tile_array[tt][left_nbr], tile_array[tt][right_nbr],
                   tile_array[tt][j], tile_array[tt+1][j],
                   n_steps_task, cfl, do_checksums, fails
#if defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
#ifdef DUMP_FAIL_DURATION
                   ,fails
#else
                   ,(tt+1)%TIMING_INTERVAL == 0
#endif
#endif
      );
    }
  }
}

void get_params(int argc, char *argv[], params_t &params)
{
  params_t defaults(params);
  for (int i=1; i < argc; ++i){
    if (strcmp(argv[i], "-ntiles") == 0){
      if (i+1 == argc){
        std::cerr << "Error: didn't specify value for ntiles\n";
        exit(EXIT_FAILURE);
      }
      params.n_tiles = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-tilesize") == 0){
      if (i+1 == argc){
        std::cerr << "Error: didn't specify value for tilesize\n";
        exit(EXIT_FAILURE);
      }
      params.tilesize = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-nsteps") == 0){
      if (i+1 == argc){
        std::cerr << "Error: didn't specify value for nsteps\n";
        exit(EXIT_FAILURE);
      }
      params.n_steps_total = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-nstepstask") == 0){
      if (i+1 == argc){
        std::cerr << "Error: didn't specify value for nstepstask\n";
        exit(EXIT_FAILURE);
      }
      params.n_steps_task = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-checksum") == 0){
      params.do_checksums = true;
    } else if (strcmp(argv[i], "-nochecksum") == 0){
      params.do_checksums = false;
    } else if (strcmp(argv[i], "-inject") == 0){
      if (i+1 == argc){
        std::cerr << "Error: didn't specify failure file for injection\n";
        exit(EXIT_FAILURE);
      }
      params.do_injection = true;
      params.failfile = std::string(argv[++i]);
    } else if ((strcmp(argv[i], "-help") == 0)
               || (strcmp(argv[i], "-h") == 0)){
      std::cout << "Options:\n"
                << "  -ntiles <int>: Number of tiles in the domain"
                << " (default " << defaults.n_tiles << ")\n"
                << "  -tilesize <int>: Number of elements per tile"
                << " (default " << defaults.tilesize << ")\n"
                << "  -nsteps <int>: Number of time steps total"
                << " (default " << defaults.n_steps_total << ")\n"
                << "  -nstepstask <int>: Number of time steps within each task"
                << " (default " << defaults.n_steps_task << ")\n"
                << "  -[no]checksum: Whether or not to use checksums"
                << " (default "
                << (defaults.do_checksums ? "checksum" : "nochecksum") << ")\n"
                << "  -inject <filename>: Failure list for failure injection"
                << " (default disabled)\n";
      exit(EXIT_SUCCESS);
    } else {
      std::cerr << "Unexpected option: " << argv[i] << "; use -h for syntax\n";
      exit(EXIT_FAILURE);
    }
  }
}

int main (int argc, char *argv[])
{
  const char *deps[] = { "system" };
  hclib::launch(deps, 1, [&]() {
    params_t params;
    get_params(argc, argv, params);

    const double cfl = 0.5;

#if defined(USE_REPLAY)
    if (!params.do_checksums){
      fprintf(stderr, "Fatal error: Replay resilience requires checksums\n");
      exit(EXIT_FAILURE);
    }
#endif
    if (params.n_tiles < 1){
      fprintf(stderr, "Fatal error: ntiles must be positive\n");
      exit(EXIT_FAILURE);
    }
    if (params.tilesize < 1){
      fprintf(stderr, "Fatal error: tilesize must be positive\n");
      exit(EXIT_FAILURE);
    }
    // we can't time-advance a tile more times than there are ghost
    // values from each neighbor
    if ((params.n_steps_task < 1) || (params.n_steps_task > params.tilesize)){
      fprintf(stderr, "Fatal error: nstepstask must not be larger than "
              "tilesize or negative\n");
      exit(EXIT_FAILURE);
    }
    // since we're doing multiple iters per task, we ask that the total
    // number of timesteps be a multiple of the number of iters per task
    if (params.n_steps_total % params.n_steps_task != 0){
      fprintf(stderr, "Fatal error: nsteps must be multiple of nstepstask\n");
      exit(EXIT_FAILURE);
    }

    // domain size
    const int_t domainsize = params.n_tiles*params.tilesize;
    // each task advances n_steps_task time steps, so need fewer rounds
    const int_t n_rounds = params.n_steps_total/params.n_steps_task;

    int num_workers = hclib::get_num_workers();

    printf("num_workers=%d\n", num_workers);
    printf("n_tiles=%d\n", params.n_tiles);
    printf("tilesize=%d\n", params.tilesize);
    printf("domainsize=%d\n", domainsize);
    printf("n_steps_task=%d\n", params.n_steps_task);
    printf("n_steps_total=%d\n", params.n_steps_total);
    printf("n_rounds=%d\n", n_rounds);
    printf("checksums=%s\n", (params.do_checksums ? "yes" : "no"));
    printf("injection=%s\n", (params.do_injection ? "yes" : "no"));
#if defined(USE_REPLAY)
    printf("resilience=replay\n");
#elif defined(USE_REPLICATION)
    printf("resilience=replication\n");
#else
    printf("resilience=none\n");
#endif

    failure_generator_t failgen;
    if (params.do_injection){
      printf("failfile=%s\n", params.failfile.c_str());
      //printf("Reading failure list\n");

      if (params.failfile.size() == 0){
        fprintf(stderr, "Error: no failure file specified for error injection\n");
        exit(EXIT_FAILURE);
      }

      FILE *pFile = fopen(params.failfile.c_str(),"r");
      if (pFile == NULL){
        fprintf(stderr, "Error: unable to open %s\n", params.failfile.c_str());
        exit(EXIT_FAILURE);
      }

      unsigned long long n_failures;
      int ret = fscanf(pFile, "%llu\n", &n_failures);
      if (ret < 1){
        fprintf(stderr, "Error: failure file does not have the right format\n");
        exit(EXIT_FAILURE);
      }
      printf("n_failures: %llu\n", n_failures);
      for (unsigned long long ff=0; ff<n_failures; ++ff){
        unsigned long long ti;
#if defined(USE_REPLICATION)
        int ri;
        int ret = fscanf(pFile, "%llu %d\n", &ti, &ri);
        if (ret < 2){
          fprintf(stderr, "Error: failure file does not have the replication format "
                  "or ended prematurely\n");
          exit(EXIT_FAILURE);
        }
        //printf("planned failure: %llu %d\n", ti, ri);
        failgen.push_failure(ti, ri);
#else
        int ret = fscanf(pFile, "%llu\n", &ti);
        if (ret < 1){
          fprintf(stderr, "Error: failure file does not have the replay format "
                  "or ended prematurely\n");
          exit(EXIT_FAILURE);
        }
        //printf("planned failure: %llu\n", ti);
        failgen.push_failure(ti);
#endif
      }
      fclose(pFile);

      //printf("Read failure list\n");
    }

    //printf("Allocating tile array\n");

    Tile_t** tile_array = (Tile_t **) malloc(sizeof(Tile_t*)*(n_rounds+1)); 
    assert(tile_array);
    for (int_t tt=0; tt < n_rounds+1; ++tt) {
      tile_array[tt] = (Tile_t *) malloc(sizeof(Tile_t)*(params.n_tiles));
      assert(tile_array[tt]);
      for (int_t j=0; j < params.n_tiles; ++j) {
#if defined(DEBUG) || defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
        tile_array[tt][j].task_idx = 1+j + tt*params.n_tiles;
#endif
#if defined(USE_REPLAY)
        tile_array[tt][j].promise
            = new hclib::resilience::replay::promise_t<DoubleVec*>(tt < n_rounds ? 3 : 1);
#elif defined(USE_REPLICATION)
        tile_array[tt][j].promise
            = new hclib::resilience::diamond::promise_t<DoubleVec*>(tt < n_rounds ? 3 : 1);
#else
        tile_array[tt][j].promise
            = new hclib::ref_count::promise_t<DoubleVec*>(tt < n_rounds ? 3 : 1);
#endif
        assert(tile_array[tt][j].promise);
      }
    }

    //printf("Allocated tile array\n");
    //printf("Initialising sine wave in domain\n");

    hclib::finish([tile_array, &params, &domainsize]() {
      for (int_t j=0; j < params.n_tiles; ++j){
        initialize_tile(tile_array[0][j], j, params.tilesize,
                        domainsize, params.n_steps_task,
                        params.do_checksums);
      }
    });

    //printf("Done Initialising; Timestepping\n");

    struct timespec begin, end;
    clock_gettime(CLOCK_MONOTONIC, &begin);
#ifdef DUMP_TIMINGS
    std::cout << "BEGINN " << begin.tv_sec << " " << begin.tv_nsec << "\n";
#endif
	
    hclib::finish([tile_array, &params, &n_rounds, &cfl, &failgen]() {
      launch_iters(tile_array, params.n_tiles, n_rounds,
                   params.n_steps_task, cfl, params.do_checksums, failgen);
    });

    clock_gettime(CLOCK_MONOTONIC, &end);
#ifdef DUMP_TIMINGS
    std::cout << "ENDN " << end.tv_sec << " " << end.tv_nsec << "\n";
#endif
    double dur = ((end.tv_sec - begin.tv_sec)*1000000000
                 +(end.tv_nsec - begin.tv_nsec))*1.0/1000000000;
    printf("The computation took %f seconds\n", dur);

#ifdef DEBUG
    // print solution
    for (int_t j=0; j < params.n_tiles; ++j) {
      DoubleVec* soln =
         (DoubleVec*)tile_array[n_rounds][j].promise->get_future()->get();
      soln->atomic_print(j);
    }
#endif

    for (int_t j=0; j < params.n_tiles; ++j) {
      tile_array[n_rounds][j].promise->get_future()->release();
    }
    for (int_t tt=0; tt < n_rounds+1; ++tt) {
      free(tile_array[tt]);
    }
    free(tile_array);
  });

  return 0;
}
