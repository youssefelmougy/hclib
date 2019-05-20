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

// This is standard 3D stencil with six neighbors.  This code allows
// mixing resilience mechanisms.

#define USE_RESILIENT_PROMISE

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
//#define TIMING_INTERVAL 32

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

#define PI 3.14159265358979323846  /* pi */

enum TASK_STATE {NON_LEAF, LEAF};

typedef int int_t;

typedef struct {
  // default values defined below; can be overridden on command line
  int_t n_tiles = 4;  // in each dimension
  int_t tilesize = 8;  // in each dimension
  int_t n_steps = 128;
  bool do_checksums = true;
  bool do_injection = false;
  float replic_frac = 0.2;
  std::string failfile;
} params_t;

typedef struct {
  int_t x;
  int_t y;
  int_t z;
} tile_idx_t;

typedef enum {
  non_resilient,
  replay_resilient,
  replication_resilient
} resilience_type_t;

class iter_resilience_strategy {
 public:
  iter_resilience_strategy(int_t totaltiles, int_t nsteps,
                           float replic_frac)
  {
    printf("resilience=iter_resilience\n");
    float replay_frac = 1.0-replic_frac;
    printf("replay_frac=%g\n", replay_frac);
    printf("replay_start_iter=0\n");
    printf("replay_start_idx=1\n");
    printf("replic_frac=%g\n", replic_frac);
    int replic_start_iter = (int)(replay_frac*nsteps);
    printf("replic_start_iter=%d\n", replic_start_iter);
    replic_start_idx = replic_start_iter*totaltiles+1;
    printf("replic_start_idx=%d\n", replic_start_idx);
  }

  resilience_type_t task_strategy(int_t task_idx) const {
    // NOTE: constraints elsewhere in the code do not allow replay tasks to follow replication tasks
    if (task_idx < replic_start_idx){
      return replay_resilient;
    } else {
      return replication_resilient;
    }
  }

 private:
  int_t replic_start_idx;
};

typedef iter_resilience_strategy resilience_strategy_t;

class deterministic_mixed_failure_log {
 public:
  void push_failure(unsigned long long task_idx, int replica_idx){
    assert(replica_idx < 3);
    if (replica_idx < 0){
      int &cnt = failures_[task_idx];
      assert(cnt <= 0);
      --cnt;
    } else { 
      int &bits = failures_[task_idx];
      assert(bits >= 0);
      bits |= 2<<replica_idx;
    }
  }

  int should_fail(unsigned long long task_idx) const {
    auto it = failures_.find(task_idx);
    if (it == failures_.end()){
      return 0;
    } else {
      return it->second;
    }
  }

  static bool instance_should_fail(int fails, int index){
    if (fails > 0){
      // replication
      int mask = 2<<index;
      return ((fails & mask) == mask);
    } else {
      // replay
      return (index < -fails);
    }
  }

 private:
  std::unordered_map<unsigned long long,
                     int/*replica_bits or neg_replay_count*/> failures_;
};

typedef deterministic_mixed_failure_log failure_generator_t;

class Tile3D : public hclib::resilience::obj,
               public std::vector<double> {
public:
  virtual ~Tile3D() {}

  Tile3D(const Tile3D &other) = delete;

  Tile3D(int_t tilesize, bool do_checksums)
    : vector((tilesize+2)*(tilesize+2)*(tilesize+2)),
      tilesize_(tilesize), fullsize_(tilesize+2), chk_(0.0),
      chkxm_(0.0), chkxp_(0.0), chkym_(0.0), chkyp_(0.0), chkzm_(0.0),
      chkzp_(0.0), do_checksums_(do_checksums) {}

  Tile3D(const Tile3D *self, const Tile3D *xm, const Tile3D *xp,
         const Tile3D *ym, const Tile3D *yp, const Tile3D *zm,
         const Tile3D *zp, bool do_checksums)
    : vector((self->fullsize_)*(self->fullsize_)*(self->fullsize_)),
      tilesize_(self->tilesize_), fullsize_(self->fullsize_), chk_(0.0),
      chkxm_(0.0), chkxp_(0.0), chkym_(0.0), chkyp_(0.0), chkzm_(0.0),
      chkzp_(0.0), do_checksums_(do_checksums)
  {
    // compute assumed checksum at the end of construction
    if (do_checksums_){
      chk_ = self->chk_ + xm->chkxp_ + xp->chkxm_ + ym->chkyp_
           + yp->chkym_ + zm->chkzp_ + zp->chkzm_;
    }

    // copy main self
    for (int_t k=1; k<=tilesize_; ++k)
      copy_z_plane(self, k, k);

    // copy z- ghosts
    copy_z_plane(zm, tilesize_, 0);
    // copy z+ ghosts
    copy_z_plane(zp, 1, tilesize_+1);

    // copy y- ghosts
    copy_y_plane(ym, tilesize_, 0);
    // copy y+ ghosts
    copy_y_plane(yp, 1, tilesize_+1);

    // copy x- ghosts
    copy_x_plane(xm, tilesize_, 0);
    // copy x+ ghosts
    copy_x_plane(xp, 1, tilesize_+1);
  }

  bool apply_stencil(int fails){
    double checksum = chk_;

    if (do_checksums_){
      // subtract off fluxes on shrinking domain at x-faces of
      // time-advanced domain (x=0 and x=tilesize_+1 were included in
      // the input domain but will not be included in the
      // time-advanced domain)
      for (int_t k=1; k<=tilesize_; ++k){
        for (int_t j=1; j<=tilesize_; ++j){
          int_t xm = findex(0,j,k);
          int_t selfm = findex(1,j,k);
          int_t selfp = findex(tilesize_,j,k);
          int_t xp = findex(tilesize_+1,j,k);
          checksum -= left_flux((*this)[xm], (*this)[selfm], cfl)
                      + right_flux((*this)[selfp], (*this)[xp], cfl);
        }
      }

      // subtract off fluxes on shrinking domain at y-faces of
      // time-advanced domain (y=0 and y=tilesize_+1 were included in
      // the input domain but will not be included in the
      // time-advanced domain)
      for (int_t k=1; k<=tilesize_; ++k){
        for (int_t i=1; i<=tilesize_; ++i){
          int_t ym = findex(i,0,k);
          int_t selfm = findex(i,1,k);
          int_t selfp = findex(i,tilesize_,k);
          int_t yp = findex(i,tilesize_+1,k);
          checksum -= left_flux((*this)[ym], (*this)[selfm], cfl)
                      + right_flux((*this)[selfp], (*this)[yp], cfl);
        }
      }

      // subtract off fluxes on shrinking domain at z-faces of
      // time-advanced domain (z=0 and z=tilesize_+1 were included in
      // the input domain but will not be included in the
      // time-advanced domain)
      for (int_t j=1; j<=tilesize_; ++j){
        for (int_t i=1; i<=tilesize_; ++i){
          int_t zm = findex(i,j,0);
          int_t selfm = findex(i,j,1);
          int_t selfp = findex(i,j,tilesize_);
          int_t zp = findex(i,j,tilesize_+1);
          checksum -= left_flux((*this)[zm], (*this)[selfm], cfl)
                      + right_flux((*this)[selfp], (*this)[zp], cfl);
        }
      }

      // subtract off xy-edge values as well as xyz-corner values that
      // were included in input domain but not the time-advanced domain
      for (int_t k=0; k<fullsize_; ++k){
        int_t xmym = findex(0,0,k);
        int_t xmyp = findex(0,tilesize_+1,k);
        int_t xpym = findex(tilesize_+1,0,k);
        int_t xpyp = findex(tilesize_+1,tilesize_+1,k);
        checksum -= (*this)[xmym] + (*this)[xmyp]
                     + (*this)[xpym] + (*this)[xpyp];
      }

      // subtract off xz-edge values that were included in input domain
      // but not the time-advanced domain
      for (int_t j=1; j<=tilesize_; ++j){
        int_t xmzm = findex(0,j,0);
        int_t xmzp = findex(0,j,tilesize_+1);
        int_t xpzm = findex(tilesize_+1,j,0);
        int_t xpzp = findex(tilesize_+1,j,tilesize_+1);
        checksum -= (*this)[xmzm] + (*this)[xmzp]
                     + (*this)[xpzm] + (*this)[xpzp];
      }

      // subtract off yz-edge values that were included in input domain
      // but not the time-advanced domain
      for (int_t i=1; i<=tilesize_; ++i){
        int_t ymzm = findex(i,0,0);
        int_t ymzp = findex(i,0,tilesize_+1);
        int_t ypzm = findex(i,tilesize_+1,0);
        int_t ypzp = findex(i,tilesize_+1,tilesize_+1);
        checksum -= (*this)[ymzm] + (*this)[ymzp]
                     + (*this)[ypzm] + (*this)[ypzp];
      }
    }

    // apply stencil and left-justify the result so we can operate in-place
    for (int_t k=0; k<tilesize_; ++k){
      for (int_t j=0; j<tilesize_; ++j){
        for (int_t i=0; i<tilesize_; ++i){
          int_t self = findex(i+1,j+1,k+1);
          int_t xm = findex(i,j+1,k+1);
          int_t xp = findex(i+2,j+1,k+1);
          int_t ym = findex(i+1,j,k+1);
          int_t yp = findex(i+1,j+2,k+1);
          int_t zm = findex(i+1,j+1,k);
          int_t zp = findex(i+1,j+1,k+2);
          (*this)[findex(i,j,k)]=stencil((*this)[self], (*this)[xm],
                                         (*this)[xp], (*this)[ym],
                                         (*this)[yp], (*this)[zm],
                                         (*this)[zp], cfl);
        }
      }
    }

    // shift the result back to its proper position
    for (int_t k=tilesize_; k>0; --k){
      for (int_t j=tilesize_; j>0; --j){
        for (int_t i=tilesize_; i>0; --i){
          (*this)[findex(i,j,k)]=(*this)[findex(i-1,j-1,k-1)];
        }
      }
    }

    // zero out stale ghost values so they won't contribute to checksum
    for (int_t k=0; k<fullsize_; ++k){
      for (int_t j=0; j<fullsize_; ++j){
        (*this)[findex(0,j,k)]=0;
        (*this)[findex(tilesize_+1,j,k)]=0;
      }
      for (int_t i=0; i<fullsize_; ++i){
        (*this)[findex(i,0,k)]=0;
        (*this)[findex(i,tilesize_+1,k)]=0;
      }
    }
    for (int_t j=0; j<fullsize_; ++j){
      for (int_t i=0; i<fullsize_; ++i){
        (*this)[findex(i,j,0)]=0;
        (*this)[findex(i,j,tilesize_+1)]=0;
      }
    }

    // inject failure, if appropriate
    if (fails){
      int idx = hclib::resilience::get_index();
      if (failure_generator_t::instance_should_fail(fails, idx))
        inject_failure(idx);
    }

    if (do_checksums_){
      // compute alternative checksum and verify
      chk_ = std::accumulate(begin(), end(), 0.0);
      bool pass = chk_within_tol(checksum);

      // if checksum matches, compute sub-checksums
      if (pass){
        chkzm_ = zchecksum(1);
        chkzp_ = zchecksum(tilesize_);
        chkym_ = ychecksum(1);
        chkyp_ = ychecksum(tilesize_);
        chkxm_ = xchecksum(1);
        chkxp_ = xchecksum(tilesize_);
      } else {
#if !defined(DUMP_TIMINGS) && !defined(DUMP_DURATION)
        //fprintf(stderr,"checksums do not match: %g %g %g\n",
        //        checksum, chk_, checksum-chk_);
#endif
      }

      // return success/fail according to checksums
      return pass;
    } else {
      // no checksums, so assume success
      return true;
    }
  }

  bool compute_checksums(double indep_checksum){
    // used after construction of initial conditions only
    // (checksums happen automatically after applying stencil)
    if (do_checksums_){
      chk_ = std::accumulate(begin(), end(), 0.0);
      bool pass = chk_within_tol(indep_checksum);

      // if checksum matches, compute sub-checksums
      if (pass){
        chkzm_ = zchecksum(1);
        chkzp_ = zchecksum(tilesize_);
        chkym_ = ychecksum(1);
        chkyp_ = ychecksum(tilesize_);
        chkxm_ = xchecksum(1);
        chkxp_ = xchecksum(tilesize_);
      } else {
#if !defined(DUMP_TIMINGS) && !defined(DUMP_DURATION)
        //fprintf(stderr,"checksums do not match: %g %g %g\n",
        //        indep_checksum, chk_, indep_checksum-chk_);
#endif
      }

      // return success/fail according to checksums
      return pass;
    } else {
      // no checksums, so assume success
      return true;
    }
  }

  bool equals(obj* obj2){
    Tile3D *o2 = (Tile3D*)obj2;
    return *this == *o2;
  }

#ifdef DEBUG
  void atomic_print(unsigned long long task_id, bool include_ghost) const {
    std::stringstream ss;
    ss << task_id << ":\n";
    if (include_ghost){
      for (int_t kz=0; kz<fullsize_; ++kz){
        for (int_t ky=0; ky<fullsize_; ++ky){
          for (int_t kx=0; kx<fullsize_; ++kx){
            int_t k = findex(kx,ky,kz);
            ss << (*this)[k] << " ";
          }
          ss << "\n";
        }
        ss << "\n";
      }
    } else {
      for (int_t kz=1; kz<=tilesize_; ++kz){
        for (int_t ky=1; ky<=tilesize_; ++ky){
          for (int_t kx=1; kx<=tilesize_; ++kx){
            int_t k = findex(kx,ky,kz);
            ss << (*this)[k] << " ";
          }
          ss << "\n";
        }
        ss << "\n";
      }
    }
    std::cout << ss.str();
  }
#endif

  double chk() const { return chk_; }

  int_t tindex(int_t kx, int_t ky, int_t kz) const {
    return (kx+1) + fullsize_*((ky+1) + fullsize_*(kz+1));
  }

private:
  int_t findex(int_t kx, int_t ky, int_t kz) const {
    return kx + fullsize_*(ky + fullsize_*kz);
  }

  void copy_z_plane(const std::vector<double> *source,
                    int_t src_kz, int_t dest_kz){ 
    for (int_t ky=1; ky<=tilesize_; ++ky){
      auto where = begin() + findex(1, ky, dest_kz);
      auto start = source->begin() + findex(1, ky, src_kz);
      auto end = start + tilesize_;
      std::copy(start, end, where);
    }
  }

  void copy_y_plane(const std::vector<double> *source,
                    int_t src_ky, int_t dest_ky){
    for (int_t kz=1; kz<=tilesize_; ++kz){
      auto where = begin() + findex(1, dest_ky, kz);
      auto start = source->begin() + findex(1, src_ky, kz);
      auto end = start + tilesize_;
      std::copy(start, end, where);
    }
  }

  void copy_x_plane(const std::vector<double> *source,
                    int_t src_kx, int_t dest_kx){
    for (int_t kz=1; kz<=tilesize_; ++kz){
      for (int_t ky=1; ky<=tilesize_; ++ky){
        auto where = begin() + findex(dest_kx, ky, kz);
        auto start = source->begin() + findex(src_kx, ky, kz);
        *where = *start;
      }
    }
  }

  double zchecksum(int_t kz) const {
    double checksum = 0.0;
    for (int_t ky=1; ky<=tilesize_; ++ky){
      auto start = begin() + findex(1, ky, kz);
      auto end = start + tilesize_;
      checksum = std::accumulate(start, end, checksum);
    }
    return checksum;
  }

  double ychecksum(int_t ky) const {
    double checksum = 0.0;
    for (int_t kz=1; kz<=tilesize_; ++kz){
      auto start = begin() + findex(1, ky, kz);
      auto end = start + tilesize_;
      checksum = std::accumulate(start, end, checksum);
    }
    return checksum;
  }

  double xchecksum(int_t kx) const {
    double checksum = 0.0;
    for (int_t kz=1; kz<=tilesize_; ++kz){
      for (int_t ky=1; ky<=tilesize_; ++ky){
        auto start = begin() + findex(kx, ky, kz);
        checksum += *start;
      }
    }
    return checksum;
  }

  bool chk_within_tol(double independent_chk) const {
    return (std::abs(chk_ - independent_chk) <= checksum_tol);
  }

  // 7-point stencil for 3D heat equation, cfl=D*dt/(dx*dx) 
  double stencil(double self, double xm, double xp, double ym,
                 double yp, double zm, double zp, double cfl) const {
    return (1.0-6.0*cfl)*self + cfl*(xm + xp + ym + yp + zm + zp);
  }

  // flux near negative edge of domain that decreases in size by 1
  double left_flux(double left, double self, double cfl) const {
    return (1.0-cfl)*left + cfl*self;
  }

  // flux near positive edge of domain that decreases in size by 1
  double right_flux(double self, double right, double cfl) const {
    return (1.0-cfl)*right + cfl*self;
  }

  void inject_failure(int idx){
    //printf("injecting failure\n");
    int_t failidx = tindex(tilesize_/2,tilesize_/2,tilesize_/2);
    double *val = &((*this)[failidx]);
    uint64_t tmp;
    std::memcpy(&tmp, val, sizeof(double));
    tmp ^= 2ULL<<(60-idx);
    std::memcpy(val, &tmp, sizeof(double));
  }

  int_t tilesize_;
  int_t fullsize_;
  double chk_;
  double chkxm_, chkxp_;
  double chkym_, chkyp_;
  double chkzm_, chkzp_;
  bool do_checksums_;

  static constexpr double checksum_tol = 1e-5;
  static constexpr double cfl = 0.5;
};

typedef struct {
  hclib::resilience::promise_t<Tile3D*>* promise;
  unsigned long long task_idx;
} Tile_t;

int check(void *args)
{
  bool succeeded = *(bool*)args;
  return (succeeded ? 1 : 0);
}

int_t index_3d(int_t kx, int_t ky, int_t kz, int_t n_tiles){
  return kx + n_tiles*(ky + n_tiles*kz);
}

double initial_conditions(int_t mx, int_t my, int_t mz, int_t domainsize)
{
  // TODO: make this something appropriate for 3D heat equation;
  // needs to be periodic due to neighbor configuration chosen
  double m = (double)mx / domainsize;
  return std::sin(2.0 * PI * m);
}

void initialize_tile(Tile_t &tile, tile_idx_t tile_idx, int_t tilesize,
                     int_t domainsize, bool do_checksums)
{
  bool *succeeded = new bool;
  hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();
  hclib::resilience::replay::async_await_check<LEAF,3>([=]{
    double checksum = 0.0;
    // compute initial conditions over tile
    Tile3D *workspace = new Tile3D(tilesize, do_checksums);
    assert(workspace);
    for (int_t kz=0; kz<tilesize; ++kz){
      int_t mz = tile_idx.z * tilesize + kz;
      for (int_t ky=0; ky<tilesize; ++ky){
        int_t my = tile_idx.y * tilesize + ky;
        for (int_t kx=0; kx<tilesize; ++kx){
          int_t mx = tile_idx.x * tilesize + kx;
          int_t k = workspace->tindex(kx,ky,kz);
          (*workspace)[k] = initial_conditions(mx,my,mz,domainsize);
          if (do_checksums)
            checksum += (*workspace)[k];
        }
      }
    }

    bool pass = workspace->compute_checksums(checksum);
    *succeeded = pass;

    // put initial conditions
    tile.promise->put(workspace);
  }, prom_res, check, (void*)(succeeded), (hclib_future_t*)NULL);

  hclib::async_await([=]{
    int res = prom_res->get_future()->get();
    if (res == 0){
      fprintf(stderr,"Fatal error: resilience failed to save you\n");
      exit(EXIT_FAILURE);
    }
    delete succeeded;
  }, prom_res->get_future());
} 

inline
bool advance_tile_body(const Tile_t &self_prev, const Tile_t &xm_prev,
                       const Tile_t &xp_prev, const Tile_t &ym_prev,
                       const Tile_t &yp_prev, const Tile_t &zm_prev,
                       const Tile_t &zp_prev, const Tile_t &self_updated,
                       double cfl, bool do_checksums, int fails
#if defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
                       ,bool do_timing
#endif
                       )
{
#if defined(DUMP_DURATION)
  struct timespec begin;
  if (do_timing)
    clock_gettime(CLOCK_MONOTONIC, &begin);
#endif

  // get previous state
  Tile3D* self = (Tile3D*)self_prev.promise->get_future()->get();
  Tile3D* xm = (Tile3D*)xm_prev.promise->get_future()->get();
  Tile3D* xp = (Tile3D*)xp_prev.promise->get_future()->get();
  Tile3D* ym = (Tile3D*)ym_prev.promise->get_future()->get();
  Tile3D* yp = (Tile3D*)yp_prev.promise->get_future()->get();
  Tile3D* zm = (Tile3D*)zm_prev.promise->get_future()->get();
  Tile3D* zp = (Tile3D*)zp_prev.promise->get_future()->get();

  // construct a new tile and copy in previous state, including ghosts
  Tile3D *workspace = new Tile3D(self, xm, xp, ym, yp, zm, zp,
                                 do_checksums);
  assert(workspace);

  // do the time-advance
  bool pass = workspace->apply_stencil(fails);

  // put the advanced values
  self_updated.promise->put(workspace);

#if defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
  if (do_timing){
    struct timespec end;
    clock_gettime(CLOCK_MONOTONIC, &end);
#if defined(DUMP_DURATION)
    double dur = ((end.tv_sec - begin.tv_sec)*1000000000
                 +(end.tv_nsec - begin.tv_nsec))*1.0/1000000000;
    printf("TASKDUR %llu %g\n", self_prev.task_idx, dur);
#else
    printf("TASKN %llu %ld %ld\n", self_prev.task_idx, end.tv_sec, end.tv_nsec);
#endif
  }
#endif

  return pass;
}

void advance_tile(Tile_t &self_prev, Tile_t &xm_prev, Tile_t &xp_prev,
                  Tile_t &ym_prev, Tile_t &yp_prev, Tile_t &zm_prev,
                  Tile_t &zp_prev, Tile_t &self_updated,
                  double cfl, bool do_checksums, int fails,
#if defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
                  bool do_timing,
#endif
                  resilience_type_t task_strategy)
{
  auto futures = new std::vector<hclib_future_t*>(7);
  (*futures)[0] = self_prev.promise->get_future();
  (*futures)[1] = xm_prev.promise->get_future();
  (*futures)[2] = xp_prev.promise->get_future();
  (*futures)[3] = ym_prev.promise->get_future();
  (*futures)[4] = yp_prev.promise->get_future();
  (*futures)[5] = zm_prev.promise->get_future();
  (*futures)[6] = zp_prev.promise->get_future();

  bool *succeeded = new bool;
  hclib::promise_t<int>* prom_res = new hclib::promise_t<int>();
  
  if (task_strategy == replay_resilient){ /* use replay */
    hclib::resilience::replay::async_await_check<LEAF,3>([=]{
      bool pass = advance_tile_body(self_prev, xm_prev, xp_prev,
                                    ym_prev, yp_prev, zm_prev, zp_prev,
                                    self_updated, cfl, do_checksums,
                                    fails
#if defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
                                   ,do_timing
#endif
                                   );
      *succeeded = pass;
    }, prom_res, check, (void*)(succeeded), futures);
  } else if (task_strategy == replication_resilient){ /* use replication */
    hclib::resilience::diamond::async_await_check<LEAF>([=]{
      advance_tile_body(self_prev, xm_prev, xp_prev, ym_prev, yp_prev,
                        zm_prev, zp_prev, self_updated, cfl,
                        false, fails  /* forces no checksums for replication! */
#if defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
                       ,do_timing
#endif
                       );
    }, prom_res, futures);
  } else {
    fprintf(stderr,"Fatal error: Task resilience strategy not implemented in advance_tile\n");
    exit(EXIT_FAILURE);
  }

  hclib::async_await([=]{
    int res = prom_res->get_future()->get();
    if (res == 0){
      fprintf(stderr,"Fatal error: resilience failed to save you\n");
      exit(EXIT_FAILURE);
    }
    delete succeeded;
    delete futures;
  }, prom_res->get_future());
} 

void launch_iters(Tile_t **tile_array, int_t n_tiles, int_t n_steps,
                  double cfl, bool do_checksums,
                  const failure_generator_t &failgen,
                  resilience_strategy_t &resil_strategy)
{
  // one task per tile per time step
  int_t totaltiles = n_tiles*n_tiles*n_tiles;
  for (int_t tt=0; tt < n_steps; ++tt){
    for (int_t jz=0; jz<n_tiles; ++jz){
      for (int_t jy=0; jy<n_tiles; ++jy){
        for (int_t jx=0; jx<n_tiles; ++jx){
          int_t j = index_3d(jx,jy,jz,n_tiles);
          int fails = failgen.should_fail(1+j + tt*totaltiles);

          // compute neighbor coordinates
          int_t xm = (n_tiles+jx-1)%n_tiles;
          int_t xp = (jx+1)%n_tiles;
          int_t ym = (n_tiles+jy-1)%n_tiles;
          int_t yp = (jy+1)%n_tiles;
          int_t zm = (n_tiles+jz-1)%n_tiles;
          int_t zp = (jz+1)%n_tiles;

          // compute neighbor ids
          int_t jxm = index_3d(xm,jy,jz,n_tiles);
          int_t jxp = index_3d(xp,jy,jz,n_tiles);
          int_t jym = index_3d(jx,ym,jz,n_tiles);
          int_t jyp = index_3d(jx,yp,jz,n_tiles);
          int_t jzm = index_3d(jx,jy,zm,n_tiles);
          int_t jzp = index_3d(jx,jy,zp,n_tiles);

          advance_tile(tile_array[tt][j], tile_array[tt][jxm],
                       tile_array[tt][jxp], tile_array[tt][jym],
                       tile_array[tt][jyp], tile_array[tt][jzm],
                       tile_array[tt][jzp], tile_array[tt+1][j],
                       cfl, do_checksums, fails,
#if defined(DUMP_TIMINGS) || defined(DUMP_DURATION)
#ifdef DUMP_FAIL_DURATION
                       fails,
#else
                       (tt+1)%TIMING_INTERVAL == 0,
#endif
#endif
                       resil_strategy.task_strategy(tile_array[tt][j].task_idx));
        }
      }
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
      params.n_steps = atoi(argv[++i]);
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
    } else if (strcmp(argv[i], "-replic") == 0){
      if (i+1 == argc){
        std::cerr << "Error: didn't specify value for replic\n";
        exit(EXIT_FAILURE);
      }
      params.replic_frac = atoi(argv[++i])/100.0;
    } else if ((strcmp(argv[i], "-help") == 0)
               || (strcmp(argv[i], "-h") == 0)){
      std::cout << "Options:\n"
                << "  -ntiles <int>: Number of tiles in the domain"
                << " (default " << defaults.n_tiles << ")\n"
                << "  -tilesize <int>: Number of elements per tile in each dimension"
                << " (default " << defaults.tilesize << ")\n"
                << "  -nsteps <int>: Number of time steps total"
                << " (default " << defaults.n_steps << ")\n"
                << "  -replic <int>: Replicate last percentage of tasks"
                << " (default " << (int)(defaults.replic_frac*100) << ")\n"
                << "  -[no]checksum: Whether or not to use checksums"
                << " (default "
                << (defaults.do_checksums ? "checksum" : "nochecksum") << ")\n"
                << "  -inject <filename>: Failure list for failure injection\n"
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

    // must be less than 1/6 for stability
    const double cfl = 0.16;

    if (!params.do_checksums){
      fprintf(stderr, "Fatal error: Replay and mixed resilience require checksums\n");
      exit(EXIT_FAILURE);
    }
    if (params.n_tiles < 1){
      fprintf(stderr, "Fatal error: ntiles must be positive\n");
      exit(EXIT_FAILURE);
    }
    if (params.tilesize < 1){
      fprintf(stderr, "Fatal error: tilesize must be positive\n");
      exit(EXIT_FAILURE);
    }

    // domain size
    const int_t domainsize = params.n_tiles*params.tilesize;
    const int_t totaltiles = params.n_tiles*params.n_tiles*params.n_tiles;
    int num_workers = hclib::get_num_workers();

    printf("num_workers=%d\n", num_workers);
    printf("n_tiles=%d\n", params.n_tiles);
    printf("totaltiles=%d\n", totaltiles);
    printf("tilesize=%d\n", params.tilesize);
    printf("domainsize=%d\n", domainsize);
    printf("n_steps=%d\n", params.n_steps);
    printf("checksums=%s\n", (params.do_checksums ? "yes" : "no"));
    printf("injection=%s\n", (params.do_injection ? "yes" : "no"));

    resilience_strategy_t strategy(totaltiles, params.n_steps, params.replic_frac);

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
        int ri;
        int ret = fscanf(pFile, "%llu %d\n", &ti, &ri);
        if (ret < 2){
          fprintf(stderr, "Error: failure file does not have mixed resilience format "
                  "or ended prematurely\n");
          exit(EXIT_FAILURE);
        }
        //printf("planned failure: %llu %d\n", ti, ri);
        failgen.push_failure(ti, ri);
      }
      fclose(pFile);

      //printf("Read failure list\n");
    }

    //printf("Allocating tile array\n");

    Tile_t** tile_array = (Tile_t **) malloc(sizeof(Tile_t*)*(params.n_steps+1));
    assert(tile_array);
    for (int_t tt=0; tt < params.n_steps+1; ++tt) {
      tile_array[tt] = (Tile_t *) malloc(sizeof(Tile_t)*(totaltiles));
      assert(tile_array[tt]);
      for (int_t j=0; j < totaltiles; ++j) {
        tile_array[tt][j].task_idx = 1+j + tt*totaltiles;
        tile_array[tt][j].promise
            = new hclib::resilience::promise_t<Tile3D*>(
                tt < params.n_steps ? 7 : 1);
        assert(tile_array[tt][j].promise);
      }
    }

    //printf("Allocated tile array\n");
    //printf("Initialising sine wave in domain\n");

    hclib::finish([tile_array, &params, &domainsize]() {
      tile_idx_t tile_idx;
      for (int_t jz=0; jz < params.n_tiles; ++jz){
        tile_idx.z = jz;
        for (int_t jy=0; jy < params.n_tiles; ++jy){
          tile_idx.y = jy;
          for (int_t jx=0; jx < params.n_tiles; ++jx){
            tile_idx.x = jx;
            int_t j = index_3d(jx,jy,jz,params.n_tiles);
            initialize_tile(tile_array[0][j], tile_idx, params.tilesize,
                            domainsize, params.do_checksums);
          }
        }
      }
    });

    //printf("Done Initialising; Timestepping\n");

    struct timespec begin, end;
    clock_gettime(CLOCK_MONOTONIC, &begin);
#ifdef DUMP_TIMINGS
    std::cout << "BEGINN " << begin.tv_sec << " " << begin.tv_nsec << "\n";
#endif
	
    hclib::finish([tile_array, &params, &cfl, &failgen, &strategy]() {
      launch_iters(tile_array, params.n_tiles, params.n_steps,
                   cfl, params.do_checksums, failgen, strategy);
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
    for (int_t j=0; j < totaltiles; ++j) {
      Tile3D* soln =
         (Tile3D*)tile_array[params.n_steps][j].promise->get_future()->get();
      soln->atomic_print(j, false);
    }
#endif

    for (int_t j=0; j < totaltiles; ++j) {
      tile_array[params.n_steps][j].promise->get_future()->release();
    }
    for (int_t tt=0; tt < params.n_steps+1; ++tt) {
      free(tile_array[tt]);
    }
    free(tile_array);
  });

  return 0;
}
