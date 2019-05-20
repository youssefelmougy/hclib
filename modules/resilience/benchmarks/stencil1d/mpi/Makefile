include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.pre.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.pre.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.post.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.post.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.post.mak

TARGETS=lw_1d_full_tile_mpi_nonres lw_1d_full_tile_mpi_replay lw_1d_full_tile_mpi_replic

FLAGS= -g -std=c++14 -Wall -O3
#MPICXX=CC $(HCLIB_CXXFLAGS)
MPICXX=mpicxx
all: $(TARGETS) targets.txt

.PHONY: targets.txt # always update
targets.txt:
	@echo "$(TARGETS)" > $@

%_mpi_nonres: %.cpp
	$(MPICXX) ${FLAGS} -DUSE_MPI $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

%_mpi_replic: %.cpp
	$(MPICXX) ${FLAGS} -DUSE_MPI -DUSE_REPLICATION $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

%_mpi_replay: %.cpp
	$(MPICXX) ${FLAGS} -DUSE_MPI -DUSE_REPLAY $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

clean:
	rm -f $(TARGETS)

