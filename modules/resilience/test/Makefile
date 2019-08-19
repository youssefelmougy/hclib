include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.pre.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.pre.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/system/inc/hclib_system.post.mak
include $(HCLIB_ROOT)/../modules/ref-count/inc/hclib_ref-count.post.mak
include $(HCLIB_ROOT)/../modules/resilience/inc/hclib_resilience.post.mak

MPIRUN=mpirun

TARGETS=resilience-d0 resilience-d1 resilience-r0 resilience-r1 resilience-a0 resilience-mix resilience-c0
MPI_TARGETS=resilience-r0-mpi resilience-d0-mpi resilience-r0-allreduce-mpi

all: $(TARGETS) $(MPI_TARGETS)

%: %.cpp
	$(CXX) -g -std=c++14 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS)

%-mpi: %-mpi.cpp
	$(MPICXX) -g -std=c++14 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -lhclib_resilience

test:
	for target in $(TARGETS); do ./$$target ; done

test_mpi:
	for target in $(MPI_TARGETS); do $(MPIRUN) -n 2 ./$$target ; done

clean:
	rm -f $(TARGETS) $(MPI_TARGETS)
