
OSHCC ?=oshcc
UPCC ?=upcc
MPICC ?=mpicc
OSHRUN ?=oshrun
SRUN ?=srun

TARGETS=histo_shmem histo_conveyor histo_upc histo_mpi3\
		ig_shmem ig_nbi_shmem ig_conveyor ig_upc ig_mpi3

all: $(TARGETS)

%_shmem: %_shmem.cpp
	$(OSHCC) -g -O3 -o $@ $^

histo_mpi3: histo_mpi3.cpp
	$(MPICC) -g -O3 -DUSE_FETCH_AND_OP -o histo_mpi3.fa $^
	$(MPICC) -g -O3 -DUSE_GET_ACCUMULATE -o histo_mpi3.gaccum $^
	$(MPICC) -g -O3 -DUSE_ACCUMULATE -o histo_mpi3.accum $^

ig_mpi3: ig_mpi3.cpp
	$(MPICC) -g -O3 -o ig_mpi3 $^

%_conveyor: %_conveyor.cpp
	$(OSHCC) -g -O3 -I$(BALE_INSTALL)/include -L$(BALE_INSTALL)/lib -o $@ $^ -lspmat -lconvey -lexstack -llibgetput

%_upc: %_upc.c
	$(UPCC) -O -o $@ $^

test:
	for target in $(TARGETS); do $(SRUN) ./$$target ; done

clean:
	rm -f $(TARGETS)

