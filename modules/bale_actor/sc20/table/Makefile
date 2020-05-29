
OSHCC ?=oshcc
OSHRUN ?=oshrun

TARGETS=histo_shmem histo_conveyor\
		ig_shmem ig_nbi_shmem ig_conveyor

all: $(TARGETS)

%_shmem: %_shmem.cpp
	$(OSHCC) -g -O3 -o $@ $^

histo_mpi3: histo_mpi3.cpp
	$(CXX) -g -O3 -DUSE_FETCH_AND_OP -o histo_mpi3.fa $^
	$(CXX) -g -O3 -DUSE_GET_ACCUMULATE -o histo_mpi3.gaccum $^
	$(CXX) -g -O3 -DUSE_ACCUMULATE -o histo_mpi3.accum $^

ig_mpi3: ig_mpi3.cpp
	$(CXX) -g -O3 -o ig_mpi3 $^

%_conveyor: %_conveyor.cpp
	$(OSHCC) -g -O3 -I$(BALE_INSTALL)/include -L$(BALE_INSTALL)/lib -o $@ $^ -lspmat -lconvey -lexstack -llibgetput

test:
	for target in $(TARGETS); do $(OSHRUN) ./$$target ; done
clean:
	rm -f $(TARGETS)
