PROJECT_CXXFLAGS=-O3 -std=c++11 -I$(HCLIB_ROOT)/include -I$(LIBXML2_INCLUDE)
PROJECT_LDFLAGS=-O3 -L$(LIBXML2_LIBS) -L$(HCLIB_ROOT)/lib
PROJECT_LDLIBS=-lhclib -lxml2
ifdef TBB_MALLOC
  PROJECT_LDFLAGS+=-L$(TBB_MALLOC)
  PROJECT_LDLIBS+=-ltbbmalloc_proxy
endif
