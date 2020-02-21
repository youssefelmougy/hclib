FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    apt-utils \
    autoconf \ 
    git \
    libboost-all-dev \
    libtool \
    pkg-config \
    vim 

ENV LOCAL=/root/local

WORKDIR /root
RUN git clone https://github.com/ofiwg/libfabric.git
WORKDIR /root/libfabric
RUN ./autogen.sh && ./configure --prefix=$LOCAL && make install

WORKDIR /root
RUN git clone https://github.com/Sandia-OpenSHMEM/SOS.git
WORKDIR /root/SOS
RUN ./autogen.sh && ./configure --prefix=$LOCAL --with-ofi=$LOCAL --enable-pmi-simple && make install

ENV CC=/root/local/bin/oshcc
ENV CXX=/root/local/bin/oshc++
ENV UPC=/root/local/bin/oshcc

WORKDIR /root
RUN git clone https://github.com/jdevinney/bale.git
WORKDIR /root/bale
RUN chmod +x ./install.sh
RUN ./install.sh -s -f

ENV BALE_INSTALL=/root/bale/build_unknown

WORKDIR /root
RUN git clone https://github.com/srirajpaul/hclib
WORKDIR /root/hclib
RUN git fetch && git checkout bale_actor
RUN ./install.sh
ENV HCLIB_ROOT=/root/hclib/hclib-install
RUN cd modules/bale_actor && make

WORKDIR /root/hclib/modules/bale_actor/test

ENV LD_LIBRARY_PATH=$LOCAL/lib:$BALE_INSTALL/lib:$HCLIB_ROOT/lib:$HCLIB_ROOT/../modules/bale_actor/lib