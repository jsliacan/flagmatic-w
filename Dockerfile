FROM sagemath/sagemath:9.5
# ubuntu:latest
# ENV DEBIAN_FRONTEND=noninteractive
# ENV TZ=Europe/Berlin

USER root

RUN apt-get update && apt-get install -y build-essential sdpa gfortran make cmake git

# Compile OpenBLAS and CSDP form scratch for best performance
RUN git clone https://github.com/xianyi/OpenBLAS.git && cd OpenBLAS && make --jobs=$(nproc) && mkdir build && cd build && cmake --parallel .. && make --jobs=$(nproc) install && cd ../.. && rm -rf OpenBlas
RUN git clone https://github.com/FordUniver/Csdp.git && cd Csdp && git checkout openblas && make --jobs=$(nproc) && make --jobs=$(nproc) install && cd .. && rm -rf Csdp

# Double precision variant of SDP is causing issues currently
#RUN git clone https://github.com/scibuilder/QD.git && cd QD && ./configure --prefix=/usr/local/qd && make --jobs=$(nproc) && make --jobs=$(nproc) install && cd / && rm -rf QD
#RUN git clone https://github.com/ericphanson/sdpa-dd.git && cd sdpa-dd && ./configure --with-qd-includedir=/usr/local/qd/include --with-qd-libdir=/usr/local/qd/lib && make --jobs=$(nproc) && cp sdpa_dd /usr/bin && cd / && rm -rf sdpa_dd
#RUN git clone https://github.com/ericphanson/sdpa-qd.git && cd sdpa-qd && ./configure --with-qd-includedir=/usr/local/qd/include --with-qd-libdir=/usr/local/qd/lib && make --jobs=$(nproc) && cp sdpa_qd /usr/bin && cd / && rm -rf sdpa_qd

USER sage

RUN sage --pip install tqdm

# ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN git clone https://github.com/FordUniver/flagmatic.git && cd flagmatic && cd pkg && sage -python setup.py install && cd ../.. && rm -rf flagmatic

ENTRYPOINT /bin/bash
