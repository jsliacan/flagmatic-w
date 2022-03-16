FROM sagemath/sagemath:7.4
# ubuntu:latest
# ENV DEBIAN_FRONTEND=noninteractive
# ENV TZ=Europe/Berlin

USER root

RUN sed -i -e 's/archive.ubuntu.com\|security.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list
RUN apt-get update && apt-get install -y build-essential sdpa gfortran libblas-dev liblapack-dev make cmake git
# RUN git clone https://github.com/xianyi/OpenBLAS.git && cd OpenBLAS && make --jobs=$(nproc) && mkdir build && cd build && cmake --parallel .. && make --jobs=$(nproc) PREFIX=/opt install && cd / && rm -rf OpenBlas
RUN git clone https://github.com/FordUniver/Csdp.git && cd Csdp && git checkout docker3 && make --jobs=$(nproc) && make --jobs=$(nproc) install && cd / && rm -rf Csdp
RUN git clone https://github.com/scibuilder/QD.git && cd QD && ./configure --prefix=/usr/local/qd && make --jobs=$(nproc) && make --jobs=$(nproc) install && cd / && rm -rf QD
RUN git clone https://github.com/ericphanson/sdpa-dd.git && cd sdpa-dd && ./configure --with-qd-includedir=/usr/local/qd/include --with-qd-libdir=/usr/local/qd/lib && make --jobs=$(nproc) && cp sdpa_dd /usr/bin && cd / && rm -rf sdpa_dd
RUN git clone https://github.com/ericphanson/sdpa-qd.git && cd sdpa-qd && ./configure --with-qd-includedir=/usr/local/qd/include --with-qd-libdir=/usr/local/qd/lib && make --jobs=$(nproc) && cp sdpa_qd /usr/bin && cd / && rm -rf sdpa_qd

USER sage

# RUN sage -i csdp

RUN sage --pip install tqdm

# RUN git clone https://github.com/jsliacan/flagmatic && cd flagmatic/pkg && sage -python setup.py install && cd ../.. && rm -rf flagmatic
RUN git clone https://github.com/FordUniver/flagmatic.git && cd flagmatic/pkg && sage -python setup.py install && cd ../.. && rm -rf flagmatic

ENTRYPOINT /bin/bash