FROM sagemath/sagemath-jupyter:7.4

USER root

RUN sed -i -e 's/archive.ubuntu.com\|security.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list
RUN apt-get update && apt-get install -y sdpa libsdpa-dev build-essential make libopenblas-dev libblas-dev liblapack-dev
RUN git clone https://github.com/scibuilder/QD.git && cd QD && ./configure --prefix=/usr/local/qd && make && make install && cd / && rm -rf QD
RUN git clone https://github.com/ericphanson/sdpa-dd.git && cd sdpa-dd && ./configure --with-qd-includedir=/usr/local/qd/include --with-qd-libdir=/usr/local/qd/lib && make && cp sdpa_dd /usr/bin && cd / && rm -rf sdpa_dd
RUN git clone https://github.com/ericphanson/sdpa-qd.git && cd sdpa-qd && ./configure --with-qd-includedir=/usr/local/qd/include --with-qd-libdir=/usr/local/qd/lib && make && cp sdpa_qd /usr/bin && cd / && rm -rf sdpa_qd
RUN git clone https://github.com/FordUniver/Csdp.git && cd Csdp && git checkout docker && make && make install && cd / && rm -rf Csdp

USER sage

# RUN sage -i csdp

RUN sage --pip install tqdm

# RUN git clone https://github.com/jsliacan/flagmatic && cd flagmatic/pkg && sage -python setup.py install && cd ../.. && rm -rf flagmatic
RUN git clone https://github.com/FordUniver/flagmatic.git && cd flagmatic/pkg && sage -python setup.py install && cd ../.. && rm -rf flagmatic

ENTRYPOINT /bin/bash