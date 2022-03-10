FROM sagemath/sagemath-jupyter:7.4

USER root

RUN sed -i -e 's/archive.ubuntu.com\|security.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list
RUN apt-get update && apt-get install -y sdpa libsdpa-dev 
RUN git clone https://github.com/scibuilder/QD.git && cd QD && ./configure --prefix=/usr/local/qd && make && make install && cd / && rm -rf QD
RUN git clone https://github.com/ericphanson/sdpa-dd.git && cd sdpa-dd && ./configure --with-qd-includedir=/usr/local/qd/include --with-qd-libdir=/usr/local/qd/lib && make && cp sdpa_dd /usr/bin && cd / && rm -rf sdpa_dd
RUN git clone https://github.com/ericphanson/sdpa-qd.git && cd sdpa-qd && ./configure --with-qd-includedir=/usr/local/qd/include --with-qd-libdir=/usr/local/qd/lib && make && cp sdpa_qd /usr/bin && cd / && rm -rf sdpa_qd

USER sage

RUN sage -i csdp

RUN sage --pip install tqdm

RUN git clone https://github.com/jsliacan/flagmatic.git && cd flagmatic/pkg && sage -python setup.py install && cd ../.. && rm -rf flagmatic
# https://github.com/FordUniver/flagmatic.git


#bc binutils bzip2 ca-certificates cliquer cmake curl ecl eclib-tools fflas-ffpack flintqs g++ g++ gcc gcc gengetopt gfan gfortran libgfortran glpk-utils gmp-ecm lcalc libatomic-ops-dev libboost-dev libbrial-dev libbrial-groebner-dev libbz2-dev libcdd-dev libcdd-tools libcliquer-dev libcurl4-openssl-dev libec-dev libecm-dev libffi-dev libflint-arb-dev libflint-dev libfreetype6-dev libgc-dev libgd-dev libgf2x-dev libgiac-dev libgivaro-dev libglpk-dev libgmp-dev libgsl-dev libiml-dev liblfunction-dev liblrcalc-dev liblzma-dev libm4rie-dev libmpc-dev libmpfi-dev libmpfr-dev libncurses5-dev libntl-dev libopenblas-dev libpari-dev libpcre3-dev libplanarity-dev libppl-dev libprimesieve-dev libpython3-dev libqhull-dev libreadline-dev librw-dev libsingular4-dev libsqlite3-dev libssl-dev libsuitesparse-dev libsymmetrica2-dev libz-dev libzmq3-dev libzn-poly-dev m4 make nauty openssl palp pari-doc pari-elldata pari-galdata pari-galpol pari-gp2c pari-seadata patch perl pkg-config planarity ppl-dev r-base-dev r-cran-lattice singular sqlite3 sympow tachyon tar tox xcas xz-utils xz-utils


# python bc binutils bzip2 ca-certificates cliquer cmake curl ecl eclib-tools fflas-ffpack flintqs g++ g++ gcc gcc gengetopt gfan gfortran libgfortran glpk-utils gmp-ecm lcalc libatomic-ops-dev libboost-dev libbrial-dev libbrial-groebner-dev libbz2-dev libcdd-dev libcdd-tools libcliquer-dev libcurl4-openssl-dev libec-dev libecm-dev libffi-dev libflint-arb-dev libflint-dev libfreetype6-dev libgc-dev libgd-dev libgf2x-dev libgiac-dev libgivaro-dev libglpk-dev libgmp-dev libgsl-dev libiml-dev liblfunction-dev liblrcalc-dev liblzma-dev libm4rie-dev libmpc-dev libmpfi-dev libmpfr-dev libncurses5-dev libntl-dev libopenblas-dev libpari-dev libpcre3-dev libplanarity-dev libppl-dev libprimesieve-dev libpython3-dev libqhull-dev libreadline-dev librw-dev libsingular4-dev libsqlite3-dev libssl-dev libsuitesparse-dev libsymmetrica2-dev libz-dev libzmq3-dev libzn-poly-dev m4 make nauty openssl palp pari-doc pari-elldata pari-galdata pari-galpol pari-gp2c pari-seadata patch perl pkg-config planarity ppl-dev r-base-dev r-cran-lattice singular sqlite3 sympow tachyon tar tox xcas xz-utils xz-utils
# RUN ln -s /usr/bin/python2 /usr/bin/python
# python3 python3 python3-distutils libbraiding-dev  libhomfly-dev

# COPY --from=build /opt/sage /opt/sage
# COPY --from=build /usr/bin/sage /usr/bin/sage
# COPY --from=build /usr/bin/sagemath /usr/bin/sagemath

# ARG SAGE_ROOT=/opt/sage

# USER root

# RUN sage -i csdp

# RUN sage --pip install tqdm

# RUN git clone https://github.com/jsliacan/flagmatic.git && cd flagmatic/pkg && sage -python setup.py install && cd ../.. && rm -rf flagmatic
# https://github.com/FordUniver/flagmatic.git

# USER root

# libsdpa-dev

# WORKDIR /home/sage/flagmatic

# COPY sdpa-dd-7.1.2 /home/sage/sdpa-dd-7.1.2
# COPY sdpa-qd-7.1.2 /home/sage/sdpa-qd-7.1.2
# COPY QD-master /home/sage/QD-master

ENTRYPOINT /bin/bash
# sage -n jupyter --ip=127.0.0.1