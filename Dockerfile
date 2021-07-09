FROM sagemath/sagemath-jupyter:7.4

RUN sage -i csdp

RUN git clone https://github.com/FordUniver/flagmatic.git

RUN cd flagmatic/pkg && sage -python setup.py install

# COPY create_SDP.sage .
WORKDIR /flagmatic

ENTRYPOINT /bin/bash
#sage -n jupyter --ip=127.0.0.1