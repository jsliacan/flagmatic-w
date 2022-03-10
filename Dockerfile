FROM sagemath/sagemath-jupyter:7.4

USER root

RUN apt-get install -y software-properties-common && add-apt-repository universe multiverse && apt-get update && apt-get install -y libsdpa-dev

RUN sage -i csdp

RUN sage --pip install tqdm

RUN git clone https://github.com/FordUniver/flagmatic.git && cd flagmatic/pkg && sage -python setup.py install && cd ../.. && rm -rf flagmatic

# USER root

# libsdpa-dev

# WORKDIR /home/sage/flagmatic

ENTRYPOINT /bin/bash
#sage -n jupyter --ip=127.0.0.1