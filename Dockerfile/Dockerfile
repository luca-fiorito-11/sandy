FROM ubuntu:18.04
MAINTAINER Jan Malec (jan.malec@gmail.com)

RUN apt update
RUN apt install -y gfortran python3 python3-setuptools\
    python3-numpy python3-dev unzip ca-certificates git

RUN git clone https://github.com/janmalec/sandy.git
RUN cd sandy && python3 setup.py install
