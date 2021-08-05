FROM debian:buster
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y locales \
    python3-pip python3-numpy python3-matplotlib python3-scipy \
    cmake build-essential git gfortran

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8

ENV LANG en_US.UTF-8 

WORKDIR /build
RUN git clone https://github.com/mxjeff/xfoil-python.git
WORKDIR /build/xfoil-python
RUN python3 setup.py install

WORKDIR /build/proply
ADD . .

RUN python3 setup.py install

WORKDIR /run
