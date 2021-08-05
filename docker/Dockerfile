FROM debian:buster
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y locales \
    python3-tz \
    python3-casacore \
    python3-pip \
    cmake build-essential git gfortran

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8

WORKDIR /build
RUN git clone https://github.com/mxjeff/xfoil-python.git
WORKDIR /build/xfoil-python
RUN pip3 install .

WORKDIR /build
RUN git clone https://github.com/elec-otago/proply.git
WORKDIR /build/proply
RUN pip3 install .

ENV LANG en_US.UTF-8 
