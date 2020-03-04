FROM pytorch/pytorch:1.4-cuda10.1-cudnn7-devel
MAINTAINER Kishwar Shafin, kishwar.shafin@gmail.com

# update and install dependencies
RUN apt-get update && \
    apt-get -y install git cmake make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev && \
    apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev liblzma-dev libhdf5-dev && \
    apt-get -y install python3-pip python3-virtualenv virtualenv && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# get HELEN
WORKDIR /opt
RUN git clone https://github.com/kishwarshafin/pepper.git && \
    cd /opt/pepper && \
    git fetch && \
    git submodule update --init && \
    git pull origin master && \
    python3 -m pip install .

WORKDIR /opt/pepper

# setup entrypoint
WORKDIR /data