FROM pytorch/pytorch:1.4-cuda10.1-cudnn7-devel
MAINTAINER Kishwar Shafin, kishwar.shafin@gmail.com

# update and install dependencies
RUN apt-get update && \
    apt-get -y install time git make wget autoconf gcc g++ && \
    apt-get -y install autoconf bzip2 lzma-dev zlib1g-dev && \
    apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev && \
    apt-get -y install liblzma-dev libhdf5-dev libncurses5-dev libncursesw5-dev && \
    apt-get -y install python3 python3-dev python3-pip && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# install cmake
WORKDIR /tmp
RUN mkdir /opt/cmake && \
    wget https://github.com/Kitware/CMake/releases/download/v3.14.4/cmake-3.14.4-Linux-x86_64.sh && \
    sh /tmp/cmake-3.14.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

# install requirements
RUN python3 -m pip install h5py tqdm numpy pyyaml

# install samtools, this is for kubernetes
WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -vxjf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    make

# get PEPPER
WORKDIR /opt
ADD https://api.github.com/repos/kishwarshafin/pepper/git/refs/heads/master version.json
RUN git clone https://github.com/kishwarshafin/pepper.git pepper/
RUN cd /opt/pepper && \
    git fetch && \
    git submodule update --init && \
    git pull origin master

WORKDIR /opt/pepper/build
RUN cmake .. -Wno-deprecated && make

WORKDIR /opt/pepper
# setup entrypoint
COPY pepper_wrapper.sh /opt/pepper_wrapper.sh

WORKDIR /data
ENTRYPOINT ["sh"]