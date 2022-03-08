# Install PEPPER locally

## Installation
We recommend using `Linux` environment to run `PEPPER`.

### Install dependencies
```bash
sudo apt-get -y install cmake make git gcc g++ autoconf bzip2 lzma-dev zlib1g-dev \
libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev \
liblzma-dev libhdf5-dev python3-pip python3-virtualenv virtualenv

git clone https://github.com/kishwarshafin/pepper.git
cd pepper
make install
. ./venv/bin/activate

pepper_variant --help
```
