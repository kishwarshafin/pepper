#!/bin/bash

# rm -rf build
mkdir build
cd build

cmake .. -Wno-deprecated && make -j 8
cd ..
# python3 main.py