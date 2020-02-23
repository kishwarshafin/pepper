#!/bin/bash

# rm -rf build
mkdir build
cd build

cmake .. -Wno-deprecated && make
cd ..
