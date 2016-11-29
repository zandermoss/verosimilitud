#!/bin/bash

CFLAGS="-I../inc"  \
LDFLAGS="-L../"     \
#CC=$(CXX) python setup.py build_ext -i
CC="clang++ -std=c++11" python setup.py build_ext -i


