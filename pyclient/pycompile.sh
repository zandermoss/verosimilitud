#! /bin/bash

CFLAGS="-I../inc"  \
LDFLAGS="-L../"     \
CC=$(CXX) python setup.py build_ext -i


