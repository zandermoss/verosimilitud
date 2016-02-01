#! /bin/bash

CFLAGS="-I../inc"  \
LDFLAGS="-L../"     \
python setup.py build_ext -i
