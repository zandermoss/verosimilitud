#! /bin/bash

CFLAGS="-I/home/pinkpig/physics/neutrino_decay/likelihood"  \
LDFLAGS="-L/home/pinkpig/physics/neutrino_decay/likelihood"     \
python setup.py build_ext -i
