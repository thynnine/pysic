#! /usr/bin/env python
#
# Testing pysic MPI parallellization
#
#===============================================================================
from ctypes import *
mpi = CDLL('libmpi.so.0', RTLD_GLOBAL)
from pysic import *
