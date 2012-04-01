#! /usr/bin/env python
import sys
import inspect

def isdebugging():
    for frame in inspect.stack():
        if frame[1].endswith("pdb.py"):
            return True
    return False

debug_on = isdebugging()

if debug_on:
    import pdb
    print "\n\n   Running with the PDB debugger   \n\n"

def bp(condition=True):
    if debug_on and condition:
        pdb.set_trace()
