#! /usr/bin/env python
#
#
#===============================================================================
from ase import Atoms
from ase.constraints import FixInternals, FixBondLengths
from ase.structure import molecule
from pysic import *
from ase.visualize import view
from ase.io import write, read
import numpy as np
from ase.optimize import BFGS
from ase.io.trajectory import PickleTrajectory
from ase.md.verlet import VelocityVerlet
from ase import units
from pysic.utility.visualization import AtomEyeViewer
from gpaw import GPAW, PW

#-------------------------------------------------------------------------------
# Prepare the system
dim = 4.0
h = Atoms('H', cell=(dim, dim, dim))
h.center()
calculator = GPAW(h=dim/8, hund=True, txt='H.paw.txt')
h.set_calculator(calculator)
e = h.get_potential_energy()
density = calculator.get_pseudo_density()
print density
print e
