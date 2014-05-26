#! /usr/bin/env python
#
# This test uses GPAW for the primary system and pysic for the secondary.
# The embedding scheme is set at mechanical embedding with hydrogen links.
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

#-------------------------------------------------------------------------------
# Prepare the system
d = 5.0 # Cell dimension
n = 8 # Number of water molecules
water = Atoms(cell=[d, d, d], pbc=False)
constraints = []

# Calculate the cell size for an uniformly filled unit cell
cell = int(np.ceil(np.power(n, 1./3.)))

for i in range(n):

    # Create a water molecule
    h2o = molecule('H2O')

    # Calculate the molecule position
    z = int(i/cell/cell)%cell
    y = int(i/cell)%cell
    x = i%cell
    position = np.array([x, y, z])*d/cell+d/cell*0.5
    center = np.array(h2o.get_center_of_mass())
    translation = position-center
    h2o.translate(translation.tolist())

    # Set parameters according to TIP3P water model
    qh = 0.417
    qo = -0.834
    bond = 0.9572
    angle = 104.52/180.0*np.pi
    h2o.set_initial_charges([qo, qh, qh])
    h2o.set_distance(0, 1, bond)
    h2o.set_distance(0, 2, bond)
    h2o.set_angle([1, 0, 2], angle)

    # Add the molecule to the system
    water += h2o

    # Fix the water molecule bond lengths and angles
    bond1 = [bond, [i*3, i*3+1]]
    bond2 = [bond, [i*3, i*3+2]]
    angle_indices = [i*3+1, i*3, i*3+2]
    angle = [angle, angle_indices]
    constraint = FixInternals(water, bonds=[bond1, bond2], angles=[angle], dihedrals=[], epsilon=1.E-3)
    constraints.append(constraint)

water.set_constraint(constraints)

#-------------------------------------------------------------------------------
# Assign a pysic calculator to the system. Use the parameters provided by the
# TIP3P model. Disregard the internal forces in a wate molecule, that is the
# hydrogens and oxygen within one molecule do not interact.
calc = Pysic()
kcalPerMoleInEv = 0.0436
A = 582.0E3*kcalPerMoleInEv
B = 595.0*kcalPerMoleInEv
epsilon = 0.0052635
kc = 1.0/(4.0*np.pi*epsilon)

pot1 = Potential('power', symbols=['O', 'O'], parameters=[A, 1, 12], cutoff=d)
pot2 = Potential('power', symbols=['O', 'O'], parameters=[-B, 1, 6], cutoff=d)

# Extract pairs that should interact with coulomb interaction
coulomb_pairs = []
for i in range(len(water)):
    for j in range(len(water)):
        if (j>i) and (int(i/3.0) is not int(j/3.0)):
            coulomb_pairs.append([i, j])

# The first potential given to the productpotential defines the targets and cutoff
coul1 = Potential('power', indices=coulomb_pairs, parameters=[1, 1, 2], cutoff=d)
coul2 = Potential('charge_pair', parameters=[kc, 1, 1])
pot3 = ProductPotential([coul1, coul2])

potential_list = [pot1, pot2, pot3]
calc.add_potential(potential_list)
water.set_calculator(calc)

#-------------------------------------------------------------------------------
# Run structure optimization
dyn = BFGS(water)
dyn.run(fmax=0.05)
visuals = AtomEyeViewer("A", water, "/home/lauri/water")
view(water)

#-------------------------------------------------------------------------------
# Run dynamics
#dyn = VelocityVerlet(water, 5*units.fs)
#visuals = AtomEyeViewer("A", water, "/home/lauri/water")
#dyn.attach(visuals.save_cfg, interval=2)
#dyn.run(400)
#visuals.view_series()

