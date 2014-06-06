#! /usr/bin/env python
#
# Simulating water with Pysic
#
#==============================================================================
from ase import Atoms
from ase.constraints import FixInternals, FixBondLengths
from ase.structure import molecule
from pysic import *
import numpy as np
from ase.optimize import BFGS
from ase.io.trajectory import PickleTrajectory
from ase.md.verlet import VelocityVerlet
from ase import units
from pysic.utility.visualization import AtomEyeViewer

#------------------------------------------------------------------------------
# The structure of the system is defined here. The initial system consists of n
# water molecules in a "lattice"

n = 2 # Number of water molecules
cell = int(np.ceil(np.power(n, 1./3.))) # Ideal cubic cell dimension
distance = 2 # Distance betweeen molecules
padding = 5 # Padding between molecules and the cell
water = Atoms()
constraints = []
max_x, max_y, max_z = (0,)*3

for i in range(n):

    # Create a water molecule
    h2o = molecule('H2O')

    # Calculate the molecule position
    z = int(i/cell/cell)%cell
    y = int(i/cell)%cell
    x = i%cell
    if x>max_x: max_x = x
    if y>max_y: max_y = y
    if z>max_z: max_z = z
    position = np.array([x, y, z])*distance+padding
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
    constraint = FixInternals(
            water,
            bonds=[bond1, bond2],
            angles=[angle],
            dihedrals=[],
            epsilon=1.E-3
            )
    constraints.append(constraint)

water.set_constraint(constraints)
cell = np.array([max_x, max_y, max_z])*distance+2*padding
water.set_cell(cell.tolist())

#------------------------------------------------------------------------------
# Assign a pysic calculator to the system. Use the parameters provided by the
# TIP3P model. Disregard the internal forces in a water molecule, that is the
# hydrogens and oxygen within one molecule do not interact.

calc = Pysic()
kcalPerMoleInEv = 0.0436
A = 582.0E3*kcalPerMoleInEv
B = 595.0*kcalPerMoleInEv
epsilon = 0.0052635
kc = 1.0/(4.0*np.pi*epsilon)
max_cutoff = np.linalg.norm(water.get_cell())

# Lennard-Jones potential in two parts
pot1 = Potential('power',symbols=['O', 'O'], parameters=[A, 1, 12], cutoff=max_cutoff)
pot2 = Potential('power', symbols=['O', 'O'], parameters=[-B, 1, 6], cutoff=max_cutoff)

# Extract pairs that should interact with coulomb interaction
coulomb_pairs = []
for i in range(len(water)):
    for j in range(len(water)):
        if (j>i) and (int(i/3.0) is not int(j/3.0)):
            coulomb_pairs.append([i, j])

print coulomb_pairs

# The first potential given to the productpotential defines the targets and cutoff
coul1 = Potential('power', indices=coulomb_pairs, parameters=[1, 1, 1], cutoff=max_cutoff)
coul2 = Potential('charge_pair', parameters=[kc, 1, 1])
pot3 = ProductPotential([coul1, coul2])

potential_list = [pot1, pot2, pot3]
calc.add_potential(potential_list)
water.set_calculator(calc)

#------------------------------------------------------------------------------
# Run dynamics
steps = 400
interval = steps/10
step = 0

def print_progress():
    global step, interval, steps
    step += interval
    print str(float(step)/steps*100)+'%'

dyn = VelocityVerlet(water, 0.4*units.fs)
visuals = AtomEyeViewer(water, "/home/lauri/water")
dyn.attach(visuals.save_cfg_frame, interval=1)
dyn.attach(print_progress, interval=interval)
dyn.run(steps)
#visuals.save_jpg_series()
visuals.view_series()

