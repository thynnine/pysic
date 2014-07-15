#! /usr/bin/env python
#
# Test dynamics in a hydrid system. The system consists of two water molecules,
# one QM the other MM. These interact with electrostatic forces.
#
#===============================================================================
from ase import Atoms
from ase.constraints import FixInternals
from ase.structure import molecule
from pysic import *
import numpy as np
from ase.md.verlet import VelocityVerlet
from ase import units
from pysic.utility.visualization import AtomEyeViewer
from gpaw import GPAW

#-------------------------------------------------------------------------------
# Prepare the system
n = 2  # Number of water molecules
distance = 5  # distance between molecules
padding = 2  # padding between molecules and cell wall
dim = int(np.ceil(np.power(n, 1./3.)))  # Ideal cubic cell dimension
water = Atoms()
constraints = []

max_x = 0
max_y = 0
max_z = 0

for i in range(n):

    # Create a water molecule
    h2o = molecule('H2O')
    h2o.set_tags([i]*3)

    # Calculate the molecule position
    z = int(i/dim/dim) % dim
    y = int(i/dim) % dim
    x = i % dim
    if x > max_x:
        max_x = x
    if y > max_y:
        max_y = y
    if z > max_z:
        max_z = z
    position = np.array([x, y, z])*distance+padding
    center = np.array(h2o.get_center_of_mass())
    translation = position-center
    h2o.translate(translation.tolist())

    # Set parameters according to TIP3P water model
    qh = 0.417
    qo = -0.834
    bond = 0.9572
    angle = 104.52/180.0*np.pi

    h2o.set_distance(0, 1, bond)
    h2o.set_distance(0, 2, bond)
    h2o.set_angle([1, 0, 2], angle)

    # Set initial charges on classical molecules
    if i != 0:
        h2o.set_initial_charges([qo, qh, qh])

    # Add the molecule to the system
    water += h2o

    # Fix the water molecule bond lengths and angles in the classical system
    if i != 0:
        bond1 = [bond, [i*3, i*3+1]]
        bond2 = [bond, [i*3, i*3+2]]
        angle_indices = [i*3+1, i*3, i*3+2]
        angle = [angle, angle_indices]
        constraint = FixInternals(water, bonds=[bond1, bond2], angles=[angle], dihedrals=[], epsilon=1.E-3)
        constraints.append(constraint)

water.set_constraint(constraints)
cell = np.array([max_x, max_y, max_z])*distance+2*padding
water.set_cell(cell.tolist())

#-------------------------------------------------------------------------------
# Prepare the subsystems
hc = HybridCalculator()

# Setup the subsystem calculators
gpaw_calc = GPAW(h=0.3, txt=None)
pysic_calc = Pysic()

PS = SubSystem("primary", tag=0, calculator=gpaw_calc)
hc.add_subsystem(PS)

SS = SubSystem("secondary", tag=1, calculator=pysic_calc)
hc.add_subsystem(SS)

#-------------------------------------------------------------------------------
# Define an embedding scheme between the subsystems
binding = Binding("primary", "secondary")
binding.set_coulomb_interaction()
hc.add_binding(binding)

#-------------------------------------------------------------------------------
# Set the calculator
water.set_calculator(hc)

#-------------------------------------------------------------------------------
# Visualize system
viewer = AtomEyeViewer(water, "/home/lauri/water")
hc.set_atoms(water)
colors = hc.get_colors()
viewer.set_colors(colors)
viewer.view()

#------------------------------------------------------------------------------
# Run some dynamics
water.set_calculator(hc)
steps = 10
interval = 1
step = 0


def print_progress():
    global step, interval, steps
    step += interval
    print str(float(step)/steps*100)+'%'

dyn = VelocityVerlet(water, 0.5*units.fs)
dyn.attach(viewer.save_cfg_frame, interval=1)
dyn.attach(print_progress, interval=interval)
dyn.run(steps)
viewer.view_series()
