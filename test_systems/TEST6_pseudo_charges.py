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
from gpaw import GPAW

#-------------------------------------------------------------------------------
# Prepare the system
n = 1 # Number of water molecules
distance = 5 # distance between molecules
padding = 2 # padding between molecules and cell wall
dim = int(np.ceil(np.power(n, 1./3.))) # Ideal cubic cell dimension
water = Atoms()

max_x = 0
max_y = 0
max_z = 0

for i in range(n):

    # Create a water molecule
    h2o = molecule('H2O')
    if i == 0:
        h2o.set_tags([0, 0, 1])

    # Calculate the molecule position
    z = int(i/dim/dim)%dim
    y = int(i/dim)%dim
    x = i%dim
    if x>max_x:
        max_x = x
    if y>max_y:
        max_y = y
    if z>max_z:
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
    h2o.set_initial_charges([qo, qh, qh])
    h2o.set_distance(0, 1, bond)
    h2o.set_distance(0, 2, bond)
    h2o.set_angle([1, 0, 2], angle)

    # Add the molecule to the system
    water += h2o

    ## Fix the water molecule bond lengths and angles in the classical system
    #if i == 0:
        #bond1 = [bond, [i*3, i*3+1]]
        #bond2 = [bond, [i*3, i*3+2]]
        #angle_indices = [i*3+1, i*3, i*3+2]
        #angle = [angle, angle_indices]
        #constraint = FixInternals(water, bonds=[bond1, bond2], angles=[angle], dihedrals=[], epsilon=1.E-3)
        #water.set_constraint(constraint)

cell = np.array([max_x, max_y, max_z])*distance+2*padding
water.set_cell(cell.tolist())
view(water)

#-------------------------------------------------------------------------------
# Prepare the hybrid calculator
hc = HybridCalculator(water)
hc.set_primary_system(tag=1)
hc.set_secondary_system(tag=0)

#-------------------------------------------------------------------------------
# Prepare the secondary pysic calculator. Use the parameters provided by the
# TIP3P model. Disregard the internal forces in a water molecule, that is the
# hydrogens and oxygen within one molecule do not interact.
pysic_calc = Pysic()
kcalPerMoleInEv = 0.0436
A = 582.0E3*kcalPerMoleInEv
B = 595.0*kcalPerMoleInEv
epsilon = 0.0052635
kc = 1.0/(4.0*np.pi*epsilon)

max_cutoff = np.linalg.norm(water.get_cell())
pot1 = Potential('power', symbols=['O', 'O'], parameters=[A, 1, 12], cutoff=max_cutoff)
pot2 = Potential('power', symbols=['O', 'O'], parameters=[-B, 1, 6], cutoff=max_cutoff)

# The first potential given to the productpotential defines the targets and cutoff
coul1 = Potential('power', indices=[], parameters=[1, 1, 1], cutoff=max_cutoff)
coul2 = Potential('charge_pair', parameters=[kc, 1, 1])
pot3 = ProductPotential([coul1, coul2])

potential_list = [pot1, pot2]
pysic_calc.add_potential(potential_list)
hc.set_secondary_calculator(pysic_calc)

#-------------------------------------------------------------------------------
# Setup the primary GPAW calculator
gpaw_calc = GPAW(h=0.5, txt=None)
hc.set_primary_calculator(gpaw_calc)

#-------------------------------------------------------------------------------
# Define an embedding scheme between the subsystems
hc.set_embedding('MEHL', 'primary', 'secondary', {})
hc.get_potential_energy()
hc.print_potential_energies()
