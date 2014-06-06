#! /usr/bin/env python
#
# This test uses GPAW for the primary system and pysic for the secondary.
# The embedding scheme is set at mechanical embedding with hydrogen links.
#
#==============================================================================
from ase import Atoms, Atom
from ase.visualize import view
from ase.structure import molecule
from ase.md.verlet import VelocityVerlet
from ase.constraints import FixInternals, FixBondLengths
from ase import units
from ase.io import write
from pysic import *

# The system consists of two water molecules, one MM other QM
water = Atoms()
constraints = []

n = 2
for i in range(n):
    h2o =  molecule('H2O')
    # Set parameters according to TIP3P water model
    qh = 0.417
    qo = -0.834
    bond = 0.9572
    angle = 104.52/180.0*np.pi
    h2o.set_initial_charges([qo, qh, qh])
    h2o.set_distance(0, 1, bond)
    h2o.set_distance(0, 2, bond)
    h2o.set_angle([1, 0, 2], angle)
    if i == 1:
        h2o.set_tags([i]*3)
        h2o.translate((1, 1, 1))
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
            epsilon=1.E-2
            )
    constraints.append(constraint)
    
water.set_constraint(constraints)
water.set_cell((10, 10, 10))
water.center()

# Setup the calculator and subsystems
hc = HybridCalculator()
ps = hc.add_subsystem("primary")
ps.set_calculator(Pysic())
ps.set_atoms(tag=0)
ss = hc.add_subsystem("secondary")
ss.set_calculator(Pysic())
ss.set_atoms(tag=1)

# Set binding options between the subsystems
binding = hc.add_binding("primary", "secondary")
binding.set_hydrogen_links((0, 3), 0.71)
binding.set_electrostatic_binding()
binding.add_binding_potential(Potential('LJ', symbols=['O','O'], cutoff=10))

# Calculate energies
water.set_calculator(hc) # This also calls hc.set_atoms(water)

# Setup viewer
viewer = AtomEyeViewer(water, "/home/lauri/water")

# Dynamics
dyn = VelocityVerlet(water, 0.4*units.fs)
dyn.attach(viewer.save_cfg_frame, 2)
dyn.run(400)
viewer.view_series()
