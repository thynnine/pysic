#! /usr/bin/env python
#
# This test uses GPAW for the primary system and pysic for the secondary.
# The embedding scheme is set at mechanical embedding with hydrogen links.
#
#===============================================================================
from ase import Atoms
from pysic import *
from pysic.utility import visualization
from gpaw import GPAW
from ase.visualize import view
from ase.io import *

#-------------------------------------------------------------------------------
# Prepare the system
h2 = Atoms('H2', [(0, 0, 0), (2, 0, 0)], charges=(1, -1))
b = 10.0
h2.set_cell((b, b, b))
#h2.set_cell([(0, b, b), (b, 0, b), (b, b, 0)])
h2.center()

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hc = HybridCalculator()

# Initialize calculator for subsystems
pysic_calc = Pysic()
potential = Potential('LJ', cutoff=4.0, symbols=['H', 'H'], parameters=[0.1, 2.5])
pysic_calc.add_potential(potential)
gpaw_calc = GPAW(h=b/16.0, txt=None)

# Define subsystems
hc.add_subsystem(SubSystem("primary", indices=0, calculator=gpaw_calc))
hc.add_subsystem(SubSystem("secondary", special_set="remaining", calculator=pysic_calc))

# Define an embedding scheme between the subsystems
# In this case the scheme is mechanical embedding with hydrogen links
binding = Binding("primary", "secondary", coulomb_interaction=True)
binding.set_hydrogen_links((0, 1), 1)
hc.add_binding(binding)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hc.get_potential_energy(h2)
hc.print_energies()
hc.view_subsystems()
