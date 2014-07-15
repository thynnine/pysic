#! /usr/bin/env python
#
# This test uses GPAW for the primary system and pysic for the secondary.
# The embedding scheme is set at mechanical embedding with hydrogen links.
#
#===============================================================================
from ase import Atoms
from pysic import *
from gpaw import GPAW
from ase.io import *

#-------------------------------------------------------------------------------
# Prepare the system
h2 = Atoms('He3', [(5, 5, 5), (6, 6, 6), (6, 6, 7)], charges=(0, 0, 1))
b = 6.0
h2.set_cell((b, b, b))
h2.center()

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hc = HybridCalculator(record_time_usage=True)

# Initialize calculator for subsystems
pysic_calc = Pysic()
potential = Potential('LJ', cutoff=4.0, symbols=['H', 'H'], parameters=[0.1, 2.5])
pysic_calc.add_potential(potential)
gpaw_calc = GPAW(h=0.4, txt=None)

# Define subsystems
hc.add_subsystem(SubSystem("primary", indices=[0, 1], calculator=gpaw_calc))
hc.add_subsystem(SubSystem("secondary", special_set="remaining", calculator=pysic_calc))

# Define an embedding scheme between the subsystems
# In this case the scheme is mechanical embedding with hydrogen links
binding = Binding("primary", "secondary")
binding.set_coulomb_interaction()
binding.set_hydrogen_links((1, 2), 1)
hc.add_binding(binding)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hc.get_potential_energy(h2)
hc.print_summary()
hc.print_charge_summary()
hc.view_subsystems()
