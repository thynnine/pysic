#! /usr/bin/env python
#
# This test uses GPAW for the primary system and pysic for the secondary.
# The embedding scheme is set at mechanical embedding with hydrogen links.
#
#===============================================================================
from ase import Atoms
from pysic import Pysic, CoulombSummation, Potential, HybridSystem
from gpaw import GPAW
from ase.visualize import view
from ase.io import *

#-------------------------------------------------------------------------------
# Prepare the system
h2 = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)])
h2.set_cell((5, 5, 5))
h2.center()

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hybrid_system = HybridSystem()
hybrid_system.set_system(h2)

# Define QM/MM regions
hybrid_system.set_primary_system([0])
hybrid_system.set_secondary_system(special_set='remaining')

# Initialize calculators
pysic_calc = Pysic()
potential = Potential('LJ', cutoff=4.0, symbols=['H', 'H'], parameters=[0.1, 2.5])
pysic_calc.add_potential(potential)
gpaw_calc = GPAW(txt=None)

# Add different calculators for the subsystems
hybrid_system.set_primary_calculator(gpaw_calc)
hybrid_system.set_secondary_calculator(pysic_calc)

#-------------------------------------------------------------------------------
# Define an embedding scheme between the subsystems
# In this case the scheme is mechanical embedding with hydrogen links
parameters = {
    'links': ((0, 1),),
    'CHL': 1,
    'epsilon': 0.0052635
}
hybrid_system.set_embedding('MEHL', 'primary', 'secondary', parameters)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hybrid_system.get_potential_energy()
hybrid_system.print_potential_energies()
hybrid_system.view_subsystems()
