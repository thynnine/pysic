#! /usr/bin/env python
#
# This test uses GPAW for the primary system and pysic for the secondary.
# The embedding scheme is set at mechanical embedding with hydrogen links.
#
#===============================================================================
from ase import Atoms
from ase.structure import molecule
from pysic import Pysic, CoulombSummation, Potential, HybridSystem
from gpaw import GPAW
from ase.visualize import view
from ase.io import write

#-------------------------------------------------------------------------------
# Prepare the system
acetone = molecule('CH3COCH3')
acetone.set_cell((10, 10, 10))
acetone.center()

# View the molecule and save the structure
view(acetone)
write('acetone.cfg', acetone)

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hybrid_calculation = HybridSystem()
hybrid_calculation.set_system(acetone)

# Define QM/MM regions
hybrid_calculation.set_primary_system([0, 1])
hybrid_calculation.set_secondary_system(special_set='remaining')

# Initialize calculators
pysic_calc = Pysic()
potential1 = Potential('LJ', cutoff=10.0, symbols=['C', 'C'], parameters=[0.1, 2.5])
potential2 = Potential('LJ', cutoff=10.0, symbols=['C', 'H'], parameters=[0.1, 2.5])
potential_list = [potential1, potential2]
pysic_calc.add_potential(potential_list)
gpaw_calc = GPAW(txt=None)

# Add different calculators for the subsystems
hybrid_calculation.set_primary_calculator(gpaw_calc)
hybrid_calculation.set_secondary_calculator(pysic_calc)

#-------------------------------------------------------------------------------
# Define an embedding scheme between the subsystems
# In this case the scheme is mechanical embedding with hydrogen links
parameters = {
    'links': ((1, 2), (1, 3)),
    'CHL': 0.71,
    'epsilon': 0.0052635
}
hybrid_calculation.set_embedding('MEHL', 'primary', 'secondary', parameters)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hybrid_calculation.get_potential_energy()
hybrid_calculation.print_potential_energies()
hybrid_calculation.view_subsystems()
