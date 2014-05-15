#! /usr/bin/env python
#
# This script tests that the hydrogen link in the MEHL scheme is working as it
# should. Namely it calculates the energy of a special system with hybrid
# division and with traditional ways and checks that they give the same results.
#
#===============================================================================
from ase import Atoms
from pysic import Pysic, CoulombSummation, Potential
from pysic import hybridcalculation
from ase.visualize import view

#-------------------------------------------------------------------------------
# Prepare the system
h2 = Atoms('H2', [(0, 0, 0), (1, 1, 1)])
h2.set_cell((2, 2, 2))

# Set periodic boundary conditions
h2.set_pbc(True)

# Use ASEs built in viewer to make sure the structure is correct:
#view(h2)

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hybrid_calculation = hybridcalculation.HybridCalculation()
hybrid_calculation.set_system(h2)

# Define QM/MM regions. You can get the indices by e.g.
# examining the the structure in ASEs viewer. The oxygen and the carbon
# connected to it are set to belong to the QM-region.
# The structure of the subsystems is stored internally within pysic, so no
# changes are made to the Atoms object.
hybrid_calculation.set_primary_system([0])
hybrid_calculation.set_secondary_system(special_set='remaining')

# Initialize calculators
pysic_calc1 = Pysic()
pysic_calc2 = Pysic()
potential = Potential('LJ', cutoff=4.0, symbols=['H', 'H'], parameters=[0.1, 2.5])
pysic_calc1.add_potential(potential)
pysic_calc2.add_potential(potential)

# Add different calculators for the subsystems
hybrid_calculation.set_primary_calculator(pysic_calc1)
hybrid_calculation.set_secondary_calculator(pysic_calc2)

#-------------------------------------------------------------------------------
# Define an embedding scheme between the subsystems
# In this case the scheme is mechanical embedding with hydrogen links
parameters = {
    'links': ((0, 1),),
    'CHL': 1
}
link = parameters['links']
hybrid_calculation.set_embedding('MEHL', 'primary', 'secondary', parameters)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hybrid_energy = hybrid_calculation.get_potential_energy()
#hybrid_calculation.view_subsystems()

# Calculate the energy of the same setup, but use only one region. In this
# special case these energies should be same.
pysic_calc3 = Pysic()
pysic_calc3.add_potential(potential)
h2.set_calculator(pysic_calc3)
real_energy = h2.get_potential_energy()

print "Energy with hybrid calculation: " + str(hybrid_energy)
print "Energy with traditional calculation: " + str(real_energy)
