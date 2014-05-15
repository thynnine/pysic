#! /usr/bin/env python
#
# Tests that the Coulomb energy calculation is done correctly for hybrid systems
# with charges
#
#===============================================================================
from ase import Atoms, Atom
from pysic import Pysic, CoulombSummation, Potential
from pysic import hybridcalculation
from ase.visualize import view

#-------------------------------------------------------------------------------
# Prepare the system
h1 = Atom('H', (0, 0, 0), charge=-1)
h2 = Atom('H', (1, 0, 0), charge=1)
h3 = Atom('H', (0, 1, 0), charge=-1)
h4 = Atom('H', (0, 0, 1), charge=1)
system = Atoms()
system.append(h1)
system.append(h2)
system.append(h3)
system.append(h4)
system.set_cell((2, 2, 2))

# Use ASEs built in viewer to make sure the structure is correct:
view(system)

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hybrid_calculation = hybridcalculation.HybridCalculation()
hybrid_calculation.set_system(system)

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
ewald = CoulombSummation()
ewald.set_parameter_value('epsilon',0.00552635)
ewald.set_parameter_value('k_cutoff',0.7)
ewald.set_parameter_value('real_cutoff',10.0) #
ewald.set_parameter_value('sigma',1.4)
pysic_calc1.set_coulomb_summation(ewald)
pysic_calc2.set_coulomb_summation(ewald)

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
pysic_calc3.set_coulomb_summation(ewald)
system.set_calculator(pysic_calc3)
real_energy = system.get_potential_energy()

print "Energy with hybrid calculation: " + str(hybrid_energy)
print "Energy with traditional calculation: " + str(real_energy)
