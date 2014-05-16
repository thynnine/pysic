#! /usr/bin/env python
#
# Tests that the Coulomb energy calculation is done correctly for hybrid systems
# with charges
#
#===============================================================================
from ase import Atoms, Atom
from pysic import Pysic, Potential, HybridSystem
from pysic.interactions.coulomb import *
from ase.visualize import view

#-------------------------------------------------------------------------------
# Prepare the system
h1 = Atom('H', (0, 0, 0), charge=-1)
h2 = Atom('H', (2, 0, 0), charge=1)
h3 = Atom('H', (0, 2, 0), charge=-1)
h4 = Atom('H', (0, 0, 2), charge=1)
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
hybrid_system = HybridSystem()
hybrid_system.set_system(system)

# Define QM/MM regions
hybrid_system.set_primary_system([0])
hybrid_system.set_secondary_system(special_set='remaining')

# Initialize calculators
pysic_calc1 = Pysic()
pysic_calc2 = Pysic()
ewald_param = estimate_ewald_parameters(1, 'high')
ewald = CoulombSummation(parameters=ewald_param)
pysic_calc1.set_coulomb_summation(ewald)
pysic_calc2.set_coulomb_summation(ewald)

# Add different calculators for the subsystems
hybrid_system.set_primary_calculator(pysic_calc1)
hybrid_system.set_secondary_calculator(pysic_calc2)

#-------------------------------------------------------------------------------
# Define an embedding scheme between the subsystems
# In this case the scheme is mechanical embedding with hydrogen links
parameters = {
    'links': ((0, 1),),
    'CHL': 1,
    'epsilon': ewald_param[3],
}
hybrid_system.set_embedding('MEHL', 'primary', 'secondary', parameters)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hybrid_energy = hybrid_system.get_potential_energy()
#hybrid_system.view_subsystems()

# Calculate the energy of the same setup, but use only one region. In this
# special case these energies should be same.
pysic_calc3 = Pysic()
pysic_calc3.set_coulomb_summation(ewald)
system.set_calculator(pysic_calc3)
real_energy = system.get_potential_energy()

# These energies differ, although the electrostatic connection energy seems to
# be correct. This means that the Either the Ewald sum really is't suitable for
# finite systems or requires better parameter setup
print "Energy with hybrid calculation: " + str(hybrid_energy)
print "Energy with traditional calculation: " + str(real_energy)

hybrid_system.view_subsystems()
hybrid_system.print_potential_energies()
