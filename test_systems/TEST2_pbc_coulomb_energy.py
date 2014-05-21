#! /usr/bin/env python
#
# Tests that the Ewald summation is done correctly for hybrid systems with
# charges and pbc.
#
#===============================================================================
from ase import Atoms, Atom
from pysic import Pysic, Potential, HybridSystem
from pysic.interactions.coulomb import *
from ase.visualize import view

#-------------------------------------------------------------------------------
# Prepare the system
h1 = Atom('H', (0, 0, 0), charge=-1.0)
h2 = Atom('H', (0.32, 0.845, 0.12), charge=2.4)
h3 = Atom('H', (1, 1.78, 1.99), charge=-0.51)
system = Atoms()
system.append(h1)
system.append(h2)
system.append(h3)
system.set_cell((2, 2, 2))

# Set periodic boundary conditions
system.set_pbc(True)

# Use ASEs built in viewer to make sure the structure is correct:
view(system)

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hybrid_system = HybridSystem()
hybrid_system.set_system(system)

# Define QM/MM regions
hybrid_system.set_primary_system([0])
hybrid_system.set_secondary_system(special_set='remaining')
print hybrid_system.get_subsystem_indices('primary')
print hybrid_system.get_subsystem_indices('secondary')

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
    'k_cutoff': ewald_param[1],
    'real_cutoff': ewald_param[0],
    'sigma': ewald_param[2]
}
hybrid_system.set_embedding('MEHL', 'primary', 'secondary', parameters)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hybrid_energy = hybrid_system.get_potential_energy()

# Calculate the energy of the same setup, but use only one region. In this
# special case these energies should be same.
pysic_calc3 = Pysic()
pysic_calc3.set_coulomb_summation(ewald)
system.set_calculator(pysic_calc3)
real_energy = system.get_potential_energy()

print "Energy with hybrid calculation: " + str(hybrid_energy)
print "Energy with traditional calculation: " + str(real_energy)

hybrid_system.view_subsystems()
hybrid_system.print_potential_energies()
