#! /usr/bin/env python
#
# Tests that the Coulomb energy calculation is done correctly for hybrid systems
# with charges
#
#===============================================================================
from ase import Atoms, Atom
from pysic import Pysic, Potential, HybridCalculator, ProductPotential
from pysic.interactions.coulomb import *
from ase.visualize import view
import numpy as np
import copy

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
hc = HybridCalculator(system)

# Define QM/MM regions
hc.set_primary_system([0])
hc.set_secondary_system(special_set='remaining')

# Initialize calculator
calc = Pysic()
epsilon = 0.0052635
kc = 1.0/(4.0*np.pi*epsilon)
coul1 = Potential('power', symbols=[['H', 'H']], parameters=[1, 1, 1], cutoff=20)
coul2 = Potential('charge_pair', parameters=[kc, 1, 1])
coulomb_potential = ProductPotential([coul1, coul2])
calc.add_potential(coulomb_potential)

# Add different calculators for the subsystems
hc.set_primary_calculator(copy.copy(calc))
hc.set_secondary_calculator(copy.copy(calc))

#-------------------------------------------------------------------------------
# Define an embedding scheme between the subsystems
# In this case the scheme is mechanical embedding with hydrogen links
parameters = {'links': ((0, 1),), 'CHL': 1}
hc.set_embedding('MEHL', 'primary', 'secondary', parameters)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hybrid_energy = hc.get_potential_energy()

# Calculate the energy of the same setup, but use only one region. In this
# special case these energies should be same.
real_system = system.copy()
real_system.set_calculator(copy.copy(calc))
real_energy = real_system.get_potential_energy()

print "Energy with hybrid calculation: " + str(hybrid_energy)
print "Energy with traditional calculation: " + str(real_energy)

hc.view_subsystems()
hc.print_potential_energies()
