#! /usr/bin/env python
#
# Tests that the Coulomb energy calculation is done correctly for hybrid systems
# with charges
#
#===============================================================================
from ase import Atoms, Atom
from pysic import *
from pysic.utility import visualization
import numpy as np
import copy
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

# Use AtomEye to view the structure:
viewer = AtomEyeViewer(system)
viewer.view()

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hc = HybridCalculator()

# Initialize calculator for subsystems
calc = Pysic()
epsilon = 0.00552635
kc = 1.0/(4.0*np.pi*epsilon)
coul1 = Potential('power', symbols=[['H', 'H']], parameters=[1, 1, 1], cutoff=20)
coul2 = Potential('charge_pair', parameters=[kc, 1, 1])
coulomb_potential = ProductPotential([coul1, coul2])
calc.add_potential(coulomb_potential)

# Define subsystems
hc.add_subsystem(SubSystem("primary", indices=0, calculator=calc))
hc.add_subsystem(SubSystem("secondary", special_set="remaining", calculator=calc))

# Define an embedding scheme between the subsystems
# In this case the scheme is mechanical embedding with hydrogen links
binding = Binding("primary", "secondary", coulomb_interaction=True)
binding.set_hydrogen_links((0, 1), 1)
hc.add_binding(binding)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hybrid_energy = hc.get_potential_energy(system)

# Calculate the energy of the same setup, but use only one region. In this
# special case these energies should be same.
real_energy = calc.get_potential_energy(system)

print "Energy with hybrid calculation: " + str(hybrid_energy)
print "Energy with traditional calculation: " + str(real_energy)

#hc.view_subsystems()
hc.print_energies()
