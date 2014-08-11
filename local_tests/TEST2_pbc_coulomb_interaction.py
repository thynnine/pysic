#! /usr/bin/env python
#
# Tests that the Ewald summation is done correctly for hybrid systems with
# charges and pbc.
#
#===============================================================================
from ase import Atoms, Atom
from pysic import *
from pysic.utility.atomeyeviewer import AtomEyeViewer
from pysic.interactions.coulomb import *

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
system.set_pbc(True)

# Use AtomEye to make sure the structure is correct:
viewer = AtomEyeViewer(system)
viewer.view()

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hc = HybridCalculator()

# Initialize calculator for subsystems
calc = Pysic()
ewald_param = estimate_ewald_parameters(1, 'high')
ewald = CoulombSummation(parameters=ewald_param)
calc.set_coulomb_summation(ewald)

# Define subsystems
hc.add_subsystem(SubSystem("primary", indices=0, calculator=calc))
hc.add_subsystem(SubSystem("secondary", indices="remaining", calculator=calc))

# Define an embedding scheme between the subsystems
# In this case the scheme is mechanical embedding with hydrogen links
interaction = Interaction("primary", "secondary")
interaction.add_hydrogen_links((0, 1), 1)
interaction.enable_coulomb_potential(sigma=ewald_param[2], k_cutoff=ewald_param[1], real_cutoff=ewald_param[0])
hc.add_interaction(interaction)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hybrid_energy = hc.get_potential_energy(system)
hybrid_forces = hc.get_forces(system)

# Calculate the energy of the same setup, but use only one region. In this
# special case these energies should be same.
real_energy = calc.get_potential_energy(system)
real_forces = calc.get_forces(system)

print "Energy with hybrid calculation: " + str(hybrid_energy)
print "Energy with traditional calculation: " + str(real_energy)
print "Forces with hybrid calculation: " + str(hybrid_forces)
print "Forces with traditional calculation: " + str(real_forces)

#hc.view_subsystems()
#hc.print_energy_summary()
#hc.print_force_summary()
