#! /usr/bin/env python
#
# This script tests that the hydrogen links in the MEHL scheme are working. We
# visually inspect the subsystems and calculate the energy of a special system
# with hybrid subsystem division and with traditional ways and checks that they
# give the same results.
#
#===============================================================================
from ase import Atoms
from pysic import *
from pysic.utility import visualization

#-------------------------------------------------------------------------------
# Prepare the system
h2 = Atoms('H2', [(0, 0, 0), (1, 1, 1)])
h2.set_cell((2, 2, 2))
h2.set_pbc(True)

# Use AtomEye to make sure the structure is correct:
viewer = AtomEyeViewer(h2)
viewer.view()

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hc = HybridCalculator()

# Define calculator for the subsystems
calc = Pysic()
potential = Potential('LJ', cutoff=4.0, symbols=['H', 'H'], parameters=[0.1, 2.5])
calc.set_potentials(potential)

# Define QM/MM regions. You can get the indices by e.g. examining the the
# structure in ASEs viewer.
primary_subsystem = SubSystem("primary", indices=0, calculator=calc)
hc.add_subsystem(primary_subsystem)

secondary_subsystem = SubSystem("secondary", special_set="remaining", calculator=calc)
hc.add_subsystem(secondary_subsystem)

#-------------------------------------------------------------------------------
# Define a binding between the subsystems
binding = Binding("primary", "secondary")
binding.set_hydrogen_links((0, 1), 1, interaction_correction=True)
hc.add_binding(binding)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hybrid_energy = hc.get_potential_energy(h2)

# Calculate the energy of the same setup, but use only one region. In this
# special case these energies should be the same.
real_energy = calc.get_potential_energy(h2)

# When periodic boundary conditions are on, the link atoms in the primary
# system will interact with each other. This energy is already included in the
# secondary system and should be corrected somehow.
print "Energy with hybrid calculation: " + str(hybrid_energy)
print "Energy with traditional calculation: " + str(real_energy)
hc.print_energies()

