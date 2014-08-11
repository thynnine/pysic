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
from pysic.utility.atomeyeviewer import AtomEyeViewer

#-------------------------------------------------------------------------------
# Prepare the system
h2 = Atoms('He2', [(4, 4, 4), (6, 6, 6)])
h2.set_cell((10, 10, 10))
h2.set_pbc(False)

# Use AtomEye to make sure the structure is correct:
viewer = AtomEyeViewer(h2)
#viewer.view()

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hc = HybridCalculator()

# Define calculator for the subsystems
calc = Pysic()
potential1 = Potential('LJ', cutoff=10, symbols=['He', 'He'], parameters=[0.1, 2.5])
potential2 = Potential('LJ', cutoff=10, symbols=['He', 'H'], parameters=[0.1, 2.5])
calc.set_potentials([potential1, potential2])

# Define QM/MM regions. You can get the indices by e.g. examining the the
# structure in ASEs viewer.
primary_subsystem = SubSystem("primary", indices=0, calculator=calc)
primary_subsystem.enable_cell_optimization(1)
hc.add_subsystem(primary_subsystem)

secondary_subsystem = SubSystem("secondary", indices=1, calculator=calc)
hc.add_subsystem(secondary_subsystem)

#-------------------------------------------------------------------------------
# Define an interaction between the subsystems
interaction = Interaction("primary", "secondary")
interaction.add_hydrogen_links((0, 1), 1)
interaction.set_link_atom_correction(True)
interaction.set_potentials(potential1)
hc.add_interaction(interaction)

#-------------------------------------------------------------------------------
# Calculate the potential energy and forces of the hybrid qm/mm system.
hybrid_energy = hc.get_potential_energy(h2)
hybrid_forces = hc.get_forces(h2)

# Calculate the potential energy and forces of the same setup, but use only one
# region. In this special case these energies and forces should be the same.
real_energy = calc.get_potential_energy(h2)
real_forces = calc.get_forces(h2)

print "Energy with hybrid calculation: " + str(hybrid_energy)
print "Energy with traditional calculation: " + str(real_energy)
print "Forces with hybrid calculation: " + str(hybrid_forces)
print "Forces with traditional calculation: " + str(real_forces)

#hc.print_energy_summary()
#hc.print_force_summary()
#hc.print_time_summary()
#hc.view_subsystems()
