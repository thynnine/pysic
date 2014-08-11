#! /usr/bin/env python
#
# Test that link atoms with differet CHL parameters work.
#
#===============================================================================
from ase import Atoms
from pysic import *
from ase.md.verlet import VelocityVerlet
from ase import units

#-------------------------------------------------------------------------------
# Prepare the system
h2 = Atoms('C2', [(0, 0, 0), (3, 0, 0)])
h2.set_cell((10, 10, 10))
h2.center()

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hc = HybridCalculator()

# Define calculator for the subsystems
calc = Pysic()
potential1 = Potential('LJ', cutoff=10, symbols=['C', 'H'], parameters=[0.1, 1])
potential2 = Potential('LJ', cutoff=10, symbols=['H', 'H'], parameters=[0.1, 1])
#potential3 = Potential('LJ', cutoff=10, symbols=['C', 'C'], parameters=[0.1, 3])
calc.set_potentials([potential1])  # , potential2])

# Define QM/MM regions. You can get the indices by e.g. examining the the
# structure in ASEs viewer.
primary_subsystem = SubSystem("primary", indices=(0), calculator=calc)
hc.add_subsystem(primary_subsystem)

secondary_subsystem = SubSystem("secondary", indices=(1), calculator=calc)
hc.add_subsystem(secondary_subsystem)

#-------------------------------------------------------------------------------
# Define an interaction between the subsystems
interaction = Interaction("primary", "secondary")
interaction.add_hydrogen_links((0, 1), 0.3)
interaction.add_hydrogen_links((0, 1), 0.4)
interaction.add_hydrogen_links((0, 1), 0.5)
interaction.add_hydrogen_links((0, 1), 0.6)
interaction.add_hydrogen_links((0, 1), 0.7)
interaction.set_link_atom_correction(False)
hc.add_interaction(interaction)

#-------------------------------------------------------------------------------
# Calculate the potential energy and forces of the hybrid qm/mm system.
h2.set_calculator(hc)
h2.get_forces()
h2.get_potential_energy()

#viewer = AtomEyeViewer(h2)
dyn = VelocityVerlet(h2, 10*units.fs)  # 5 fs time step.
#dyn.attach(viewer.save_cfg_frame)
dyn.run(20)
#viewer.view_series()

#hc.print_energy_summary()
#hc.print_force_summary()
#hc.print_time_summary()
hc.view_subsystems()
