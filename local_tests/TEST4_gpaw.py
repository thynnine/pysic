#! /usr/bin/env python
#
# This test uses GPAW for the primary system and pysic for the secondary.
# The embedding scheme is set at mechanical embedding with hydrogen links.
#
#===============================================================================
from ase import Atoms
from pysic import *
from gpaw import GPAW
from ase.io import *

#-------------------------------------------------------------------------------
# Prepare the system
h2 = Atoms('He3', [(5, 5, 5), (6, 6, 6), (6, 6, 7)], charges=(0, 0, 1))
b = 6.0
h2.set_cell((b, b, b))
h2.center()

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hc = HybridCalculator()

# Initialize calculator for subsystems
pysic_calc = Pysic()
potential1 = Potential('LJ', cutoff=4.0, symbols=['He', 'He'], parameters=[0.1, 1])
potential2 = Potential('LJ', cutoff=4.0, symbols=['He', 'H'], parameters=[0.1, 1])
pysic_calc.set_potentials([potential1, potential2])
gpaw_calc = GPAW(h=0.4, txt=None)

# Define subsystems
PS = SubSystem("primary", indices=[0, 1], calculator=gpaw_calc)
PS.enable_charge_calculation(division="Bader", source="all-electron", gridrefinement=4)
hc.add_subsystem(PS)
hc.add_subsystem(SubSystem("secondary", indices="remaining", calculator=pysic_calc))

# Define Interaction
interaction = Interaction("primary", "secondary")
interaction.enable_coulomb_potential()
interaction.add_hydrogen_links((1, 2), 0.5)
interaction.set_link_atom_correction(True)
interaction.set_potentials(potential1)
hc.add_interaction(interaction)

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hc.get_potential_energy(h2)
hc.get_forces()
hc.print_interaction_charge_summary()
hc.print_energy_summary()
hc.print_force_summary()
hc.print_time_summary()
hc.view_subsystems()
