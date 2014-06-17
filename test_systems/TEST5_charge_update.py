#! /usr/bin/env python
#
# Test that the charges in the primary system are calculated from the
# pseudo-density provided by the DFT-calculator.
#
#===============================================================================
from ase import Atoms
from pysic import *
from pysic.utility import visualization
from gpaw import GPAW

#-------------------------------------------------------------------------------
# Prepare the system
h2 = Atoms('H2', [(0, 0, 0), (0.7, 0, 0)], charges=(0, 1))
b = 3
h2.set_cell((b, b, b))
h2.center()

# Use AtomEye to view the structure:
viewer = AtomEyeViewer(h2)
viewer.view()

#-------------------------------------------------------------------------------
# Setup a hybrid calculation environment
hc = HybridCalculator()

# Initialize calculators
pysic_calc = Pysic()
potential = Potential('LJ', cutoff=4.0, symbols=['H', 'H'], parameters=[0.1, 2.5])
pysic_calc.add_potential(potential)
gpaw_calc = GPAW(h=b/8.0, txt=None)

# Define subsystems
hc.add_subsystem("primary", indices=0, calculator=gpaw_calc)
hc.add_subsystem("secondary", special_set="remaining", calculator=pysic_calc)

# Define an embedding scheme between the subsystems
binding = hc.add_binding("primary", "secondary")
binding.set_hydrogen_links((0, 1), 1)
binding.set_electrostatic_binding()

#-------------------------------------------------------------------------------
# Calculate the potential energy of the hybrid qm/mm system.
hc.get_potential_energy(h2)
hc.print_energies()
