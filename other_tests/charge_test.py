#! /usr/bin/env python
#
# Simulating a NaCl-H2O solution. The NaCl subsystem employs a GPAW calculator,
# and the H2O subsystem uses the Pysic calculator. The two subsystems are
# electrically bound together.
#
#===============================================================================
from pysic import *
import numpy as np
from ase import Atoms
from pysic.utility.visualization import AtomEyeViewer
from pysic.utility.make_solution import *
from gpaw import GPAW, mpi
rank = mpi.world.rank
n_ranks = mpi.world.size

#-------------------------------------------------------------------------------
# Create a hydrogen solute
a = 5.64
solute = Atoms('Na4Cl4', positions=[(0, 0, a), (a, 0, 0), (0, a, 0), (a, a, a), (a, a, 0), (a, 0, a), (0, a, a), (0, 0, 0)])
solute.center()
solute.set_tags([0] * len(solute))
#solute.rotate("z", np.pi/3, center="COP")

# Create the hydrogen used as solvent
solvent = Atoms('H')
solvent.set_tags([1])
qh = 1
solvent.set_initial_charges([qh])

# Make the solution
d = 2
l = 9
nacl_solution, n_solvents = make_solution(solute, solvent, (d, d, d), (l, l, l), 1, 2)
n_solute = int(len(solute))

if rank == 0:
    viewer = AtomEyeViewer(nacl_solution, "/home/lauri/Dropbox/SIN", "nacl_solution27.6")
    viewer.view()

#-------------------------------------------------------------------------------
# Setup the primary subsystem. A GPAW calculator is used.
hc = HybridCalculator(record_time_usage=True)
gpaw_calc = GPAW(h=0.25, xc="PBE", convergence={'eigenstates': 1E-8}, txt=None)
primary_subsystem = SubSystem("primary", tag=0, calculator=gpaw_calc)
#primary_subsystem.enable_cell_size_optimization(3.5)
hc.add_subsystem(primary_subsystem)

#-------------------------------------------------------------------------------
# Setup the secondary subsystem. The water is modeled with Pysic. The water is modeled according to the TIP3P
# model. We disregard the internal forces in a water molecule, that is the
# hydrogens and oxygen within one molecule do not interact.

pysic_calc = Pysic()
kcalPerMoleInEv = 0.0436
A = 582.0E3 * kcalPerMoleInEv
B = 595.0 * kcalPerMoleInEv
epsilon = 0.0052635
kc = 1.0 / (4.0*np.pi*epsilon)
max_cutoff = np.linalg.norm(nacl_solution.get_cell())
real_cutoff = 6

# Lennard-Jones potential in two parts
pot1 = Potential('power', symbols=['O', 'O'], parameters=[A, 1, 12], cutoff=max_cutoff)
pot2 = Potential('power', symbols=['O', 'O'], parameters=[-B, 1, 6], cutoff=max_cutoff)

coulomb_pairs = []
for i in range(n_solvents*3):
    for j in range(n_solvents*3):
        if (j > i) and (int(i/3.0) is not int(j/3.0)):
            coulomb_pairs.append([i, j])

# Coulomb potential in two parts, combined with ProductPotential. The first
# potential given to the productpotential defines the targets and cutoff
coul1 = Potential('power', indices=coulomb_pairs, parameters=[1, 1, 1], cutoff=max_cutoff)
coul2 = Potential('charge_pair', parameters=[kc, 1, 1])
pot3 = ProductPotential([coul1, coul2])

potential_list = [pot1, pot2, pot3]
pysic_calc.add_potential(potential_list)

secondary_subsystem = SubSystem("secondary", tag=1, calculator=pysic_calc)
hc.add_subsystem(secondary_subsystem)

#-------------------------------------------------------------------------------
# Define interactions between subsystems
binding = Binding("primary", "secondary", coulomb_interaction=True)
hc.add_binding(binding)

#-------------------------------------------------------------------------------
# Run dynamics
nacl_solution.set_calculator(hc)
nacl_solution.get_potential_energy()

if rank == 0:
    hc.print_charge_summary()
    hc.print_summary()
