#! /usr/bin/env python
#
# Simulating a NaCl-H2O solution. The NaCl subsystem employs a GPAW calculator,
# and the H2O subsystem uses the Pysic calculator. The two subsystems are
# electrically bound together.
#
#===============================================================================
from ase import Atoms
from ase.lattice.spacegroup import crystal
from ase.structure import molecule
from pysic import *
import numpy as np
from ase.md.verlet import VelocityVerlet
from ase.constraints import FixInternals, FixBondLengths
from ase import units
from pysic.utility.visualization import AtomEyeViewer
from pysic.utility.make_solution import *
from gpaw import GPAW
from ase.visualize import view
from ase.parallel import rank
from ase.parallel import parprint
import time

#-------------------------------------------------------------------------------
# Create a NaCl lattice used as solute
a = 5.64
nacl = crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
               cellpar=[a, a, a, 90, 90, 90])
nacl.set_pbc(False)
nacl.set_cell((a*2, a*2, a*2))
nacl.center()
nacl.set_tags([0]*len(nacl))

# Create the water molecule used as solvent
h2o = molecule('H2O')
h2o.set_tags([1]*len(h2o))
qh = 0.417
qo = -0.834
bond = 0.9572
angle = 104.52/180.0*np.pi

h2o.set_distance(0, 1, bond)
h2o.set_distance(0, 2, bond)
h2o.set_angle([1, 0, 2], angle)
h2o.set_initial_charges([qo, qh, qh])

# Make the water solution
d = 2
l = 6
nacl_solution, n_solvents = make_solution(nacl, h2o, (d, d, d), (l, l, l), 2, 2)
n_solute = int(len(nacl))

# Add the constraints
constraints = []
for i in range(n_solvents):
    a = n_solute+i*3
    b = n_solute+i*3+1
    c = n_solute+i*3+2
    bond1 = [bond, [a, b]]
    bond2 = [bond, [a, c]]
    angle_indices = [b, a, c]
    c_angle = [angle, angle_indices]
    constraint = FixInternals(nacl_solution, bonds=[bond1, bond2], angles=[c_angle], dihedrals=[], epsilon=1.E-4)
    constraints.append(constraint)

nacl_solution.set_constraint(constraints)

if rank == 0:
    viewer = AtomEyeViewer(nacl_solution, "/home/lauri/Dropbox/SIN")
    #viewer.view()

#-------------------------------------------------------------------------------
# Setup the primary subsystem
hc = HybridCalculator()
primary_subsystem = SubSystem("primary", tag=0, calculator=GPAW(h=0.4, txt=None))
hc.add_subsystem(primary_subsystem)

#-------------------------------------------------------------------------------
# Setup the secondary subsystem

# Assign a pysic calculator to the system. Use the parameters provided by the
# TIP3P model. Disregard the internal forces in a water molecule, that is the
# hydrogens and oxygen within one molecule do not interact.

pysic_calc = Pysic()
kcalPerMoleInEv = 0.0436
A = 582.0E3*kcalPerMoleInEv
B = 595.0*kcalPerMoleInEv
epsilon = 0.0052635
kc = 1.0/(4.0*np.pi*epsilon)
max_cutoff = np.linalg.norm(nacl_solution.get_cell())

# Lennard-Jones potential in two parts
pot1 = Potential('power',symbols=['O', 'O'], parameters=[A, 1, 12], cutoff=max_cutoff)
pot2 = Potential('power', symbols=['O', 'O'], parameters=[-B, 1, 6], cutoff=max_cutoff)

coulomb_pairs = []
for i in range(n_solvents*3):
    for j in range(n_solvents*3):
        if (j>i) and (int(i/3.0) is not int(j/3.0)):
            coulomb_pairs.append([i, j])

# The first potential given to the productpotential defines the targets and cutoff
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
steps = 100
interval = 1
step = 0
start = time.time()
end = time.time()

def print_progress():
    global step, interval, steps, start, end
    if step == 0:
        start = time.time()
        end = time.time()
    else:
        start = end
        end = time.time()
        elapsed = end-start
        print "Time remaining: "+str((steps-step)*elapsed/60.0)+" minutes"

    step += interval
    print str(float(step)/steps*100)+'%'

dyn = VelocityVerlet(nacl_solution, 1*units.fs)

if rank == 0:
    dyn.attach(viewer.save_cfg_frame, interval=1)
    dyn.attach(print_progress, interval=interval)

dyn.run(steps)

if rank == 0:
    viewer.set_colors(hc.get_colors())
    viewer.view_series()
