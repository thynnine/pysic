#! /usr/bin/env python
import pysic
import pysic_fortran as pf
import ase
import numpy as np
import debug as d
import math
"""
system = ase.Atoms('H2',[[0,0,0],[2,0,0]])
system.set_cell([[0.5,0.5,0.0],
                 [0.0,1.0,0.0],
                 [0.0,0.0,1.0]])
system.set_charges([2,0])
calc = pysic.Pysic()
#pot = pysic.Potential('LJ',cutoff=10.0,symbols=['H','H'],parameters=[1.0,1.0])
pot = pysic.Potential('exponential',cutoff=10.0,symbols=['H','H'])
pot.set_parameter_value('epsilon',1.0)
pot.set_parameter_value('zeta',0.2)
pot.set_parameter_value('xi1',1.0)
pot.set_parameter_value('xi2',1.0)
pot.set_parameter_value('Rmax1',-0.1)
pot.set_parameter_value('Rmax2',-0.1)
pot.set_parameter_value('Rmin1',0.2)
pot.set_parameter_value('Rmin2',0.2)
pot.set_parameter_value('Qmax1',4.0)
pot.set_parameter_value('Qmax2',4.0)
pot.set_parameter_value('Qmin1',-4.0)
pot.set_parameter_value('Qmin2',-4.0)

system.set_calculator(calc)
calc.add_potential(pot)

print ""
print "initial charges: \n", system.get_charges()
print "energy: \n", system.get_potential_energy()
print "forces: \n", system.get_forces()
print "numeric forces: \n", np.array( [ calc.get_numerical_energy_gradient(0),  calc.get_numerical_energy_gradient(1) ] )
print "e-negativities: \n", calc.get_electronegativities(system)
print "numeric e-negativities: \n", np.array( [ calc.get_numerical_electronegativity(0), calc.get_numerical_electronegativity(1) ] )
print "e-negativity differences: \n", calc.get_electronegativity_differences(system)
print ""

d.bp()


chrg = pysic.ChargeRelaxation('dynamic')
chrg.set_calculator(calc)
chrg.set_parameter_value('timestep',1.0)
chrg.set_parameter_value('n_steps',1000)
chrg.set_parameter_value('friction',0.01)
chrg.set_parameter_value('tolerance',0.000001)
chrg.set_parameter_value('inertia',1.0)
chrg.set_atoms(system)
chrg.charge_relaxation()

print ""
print "equilibrated charges: \n", system.get_charges()
print "energy: \n", system.get_potential_energy()
print "forces: \n", system.get_forces()
print "numeric forces: \n", np.array( [ calc.get_numerical_energy_gradient(0),  calc.get_numerical_energy_gradient(1) ] )
print "e-negativities: \n", calc.get_electronegativities(system)
print "numeric e-negativities: \n", np.array( [ calc.get_numerical_electronegativity(0), calc.get_numerical_electronegativity(1) ] )
print "e-negativity differences: \n", calc.get_electronegativity_differences(system)
print ""

#pysic.Pysic.core.view_fortran()
print system.get_cell()
print system.get_reciprocal_cell()
print np.linalg.inv(system.get_cell()).transpose()
"""


"""
system = ase.Atoms('Na4Cl4',[[0.5,0.5,0.5],
                             [1.5,1.5,0.5],
                             [1.5,0.5,1.5],
                             [0.5,1.5,1.5],
                             [1.5,0.5,0.5],
                             [0.5,1.5,0.5],
                             [0.5,0.5,1.5],
                             [1.5,1.5,1.5]])

system.set_charges([2,1,1,1,-1,-1,-1,-1])
system.set_cell([2,2,2])
system.set_pbc([True,True,True])

calc = pysic.Pysic()
system.set_calculator(calc)

ewald = pysic.CoulombSummation()
ewald.set_parameter_value('real_cutoff',6.0)
ewald.set_parameter_value('k_cutoff',6.0)
ewald.set_parameter_value('sigma',0.3)
ewald.set_parameter_value('epsilon',1.0/(math.pi*4.0))
calc.set_coulomb_summation(ewald)


calc.set_core()


print "initial charges: \n", system.get_charges()
print "energy: \n", system.get_potential_energy()
print "forces: \n", system.get_forces()[0:2]
print "numeric forces: \n", np.array( [ calc.get_numerical_energy_gradient(0),  calc.get_numerical_energy_gradient(1) ] )
print "e-negativities: \n", calc.get_electronegativities(system)[0:2]
print "numeric e-negativities: \n", np.array( [ calc.get_numerical_electronegativity(0), calc.get_numerical_electronegativity(1) ] )
print "e-negativity differences: \n", calc.get_electronegativity_differences(system)
print ""
"""

"""
for sigma in [0.1, 0.3, 0.9]:
    for real in [2.0, 6.0, 10.0]:
        for reci in [2.0, 4.0, 6.0, 8.0]:
            
            ewald = pysic.CoulombSummation()
            ewald.set_parameter_value('real_cutoff',real)
            ewald.set_parameter_value('k_cutoff',reci)
            ewald.set_parameter_value('sigma',sigma)
            ewald.set_parameter_value('epsilon',1.0/(math.pi*4.0))
            calc.set_coulomb_summation(ewald)
            calc.set_core()
            #ewald_ene = pf.pysic_interface.get_ewald_energy(real,[reci,reci,reci],sigma,1.0/(math.pi*4.0))
            ewald_ene = system.get_potential_energy()
            print "{sg} {di} {re} madelung constant = {mc} (error = {er})".format(sg=str(sigma),
                                                                              di=str(real),
                                                                              re=str(reci),
                                                                              mc=str(ewald_ene/4),
                                                                              er=str(ewald_ene/4+1.747564594633))
"""            

"""
#pot = pysic.Potential('force',parameters=[1.0,2.0,3.0],symbols=[['Na'],['Cl']])
#calc.add_potential(pot)

system.set_charges([1,1,1,1,-1,-1,-1,-1])
ewald_ene = system.get_potential_energy()
print "accuracy: \n", ewald_ene/4 +1.747564594633

#print ewald_ene
system[0].x += 0.1
system[0].y += 0.2
system[0].z += 0.3
#ewald_ene = system.get_potential_energy()
#print ewald_ene

print "forces: \n", system.get_forces()[0:2]
print "numeric forces: \n", np.array( [ calc.get_numerical_energy_gradient(0),  calc.get_numerical_energy_gradient(1) ] )

#pysic.Pysic.core.view_fortran()
"""

import ase.calculators.neighborlist as nbl
import pysic_utility as pu
import interface as s
import time

atoms = s.system()[0:8000]
atoms.pbc = [True,True,True]
#atoms.set_cell([[30,20,2],[5,35,-15],[0,13,30]])


"""
atoms = ase.Atoms('Na4Cl4',[[0.5,0.5,0.5],
                             [1.5,1.5,0.5],
                             [1.5,0.5,1.5],
                             [0.5,1.5,1.5],
                             [1.5,0.5,0.5],
                             [0.5,1.5,0.5],
                             [0.5,0.5,1.5],
                             [1.5,1.5,1.5]])

atoms.set_cell([2,2,2])
atoms.set_pbc([True,True,True])
"""

natoms = atoms.get_number_of_atoms()
print natoms, "atoms"
nbors1 = nbl.NeighborList([2.0]*natoms,skin=0.5,bothways=True,self_interaction=False)
t1 = time.time()
print "building list 1"
nbors1.build(atoms)
t2 = time.time()
print "done", t2-t1

nbors2 = pu.FastNeighborList([4.0]*natoms,skin=0.5)

t3 = time.time()
print "building list 2"
nbors2.build(atoms)
t4 = time.time()
print "done", t4-t3


subs = nbors2.cell.subcells
print "splits \n"
print len(subs), len(subs[0]), len(subs[0][0])
print "atom 0 subcell \n"
print nbors2.cell.atom_subcell_indices[0]
print "subcell 0,0,0 matrix\n"
print nbors2.cell.subcells[0][0][0].matrix
print "subcell 0,0,0 atoms\n"
print nbors2.cell.subcells[0][0][0].atoms
#print "subcell 0,0,0 offsets\n"
#print nbors2.cell.subcells[0][0][0].offsets
#print "subcell 0,0,0 include\n"
#print nbors2.cell.subcells[0][0][0].include_nbor
print "subcell 0,0,0 neighbor cells, indices, include, offsets \n"
for i in range(-1,2):
    for j in range(-1,2):
        for k in range(-1,2):
            print i,j,k, \
                nbors2.cell.subcells[0][0][0].neighbors[i][j][k].subcell_indices, \
                nbors2.cell.subcells[0][0][0].include_nbor[i][j][k], \
                nbors2.cell.subcells[0][0][0].offsets[i][j][k]


print "atom 0 neighbors\n"
print sorted(nbors2.get_neighbors(0)[0])
print sorted(nbors1.get_neighbors(0)[0])
ds = []
"""
for n in sorted(nbors2.get_neighbors(0)[0]):
    i = 0
    off = [0,0,0]
    for index in range(len(nbors2.get_neighbors(0)[0])):
        #print i, n, index, nbors2.get_neighbors(0)[0][index]
        if nbors2.get_neighbors(0)[0][index] == n:
            i = index
            off = np.array(nbors2.get_neighbors(0)[1][i])
    print n,i,off,nbors2.cell.get_distance(atoms[0],atoms[n],off)
"""
print "\n"
#print sorted(nbors2.get_neighbors(224)[0])
#print sorted(nbors1.get_neighbors(224)[0])
print "last atom neighbors\n"
print sorted(nbors2.get_neighbors(natoms-1)[0])
print sorted(nbors1.get_neighbors(natoms-1)[0])

print "runtimes", t2-t1, t4-t3









