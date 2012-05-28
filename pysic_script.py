#! /usr/bin/env python
import pysic
import pysic_fortran as pf
import ase
import numpy as np
import debug as d
import math
import time

system = ase.Atoms('CNOH',[[0.01,0.02,0.03], [2.03,-0.02,0.1], [2.3,2.5,2.6], [1.8,2.1,0.02]])
#system = ase.Atoms('HHeLi',[[0.1,0.2,0.3], [2.1,-0.2,0.01], [0,2,0]])
#system = ase.Atoms('HLiHe',[[0.0,0.0,0.0], [2.0,-0.0,0.0], [0,2,0]])
#system = ase.Atoms('H3',[[0.0,0.0,0.0], [2.0,-0.0,0.0], [0,2,0]])
system.set_cell([[30.5,0.5,0.0],
                 [0.0,31.0,0.0],
                 [0.0,0.0,31.0]])
#system.set_charges([2,0,0])
system.set_charges([2,0,0,0])
system.set_pbc([True,True,True])
calc = pysic.Pysic()
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
#pot = pysic.Potential('LJ',cutoff=6.0,symbols=['H','H'],parameters=[1.0,1.0])
#pot = pysic.Potential('Buckingham',cutoff=10.0,symbols=['H','He'],parameters=[1.0,1.0,1.0])
#pot = pysic.Potential('constant',symbols=[['H'],['He'],['Li']],parameters=[1.0])
#pot = pysic.Potential('bond_bend', cutoff=12.5, symbols=['He','H','He'])
pot = pysic.Potential('bond_bend', cutoff=2.2, symbols=[['He','H','Li']])
#pot = pysic.Potential('bond_bend', cutoff=12.5, symbols=['H','H','H'])
pot = pysic.Potential('dihedral', cutoff=13.2, symbols=['H','O','N','C'])
pot.set_parameter_value('k',1.0)
pot.set_parameter_value('theta_0',0.0)
pot.set_cutoff_margin(0.4)

bonds = pysic.BondOrderParameters('tersoff')
bonds.set_cutoff(12.2)
bonds.set_cutoff_margin(0.1)
bonds.set_symbols([['H', 'He', 'Li']])
bonds.set_parameter_value('beta', 1.0)
bonds.set_parameter_value('eta', 3.0)
bonds.set_parameter_value('mu', 1.0)
bonds.set_parameter_value('a', 1.0)
bonds.set_parameter_value('c', 1.0)
bonds.set_parameter_value('d', 1.0)
bonds.set_parameter_value('h', 1.0)

bonds3 = pysic.BondOrderParameters('tersoff')
bonds3.set_cutoff(12.2)
bonds3.set_cutoff_margin(0.1)
bonds3.set_symbols([['H','Li','He']])
bonds3.set_parameter_value('beta', 1.0)
bonds3.set_parameter_value('eta', 3.0)
bonds3.set_parameter_value('mu', 1.0)
bonds3.set_parameter_value('a', 1.0)
bonds3.set_parameter_value('c', 2.0)
bonds3.set_parameter_value('d', 1.0)
bonds3.set_parameter_value('h', 1.0)

bonds2 = pysic.BondOrderParameters('neighbors')
bonds2.set_cutoff(3.2)
bonds2.set_cutoff_margin(0.4)
bonds2.set_symbols([['H', 'H']])

if False:
    bonds = pysic.BondOrderParameters('triplet')
    bonds.set_cutoff(12.2)
    bonds.set_cutoff_margin(0.1)
    bonds.set_symbols([['H','He','Li']])
    bonds3 = pysic.BondOrderParameters('triplet')
    bonds3.set_cutoff(12.2)
    bonds3.set_cutoff_margin(0.1)
    bonds3.set_symbols([['H','Li','He']])

crd = pysic.Coordinator( [bonds,bonds3] )
#pot.set_coordinator(crd)

system.set_calculator(calc)
calc.add_potential(pot)

#print pot
#print system.get_pbc()

print ""
print "initial charges: \n", system.get_charges()
calc.set_core()
t0 = time.time()
print "energy: \n", system.get_potential_energy()
t1 = time.time()
print "forces: \n", system.get_forces()
t2 = time.time()
if False:
    print "numeric forces: \n", np.array( [ calc.get_numerical_energy_gradient(0),  
                                       calc.get_numerical_energy_gradient(1),  
                                       calc.get_numerical_energy_gradient(2) ] )
else:
    print "numeric forces: \n", np.array( [ calc.get_numerical_energy_gradient(0),  
                                       calc.get_numerical_energy_gradient(1),  
                                       calc.get_numerical_energy_gradient(2),  
                                       calc.get_numerical_energy_gradient(3) ] )
t3 = time.time()
print "e-negativities: \n", calc.get_electronegativities(system)
t4 = time.time()
print "numeric e-negativities: \n", np.array( [ calc.get_numerical_electronegativity(0), 
                                               calc.get_numerical_electronegativity(1), 
                                               calc.get_numerical_electronegativity(2) ] )
print "e-negativity differences: \n", calc.get_electronegativity_differences(system)
print ""

if False:
    crd.calculate_bond_order_factors()
    print "bond orders: \n", crd.get_bond_order_factors()
    print "bond order gradients 0: \n", crd.get_bond_order_gradients_of_factor(0)
    numeric_grad_0 = np.array( [ calc.get_numerical_bond_order_gradient(crd,0,0),
                                calc.get_numerical_bond_order_gradient(crd,0,1), 
                                calc.get_numerical_bond_order_gradient(crd,0,2)] )
    print "numeric bond order gradients 0: \n", numeric_grad_0
    print "bond order gradients 1: \n", crd.get_bond_order_gradients_of_factor(1)   
    numeric_grad_1 = np.array( [ calc.get_numerical_bond_order_gradient(crd,1,0),
                                calc.get_numerical_bond_order_gradient(crd,1,1), 
                                calc.get_numerical_bond_order_gradient(crd,1,2)] )
    print "numeric bond order gradients 1: \n", numeric_grad_1
    print "bond order gradients 2: \n", crd.get_bond_order_gradients_of_factor(2) 
    numeric_grad_2 = np.array( [ calc.get_numerical_bond_order_gradient(crd,2,0),
                                calc.get_numerical_bond_order_gradient(crd,2,1), 
                                calc.get_numerical_bond_order_gradient(crd,2,2)] )
    print "numeric bond order gradients 2: \n", numeric_grad_2

    print "bond order gradient sum: \n", -crd.get_bond_order_gradients_of_factor(0)-\
        crd.get_bond_order_gradients_of_factor(1)-crd.get_bond_order_gradients_of_factor(2)

    print "numeric bond order gradient sum: \n", -numeric_grad_0-numeric_grad_1-numeric_grad_2

print "timing (E, F, chi) : ", t1-t0, t2-t1, t4-t3


print calc.get_neighbor_lists().get_neighbors(0)
print calc.get_neighbor_lists().get_neighbors(1)
print calc.get_neighbor_lists().get_neighbors(2)

quit()
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
    """


"""
    
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


"""
    import ase.calculators.neighborlist as nbl
    import pysic_fortran as pf
    import interface as s
    import time
    import pysic_utility as pu
    
    n_atoms = 75712
    system = s.system()[0:n_atoms]
    #print system.get_cell()
    #system.set_cell(atoms.get_cell())
    #system.set_cell([[41,23,2],[5,45,-15],[0,13,37]])
    #system.set_cell([20,20,20])
    system.set_pbc([True,True,True])
    
    calc = pysic.Pysic()
    system.set_calculator(calc)
    
    
    masses = np.array( system.get_masses() )
    charges = np.array( system.get_charges() )
    positions = np.array( system.get_positions() ).transpose()
    momenta = np.array( system.get_momenta() ).transpose()
    tags = np.array( system.get_tags() )
    elements = system.get_chemical_symbols()
    
    for index in range(len(elements)):
    elements[index] = pu.str2ints(elements[index],2)
    
    elements = np.array( elements ).transpose()
    
    pf.pysic_interface.create_atoms(masses,
    charges,
    positions,
    momenta,
    tags,
    elements)
    calc.update_core_supercell()
    pf.pysic_interface.distribute_mpi(system.get_number_of_atoms())
    
    cutoffs = [6.0]*n_atoms
    
    print ""
    print "building ASE list"
    t1 = time.time()
    #nbors1 = nbl.NeighborList([4.0]*n_atoms,skin=1.0,bothways=True,self_interaction=False)
    #nbors1.build(system)
    t2 = time.time()
    print "done ", t2-t1
    print "building Fortran list"
    t3 = time.time()
    pf.pysic_interface.generate_neighbor_lists(cutoffs)
    t4 = time.time()
    print "done ", t4-t3
    
    """    

"""
    import random
    random.seed(5801)
    for i in range(n_atoms):
    rnd = random.randint(0,n_atoms/5)
    ase = sorted(nbors1.get_neighbors(i)[0])
    n_nbs = pf.pysic_interface.get_number_of_neighbors_of_atom(i)
    fortran = sorted(pf.pysic_interface.get_neighbor_list_of_atom(i,n_nbs)[0])
    
    if rnd == 0 or ase != fortran:
    print "atom ",i
    print "ASE:    ", ase
    print "Fortran:", fortran
    """

#print """The run times were (in seconds)
#  ASE:     {a}
#  Python:  {p}
#  Fortran: {f}""".format(a=str(t2-t1),p=str(t3-t2),f=str(t4-t3))





"""
import ase.calculators.neighborlist as nbl
import pysic_fortran as pf
import interface as s
import time

system = s.system()
system.set_pbc([True,True,True])

calc = pysic.Pysic()
pot = pysic.Potential('LJ')
pot.set_symbols([['Si','Si'],
                 ['Si','O'],
                 ['Hf','Hf'],
                 ['Hf','O'],
                 ['O','O']])
pot.set_cutoff(8.0)
pot.set_cutoff_margin(2.0)
pot.set_parameter_value('epsilon',1.0)
pot.set_parameter_value('sigma',3.0)
calc.add_potential(pot)
system.set_calculator(calc)

print "set up"
t0 = time.time()
calc.set_core()
print "neighbor lists ready"
t1 = time.time()
print system.get_potential_energy()
t2 = time.time()
print "time: ", t1-t0, t2-t1, t2-t0

"""
