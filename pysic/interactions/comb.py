#  COMB  POTENTIAL implementation - COMB10 parameters
#
# Teemu Hynninen, Tiziana Musso

import pysic
from pysic.core import *
from pysic.utility.error import InvalidPotentialError
from pysic.interactions.local import Potential, ProductPotential
from pysic.interactions.compound import CompoundPotential
from pysic.interactions.bondorder import Coordinator, BondOrderParameters
from pysic.interactions.coulomb import CoulombSummation, estimate_ewald_parameters
import copy
import math

class CombPotential(CompoundPotential):
    """Class representing a COMB potential.
    
    Use::
    
     import pysic
     import pysic.interactions.comb as comb_sio
    
     calc = pysic.Pysic()
     pot = comb_sio.CombPotential()     
     pot.set_calculator(calc, True)
     
    Since the potential may also define CoulombSummation (which it by default
    does), you should always pass the calculator to the potential in addition
    to passing the potential to the calculator.
    
    The potential can be modified through the exclusion list. This
    defines which chuncks of the potential are incorporated at any time.
    For instance::
    
     pot.exclude('si_self')
     calc.set_potentials(pot)
     
    will remove the self energy contribution from Si atoms. If you modify
    the potential, you need to confirm the changes by passing the potential
    to the calculator. This is because the potential is built in the calculator
    as it is passed. Therefore any changes made after that to the potential will not
    propagate to the calculator. 
    """

    def __init__(self,excludes=['direct_coulomb','long_coulomb']):
        super(CombPotential,self).__init__(n_params=0,
                                           n_targets=3)
        self.set_potential_type('comb')
        self.description = "Comb trial."
        self.excludes = excludes
        self.calc = None
        self.possible_excludes = ['siosi_bend', # 0
                                  'osio_bend', # 1
                                  'sisi_attractive', # 2
                                  'oo_attractive', # 3
                                  'sio_attractive', # 4
                                  'sisi_bondorder', # 5
                                  'oo_bondorder', # 6
                                  'sio_bondorder', # 7
                                  'sisi_repulsive', # 8
                                  'oo_repulsive', # 9
                                  'sio_repulsive', # 10
                                  'si_self', # 11
                                  'o_self', # 12
                                  'sisi_penalty', #13
                                  'oo_penalty', #14
                                  'sio_penalty', #15
                                  'sisi_coulomb', #16
                                  'oo_coulomb', #17
                                  'sio_coulomb', #18
                                  'direct_coulomb', #19
                                  'long_coulomb' #20
                                  ] 

    def set_calculator(self, calc, reciprocal = False):
        """Attaches the calculator to the potential. This is needed so that a CoulombSummation
        can be set. If the argument ``reciprocal=True`` is given, also the command
        ``calc.set_potentials(self)`` is run.
        """
        self.calc = calc
        if reciprocal:
            calc.set_potentials(self)

    def check_exclude(self, ex):
        """Returns True, if ``ex`` is a valid keyword of the list of components
        for this potential.
        """
        if ex not in self.possible_excludes:
            print " **** unknown keyword ****", ex

    def toggle_exclude(self, ex):
        """Changes the state of the given component from included to excluded, or from excluded to included.
        """
        self.check_exclude(ex)
        if ex in self.excludes:
            self.excludes.remove(ex)
        else:
            self.excludes.append(ex)
            #self.define_elements()
    
    def exclude(self, ex):
        """Excludes the given component from the potential.
        """
        self.check_exclude(ex)
        if ex in self.excludes:
            pass
        else:
            self.excludes.append(ex)
            #self.define_elements()

    def include(self, ex):
        """Includes the given component in the potential.
        """
        self.check_exclude(ex)
        if ex in self.excludes:
            self.excludes.remove(ex)
        else:
            pass
        #self.define_elements()

    def list_excludes(self):
        """Lists the excluded components.
        """
        return self.excludes

    def list_possible_excludes(self):
        """Lists the keywords of the list of components of the potential.
        """
        return self.possible_excludes

    def is_included(self, ex):
        """Returns ``True`` if the given component is included in the potential.
        """
        self.check_exclude(ex)
        if ex in self.excludes:
            return False
        else:
            return True
            
    def is_excluded(self, ex):
        """Returns ``True`` if the given component is excluded from the potential.
        """
        return not self.is_included(ex)
        

    def include_all(self):
        """Includes all components in the potential, except long range Coulomb interaction.
        """
        self.excludes = ['direct_coulomb','long_coulomb']
        #self.define_elements()

    def exclude_all(self):
        """Excludes everything.
        """
        self.excludes = copy.deepcopy(self.possible_excludes)
        #self.define_elements()
        
    def set_ewald(self,calc,cheat):
        """Switches on Ewald summation in ``calc``. If ``cheat == True``, a short ranged screened
        interaction is used. For internal use.
        """
        ewald = pysic.CoulombSummation()
        if cheat:
            #print "Ewald cheat"
            real_cutoff = 12.0
            k_cutoff = 0.0
            sigma = 3.535533906
            epsilon = .00552635
        else:
            #print "Ewald real"
            real_cutoff, k_cutoff, sigma, epsilon = estimate_ewald_parameters(12.0)
        ewald.set_parameter_value('real_cutoff',real_cutoff)
        ewald.set_parameter_value('k_cutoff',k_cutoff)
        ewald.set_parameter_value('sigma',sigma)
        ewald.set_parameter_value('epsilon',epsilon)
        calc.set_coulomb_summation(ewald)


    def define_elements(self):
        """Creates the components of the potential.
        """
        #print "\nBuilding Comb potential:"
        self.pieces = []

        # Si-O-Si bond bending
        if self.possible_excludes[0] not in self.excludes:
            #print "including Si-O-Si bond bending"
            bond_bend_pot = pysic.Potential('bond_bend')
            bond_bend_pot.set_parameter_value('epsilon',2.6)
            bond_bend_pot.set_parameter_value('theta_0',2.508561734)
            bond_bend_pot.set_parameter_value('n',1)
            bond_bend_pot.set_parameter_value('m',2)
            bond_bend_pot.set_symbols(['Si','O','Si'])
            bond_bend_pot.set_cutoff(3.05)
            bond_bend_pot.set_cutoff_margin(0.5) #is cutoff_margin correct?

            self.pieces.append(bond_bend_pot)

        # O-Si-O bond bending
        if self.possible_excludes[1] not in self.excludes:
            #print "including O-Si-O bond bending"
            bond_bend_pot2 = pysic.Potential('bond_bend')
            bond_bend_pot2.set_parameter_value('epsilon',0.3122)
            bond_bend_pot2.set_parameter_value('theta_0',1.910611932)
            bond_bend_pot2.set_parameter_value('n',1)
            bond_bend_pot2.set_parameter_value('m',2)
            bond_bend_pot2.set_symbols(['O','Si','O'])
            bond_bend_pot2.set_cutoff(3.05)
            bond_bend_pot2.set_cutoff_margin(0.5)

            self.pieces.append(bond_bend_pot2)
            
        # attractive Si-Si
        if self.possible_excludes[2] not in self.excludes:
            #print "including Si-Si attraction"
            B_ISiSi = pysic.Potential('charge_exp')
            B_ISiSi.set_parameter_value('epsilon',471.18)  #B_Si 
            B_ISiSi.set_parameter_value('Rmax1',-1.658949)
            B_ISiSi.set_parameter_value('Rmin1',1.651725)
            B_ISiSi.set_parameter_value('Qmax1',4)
            B_ISiSi.set_parameter_value('Qmin1',-4)
            B_ISiSi.set_parameter_value('Rmax2',-1.658949)
            B_ISiSi.set_parameter_value('Rmin2',1.651725)
            B_ISiSi.set_parameter_value('Qmax2',4)
            B_ISiSi.set_parameter_value('Qmin2',-4)
            B_ISiSi.set_parameter_value('xi1',1.7322)  #alpha_Si
            B_ISiSi.set_parameter_value('xi2',1.7322)

            B_IISiSi = pysic.Potential('charge_abs')
            B_IISiSi.set_parameter_value('a1',1) #aB
            B_IISiSi.set_parameter_value('b1',-9.536743164062500e-07) #-abs(bB)**nB = -9.53674E-7
            B_IISiSi.set_parameter_value('Q1',0) #QO
            B_IISiSi.set_parameter_value('n1',10) #nB
            B_IISiSi.set_parameter_value('a2',1)
            B_IISiSi.set_parameter_value('b2',-9.536743164062500e-07) 
            B_IISiSi.set_parameter_value('Q2',0) #QO
            B_IISiSi.set_parameter_value('n2',10)

            exp_potBSiSi = pysic.Potential('exponential')
            exp_potBSiSi.set_parameter_value('epsilon',-1)
            exp_potBSiSi.set_parameter_value('zeta',1.7322)  #alpha

            prodBSiSi = pysic.ProductPotential([B_ISiSi,B_IISiSi,exp_potBSiSi])

            prodBSiSi.set_cutoff(3.0) # was 3.0
            prodBSiSi.set_cutoff_margin(0.2) # was 0.2
            prodBSiSi.set_symbols(['Si','Si'])

            self.pieces.append(prodBSiSi)

            # Si-Si bond order factor
            if self.possible_excludes[5] not in self.excludes:
                #print "including Si-Si bond order factor"

                bijSiSi1 = pysic.BondOrderParameters('tersoff', symbols=['Si','Si','O'])
                bijSiSi1.set_parameter_value('beta',1.0999E-6)  #all Si parameters
                bijSiSi1.set_parameter_value('eta',0.78734)
                bijSiSi1.set_parameter_value('mu',3)
                bijSiSi1.set_parameter_value('a',1.7322)
                bijSiSi1.set_parameter_value('c',100390)
                bijSiSi1.set_parameter_value('d',16.218)
                bijSiSi1.set_parameter_value('h',-0.59826)
                bijSiSi1.set_cutoff(3.05)
                bijSiSi1.set_cutoff_margin(0.5)

                bijSiSi2 = pysic.BondOrderParameters('tersoff', symbols=['Si','Si','Si'])
                bijSiSi2.set_parameter_value('beta',1.0999E-6)  #all Si parameters
                bijSiSi2.set_parameter_value('eta',0.78734)
                bijSiSi2.set_parameter_value('mu',3)
                bijSiSi2.set_parameter_value('a',1.7322)
                bijSiSi2.set_parameter_value('c',100390)
                bijSiSi2.set_parameter_value('d',16.218)
                bijSiSi2.set_parameter_value('h',-0.59826)
                bijSiSi2.set_cutoff(3.0)
                bijSiSi2.set_cutoff_margin(0.2)

                coordSiSi = pysic.Coordinator([bijSiSi1,bijSiSi2])

                prodBSiSi.set_coordinator(coordSiSi)


        # attractive O-O
        if self.possible_excludes[3] not in self.excludes:
            #print "including O-O attraction"

            B_IOO = pysic.Potential('charge_exp')
            B_IOO.set_parameter_value('epsilon',260.8931)  #B_O
            B_IOO.set_parameter_value('Rmax1',-0.00112)
            B_IOO.set_parameter_value('Rmin1',0.00148)
            B_IOO.set_parameter_value('Qmax1',5.5046)
            B_IOO.set_parameter_value('Qmin1',-1.8349)
            B_IOO.set_parameter_value('Rmax2',-0.00112)
            B_IOO.set_parameter_value('Rmin2',0.00148)
            B_IOO.set_parameter_value('Qmax2',5.5046)
            B_IOO.set_parameter_value('Qmin2',-1.8349)
            B_IOO.set_parameter_value('xi1',2.68)  #alpha_O
            B_IOO.set_parameter_value('xi2',2.68)

            B_IIOO = pysic.Potential('charge_abs')
            B_IIOO.set_parameter_value('a1',1.00098)#1.071772002265933) #aB
            B_IIOO.set_parameter_value('b1',-.00000226)#-2.419566797838410e-06) #-abs(bB)**nB = -2.41926E-6
            B_IIOO.set_parameter_value('Q1',1.83485) #QO
            B_IIOO.set_parameter_value('n1',10) #nB
            B_IIOO.set_parameter_value('a2',1.00098)#1.071772002265933)
            B_IIOO.set_parameter_value('b2',-.00000226)#-2.419566797838410e-06) 
            B_IIOO.set_parameter_value('Q2',1.83485) #QO
            B_IIOO.set_parameter_value('n2',10)
            
            exp_potBOO = pysic.Potential('exponential')
            exp_potBOO.set_parameter_value('epsilon',-1)
            exp_potBOO.set_parameter_value('zeta',2.68) #alpha

            prodBOO = pysic.ProductPotential([B_IOO,B_IIOO,exp_potBOO])
            
            prodBOO.set_cutoff(3.0)
            prodBOO.set_cutoff_margin(0.4)
            prodBOO.set_symbols(['O','O'])

            self.pieces.append(prodBOO)
            
            # O-O bond order factor
            if self.possible_excludes[6] not in self.excludes:
                #print "including O-O bond order factor"

                bijOO1 = pysic.BondOrderParameters('tersoff', symbols=['O','O','Si'])
                bijOO1.set_parameter_value('beta',2.0) #all O parameters
                bijOO1.set_parameter_value('eta',1)
                bijOO1.set_parameter_value('mu',1)
                bijOO1.set_parameter_value('a',2.68)
                bijOO1.set_parameter_value('c',6.6)
                bijOO1.set_parameter_value('d',1)
                bijOO1.set_parameter_value('h',-0.229)
                bijOO1.set_cutoff(3.05)
                bijOO1.set_cutoff_margin(0.5)

                bijOO2 = pysic.BondOrderParameters('tersoff', symbols=['O','O','O'])
                bijOO2.set_parameter_value('beta',2.0) #all O parameters
                bijOO2.set_parameter_value('eta',1)
                bijOO2.set_parameter_value('mu',1)
                bijOO2.set_parameter_value('a',2.68)
                bijOO2.set_parameter_value('c',6.6)
                bijOO2.set_parameter_value('d',1)
                bijOO2.set_parameter_value('h',-0.229)
                bijOO2.set_cutoff(3.0)
                bijOO2.set_cutoff_margin(0.4)

                coordOO = pysic.Coordinator([bijOO1,bijOO2])

                prodBOO.set_coordinator(coordOO)



        # attractive Si-O
        if self.possible_excludes[4] not in self.excludes:
            #print "including Si-O attraction"

            B_ISiO = pysic.Potential('charge_exp')
            B_ISiO.set_parameter_value('epsilon',350.6103405)  # sqrt(B_Si*B_O)
            B_ISiO.set_parameter_value('Rmax1',-1.658949)
            B_ISiO.set_parameter_value('Rmin1',1.651725)
            B_ISiO.set_parameter_value('Qmax1',4)
            B_ISiO.set_parameter_value('Qmin1',-4)
            B_ISiO.set_parameter_value('Rmax2',-0.00112)
            B_ISiO.set_parameter_value('Rmin2',0.00148)
            B_ISiO.set_parameter_value('Qmax2',5.5046)
            B_ISiO.set_parameter_value('Qmin2',-1.8349)
            B_ISiO.set_parameter_value('xi1',1.7322)  #alpha_Si
            B_ISiO.set_parameter_value('xi2',2.68)  #alpha O

            B_IISiO = pysic.Potential('charge_abs')
            B_IISiO.set_parameter_value('a1',1) #aB_Si
            B_IISiO.set_parameter_value('b1', -9.536743164062500e-07) #-abs(bB)**nB for Si
            B_IISiO.set_parameter_value('Q1',0) #QO of Si
            B_IISiO.set_parameter_value('n1',10) #nB
            B_IISiO.set_parameter_value('a2',1.00098)#1.071772002265933) #aB_O
            B_IISiO.set_parameter_value('b2',-.00000226)#-2.419566797838410e-06) #-abs(bB)**nB for O
            B_IISiO.set_parameter_value('Q2',1.83485) #QO of O
            B_IISiO.set_parameter_value('n2',10)
            
            exp_potBSiO = pysic.Potential('exponential')
            exp_potBSiO.set_parameter_value('epsilon',-1) 
            exp_potBSiO.set_parameter_value('zeta',2.2061)  #alpha_SiO = (alpha_Si+alpha_O)/2

            prodBSiO = pysic.ProductPotential([B_ISiO,B_IISiO,exp_potBSiO])
            
            prodBSiO.set_cutoff(3.05)        # sqrt
            prodBSiO.set_cutoff_margin(0.5)#0.3018)  # Rs_both is sqrt(Rs1*Rs2) 
            prodBSiO.set_symbols(['Si','O'])

            self.pieces.append(prodBSiO)
            
            # Si-O bond order factor
            if self.possible_excludes[7] not in self.excludes:
                #print "including Si-O bond order factor"

                bijOSi1 = pysic.BondOrderParameters('tersoff', symbols=['O','Si','Si'])
                bijOSi1.set_parameter_value('beta',2)  #O
                bijOSi1.set_parameter_value('eta',1) #O
                bijOSi1.set_parameter_value('mu',1) #O
                bijOSi1.set_parameter_value('a',2.2061) #media between Si and O parameters
                bijOSi1.set_parameter_value('c',6.6) # O
                bijOSi1.set_parameter_value('d',1) #O
                bijOSi1.set_parameter_value('h',-0.229) #O
                bijOSi1.set_cutoff(3.05)
                bijOSi1.set_cutoff_margin(0.5) #for O

                bijOSi2 = pysic.BondOrderParameters('tersoff', symbols=['O','Si','O'])
                bijOSi2.set_parameter_value('beta',2)  #O
                bijOSi2.set_parameter_value('eta',1) #O
                bijOSi2.set_parameter_value('mu',1) #O
                bijOSi2.set_parameter_value('a',2.2061) #media between Si and O parameters
                bijOSi2.set_parameter_value('c',6.6) # O
                bijOSi2.set_parameter_value('d',1) #O
                bijOSi2.set_parameter_value('h',-0.229) #O
                bijOSi2.set_cutoff(3.05)
                bijOSi2.set_cutoff_margin(0.5) #for O

                bijSiO = pysic.BondOrderParameters('tersoff', symbols=[['Si','O','O'],['Si','O','Si']])
                bijSiO.set_parameter_value('beta',1.0999E-6)  #Si
                bijSiO.set_parameter_value('eta',0.78734) #Si
                bijSiO.set_parameter_value('mu',3) #Si
                bijSiO.set_parameter_value('a',2.2061)   #media between Si and O parameters
                bijSiO.set_parameter_value('c',100390) #Si
                bijSiO.set_parameter_value('d',16.218) #Si
                bijSiO.set_parameter_value('h',-0.59826) #Si
                bijSiO.set_cutoff(3.05)
                bijSiO.set_cutoff_margin(0.5) #for Si

                coordSiO = pysic.Coordinator([bijOSi1,bijOSi2, bijSiO])

                prodBSiO.set_coordinator(coordSiO)

            
        # repulsive Si-Si
        if self.possible_excludes[8] not in self.excludes:
            #print "including Si-Si repulsion"

            ASiSi = pysic.Potential('charge_exp')
            ASiSi.set_parameter_value('epsilon',1830.8)  #A_Si
            ASiSi.set_parameter_value('Rmax1',-1.658949)
            ASiSi.set_parameter_value('Rmin1',1.651725)
            ASiSi.set_parameter_value('Qmax1',4)
            ASiSi.set_parameter_value('Qmin1',-4)
            ASiSi.set_parameter_value('Rmax2',-1.658949)
            ASiSi.set_parameter_value('Rmin2',1.651725)
            ASiSi.set_parameter_value('Qmax2',4)
            ASiSi.set_parameter_value('Qmin2',-4)
            ASiSi.set_parameter_value('xi1',2.4799) #lambda
            ASiSi.set_parameter_value('xi2',2.4799)

            exp_potSi = pysic.Potential('exponential')
            exp_potSi.set_parameter_value('epsilon',1)
            exp_potSi.set_parameter_value('zeta', 2.4799)  #lambda_Si
            
            prodASiSi = pysic.ProductPotential([ASiSi,exp_potSi])

            prodASiSi.set_cutoff(3.0) # was 3.0
            prodASiSi.set_cutoff_margin(0.2) # was 0.2
            prodASiSi.set_symbols(['Si','Si'])

            self.pieces.append(prodASiSi)

        # repulsive O-O
        if self.possible_excludes[9] not in self.excludes:
            #print "including O-O repulsion"

            AOO = pysic.Potential('charge_exp')
            AOO.set_parameter_value('epsilon',3326.699)  #A_O
            AOO.set_parameter_value('Rmax1',-0.00112)
            AOO.set_parameter_value('Rmin1',0.00148)
            AOO.set_parameter_value('Qmax1',5.5046)
            AOO.set_parameter_value('Qmin1',-1.8349)
            AOO.set_parameter_value('Rmax2',-0.00112)
            AOO.set_parameter_value('Rmin2',0.00148)
            AOO.set_parameter_value('Qmax2',5.5046)
            AOO.set_parameter_value('Qmin2',-1.8349)
            AOO.set_parameter_value('xi1',5.36) #lambda
            AOO.set_parameter_value('xi2',5.36)

            exp_potO = pysic.Potential('exponential')
            exp_potO.set_parameter_value('epsilon',1)
            exp_potO.set_parameter_value('zeta',5.36)  #lambda_O

            prodAOO = pysic.ProductPotential([AOO,exp_potO])
 
            prodAOO.set_cutoff(3.0)
            prodAOO.set_cutoff_margin(0.4)
            prodAOO.set_symbols(['O','O'])

            self.pieces.append(prodAOO)

            
        # repulsive Si-O
        if self.possible_excludes[10] not in self.excludes:
            #print "including Si-O repulsion"

            ASiO = pysic.Potential('charge_exp')
            ASiO.set_parameter_value('epsilon',2467.894659826) #sqrt(A_Si*A_O)
            ASiO.set_parameter_value('Rmax1',-1.658949)
            ASiO.set_parameter_value('Rmin1',1.651725)
            ASiO.set_parameter_value('Qmax1',4)
            ASiO.set_parameter_value('Qmin1',-4)
            ASiO.set_parameter_value('Rmax2',-0.00112)
            ASiO.set_parameter_value('Rmin2',0.00148)
            ASiO.set_parameter_value('Qmax2',5.5046)
            ASiO.set_parameter_value('Qmin2',-1.8349)
            ASiO.set_parameter_value('xi1',2.4799) #lambda
            ASiO.set_parameter_value('xi2',5.36)
            
            exp_potSiO = pysic.Potential('exponential')
            exp_potSiO.set_parameter_value('epsilon',1)
            exp_potSiO.set_parameter_value('zeta', 3.91995)  # (lambda_Si + lambda_O)/2

            prodASiO = pysic.ProductPotential([ASiO,exp_potSiO])

            prodASiO.set_cutoff(3.05)
            prodASiO.set_cutoff_margin(0.5)#0.3018)  # Rs_both is sqrt(Rs1*Rs2) 
            prodASiO.set_symbols(['Si','O'])

            self.pieces.append(prodASiO)


        # charge self-energy Si
        if self.possible_excludes[11] not in self.excludes:
            #print "including Si self energy"

            chiSiJ = pysic.Potential('charge_self')
            chiSiJ.set_parameter_value('epsilon', 3.625144859)      #J, X = 0 and K = 0 for Si
            chiSiJ.set_parameter_value('n',2)

            chiSiL = pysic.Potential('charge_self')
            chiSiL.set_parameter_value('epsilon', 0.087067714)       # L
            chiSiL.set_parameter_value('n',4)
            chiSiJ.set_symbols(['Si'])
            chiSiL.set_symbols(['Si'])
            
            #selfpot offset
            selfpotsi = pysic.Potential('charge_self')
            selfpotsi.set_parameter_value('epsilon',-0.000413097)
            selfpotsi.set_parameter_value('n',2)
            selfpotsi.set_symbols(['Si'])

            self.pieces.append([chiSiJ])
            self.pieces.append([chiSiL])
            #self.pieces.append([selfpotsi])

            
        # charge self-energy O
        if self.possible_excludes[12] not in self.excludes:
            #print "including O self energy"

            chiOX = pysic.Potential('charge_self')
            chiOX.set_parameter_value('epsilon', 5.63441383)     #X
            chiOX.set_parameter_value('n',1)
            chiOJ = pysic.Potential('charge_self')
            chiOJ.set_parameter_value('epsilon', 7.689598017)   #J
            chiOJ.set_parameter_value('n',2)
            chiOK = pysic.Potential('charge_self')
            chiOK.set_parameter_value('epsilon', 4.51426991)    #K
            chiOK.set_parameter_value('n', 3)
            chiOL = pysic.Potential('charge_self')
            chiOL.set_parameter_value('epsilon', 1.330079082)      #L
            chiOL.set_parameter_value('n', 4)
            chiOX.set_symbols(['O'])
            chiOJ.set_symbols(['O'])
            chiOK.set_symbols(['O'])
            chiOL.set_symbols(['O'])

            #selfpot contribution
            selfpoto = pysic.Potential('charge_self')
            selfpoto.set_parameter_value('epsilon',-0.000413097)
            selfpoto.set_parameter_value('n',2)
            selfpoto.set_symbols(['O'])

            self.pieces.append(chiOX)
            self.pieces.append(chiOJ)
            self.pieces.append(chiOK)
            self.pieces.append(chiOL)
            #self.pieces.append(selfpoto)

        pow1 = pysic.Potential('power')           #just multiplies by r
        pow1.set_parameter_value('epsilon',1)
        pow1.set_parameter_value('a',1)
        pow1.set_parameter_value('n',-1)
        pow2 = pysic.Potential('power')           #just multiplies by (1/r^5)
        pow2.set_parameter_value('epsilon',1)
        pow2.set_parameter_value('a',1)
        pow2.set_parameter_value('n',5) 


        chp1 = pysic.Potential('charge_pair')
        chp1.set_parameter_value('epsilon',-8.362034089934842e-07)  #is (5/rc^6 * rho1_Si)
        chp1.set_parameter_value('n1',0)
        chp1.set_parameter_value('n2',1)          # multiplies by q_j 

        chp4 = pysic.Potential('charge_pair')
        chp4.set_parameter_value('epsilon', -0.499378)  #is (rho1_Si)
        chp4.set_parameter_value('n1',0)
        chp4.set_parameter_value('n2',1)          # multiplies by q_j 

        chp1b = pysic.Potential('charge_pair')
        chp1b.set_parameter_value('epsilon', 5.023320620606138e-06)  #is (5/rc^6 * rho2_Si)
        chp1b.set_parameter_value('n1',0)
        chp1b.set_parameter_value('n2',2)          # multiplies by q_j^2 

        chp4b = pysic.Potential('charge_pair')
        chp4b.set_parameter_value('epsilon', 2.999911)  #is (rho2_Si)
        chp4b.set_parameter_value('n1',0)
        chp4b.set_parameter_value('n2',2)          # multiplies by q_j^2 

        chp1c = pysic.Potential('charge_pair')
        chp1c.set_parameter_value('epsilon',-8.362034089934842e-07)  #is (5/rc^6 * rho1_Si)
        chp1c.set_parameter_value('n1',1)
        chp1c.set_parameter_value('n2',0)          # multiplies by q_i

        chp4c = pysic.Potential('charge_pair')
        chp4c.set_parameter_value('epsilon', -0.499378)  #is (rho1_Si)
        chp4c.set_parameter_value('n1',1)
        chp4c.set_parameter_value('n2',0)          # multiplies by q_i 

        chp1d = pysic.Potential('charge_pair')
        chp1d.set_parameter_value('epsilon', 5.023320620606138e-06)  #is (5/rc^6 * rho2_Si)
        chp1d.set_parameter_value('n1',2)
        chp1d.set_parameter_value('n2',0)          # multiplies by q_i^2 

        chp4d = pysic.Potential('charge_pair')
        chp4d.set_parameter_value('epsilon', 2.999911)  #is (rho2_Si)
        chp4d.set_parameter_value('n1',2)
        chp4d.set_parameter_value('n2',0)          # multiplies by q_i^2 

        
        # penalty function Si-Si
        if self.possible_excludes[13] not in self.excludes:
            #print "including Si-Si penalty"

            prod1 = pysic.ProductPotential([chp1, pow1])
            prod1.set_symbols(['Si','Si'])
            prod1.set_cutoff(12)   
            #prod1.set_cutoff_margin(0.2)
            self.pieces.append(prod1)

            chp2 = pysic.Potential('charge_pair')
            chp2.set_parameter_value('epsilon',1.003444090792181e-05)  #is -(5/rc^5 * rho1_Si)
            chp2.set_parameter_value('n1',0)
            chp2.set_parameter_value('n2',1)          # multiplies by q_j 
            chp2.set_symbols(['Si','Si'])
            chp2.set_cutoff(12)
            #chp2.set_cutoff_margin(0.2)
            self.pieces.append(chp2)

            chp3 = pysic.Potential('charge_pair')
            chp3.set_parameter_value('epsilon', 2.006888181584362e-06)  #is -(1/rc^5 * rho1_Si)
            chp3.set_parameter_value('n1',0)
            chp3.set_parameter_value('n2',1)          # multiplies by q_j 
            chp3.set_symbols(['Si','Si'])
            chp3.set_cutoff(12)
            #chp3.set_cutoff_margin(0.2)
            self.pieces.append(chp3)

            chp4.set_symbols(['Si','Si'])
            chp4.set_cutoff(12)
            #chp4.set_cutoff_margin(0.2)
            prod2 = pysic.ProductPotential([chp4, pow2])
            prod2.set_symbols(['Si','Si'])
            prod2.set_cutoff(12)
            #prod2.set_cutoff_margin(0.2)
            self.pieces.append(prod2)

            #Second term------------------------------------------------------------------------------------------
            prod1b = pysic.ProductPotential([chp1b, pow1])
            prod1b.set_symbols(['Si','Si'])
            prod1b.set_cutoff(12)
            #prod1b.set_cutoff_margin(0.2)
            self.pieces.append(prod1b)

            chp2b = pysic.Potential('charge_pair')
            chp2b.set_parameter_value('epsilon',-6.027984744727366e-05)  #is -(5/rc^5 * rho2_Si)
            chp2b.set_parameter_value('n1',0)
            chp2b.set_parameter_value('n2',2)          # multiplies by q_j^2 
            chp2b.set_symbols(['Si','Si'])
            chp2b.set_cutoff(12)
            #chp2b.set_cutoff_margin(0.2)
            self.pieces.append(chp2b)

            chp3b = pysic.Potential('charge_pair')
            chp3b.set_parameter_value('epsilon',-1.205596948945473e-05)  #is -(1/rc^5 * rho2_Si)
            chp3b.set_parameter_value('n1',0)
            chp3b.set_parameter_value('n2',2)          # multiplies by q_j^2 
            chp3b.set_symbols(['Si','Si'])
            chp3b.set_cutoff(12)
            #chp3b.set_cutoff_margin(0.2)
            self.pieces.append(chp3b)

            chp4b.set_symbols(['Si','Si'])
            chp4b.set_cutoff(12)
            chp4b.set_cutoff_margin(0.2)
            prod2b = pysic.ProductPotential([chp4b, pow2])
            prod2b.set_symbols(['Si','Si'])
            prod2b.set_cutoff(12)
            #prod2b.set_cutoff_margin(0.2)
            self.pieces.append(prod2b)

            #Third term---------------------------------------------------------------------------------------
            prod1c = pysic.ProductPotential([chp1c, pow1])
            prod1c.set_symbols(['Si','Si'])
            prod1c.set_cutoff(12)
            #prod1c.set_cutoff_margin(0.2)
            self.pieces.append(prod1c)

            chp2c = pysic.Potential('charge_pair')
            chp2c.set_parameter_value('epsilon',1.003444090792181e-05)  #is -(5/rc^5 * rho1_Si)
            chp2c.set_parameter_value('n1',1)
            chp2c.set_parameter_value('n2',0)          # multiplies by q_i
            chp2c.set_symbols(['Si','Si'])
            chp2c.set_cutoff(12)
            #chp2c.set_cutoff_margin(0.2)
            self.pieces.append(chp2c)

            chp3c = pysic.Potential('charge_pair')
            chp3c.set_parameter_value('epsilon', 2.006888181584362e-06)  #is -(1/rc^5 * rho1_Si)
            chp3c.set_parameter_value('n1',1)
            chp3c.set_parameter_value('n2',0)          # multiplies by q_i
            chp3c.set_symbols(['Si','Si'])
            chp3c.set_cutoff(12)
            #chp3c.set_cutoff_margin(0.2)
            self.pieces.append(chp3c)

            chp4c.set_symbols(['Si','Si'])
            chp4c.set_cutoff(12)
            #chp4c.set_cutoff_margin(0.2)
            prod2c = pysic.ProductPotential([chp4c, pow2])
            prod2c.set_symbols(['Si','Si'])
            prod2c.set_cutoff(12)
            #prod2c.set_cutoff_margin(0.2)
            self.pieces.append(prod2c)

            #Forth term-------------------------------------------------------------------------------------
            prod1d = pysic.ProductPotential([chp1d, pow1])
            prod1d.set_symbols(['Si','Si'])
            prod1d.set_cutoff(12)
            #prod1d.set_cutoff_margin(0.2)
            self.pieces.append(prod1d)

            chp2d = pysic.Potential('charge_pair')
            chp2d.set_parameter_value('epsilon',-6.027984744727366e-05)  #is -(5/rc^5 * rho2_Si)
            chp2d.set_parameter_value('n1',2)
            chp2d.set_parameter_value('n2',0)          # multiplies by q_i^2 
            chp2d.set_symbols(['Si','Si'])
            chp2d.set_cutoff(12)
            #chp2d.set_cutoff_margin(0.2)
            self.pieces.append(chp2d)

            chp3d = pysic.Potential('charge_pair')
            chp3d.set_parameter_value('epsilon',-1.205596948945473e-05)  #is -(1/rc^5 * rho2_Si)
            chp3d.set_parameter_value('n1',2)
            chp3d.set_parameter_value('n2',0)          # multiplies by q_i^2 
            chp3d.set_symbols(['Si','Si'])
            chp3d.set_cutoff(12)
            #chp3d.set_cutoff_margin(0.2)
            self.pieces.append(chp3d)

            chp4d.set_symbols(['Si','Si'])
            chp4d.set_cutoff(12)
            #chp4d.set_cutoff_margin(0.2)
            prod2d = pysic.ProductPotential([chp4d, pow2])
            prod2d.set_symbols(['Si','Si'])
            prod2d.set_cutoff(12)
            #prod2d.set_cutoff_margin(0.2)
            self.pieces.append(prod2d)



        chp1O = pysic.Potential('charge_pair')
        chp1O.set_parameter_value('epsilon',-6.567367742091049e-06)  #is (5/rc^6 * rho1_O)
        chp1O.set_parameter_value('n1',0)
        chp1O.set_parameter_value('n2',1)          # multiplies by q_j 

        chp4O = pysic.Potential('charge_pair')
        chp4O.set_parameter_value('epsilon', -3.9220110000000)  #is (rho1_O)
        chp4O.set_parameter_value('n1',0)
        chp4O.set_parameter_value('n2',1)          # multiplies by q_j 
 
        chp1bO = pysic.Potential('charge_pair')
        chp1bO.set_parameter_value('epsilon', 1.626073682913237e-06)  #is (5/rc^6 * rho2_O)
        chp1bO.set_parameter_value('n1',0)
        chp1bO.set_parameter_value('n2',2)          # multiplies by q_j^2 

        chp4bO = pysic.Potential('charge_pair')
        chp4bO.set_parameter_value('epsilon',0.971086)  #is (rho2_O)
        chp4bO.set_parameter_value('n1',0)
        chp4bO.set_parameter_value('n2',2)          # multiplies by q_j^2 

        chp1cO = pysic.Potential('charge_pair')
        chp1cO.set_parameter_value('epsilon',-6.567367742091049e-06)  #is (5/rc^6 * rho1_O)
        chp1cO.set_parameter_value('n1',1)
        chp1cO.set_parameter_value('n2',0)          # multiplies by q_i

        chp4cO = pysic.Potential('charge_pair')
        chp4cO.set_parameter_value('epsilon', -3.922011000)  # is (rho1_O)
        chp4cO.set_parameter_value('n1',1)
        chp4cO.set_parameter_value('n2',0)                   # multiplies by q_i 

        chp1dO = pysic.Potential('charge_pair')
        chp1dO.set_parameter_value('epsilon', 1.626073682913237e-06)  #is (5/rc^6 * rho2_O)
        chp1dO.set_parameter_value('n1',2)
        chp1dO.set_parameter_value('n2',0)          # multiplies by q_i^2 

        chp4dO = pysic.Potential('charge_pair')
        chp4dO.set_parameter_value('epsilon', 0.971086)  #is (rho2_O)
        chp4dO.set_parameter_value('n1',2)
        chp4dO.set_parameter_value('n2',0)          # multiplies by q_i^2 


        # penalty function O-O
        if self.possible_excludes[14] not in self.excludes: 
            #print "including O-O penalty"

            #First term--------------------------------------------------------------------------------------
            prod1O = pysic.ProductPotential([chp1O, pow1])
            prod1O.set_symbols(['O','O'])
            prod1O.set_cutoff(12)
            #prod1O.set_cutoff_margin(0.2)
            self.pieces.append(prod1O)

            chp2O = pysic.Potential('charge_pair')
            chp2O.set_parameter_value('epsilon',7.880841290509259e-05)  #is -(5/rc^5 * rho1_O)
            chp2O.set_parameter_value('n1',0)
            chp2O.set_parameter_value('n2',1)          # multiplies by q_j 
            chp2O.set_symbols(['O','O'])
            chp2O.set_cutoff(12)
            #chp2O.set_cutoff_margin(0.2)
            self.pieces.append(chp2O)

            chp3O = pysic.Potential('charge_pair')
            chp3O.set_parameter_value('epsilon',1.576168258101852e-05)  #is -(1/rc^5 * rho1_O)
            chp3O.set_parameter_value('n1',0)
            chp3O.set_parameter_value('n2',1)          # multiplies by q_j 
            chp3O.set_symbols(['O','O'])
            chp3O.set_cutoff(12)
            #chp3O.set_cutoff_margin(0.2)
            self.pieces.append(chp3O)

            prod2O = pysic.ProductPotential([chp4O, pow2])
            prod2O.set_symbols(['O','O'])
            prod2O.set_cutoff(12)
            #prod2O.set_cutoff_margin(0.2)
            self.pieces.append(prod2O)


            #Second term-----------------------------------------------------------------------------------------
            prod1bO = pysic.ProductPotential([chp1bO, pow1])
            prod1bO.set_symbols(['O','O'])
            prod1bO.set_cutoff(12)
            #prod1bO.set_cutoff_margin(0.2)
            self.pieces.append(prod1bO)
 
            chp2bO = pysic.Potential('charge_pair')
            chp2bO.set_parameter_value('epsilon',-1.951288419495885e-05)  #is -(5/rc^5 * rho2_O)
            chp2bO.set_parameter_value('n1',0)
            chp2bO.set_parameter_value('n2',2)          # multiplies by q_j^2 
            chp2bO.set_symbols(['O','O'])
            chp2bO.set_cutoff(12)
            #chp2bO.set_cutoff_margin(0.2)
            self.pieces.append(chp2bO)

            chp3bO = pysic.Potential('charge_pair')
            chp3bO.set_parameter_value('epsilon',-3.902576838991770e-06)  #is -(1/rc^5 * rho2_O)
            chp3bO.set_parameter_value('n1',0)
            chp3bO.set_parameter_value('n2',2)          # multiplies by q_j^2 
            chp3bO.set_symbols(['O','O'])
            chp3bO.set_cutoff(12)
            #chp3bO.set_cutoff_margin(0.2)
            self.pieces.append(chp3bO)

            prod2bO = pysic.ProductPotential([chp4bO, pow2])
            prod2bO.set_symbols(['O','O'])
            prod2bO.set_cutoff(12)
            #prod2bO.set_cutoff_margin(0.2)
            self.pieces.append(prod2bO)

            #Third term---------------------------------------------------------------------------------------
            prod1cO = pysic.ProductPotential([chp1cO, pow1])
            prod1cO.set_symbols(['O','O'])
            prod1cO.set_cutoff(12)
            #prod1cO.set_cutoff_margin(0.2)
            self.pieces.append(prod1cO)

            chp2cO = pysic.Potential('charge_pair')
            chp2cO.set_parameter_value('epsilon',7.880841290509259e-05)  #is -(5/rc^5 * rho1_O)
            chp2cO.set_parameter_value('n1',1)
            chp2cO.set_parameter_value('n2',0)          # multiplies by q_i
            chp2cO.set_symbols(['O','O'])
            chp2cO.set_cutoff(12)
            #chp2cO.set_cutoff_margin(0.2)
            self.pieces.append(chp2cO)

            chp3cO = pysic.Potential('charge_pair')
            chp3cO.set_parameter_value('epsilon', 1.576168258101852e-05)  #is -(1/rc^5 * rho1_O)
            chp3cO.set_parameter_value('n1',1)
            chp3cO.set_parameter_value('n2',0)          # multiplies by q_i
            chp3cO.set_symbols(['O','O'])
            chp3cO.set_cutoff(12)
            #chp3cO.set_cutoff_margin(0.2)
            self.pieces.append(chp3cO)

            prod2cO = pysic.ProductPotential([chp4cO, pow2])
            prod2cO.set_symbols(['O','O'])
            prod2cO.set_cutoff(12)
            #prod2cO.set_cutoff_margin(0.2)
            self.pieces.append(prod2cO)

            #Forth term------------------------------------------------------------------------------------------
            prod1dO = pysic.ProductPotential([chp1dO, pow1])
            prod1dO.set_symbols(['O','O'])
            prod1dO.set_cutoff(12)
            #prod1dO.set_cutoff_margin(0.2)
            self.pieces.append(prod1dO)

            chp2dO = pysic.Potential('charge_pair')
            chp2dO.set_parameter_value('epsilon', -1.951288419495885e-05)  #is -(5/rc^5 * rho2_O)
            chp2dO.set_parameter_value('n1',2)
            chp2dO.set_parameter_value('n2',0)          # multiplies by q_i^2 
            chp2dO.set_symbols(['O','O'])
            chp2dO.set_cutoff(12)
            #chp2dO.set_cutoff_margin(0.2)
            self.pieces.append(chp2dO)

            chp3dO = pysic.Potential('charge_pair')
            chp3dO.set_parameter_value('epsilon',-3.902576838991770e-06)  #is -(1/rc^5 * rho2_O)
            chp3dO.set_parameter_value('n1',2)
            chp3dO.set_parameter_value('n2',0)          # multiplies by q_i^2 
            chp3dO.set_symbols(['O','O'])
            chp3dO.set_cutoff(12)
            #chp3dO.set_cutoff_margin(0.2)
            self.pieces.append(chp3dO)

            prod2dO = pysic.ProductPotential([chp4dO, pow2])
            prod2dO.set_symbols(['O','O'])
            prod2dO.set_cutoff(12)
            #prod2dO.set_cutoff_margin(0.2)
            self.pieces.append(prod2dO)


        # penalty function Si-O
        if self.possible_excludes[15] not in self.excludes:
            #print "including Si-O penalty"

            #First term--------------------------------------------------------------------------------------
            prod1mix = pysic.ProductPotential([chp1, pow1])
            prod1mix.set_symbols(['Si','O'])
            prod1mix.set_cutoff(12)
            #prod1mix.set_cutoff_margin(0.2)
            self.pieces.append(prod1mix)

            chp2mix = pysic.Potential('charge_pair')
            chp2mix.set_parameter_value('epsilon',1.003444090792181e-05)  #is -(5/rc^5 * rho1_Si)
            chp2mix.set_parameter_value('n1',0)
            chp2mix.set_parameter_value('n2',1)          # multiplies by q_j 
            chp2mix.set_symbols(['Si','O'])
            chp2mix.set_cutoff(12)
            #chp2mix.set_cutoff_margin(0.2)
            self.pieces.append(chp2mix)

            chp3mix = pysic.Potential('charge_pair')
            chp3mix.set_parameter_value('epsilon', 2.006888181584362e-06)  #is -(1/rc^5 * rho1_Si)
            chp3mix.set_parameter_value('n1',0)
            chp3mix.set_parameter_value('n2',1)          # multiplies by q_j 
            chp3mix.set_symbols(['Si','O'])
            chp3mix.set_cutoff(12)
            #chp3mix.set_cutoff_margin(0.2)
            self.pieces.append(chp3mix)

            prod2mix = pysic.ProductPotential([chp4, pow2])
            prod2mix.set_symbols(['Si','O'])
            prod2mix.set_cutoff(12)
            #prod2mix.set_cutoff_margin(0.2)
            self.pieces.append(prod2mix)

            #Second term--------------------------------------------------------------------------------------
            prod1bmix = pysic.ProductPotential([chp1b, pow1])
            prod1bmix.set_symbols(['Si','O'])
            prod1bmix.set_cutoff(12)
            #prod1bmix.set_cutoff_margin(0.2)
            self.pieces.append(prod1bmix)

            chp2bmix = pysic.Potential('charge_pair')
            chp2bmix.set_parameter_value('epsilon',-6.027984744727366e-05)  #is -(5/rc^5 * rho2_Si)
            chp2bmix.set_parameter_value('n1',0)
            chp2bmix.set_parameter_value('n2',2)          # multiplies by q_j^2 
            chp2bmix.set_symbols(['Si','O'])
            chp2bmix.set_cutoff(12)
            #chp2bmix.set_cutoff_margin(0.2)
            self.pieces.append(chp2bmix)

            chp3bmix = pysic.Potential('charge_pair')
            chp3bmix.set_parameter_value('epsilon',-1.205596948945473e-05)  #is -(1/rc^5 * rho2_Si)
            chp3bmix.set_parameter_value('n1',0)
            chp3bmix.set_parameter_value('n2',2)          # multiplies by q_j^2 
            chp3bmix.set_symbols(['Si','O'])
            chp3bmix.set_cutoff(12)
            #chp3bmix.set_cutoff_margin(0.2)
            self.pieces.append(chp3bmix)

            prod2bmix = pysic.ProductPotential([chp4b, pow2])
            prod2bmix.set_symbols(['Si','O'])
            prod2bmix.set_cutoff(12)
            #prod2bmix.set_cutoff_margin(0.2)
            self.pieces.append(prod2bmix)

            #Third term---------------------------------------------------------------------------------------------
            prod1cmix = pysic.ProductPotential([chp1cO, pow1])
            prod1cmix.set_symbols(['Si','O'])
            prod1cmix.set_cutoff(12)
            #prod1cmix.set_cutoff_margin(0.2)
            self.pieces.append(prod1cmix)

            chp2cmix = pysic.Potential('charge_pair')
            chp2cmix.set_parameter_value('epsilon',7.880841290509259e-05)  #is -(5/rc^5 * rho1_O)
            chp2cmix.set_parameter_value('n1',1)
            chp2cmix.set_parameter_value('n2',0)          # multiplies by q_i
            chp2cmix.set_symbols(['Si','O'])
            chp2cmix.set_cutoff(12)
            #chp2cmix.set_cutoff_margin(0.2)
            self.pieces.append(chp2cmix)

            chp3cmix = pysic.Potential('charge_pair')
            chp3cmix.set_parameter_value('epsilon', 1.576168258101852e-05)  #is -(1/rc^5 * rho1_O)
            chp3cmix.set_parameter_value('n1',1)
            chp3cmix.set_parameter_value('n2',0)          # multiplies by q_i
            chp3cmix.set_symbols(['Si','O'])
            chp3cmix.set_cutoff(12)
            #chp3cmix.set_cutoff_margin(0.2)
            self.pieces.append(chp3cmix)

            prod2cmix = pysic.ProductPotential([chp4cO, pow2])
            prod2cmix.set_symbols(['Si','O'])
            prod2cmix.set_cutoff(12)
            #prod2cmix.set_cutoff_margin(0.2)
            self.pieces.append(prod2cmix)

            #Forth term------------------------------------------------------------------------------------------
            prod1dmix = pysic.ProductPotential([chp1dO, pow1])
            prod1dmix.set_symbols(['Si','O'])
            prod1dmix.set_cutoff(12)
            #prod1dmix.set_cutoff_margin(0.2)
            self.pieces.append(prod1dmix)

            chp2dmix = pysic.Potential('charge_pair')
            chp2dmix.set_parameter_value('epsilon', -1.951288419495885e-05)  #is -(5/rc^5 * rho2_O)
            chp2dmix.set_parameter_value('n1',2)
            chp2dmix.set_parameter_value('n2',0)          # multiplies by q_i^2 
            chp2dmix.set_symbols(['Si','O'])
            chp2dmix.set_cutoff(12)
            #chp2dmix.set_cutoff_margin(0.2)
            self.pieces.append(chp2dmix)

            chp3dmix = pysic.Potential('charge_pair')
            chp3dmix.set_parameter_value('epsilon',-3.902576838991770e-06)  #is -(1/rc^5 * rho2_O)
            chp3dmix.set_parameter_value('n1',2)
            chp3dmix.set_parameter_value('n2',0)          # multiplies by q_i^2 
            chp3dmix.set_symbols(['Si','O'])
            chp3dmix.set_cutoff(12)
            #chp3dmix.set_cutoff_margin(0.2)
            self.pieces.append(chp3dmix)

            prod2dmix = pysic.ProductPotential([chp4dO, pow2])
            prod2dmix.set_symbols(['Si','O'])
            prod2dmix.set_cutoff(12)
            #prod2dmix.set_cutoff_margin(0.2)
            self.pieces.append(prod2dmix)

        charged = pysic.Potential('charge_pair')
        charged.set_parameter_value('epsilon',1)
        charged.set_parameter_value('n1',1)
        charged.set_parameter_value('n2',1)

        decay = pysic.Potential('power')
        decay.set_parameter_value('epsilon', 14.3996)
        decay.set_parameter_value('a',1)
        decay.set_parameter_value('n',1)

        si_zeta = 0.772871
        o_zeta  = 2.243072

        si_exp_decay = pysic.Potential('exponential')
        si_exp_decay.set_parameter_value('epsilon',-1)
        si_exp_decay.set_parameter_value('zeta',2*si_zeta)

        o_exp_decay = pysic.Potential('exponential')
        o_exp_decay.set_parameter_value('epsilon',-1)
        o_exp_decay.set_parameter_value('zeta',2*o_zeta)


        # coulomb Si-Si
        if self.possible_excludes[16] not in self.excludes:
            #print "including Si-Si coulomb"
            
            poly_sisi1 = pysic.Potential('shift_power')
            poly_sisi1.set_parameter_value('epsilon',11.0/8.0*si_zeta)
            poly_sisi1.set_parameter_value('r1',0)
            poly_sisi1.set_parameter_value('r2',1)
            poly_sisi1.set_parameter_value('n',1)

            poly_sisi2 = pysic.Potential('shift_power')
            poly_sisi2.set_parameter_value('epsilon',3.0/4.0*si_zeta*si_zeta)
            poly_sisi2.set_parameter_value('r1',0)
            poly_sisi2.set_parameter_value('r2',1)
            poly_sisi2.set_parameter_value('n',2)

            poly_sisi3 = pysic.Potential('shift_power')
            poly_sisi3.set_parameter_value('epsilon',1.0/6.0*si_zeta*si_zeta*si_zeta)
            poly_sisi3.set_parameter_value('r1',0)
            poly_sisi3.set_parameter_value('r2',1)
            poly_sisi3.set_parameter_value('n',3)
            
            
            # offsets
            sisi_offset1 = pysic.Potential('charge_pair')
            sisi_offset1.set_parameter_value('epsilon',-0.000919025)
            sisi_offset1.set_parameter_value('n1',1)
            sisi_offset1.set_parameter_value('n2',1)
            
            sisi_offset2 = pysic.Potential('shift_power')
            sisi_offset2.set_parameter_value('epsilon',1.0)
            sisi_offset2.set_parameter_value('r1',0)
            sisi_offset2.set_parameter_value('r2',1)
            sisi_offset2.set_parameter_value('n',1)
                        
            sisi_offsetB = pysic.Potential('charge_pair')
            sisi_offsetB.set_parameter_value('epsilon',0.0118523)
            sisi_offsetB.set_parameter_value('n1',1)
            sisi_offsetB.set_parameter_value('n2',1)
                        
            

            coul_sisi = []
            if self.possible_excludes[19] not in self.excludes:
                #print "including Si-Si direct coulomb"
                coul_sisi.append( pysic.ProductPotential([charged,decay]) )
            coul_sisi.append( pysic.ProductPotential([charged,decay,si_exp_decay]) )
            coul_sisi.append( pysic.ProductPotential([charged,decay,si_exp_decay,poly_sisi1]) )
            coul_sisi.append( pysic.ProductPotential([charged,decay,si_exp_decay,poly_sisi2]) )
            coul_sisi.append( pysic.ProductPotential([charged,decay,si_exp_decay,poly_sisi3]) )

            #coul_sisi.append( pysic.ProductPotential([sisi_offset1,sisi_offset2]) )
            #coul_sisi.append( sisi_offsetB )

            for c_sisi in coul_sisi:
                c_sisi.set_symbols(['Si','Si'])
                c_sisi.set_cutoff(12.0)
                self.pieces.append(c_sisi)

            
        # coulomb O-O
        if self.possible_excludes[17] not in self.excludes:
            #print "including O-O coulomb"

            poly_oo1 = pysic.Potential('shift_power')
            poly_oo1.set_parameter_value('epsilon',11.0/8.0*o_zeta)
            poly_oo1.set_parameter_value('r1',0)
            poly_oo1.set_parameter_value('r2',1)
            poly_oo1.set_parameter_value('n',1)

            poly_oo2 = pysic.Potential('shift_power')
            poly_oo2.set_parameter_value('epsilon',3.0/4.0*o_zeta*o_zeta)
            poly_oo2.set_parameter_value('r1',0)
            poly_oo2.set_parameter_value('r2',1)
            poly_oo2.set_parameter_value('n',2)

            poly_oo3 = pysic.Potential('shift_power')
            poly_oo3.set_parameter_value('epsilon',1.0/6.0*o_zeta*o_zeta*o_zeta)
            poly_oo3.set_parameter_value('r1',0)
            poly_oo3.set_parameter_value('r2',1)
            poly_oo3.set_parameter_value('n',3)


            # offsets
            oo_offset1 = pysic.Potential('charge_pair')
            oo_offset1.set_parameter_value('epsilon',-0.000922181)
            oo_offset1.set_parameter_value('n1',1)
            oo_offset1.set_parameter_value('n2',1)
            
            oo_offset2 = pysic.Potential('shift_power')
            oo_offset2.set_parameter_value('epsilon',1.0)
            oo_offset2.set_parameter_value('r1',0)
            oo_offset2.set_parameter_value('r2',1)
            oo_offset2.set_parameter_value('n',1)
                        
            oo_offsetB = pysic.Potential('charge_pair')
            oo_offsetB.set_parameter_value('epsilon',0.0118924)
            oo_offsetB.set_parameter_value('n1',1)
            oo_offsetB.set_parameter_value('n2',1)


            coul_oo = []
            if self.possible_excludes[19] not in self.excludes:
                #print "including O-O direct coulomb"
                coul_oo.append( pysic.ProductPotential([charged,decay]) )
            coul_oo.append( pysic.ProductPotential([charged,decay,o_exp_decay]) )
            coul_oo.append( pysic.ProductPotential([charged,decay,o_exp_decay,poly_oo1]) )
            coul_oo.append( pysic.ProductPotential([charged,decay,o_exp_decay,poly_oo2]) )
            coul_oo.append( pysic.ProductPotential([charged,decay,o_exp_decay,poly_oo3]) )

            #coul_oo.append( pysic.ProductPotential([oo_offset1,oo_offset2]) )
            #coul_oo.append( oo_offsetB )
            
            for c_oo in coul_oo:
                c_oo.set_symbols(['O','O'])
                c_oo.set_cutoff(12.0)
                self.pieces.append(c_oo)

        # coulomb Si-O
        if self.possible_excludes[18] not in self.excludes:
            #print "including Si-O coulomb"

            kappa = (si_zeta*si_zeta + o_zeta*o_zeta)/(si_zeta*si_zeta - o_zeta*o_zeta)
            
            poly_sio1 = pysic.Potential('shift_power')
            poly_sio1.set_parameter_value('epsilon',(1-kappa)*(1-kappa)*1.0/4.0*si_zeta)
            poly_sio1.set_parameter_value('r1',0)
            poly_sio1.set_parameter_value('r2',1)
            poly_sio1.set_parameter_value('n',1)

            poly_sio2 = pysic.Potential('shift_power')
            poly_sio2.set_parameter_value('epsilon',(1+kappa)*(1+kappa)*1.0/4.0*o_zeta)
            poly_sio2.set_parameter_value('r1',0)
            poly_sio2.set_parameter_value('r2',1)
            poly_sio2.set_parameter_value('n',1)

            scaled_si_exp_decay = pysic.Potential('exponential')
            scaled_si_exp_decay.set_parameter_value('epsilon',-(1-kappa)*(1-kappa)*(2+kappa)/4.0)
            scaled_si_exp_decay.set_parameter_value('zeta',2*si_zeta)
        
            scaled_o_exp_decay = pysic.Potential('exponential')
            scaled_o_exp_decay.set_parameter_value('epsilon',-(1+kappa)*(1+kappa)*(2-kappa)/4.0)
            scaled_o_exp_decay.set_parameter_value('zeta',2*o_zeta)



            # offsets
            sio_offset1 = pysic.Potential('charge_pair')
            sio_offset1.set_parameter_value('epsilon',-0.000921978)
            sio_offset1.set_parameter_value('n1',1)
            sio_offset1.set_parameter_value('n2',1)
            
            sio_offset2 = pysic.Potential('shift_power')
            sio_offset2.set_parameter_value('epsilon',1.0)
            sio_offset2.set_parameter_value('r1',0)
            sio_offset2.set_parameter_value('r2',1)
            sio_offset2.set_parameter_value('n',1)
                        
            sio_offsetB = pysic.Potential('charge_pair')
            sio_offsetB.set_parameter_value('epsilon',0.0118897)
            sio_offsetB.set_parameter_value('n1',1)
            sio_offsetB.set_parameter_value('n2',1)


            coul_sio = []
            if self.possible_excludes[19] not in self.excludes:
                #print "including Si-O direct coulomb"
                coul_sio.append( pysic.ProductPotential([charged,decay]) )
            coul_sio.append( pysic.ProductPotential([charged,decay,scaled_si_exp_decay]) )
            coul_sio.append( pysic.ProductPotential([charged,decay,scaled_o_exp_decay]) )
            coul_sio.append( pysic.ProductPotential([charged,decay,poly_sio1,si_exp_decay]) )
            coul_sio.append( pysic.ProductPotential([charged,decay,poly_sio2,o_exp_decay]) )

            #coul_sio.append( pysic.ProductPotential([sio_offset1,sio_offset2]) )
            #coul_sio.append( sio_offsetB )

            for c_sio in coul_sio:
                c_sio.set_symbols(['Si','O'])
                c_sio.set_cutoff(12.0)
                self.pieces.append(c_sio)
                
        if self.possible_excludes[19] in self.excludes:
            if self.calc == None:
                raise Exception("A Pysic calculator needed for Ewald. Use set_calculator(calc).")
            if self.possible_excludes[20] in self.excludes:
                #print "including Ewald without k-summation"
                self.set_ewald(self.calc,cheat=True)
            else:
                #print "including Ewald with k-summation"
                self.set_ewald(self.calc,cheat=False)
                
        #print "\n"
