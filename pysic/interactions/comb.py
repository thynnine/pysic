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
        self.description = "Comb potential."
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
    
        
        self.names_of_params = []
        self.params = {}
        self.n_params = 0


        newparam = 'coulomb cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 12.0
        self.n_params += 1
        newparam = 'coulomb sigma'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.535533906
        self.n_params += 1
        newparam = 'coulomb epsilon'
        self.names_of_params.append(newparam)
        self.params[newparam] = .00552635
        self.n_params += 1

        newparam = 'si-o-si bend epsilon'
        self.names_of_params.append(newparam)
        self.params[newparam] = 2.6
        self.n_params += 1
        newparam = 'si-o-si bend theta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 2.508561734
        self.n_params += 1
        newparam = 'si-o-si bend cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.05
        self.n_params += 1
        newparam = 'si-o-si bend cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.5
        self.n_params += 1
        
        newparam = 'o-si-o bend epsilon'
        self.names_of_params.append('osio bend epsilon')
        self.params[newparam] = 0.3122
        self.n_params += 1
        newparam = 'o-si-o bend theta'
        self.names_of_params.append('osio bend theta')
        self.params[newparam] = 1.910611932
        self.n_params += 1
        newparam = 'o-si-o bend cutoff'
        self.names_of_params.append('osio bend cutoff')
        self.params[newparam] = 3.05
        self.n_params += 1
        newparam = 'o-si-o bend cutoff margin'
        self.names_of_params.append('osio bend cutoff margin')
        self.params[newparam] = 0.5
        self.n_params += 1

        newparam = 'si attractive epsilon'
        self.names_of_params.append(newparam)
        self.params[newparam] = 471.18
        self.n_params += 1
        newparam = 'si rmax'
        self.names_of_params.append(newparam)
        self.params[newparam] = -1.658949
        self.n_params += 1
        newparam = 'si rmin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1.651725
        self.n_params += 1
        newparam = 'si qmax'
        self.names_of_params.append(newparam)
        self.params[newparam] = 4
        self.n_params += 1
        newparam = 'si qmin'
        self.names_of_params.append(newparam)
        self.params[newparam] = -4
        self.n_params += 1
        newparam = 'si attractive alpha'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1.7322
        self.n_params += 1
        newparam = 'si attractive a'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1.0
        self.n_params += 1
        newparam = 'si attractive b'
        self.names_of_params.append(newparam)
        self.params[newparam] = -9.536743164062500e-07
        self.n_params += 1
        newparam = 'si attractive q'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.0
        self.n_params += 1
        newparam = 'si attractive n'
        self.names_of_params.append(newparam)
        self.params[newparam] = 10
        self.n_params += 1
        newparam = 'si attractive zeta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1.7322
        self.n_params += 1
        newparam = 'si attractive cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.0
        self.n_params += 1
        newparam = 'si attractive cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.2
        self.n_params += 1
    
        newparam = 'si bondorder beta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1.0999E-6
        self.n_params += 1
        newparam = 'si bondorder eta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.78734
        self.n_params += 1
        newparam = 'si bondorder mu'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3
        self.n_params += 1
        newparam = 'si bondorder a'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1.7322
        self.n_params += 1
        newparam = 'si bondorder c'
        self.names_of_params.append(newparam)
        self.params[newparam] = 100390
        self.n_params += 1
        newparam = 'si bondorder d'
        self.names_of_params.append(newparam)
        self.params[newparam] = 16.218
        self.n_params += 1
        newparam = 'si bondorder h'
        self.names_of_params.append(newparam)
        self.params[newparam] = -0.59826
        self.n_params += 1
        newparam = 'si-si-o bondorder cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.05
        self.n_params += 1
        newparam = 'si-si-o bondorder cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.5
        self.n_params += 1
        newparam = 'si-si-si bondorder cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.0
        self.n_params += 1
        newparam = 'si-si-si bondorder cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.2
        self.n_params += 1
        
        newparam = 'o attractive epsilon'
        self.names_of_params.append(newparam)
        self.params[newparam] = 260.8931
        self.n_params += 1
        newparam = 'o rmax'
        self.names_of_params.append(newparam)
        self.params[newparam] = -0.00112
        self.n_params += 1
        newparam = 'o rmin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.00148
        self.n_params += 1
        newparam = 'o qmax'
        self.names_of_params.append(newparam)
        self.params[newparam] = 5.5046
        self.n_params += 1
        newparam = 'o qmin'
        self.names_of_params.append(newparam)
        self.params[newparam] = -1.8349
        self.n_params += 1
        newparam = 'o attractive alpha'
        self.names_of_params.append(newparam)
        self.params[newparam] = 2.68
        self.n_params += 1
        newparam = 'o attractive a'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1.00098
        self.n_params += 1
        newparam = 'o attractive b'
        self.names_of_params.append(newparam)
        self.params[newparam] = -.00000226
        self.n_params += 1
        newparam = 'o attractive q'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1.83485
        self.n_params += 1
        newparam = 'o attractive n'
        self.names_of_params.append(newparam)
        self.params[newparam] = 10
        self.n_params += 1
        newparam = 'o attractive zeta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 2.68
        self.n_params += 1
        newparam = 'o attractive cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.0
        self.n_params += 1
        newparam = 'o attractive cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.4
        self.n_params += 1
    
        newparam = 'o bondorder beta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 2.0
        self.n_params += 1
        newparam = 'o bondorder eta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1.0
        self.n_params += 1
        newparam = 'o bondorder mu'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1
        self.n_params += 1
        newparam = 'o bondorder a'
        self.names_of_params.append(newparam)
        self.params[newparam] = 2.68
        self.n_params += 1
        newparam = 'o bondorder c'
        self.names_of_params.append(newparam)
        self.params[newparam] = 6.6
        self.n_params += 1
        newparam = 'o bondorder d'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1
        self.n_params += 1
        newparam = 'o bondorder h'
        self.names_of_params.append(newparam)
        self.params[newparam] = -0.229
        self.n_params += 1
        newparam = 'o-o-si bondorder cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.05
        self.n_params += 1
        newparam = 'o-o-si bondorder cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.5
        self.n_params += 1
        newparam = 'o-o-o bondorder cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.0
        self.n_params += 1
        newparam = 'o-o-o bondorder cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.4
        self.n_params += 1
        
        newparam = 'si-o attractive epsilon'
        self.names_of_params.append(newparam)
        self.params[newparam] = math.sqrt( \
            self.params['si attractive epsilon'] * \
            self.params['o attractive epsilon'] )
        self.n_params += 1
        newparam = 'si-o attractive zeta'
        self.names_of_params.append(newparam)
        self.params[newparam] = ( \
            self.params['si attractive zeta'] + \
            self.params['o attractive zeta'] ) * 0.5
        self.n_params += 1
        newparam = 'si-o attractive cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.05
        self.n_params += 1
        newparam = 'si-o attractive cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.5
        self.n_params += 1
    
        newparam = 'si-o bondorder a'
        self.names_of_params.append(newparam)
        self.params[newparam] = ( \
            self.params['si bondorder a'] + \
            self.params['o bondorder a'] ) * 0.5
        self.n_params += 1
        newparam = 'o-si-si bondorder cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.05
        self.n_params += 1
        newparam = 'o-si-si bondorder cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.5
        self.n_params += 1
        newparam = 'o-si-o bondorder cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.05
        self.n_params += 1
        newparam = 'o-si-o bondorder cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.5
        self.n_params += 1
        newparam = 'si-o-si bondorder cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.05
        self.n_params += 1
        newparam = 'si-o-si bondorder cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.5
        self.n_params += 1
        newparam = 'si-o-o bondorder cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.05
        self.n_params += 1
        newparam = 'si-o-o bondorder cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.5
        self.n_params += 1
    
        newparam = 'si repulsive epsilon'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1830.8
        self.n_params += 1
        newparam = 'si repulsive lambda'
        self.names_of_params.append(newparam)
        self.params[newparam] = 2.4799
        self.n_params += 1
        newparam = 'si repulsive zeta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 2.4799
        self.n_params += 1
        newparam = 'si repulsive cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.0
        self.n_params += 1
        newparam = 'si repulsive cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.2
        self.n_params += 1
    
        newparam = 'o repulsive epsilon'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3326.699
        self.n_params += 1
        newparam = 'o repulsive lambda'
        self.names_of_params.append(newparam)
        self.params[newparam] = 5.36
        self.n_params += 1
        newparam = 'o repulsive zeta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 5.36
        self.n_params += 1
        newparam = 'o repulsive cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.0
        self.n_params += 1
        newparam = 'o repulsive cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.4
        self.n_params += 1
    
        newparam = 'si-o repulsive epsilon'
        self.names_of_params.append(newparam)
        self.params[newparam] = math.sqrt( \
            self.params['si repulsive epsilon'] * \
            self.params['o repulsive epsilon'] )
        self.n_params += 1
        newparam = 'si-o repulsive zeta'
        self.names_of_params.append(newparam)
        self.params[newparam] = ( \
            self.params['si repulsive zeta'] + \
            self.params['o repulsive zeta'] ) * 0.5
        self.n_params += 1
        newparam = 'si-o repulsive cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.05
        self.n_params += 1
        newparam = 'si-o repulsive cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.5
        self.n_params += 1

        newparam = 'si charge c1'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.0
        self.n_params += 1
        newparam = 'si charge c2'
        self.names_of_params.append(newparam)
        self.params[newparam] = 3.625144859
        self.n_params += 1
        newparam = 'si charge c3'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.0
        self.n_params += 1
        newparam = 'si charge c4'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.087067714
        self.n_params += 1

        newparam = 'o charge c1'
        self.names_of_params.append(newparam)
        self.params[newparam] = 5.63441383
        self.n_params += 1
        newparam = 'o charge c2'
        self.names_of_params.append(newparam)
        self.params[newparam] = 7.689598017
        self.n_params += 1
        newparam = 'o charge c3'
        self.names_of_params.append(newparam)
        self.params[newparam] = 4.51426991
        self.n_params += 1
        newparam = 'o charge c4'
        self.names_of_params.append(newparam)
        self.params[newparam] = 1.330079082
        self.n_params += 1
    
        newparam = 'si field rho1'
        self.names_of_params.append(newparam)
        self.params[newparam] = -0.499378
        self.n_params += 1
        newparam = 'si field rho2'
        self.names_of_params.append(newparam)
        self.params[newparam] = 2.999911
        self.n_params += 1
        newparam = 'o field rho1'
        self.names_of_params.append(newparam)
        self.params[newparam] = -3.922011000
        self.n_params += 1
        newparam = 'o field rho2'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.971086
        self.n_params += 1
        newparam = 'field cutoff'
        self.names_of_params.append(newparam)
        self.params[newparam] = 12.0
        self.n_params += 1
        newparam = 'field cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.1
        self.n_params += 1
        
        newparam = 'si slater zeta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.772871
        self.n_params += 1
        newparam = 'o slater zeta'
        self.names_of_params.append(newparam)
        self.params[newparam] = 2.243072
        self.n_params += 1
        newparam = 'slater cutoff margin'
        self.names_of_params.append(newparam)
        self.params[newparam] = 0.1
        self.n_params += 1
    

    def get_parameter_value(self,param_name):
        return self.params[param_name]

    def set_parameter_value(self,param_name,value):
        self.params[param_name] = value


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
            real_cutoff = self.params['coulomb cutoff']
            k_cutoff = 0.0
            sigma = self.params['coulomb sigma']
            epsilon = self.params['coulomb epsilon']
        else:
            #print "Ewald real"
            real_cutoff, k_cutoff, sigma, epsilon = estimate_ewald_parameters(self.params['coulomb cutoff'])
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
            bond_bend_pot.set_parameter_value('epsilon', self.params['si-o-si bend epsilon'])
            bond_bend_pot.set_parameter_value('theta_0', self.params['si-o-si bend theta'])
            bond_bend_pot.set_parameter_value('n',1)
            bond_bend_pot.set_parameter_value('m',2)
            bond_bend_pot.set_symbols(['Si','O','Si'])
            bond_bend_pot.set_cutoff(self.params['si-o-si bend cutoff'])
            bond_bend_pot.set_cutoff_margin(self.params['si-o-si bend cutoff margin'])

            self.pieces.append(bond_bend_pot)

        # O-Si-O bond bending
        if self.possible_excludes[1] not in self.excludes:
            #print "including O-Si-O bond bending"
            bond_bend_pot2 = pysic.Potential('bond_bend')
            bond_bend_pot2.set_parameter_value('epsilon', self.params['o-si-o bend epsilon'])
            bond_bend_pot2.set_parameter_value('theta_0', self.params['o-si-o bend theta'])
            bond_bend_pot2.set_parameter_value('n',1)
            bond_bend_pot2.set_parameter_value('m',2)
            bond_bend_pot2.set_symbols(['O','Si','O'])
            bond_bend_pot2.set_cutoff(self.params['o-si-o bend cutoff'])
            bond_bend_pot2.set_cutoff_margin(self.params['o-si-o bend cutoff margin'])

            self.pieces.append(bond_bend_pot2)
            
        # attractive Si-Si
        if self.possible_excludes[2] not in self.excludes:
            #print "including Si-Si attraction"
            B_ISiSi = pysic.Potential('charge_exp')
            B_ISiSi.set_parameter_value('epsilon', self.params['si attractive epsilon'])  #B_Si
            B_ISiSi.set_parameter_value('Rmax1', self.params['si rmax'])
            B_ISiSi.set_parameter_value('Rmin1', self.params['si rmin'])
            B_ISiSi.set_parameter_value('Qmax1', self.params['si qmax'])
            B_ISiSi.set_parameter_value('Qmin1', self.params['si qmin'])
            B_ISiSi.set_parameter_value('Rmax2', self.params['si rmax'])
            B_ISiSi.set_parameter_value('Rmin2', self.params['si rmin'])
            B_ISiSi.set_parameter_value('Qmax2', self.params['si qmax'])
            B_ISiSi.set_parameter_value('Qmin2', self.params['si qmin'])
            B_ISiSi.set_parameter_value('xi1', self.params['si attractive alpha'])  #alpha_Si
            B_ISiSi.set_parameter_value('xi2', self.params['si attractive alpha'])

            B_IISiSi = pysic.Potential('charge_abs')
            B_IISiSi.set_parameter_value('a1', self.params['si attractive a']) #aB
            B_IISiSi.set_parameter_value('b1', self.params['si attractive b']) #-abs(bB)**nB = -9.53674E-7
            B_IISiSi.set_parameter_value('Q1', self.params['si attractive q']) #QO
            B_IISiSi.set_parameter_value('n1', self.params['si attractive n']) #nB
            B_IISiSi.set_parameter_value('a2', self.params['si attractive a'])
            B_IISiSi.set_parameter_value('b2', self.params['si attractive b'])
            B_IISiSi.set_parameter_value('Q2', self.params['si attractive q']) #QO
            B_IISiSi.set_parameter_value('n2', self.params['si attractive n'])

            exp_potBSiSi = pysic.Potential('exponential')
            exp_potBSiSi.set_parameter_value('epsilon',-1)
            exp_potBSiSi.set_parameter_value('zeta', self.params['si attractive zeta'])  #alpha

            prodBSiSi = pysic.ProductPotential([B_ISiSi,B_IISiSi,exp_potBSiSi])

            prodBSiSi.set_cutoff(self.params['si attractive cutoff'])
            prodBSiSi.set_cutoff_margin(self.params['si attractive cutoff margin'])
            prodBSiSi.set_symbols(['Si','Si'])

            self.pieces.append(prodBSiSi)

            # Si-Si bond order factor
            if self.possible_excludes[5] not in self.excludes:
                #print "including Si-Si bond order factor"

                bijSiSi1 = pysic.BondOrderParameters('tersoff', symbols=['Si','Si','O'])
                bijSiSi1.set_parameter_value('beta', self.params['si bondorder beta'])  #all Si parameters
                bijSiSi1.set_parameter_value('eta', self.params['si bondorder eta'])
                bijSiSi1.set_parameter_value('mu', self.params['si bondorder mu'])
                bijSiSi1.set_parameter_value('a', self.params['si bondorder a'])
                bijSiSi1.set_parameter_value('c', self.params['si bondorder c'])
                bijSiSi1.set_parameter_value('d', self.params['si bondorder d'])
                bijSiSi1.set_parameter_value('h', self.params['si bondorder h'])
                bijSiSi1.set_cutoff(self.params['si-si-o bondorder cutoff'])
                bijSiSi1.set_cutoff_margin(self.params['si-si-o bondorder cutoff margin'])

                bijSiSi2 = pysic.BondOrderParameters('tersoff', symbols=['Si','Si','Si'])
                bijSiSi2.set_parameter_value('beta', self.params['si bondorder beta'])  #all Si parameters
                bijSiSi2.set_parameter_value('eta', self.params['si bondorder eta'])
                bijSiSi2.set_parameter_value('mu', self.params['si bondorder mu'])
                bijSiSi2.set_parameter_value('a', self.params['si bondorder a'])
                bijSiSi2.set_parameter_value('c', self.params['si bondorder c'])
                bijSiSi2.set_parameter_value('d', self.params['si bondorder d'])
                bijSiSi2.set_parameter_value('h', self.params['si bondorder h'])
                bijSiSi2.set_cutoff(self.params['si-si-si bondorder cutoff'])
                bijSiSi2.set_cutoff_margin(self.params['si-si-si bondorder cutoff margin'])

                coordSiSi = pysic.Coordinator([bijSiSi1,bijSiSi2])

                prodBSiSi.set_coordinator(coordSiSi)


        # attractive O-O
        if self.possible_excludes[3] not in self.excludes:
            #print "including O-O attraction"

            B_IOO = pysic.Potential('charge_exp')
            B_IOO.set_parameter_value('epsilon', self.params['o attractive epsilon'])  #B_O
            B_IOO.set_parameter_value('Rmax1', self.params['o rmax'])
            B_IOO.set_parameter_value('Rmin1', self.params['o rmin'])
            B_IOO.set_parameter_value('Qmax1', self.params['o qmax'])
            B_IOO.set_parameter_value('Qmin1', self.params['o qmin'])
            B_IOO.set_parameter_value('Rmax2', self.params['o rmax'])
            B_IOO.set_parameter_value('Rmin2', self.params['o rmin'])
            B_IOO.set_parameter_value('Qmax2', self.params['o qmax'])
            B_IOO.set_parameter_value('Qmin2', self.params['o qmin'])
            B_IOO.set_parameter_value('xi1', self.params['o attractive alpha'])  #alpha_O
            B_IOO.set_parameter_value('xi2', self.params['o attractive alpha'])

            B_IIOO = pysic.Potential('charge_abs')
            B_IIOO.set_parameter_value('a1', self.params['o attractive a']) #aB
            B_IIOO.set_parameter_value('b1', self.params['o attractive b'])# #-abs(bB)**nB = -2.41926E-6
            B_IIOO.set_parameter_value('Q1', self.params['o attractive q']) #QO
            B_IIOO.set_parameter_value('n1', self.params['o attractive n']) #nB
            B_IIOO.set_parameter_value('a2', self.params['o attractive a'])#
            B_IIOO.set_parameter_value('b2', self.params['o attractive b'])#
            B_IIOO.set_parameter_value('Q2', self.params['o attractive q']) #QO
            B_IIOO.set_parameter_value('n2', self.params['o attractive n'])
            
            exp_potBOO = pysic.Potential('exponential')
            exp_potBOO.set_parameter_value('epsilon',-1)
            exp_potBOO.set_parameter_value('zeta', self.params['o attractive zeta']) #alpha

            prodBOO = pysic.ProductPotential([B_IOO,B_IIOO,exp_potBOO])
            
            prodBOO.set_cutoff(self.params['o attractive cutoff'])
            prodBOO.set_cutoff_margin(self.params['o attractive cutoff margin'])
            prodBOO.set_symbols(['O','O'])

            self.pieces.append(prodBOO)
            
            # O-O bond order factor
            if self.possible_excludes[6] not in self.excludes:
                #print "including O-O bond order factor"

                bijOO1 = pysic.BondOrderParameters('tersoff', symbols=['O','O','Si'])
                bijOO1.set_parameter_value('beta', self.params['o bondorder beta']) #all O parameters
                bijOO1.set_parameter_value('eta', self.params['o bondorder eta'])
                bijOO1.set_parameter_value('mu', self.params['o bondorder mu'])
                bijOO1.set_parameter_value('a', self.params['o bondorder a'])
                bijOO1.set_parameter_value('c', self.params['o bondorder c'])
                bijOO1.set_parameter_value('d', self.params['o bondorder d'])
                bijOO1.set_parameter_value('h', self.params['o bondorder h'])
                bijOO1.set_cutoff(self.params['o-o-si bondorder cutoff'])
                bijOO1.set_cutoff_margin(self.params['o-o-si bondorder cutoff margin'])

                bijOO2 = pysic.BondOrderParameters('tersoff', symbols=['O','O','O'])
                bijOO2.set_parameter_value('beta', self.params['o bondorder beta']) #all O parameters
                bijOO2.set_parameter_value('eta', self.params['o bondorder eta'])
                bijOO2.set_parameter_value('mu', self.params['o bondorder mu'])
                bijOO2.set_parameter_value('a', self.params['o bondorder a'])
                bijOO2.set_parameter_value('c', self.params['o bondorder c'])
                bijOO2.set_parameter_value('d', self.params['o bondorder d'])
                bijOO2.set_parameter_value('h', self.params['o bondorder h'])
                bijOO2.set_cutoff(self.params['o-o-o bondorder cutoff'])
                bijOO2.set_cutoff_margin(self.params['o-o-o bondorder cutoff margin'])

                coordOO = pysic.Coordinator([bijOO1,bijOO2])

                prodBOO.set_coordinator(coordOO)



        # attractive Si-O
        if self.possible_excludes[4] not in self.excludes:
            #print "including Si-O attraction"

            B_ISiO = pysic.Potential('charge_exp')
            B_ISiO.set_parameter_value('epsilon', self.params['si-o attractive epsilon'])  # sqrt(B_Si*B_O)
            B_ISiO.set_parameter_value('Rmax1', self.params['si rmax'])
            B_ISiO.set_parameter_value('Rmin1', self.params['si rmin'])
            B_ISiO.set_parameter_value('Qmax1', self.params['si qmax'])
            B_ISiO.set_parameter_value('Qmin1', self.params['si qmin'])
            B_ISiO.set_parameter_value('Rmax2', self.params['o rmax'])
            B_ISiO.set_parameter_value('Rmin2', self.params['o rmin'])
            B_ISiO.set_parameter_value('Qmax2', self.params['o qmax'])
            B_ISiO.set_parameter_value('Qmin2', self.params['o qmin'])
            B_ISiO.set_parameter_value('xi1', self.params['si attractive alpha'])  #alpha_Si
            B_ISiO.set_parameter_value('xi2', self.params['o attractive alpha'])  #alpha O

            B_IISiO = pysic.Potential('charge_abs')
            B_IISiO.set_parameter_value('a1', self.params['si attractive a']) #aB_Si
            B_IISiO.set_parameter_value('b1', self.params['si attractive b']) #-abs(bB)**nB for Si
            B_IISiO.set_parameter_value('Q1', self.params['si attractive q']) #QO of Si
            B_IISiO.set_parameter_value('n1', self.params['si attractive n']) #nB
            B_IISiO.set_parameter_value('a2', self.params['o attractive a'])# #aB_O
            B_IISiO.set_parameter_value('b2', self.params['o attractive b'])# #-abs(bB)**nB for O
            B_IISiO.set_parameter_value('Q2', self.params['o attractive q']) #QO of O
            B_IISiO.set_parameter_value('n2', self.params['o attractive n'])
            
            exp_potBSiO = pysic.Potential('exponential')
            exp_potBSiO.set_parameter_value('epsilon',-1) 
            exp_potBSiO.set_parameter_value('zeta', self.params['si-o attractive zeta'])  #alpha_SiO = (alpha_Si+alpha_O)/2

            prodBSiO = pysic.ProductPotential([B_ISiO,B_IISiO,exp_potBSiO])
            
            prodBSiO.set_cutoff(self.params['si-o attractive cutoff'])
            prodBSiO.set_cutoff_margin(self.params['si-o attractive cutoff margin'])
            prodBSiO.set_symbols(['Si','O'])

            self.pieces.append(prodBSiO)
            
            # Si-O bond order factor
            if self.possible_excludes[7] not in self.excludes:
                #print "including Si-O bond order factor"

                bijOSi1 = pysic.BondOrderParameters('tersoff', symbols=['O','Si','Si'])
                bijOSi1.set_parameter_value('beta', self.params['o bondorder beta'])  #O
                bijOSi1.set_parameter_value('eta', self.params['o bondorder eta']) #O
                bijOSi1.set_parameter_value('mu', self.params['o bondorder mu']) #O
                bijOSi1.set_parameter_value('a', self.params['si-o bondorder a']) #media between Si and O parameters
                bijOSi1.set_parameter_value('c', self.params['o bondorder c']) # O
                bijOSi1.set_parameter_value('d', self.params['o bondorder d']) #O
                bijOSi1.set_parameter_value('h', self.params['o bondorder h']) #O
                bijOSi1.set_cutoff(self.params['o-si-si bondorder cutoff'])
                bijOSi1.set_cutoff_margin(self.params['o-si-si bondorder cutoff margin'])

                bijOSi2 = pysic.BondOrderParameters('tersoff', symbols=['O','Si','O'])
                bijOSi2.set_parameter_value('beta', self.params['o bondorder beta'])  #O
                bijOSi2.set_parameter_value('eta', self.params['o bondorder eta']) #O
                bijOSi2.set_parameter_value('mu', self.params['o bondorder mu']) #O
                bijOSi2.set_parameter_value('a', self.params['si-o bondorder a']) #media between Si and O parameters
                bijOSi2.set_parameter_value('c', self.params['o bondorder c']) # O
                bijOSi2.set_parameter_value('d', self.params['o bondorder d']) #O
                bijOSi2.set_parameter_value('h', self.params['o bondorder h']) #O
                bijOSi2.set_cutoff(self.params['o-si-o bondorder cutoff'])
                bijOSi2.set_cutoff_margin(self.params['o-si-o bondorder cutoff margin'])

                bijSiO1 = pysic.BondOrderParameters('tersoff', symbols=['Si','O','O'])
                bijSiO1.set_parameter_value('beta', self.params['si bondorder beta'])  #Si
                bijSiO1.set_parameter_value('eta', self.params['si bondorder eta']) #Si
                bijSiO1.set_parameter_value('mu', self.params['si bondorder mu']) #Si
                bijSiO1.set_parameter_value('a', self.params['si-o bondorder a'])   #media between Si and O parameters
                bijSiO1.set_parameter_value('c', self.params['si bondorder c']) #Si
                bijSiO1.set_parameter_value('d', self.params['si bondorder d']) #Si
                bijSiO1.set_parameter_value('h', self.params['si bondorder h']) #Si
                bijSiO1.set_cutoff(self.params['si-o-o bondorder cutoff'])
                bijSiO1.set_cutoff_margin(self.params['si-o-o bondorder cutoff margin'])

                bijSiO2 = pysic.BondOrderParameters('tersoff', symbols=['Si','O','Si'])
                bijSiO2.set_parameter_value('beta', self.params['si bondorder beta'])  #Si
                bijSiO2.set_parameter_value('eta', self.params['si bondorder eta']) #Si
                bijSiO2.set_parameter_value('mu', self.params['si bondorder mu']) #Si
                bijSiO2.set_parameter_value('a', self.params['si-o bondorder a'])   #media between Si and O parameters
                bijSiO2.set_parameter_value('c', self.params['si bondorder c']) #Si
                bijSiO2.set_parameter_value('d', self.params['si bondorder d']) #Si
                bijSiO2.set_parameter_value('h', self.params['si bondorder h']) #Si
                bijSiO2.set_cutoff(self.params['si-o-si bondorder cutoff'])
                bijSiO2.set_cutoff_margin(self.params['si-o-si bondorder cutoff margin'])
                
                coordSiO = pysic.Coordinator([bijOSi1,bijOSi2, bijSiO1, bijSiO2])

                prodBSiO.set_coordinator(coordSiO)

            
        # repulsive Si-Si
        if self.possible_excludes[8] not in self.excludes:
            #print "including Si-Si repulsion"

            ASiSi = pysic.Potential('charge_exp')
            ASiSi.set_parameter_value('epsilon', self.params['si repulsive epsilon'])  #A_Si
            ASiSi.set_parameter_value('Rmax1', self.params['si rmax'])
            ASiSi.set_parameter_value('Rmin1', self.params['si rmin'])
            ASiSi.set_parameter_value('Qmax1', self.params['si qmax'])
            ASiSi.set_parameter_value('Qmin1', self.params['si qmin'])
            ASiSi.set_parameter_value('Rmax2', self.params['si rmax'])
            ASiSi.set_parameter_value('Rmin2', self.params['si rmin'])
            ASiSi.set_parameter_value('Qmax2', self.params['si qmax'])
            ASiSi.set_parameter_value('Qmin2', self.params['si qmin'])
            ASiSi.set_parameter_value('xi1', self.params['si repulsive lambda']) #lambda
            ASiSi.set_parameter_value('xi2', self.params['si repulsive lambda'])

            exp_potSi = pysic.Potential('exponential')
            exp_potSi.set_parameter_value('epsilon',1)
            exp_potSi.set_parameter_value('zeta',  self.params['si repulsive zeta'])  #lambda_Si
            
            prodASiSi = pysic.ProductPotential([ASiSi,exp_potSi])

            prodASiSi.set_cutoff(self.params['si repulsive cutoff'])
            prodASiSi.set_cutoff_margin(self.params['si repulsive cutoff margin'])
            prodASiSi.set_symbols(['Si','Si'])

            self.pieces.append(prodASiSi)

        # repulsive O-O
        if self.possible_excludes[9] not in self.excludes:
            #print "including O-O repulsion"

            AOO = pysic.Potential('charge_exp')
            AOO.set_parameter_value('epsilon', self.params['o repulsive epsilon'])  #A_O
            AOO.set_parameter_value('Rmax1', self.params['o rmax'])
            AOO.set_parameter_value('Rmin1', self.params['o rmin'])
            AOO.set_parameter_value('Qmax1', self.params['o qmax'])
            AOO.set_parameter_value('Qmin1', self.params['o qmin'])
            AOO.set_parameter_value('Rmax2', self.params['o rmax'])
            AOO.set_parameter_value('Rmin2', self.params['o rmin'])
            AOO.set_parameter_value('Qmax2', self.params['o qmax'])
            AOO.set_parameter_value('Qmin2', self.params['o qmin'])
            AOO.set_parameter_value('xi1', self.params['o repulsive lambda']) #lambda
            AOO.set_parameter_value('xi2', self.params['o repulsive lambda'])

            exp_potO = pysic.Potential('exponential')
            exp_potO.set_parameter_value('epsilon',1)
            exp_potO.set_parameter_value('zeta', self.params['o repulsive zeta'])  #lambda_O

            prodAOO = pysic.ProductPotential([AOO,exp_potO])
 
            prodAOO.set_cutoff(self.params['o repulsive cutoff'])
            prodAOO.set_cutoff_margin(self.params['o repulsive cutoff margin'])
            prodAOO.set_symbols(['O','O'])

            self.pieces.append(prodAOO)

            
        # repulsive Si-O
        if self.possible_excludes[10] not in self.excludes:
            #print "including Si-O repulsion"

            ASiO = pysic.Potential('charge_exp')
            ASiO.set_parameter_value('epsilon', self.params['si-o repulsive epsilon']) #sqrt(A_Si*A_O)
            ASiO.set_parameter_value('Rmax1', self.params['si rmax'])
            ASiO.set_parameter_value('Rmin1', self.params['si rmin'])
            ASiO.set_parameter_value('Qmax1', self.params['si qmax'])
            ASiO.set_parameter_value('Qmin1', self.params['si qmin'])
            ASiO.set_parameter_value('Rmax2', self.params['o rmax'])
            ASiO.set_parameter_value('Rmin2', self.params['o rmin'])
            ASiO.set_parameter_value('Qmax2', self.params['o qmax'])
            ASiO.set_parameter_value('Qmin2', self.params['o qmin'])
            ASiO.set_parameter_value('xi1', self.params['si repulsive lambda']) #lambda
            ASiO.set_parameter_value('xi2', self.params['o repulsive lambda'])
            
            exp_potSiO = pysic.Potential('exponential')
            exp_potSiO.set_parameter_value('epsilon',1)
            exp_potSiO.set_parameter_value('zeta', self.params['si-o repulsive zeta'])  # (lambda_Si + lambda_O)/2

            prodASiO = pysic.ProductPotential([ASiO,exp_potSiO])

            prodASiO.set_cutoff(self.params['si-o repulsive cutoff'])
            prodASiO.set_cutoff_margin(self.params['si-o repulsive cutoff margin'])
            prodASiO.set_symbols(['Si','O'])

            self.pieces.append(prodASiO)


        # charge self-energy Si
        if self.possible_excludes[11] not in self.excludes:
            #print "including Si self energy"

            chiSiX = pysic.Potential('charge_self')
            chiSiX.set_parameter_value('epsilon', self.params['si charge c1'])
            chiSiX.set_parameter_value('n',1)
            chiSiX.set_symbols(['Si'])
            chiSiJ = pysic.Potential('charge_self')
            chiSiJ.set_parameter_value('epsilon', self.params['si charge c2'])
            chiSiJ.set_parameter_value('n',2)
            chiSiJ.set_symbols(['Si'])
            chiSiK = pysic.Potential('charge_self')
            chiSiK.set_parameter_value('epsilon', self.params['si charge c3'])
            chiSiK.set_parameter_value('n',3)
            chiSiK.set_symbols(['Si'])
            chiSiL = pysic.Potential('charge_self')
            chiSiL.set_parameter_value('epsilon', self.params['si charge c4'])
            chiSiL.set_parameter_value('n',4)
            chiSiL.set_symbols(['Si'])
            
            self.pieces.append([chiSiX])
            self.pieces.append([chiSiJ])
            self.pieces.append([chiSiK])
            self.pieces.append([chiSiL])
            
        # charge self-energy O
        if self.possible_excludes[12] not in self.excludes:
            #print "including O self energy"

            chiOX = pysic.Potential('charge_self')
            chiOX.set_parameter_value('epsilon', self.params['o charge c1'])
            chiOX.set_parameter_value('n',1)
            chiOX.set_symbols(['O'])
            chiOJ = pysic.Potential('charge_self')
            chiOJ.set_parameter_value('epsilon', self.params['o charge c2'])
            chiOJ.set_parameter_value('n',2)
            chiOJ.set_symbols(['O'])
            chiOK = pysic.Potential('charge_self')
            chiOK.set_parameter_value('epsilon', self.params['o charge c3'])
            chiOK.set_parameter_value('n', 3)
            chiOK.set_symbols(['O'])
            chiOL = pysic.Potential('charge_self')
            chiOL.set_parameter_value('epsilon', self.params['o charge c4'])
            chiOL.set_parameter_value('n', 4)
            chiOL.set_symbols(['O'])

            self.pieces.append(chiOX)
            self.pieces.append(chiOJ)
            self.pieces.append(chiOK)
            self.pieces.append(chiOL)

        pow1 = pysic.Potential('power') # just multiplies by r
        pow1.set_parameter_value('epsilon',1)
        pow1.set_parameter_value('a',1)
        pow1.set_parameter_value('n',-1)
        
        pow2 = pysic.Potential('power') # just multiplies by (1/r^5)
        pow2.set_parameter_value('epsilon',1)
        pow2.set_parameter_value('a',1)
        pow2.set_parameter_value('n',5) 

        chp4 = pysic.Potential('charge_pair')
        chp4.set_parameter_value('epsilon', self.params['si field rho1']) # is (rho1_Si)
        chp4.set_parameter_value('n1',0)
        chp4.set_parameter_value('n2',1) # multiplies by q_j

        chp4b = pysic.Potential('charge_pair')
        chp4b.set_parameter_value('epsilon', self.params['si field rho2']) # is (rho2_Si)
        chp4b.set_parameter_value('n1',0)
        chp4b.set_parameter_value('n2',2) # multiplies by q_j^2

        chp4c = pysic.Potential('charge_pair')
        chp4c.set_parameter_value('epsilon', self.params['si field rho1'])
        chp4c.set_parameter_value('n1',1)
        chp4c.set_parameter_value('n2',0) # multiplies by q_i

        chp4d = pysic.Potential('charge_pair')
        chp4d.set_parameter_value('epsilon', self.params['si field rho2'])
        chp4d.set_parameter_value('n1',2)
        chp4d.set_parameter_value('n2',0) # multiplies by q_i^2

        
        # penalty function Si-Si
        if self.possible_excludes[13] not in self.excludes:
            #print "including Si-Si penalty"

            chp4.set_symbols(['Si','Si'])
            prod2 = pysic.ProductPotential([chp4, pow2])
            prod2.set_symbols(['Si','Si'])
            prod2.set_cutoff(self.params['field cutoff'])
            prod2.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2)

            chp4b.set_symbols(['Si','Si'])
            prod2b = pysic.ProductPotential([chp4b, pow2])
            prod2b.set_symbols(['Si','Si'])
            prod2b.set_cutoff(self.params['field cutoff'])
            prod2b.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2b)

            chp4c.set_symbols(['Si','Si'])
            prod2c = pysic.ProductPotential([chp4c, pow2])
            prod2c.set_symbols(['Si','Si'])
            prod2c.set_cutoff(self.params['field cutoff'])
            prod2c.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2c)

            chp4d.set_symbols(['Si','Si'])
            prod2d = pysic.ProductPotential([chp4d, pow2])
            prod2d.set_symbols(['Si','Si'])
            prod2d.set_cutoff(self.params['field cutoff'])
            prod2d.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2d)


        chp4O = pysic.Potential('charge_pair')
        chp4O.set_parameter_value('epsilon', self.params['o field rho1'])
        chp4O.set_parameter_value('n1',0)
        chp4O.set_parameter_value('n2',1)          # multiplies by q_j

        chp4bO = pysic.Potential('charge_pair')
        chp4bO.set_parameter_value('epsilon', self.params['o field rho2'])
        chp4bO.set_parameter_value('n1',0)
        chp4bO.set_parameter_value('n2',2)          # multiplies by q_j^2

        chp4cO = pysic.Potential('charge_pair')
        chp4cO.set_parameter_value('epsilon', self.params['o field rho1'])
        chp4cO.set_parameter_value('n1',1)
        chp4cO.set_parameter_value('n2',0)          # multiplies by q_i

        chp4dO = pysic.Potential('charge_pair')
        chp4dO.set_parameter_value('epsilon', self.params['o field rho2'])
        chp4dO.set_parameter_value('n1',2)
        chp4dO.set_parameter_value('n2',0)          # multiplies by q_i^2 


        # penalty function O-O
        if self.possible_excludes[14] not in self.excludes: 
            #print "including O-O penalty"

            prod2O = pysic.ProductPotential([chp4O, pow2])
            prod2O.set_symbols(['O','O'])
            prod2O.set_cutoff(self.params['field cutoff'])
            prod2O.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2O)

            prod2bO = pysic.ProductPotential([chp4bO, pow2])
            prod2bO.set_symbols(['O','O'])
            prod2bO.set_cutoff(self.params['field cutoff'])
            prod2bO.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2bO)

            prod2cO = pysic.ProductPotential([chp4cO, pow2])
            prod2cO.set_symbols(['O','O'])
            prod2cO.set_cutoff(self.params['field cutoff'])
            prod2cO.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2cO)

            prod2dO = pysic.ProductPotential([chp4dO, pow2])
            prod2dO.set_symbols(['O','O'])
            prod2dO.set_cutoff(self.params['field cutoff'])
            prod2dO.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2dO)


        # penalty function Si-O
        if self.possible_excludes[15] not in self.excludes:
            #print "including Si-O penalty"

            prod2mix = pysic.ProductPotential([chp4, pow2])
            prod2mix.set_symbols(['Si','O'])
            prod2mix.set_cutoff(self.params['field cutoff'])
            prod2mix.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2mix)

            prod2bmix = pysic.ProductPotential([chp4b, pow2])
            prod2bmix.set_symbols(['Si','O'])
            prod2bmix.set_cutoff(self.params['field cutoff'])
            prod2bmix.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2bmix)

            prod2cmix = pysic.ProductPotential([chp4cO, pow2])
            prod2cmix.set_symbols(['Si','O'])
            prod2cmix.set_cutoff(self.params['field cutoff'])
            prod2cmix.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2cmix)

            prod2dmix = pysic.ProductPotential([chp4dO, pow2])
            prod2dmix.set_symbols(['Si','O'])
            prod2dmix.set_cutoff(self.params['field cutoff'])
            prod2dmix.set_cutoff_margin(self.params['field cutoff margin'])
            self.pieces.append(prod2dmix)

        charged = pysic.Potential('charge_pair')
        charged.set_parameter_value('epsilon',1)
        charged.set_parameter_value('n1',1)
        charged.set_parameter_value('n2',1)

        decay = pysic.Potential('power')
        decay.set_parameter_value('epsilon', 14.3996)
        decay.set_parameter_value('a',1)
        decay.set_parameter_value('n',1)

        si_zeta = self.params['si slater zeta']
        o_zeta  = self.params['o slater zeta']

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
            

            coul_sisi = []
            if self.possible_excludes[19] not in self.excludes:
                #print "including Si-Si direct coulomb"
                coul_sisi.append( pysic.ProductPotential([charged,decay]) )
            coul_sisi.append( pysic.ProductPotential([charged,decay,si_exp_decay]) )
            coul_sisi.append( pysic.ProductPotential([charged,decay,si_exp_decay,poly_sisi1]) )
            coul_sisi.append( pysic.ProductPotential([charged,decay,si_exp_decay,poly_sisi2]) )
            coul_sisi.append( pysic.ProductPotential([charged,decay,si_exp_decay,poly_sisi3]) )

            for c_sisi in coul_sisi:
                c_sisi.set_symbols(['Si','Si'])
                c_sisi.set_cutoff(self.params['coulomb cutoff'])
                c_sisi.set_cutoff_margin(self.params['slater cutoff margin'])
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

            coul_oo = []
            if self.possible_excludes[19] not in self.excludes:
                #print "including O-O direct coulomb"
                coul_oo.append( pysic.ProductPotential([charged,decay]) )
            coul_oo.append( pysic.ProductPotential([charged,decay,o_exp_decay]) )
            coul_oo.append( pysic.ProductPotential([charged,decay,o_exp_decay,poly_oo1]) )
            coul_oo.append( pysic.ProductPotential([charged,decay,o_exp_decay,poly_oo2]) )
            coul_oo.append( pysic.ProductPotential([charged,decay,o_exp_decay,poly_oo3]) )
            
            for c_oo in coul_oo:
                c_oo.set_symbols(['O','O'])
                c_oo.set_cutoff(self.params['coulomb cutoff'])
                c_oo.set_cutoff(self.params['slater cutoff margin'])
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

            coul_sio = []
            if self.possible_excludes[19] not in self.excludes:
                #print "including Si-O direct coulomb"
                coul_sio.append( pysic.ProductPotential([charged,decay]) )
            coul_sio.append( pysic.ProductPotential([charged,decay,scaled_si_exp_decay]) )
            coul_sio.append( pysic.ProductPotential([charged,decay,scaled_o_exp_decay]) )
            coul_sio.append( pysic.ProductPotential([charged,decay,poly_sio1,si_exp_decay]) )
            coul_sio.append( pysic.ProductPotential([charged,decay,poly_sio2,o_exp_decay]) )

            for c_sio in coul_sio:
                c_sio.set_symbols(['Si','O'])
                c_sio.set_cutoff(self.params['coulomb cutoff'])
                c_sio.set_cutoff_margin(self.params['slater cutoff margin'])
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
