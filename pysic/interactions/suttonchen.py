#! /usr/bin/env python

from pysic.core import *
from pysic.utility.error import InvalidPotentialError
from pysic.interactions.local import Potential, ProductPotential
from pysic.interactions.compound import CompoundPotential
from pysic.interactions.bondorder import Coordinator, BondOrderParameters
from pysic.interactions.coulomb import CoulombSummation
import copy

class SuttonChenPotential(CompoundPotential):
    """Class representing the Sutton-Chen potential as a :class:`~pysic.interactions.local.CompoundPotential` object.
        
        Parameters:
        
        symbols: list of string
            the chemical symbols (elements) on which the potential acts
        tags: integer
            atoms with specific tags on which the potential acts
        indices: list of integers
            atoms with specific indices on which the potential acts
        parameters: list of doubles
            a list of parameters for characterizing the potential; their meaning depends on the type of potential
        cutoff: double
            the maximum atomic separation at which the potential is applied
        cutoff_margin: double
            the margin in which the potential is smoothly truncated to zero
    """

    def __init__(self,
                 symbols=None,
                 tags=None,
                 indices=None,
                 parameters=None,
                 cutoff=0.0,
                 cutoff_margin=0.0):
        super(SuttonChenPotential, self).__init__(n_params=5,
                                                 n_targets=2,
                                                 symbols=symbols,
                                                 tags=tags,
                                                 indices=indices,
                                                 parameters=parameters,
                                                 cutoff=cutoff,
                                                 cutoff_margin=cutoff_margin)
        self.density_symbols = None
        self.set_potential_type('suttonchen')
        self.names_of_params = ['epsilon', 'a', 'c', 'n', 'm']
        self.descriptions_of_params = ['energy scale constant',
                                       'lenght scale or lattice constant',
                                       'density coefficient',
                                       'power decay exponent',
                                       'density decay exponent']
        self.description = """
A potential of type
            
 U = epsilon [ sum_ij (a/r)^n - c sum_i sqrt(rho) ],
            
where
            
 rho = sum_j (a/r)^m.
            
and r is the interatomic distance.
"""

    def set_density_symbols(self,symbols):
        self.density_symbols = symbols

    def define_elements(self):        

        # initialize the collection of elemental potentials
        self.pieces = []
        
        # create the pairwise power decay potential
        SC_pairwise_potential = Potential('power',
                                                cutoff=self.cutoff,
                                                cutoff_margin=self.cutoff_margin,
                                                symbols=self.symbols,
                                                tags=self.tags,
                                                indices=self.indices)
        SC_pairwise_potential.set_parameter_value('epsilon',self.parameters[0])
        SC_pairwise_potential.set_parameter_value('a',self.parameters[1])
        SC_pairwise_potential.set_parameter_value('n',self.parameters[3])

        # add the potential to the list of pieces
        self.pieces.append(SC_pairwise_potential)


        # get the possible targets of the potential:
        # all different symbols, tags, indices
        all_symbols = []
        symbol_pairs = []
        try:
            for symblist in self.symbols:
                for symb in symblist:
                    if all_symbols.count(symb) == 0:
                        all_symbols.append(symb)
                symbol_pairs.append(symblist)
                if symblist[0] != symblist[1] and len(symblist) == 2:
                    symbol_pairs.append([symblist[1], symblist[0]])
        except:
            pass
             
        try:
            all_tags = []
            tag_pairs = []
            for taglist in self.tags:
                for tagi in taglist:
                    if all_tags.count(tagi) == 0:
                        all_tags.append(tagi)
                tag_pairs.append(taglist)
                if taglist[0] != taglist[1] and len(taglist) == 2:
                    tag_pairs.append([taglist[1], taglist[0]])
        except:
            pass

        try:
            all_indices = []
            index_pairs = []
            for indlist in self.indices:
                for ind in indlist:
                    if all_indices.count(ind) == 0:
                        all_indices.append(ind)
                ind_pairs.append(indlist)
                if indlist[0] != indlist[1] and len(indlist) == 2:
                    ind_pairs.append([indlist[1], indlist[0]])
        except:
            pass


        # set up potential sets with the correct set of symbols/tags/indices to avoid excess loops
        lists = [all_symbols, all_tags, all_indices]
        for type in [0,1,2]: # loop over symbols, tags, indices
            for symbol in lists[type]:

                # create the density potential as a combination of a unit potential and a bond order term
                SC_density_potential = Potential('constant',
                                                 parameters = [self.parameters[0]] )
    
                # check that the density term should be created for this symbol
                if self.density_symbols is None or self.density_symbols.count(symbol) > 0:
                
                    # create the bond factor
                    SC_density_scaler = BondOrderParameters('sqrt_scale')
                    SC_density_scaler.set_parameter_value('epsilon',-self.parameters[2])

                    # pick the pairs where this symbol is the leader
                    right_pairs = []
                    for pair in symbol_pairs:
                        if pair[0] == symbol:
                            right_pairs.append(pair)
                    
                    SC_density_factor = BondOrderParameters('power_bond',
                                                            cutoff=self.cutoff,
                                                            cutoff_margin=self.cutoff_margin)
                    SC_density_factor.set_parameter_value('a',self.parameters[1])
                    SC_density_factor.set_parameter_value('n',self.parameters[4])

                    if type == 0:
                        SC_density_potential.set_symbols([[symbol]])
                        SC_density_scaler.set_symbols([[symbol]])
                        SC_density_factor.set_symbols(right_pairs)
                    elif type == 1:
                        SC_density_potential.set_tags([[symbol]])
                        SC_density_scaler.set_tags([[symbol]])
                        SC_density_factor.set_tags(right_pairs)
                    else:
                        SC_density_potential.set_indices([[symbol]])
                        SC_density_scaler.set_indices([[symbol]])
                        SC_density_factor.set_indices(right_pairs)
                
                    SC_density_coordinator = Coordinator([SC_density_scaler,SC_density_factor])
                    SC_density_potential.set_coordinator(SC_density_coordinator)
                        
                    # add the potential to the list of pieces
                    self.pieces.append(SC_density_potential)

    def describe(self):
        super(SuttonChenPotential,self).describe()
        if not self.density_symbols is None:
            message = "density target symbols: "
            for ele in self.density_symbols:
                message += " '{0}' ".format(str(ele))

            print message
