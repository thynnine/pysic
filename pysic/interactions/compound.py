#! /usr/bin/env python

from pysic.core import *
from pysic.utility.error import InvalidPotentialError
from pysic.interactions.local import Potential, ProductPotential
import copy

class CompoundPotential(Potential):
    """Class representing an interaction constructed of several :class:`~pysic.interactions.local.Potential` and :class:`~pysic.interactions.local.ProductPotential` objects.
        
        This class implements the framework for wrapping complicated potentials as
        Python objects.
        The class itself implements an empty potential. To utilize the compound potential
        for real calculations, subclasses should be used for defining the actual contents
        of the potential.
        
        Parameters:
        
        n_targets: integer
            number of targets the potential acts on
        n_params: integer
            number of parameters needed for describing the potential
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
                 n_targets,
                 n_params,
                 symbols=None,
                 tags=None,
                 indices=None,
                 parameters=None,
                 cutoff=0.0,
                 cutoff_margin=0.0):
        self.n_targets = n_targets
        self.n_params = n_params
        self.symbols = None
        self.tags = None
        self.indices = None
        self.cutoff = cutoff
        self.cutoff_margin = 0.0
        if parameters is None:
            self.parameters = [0.0]*n_params
        else:
            self.parameters = parameters
        self.set_cutoff_margin(cutoff_margin)
        self.set_symbols(symbols)
        self.set_tags(tags)
        self.set_indices(indices)
        self.potential_type = None
        self.names_of_params = ["anonymous"]*n_params
        self.description_of_params = ["description missing"]*n_params
        self.description = "This compound potential has no description."
        self.pieces = []

    
    def __eq__(self,other):
        try:
            if self.potential_type != other.potential_type:
                return False
            if self.symbols != other.symbols:
                return False
            if self.tags != other.tags:
                return False
            if self.indices != other.indices:
                return False
            if self.parameters != other.parameters:
                return False
            if self.cutoff != other.cutoff:
                return False
            if len(self.pieces) != len(other.pieces):
                return False
            for p1, p2 in zip(self.pieces, other.pieces):
                if(p1 != p2):
                    return False
        except:
            return False

        return True

    def __ne__(self,other):
        return not self.__eq__(other)

    def get_number_of_parameters(self):
        """Returns the number of parameters the potential accepts.
        """
        return self.n_params

    def set_potential_type(self,type):
        """Sets the name of the potential.
        
        This can be used by subclasses to set the name of the potential.
        Since normal potentials are distinguished by their names (keywords)
        compound potentials should also have a name for consistency.
        
        Parameters:
        
        type: string
            The name of the potential
        """
        self.potential_type = type

    def get_elements(self):
        """Returns the elemental potentials this compund potential consists of.
            """
        return self.pieces
    
    def build(self,calculator):
        """Constructs the potential for the calculator as a collection of elemental potentials.
        
        The method adds all the elemental potentials it consists of one by one through
        the :meth:`pysic.calculator.Pysic.add_potential` method.
            
        Parameters:
            
        calculator: :class:`~pysic.calculator.Pysic` object
            the Pysic calculator to which the potential is added
        """
        
        self.define_elements()
        for potential in self.pieces:
            calculator.add_potential(potential)

    def remove(self,calculator):
        """Removes the elemental potentials this compound contains from the calculator.
        
        The method removes all the elemental potentials it consists of one by one through
        the :meth:`pysic.calculator.Pysic.remove_potential` method.
        If a match is not found in the set of potentials in the calculator, a message will be 
        printed but the removal process continues.
        
        Parameters:
                        
        calculator: :class:`~pysic.calculator.Pysic` object
            the Pysic calculator from which the potential is removed
        """

        self.define_elements()
        for pot in self.pieces:
            try:
                calculator.remove_potential(pot)
            except:
                print "Could not remove a potential "+pot.get_potential_type+" from the calculator."


    def define_elements(self):
        """Fills the compound potential with the elemental potentials it is made of.
        
        This method just stores an empty list. Any subclass to be used as a compound potential
        should re-implement this method.
        """
        self.pieces = []


    def describe(self):
        """Prints a description of the potential on-screen.
        """
    
        message = """
potential '{pot}'
    
{n_targ}-body interaction

{descr}
parameters ({n_par}):
""".format(pot=self.potential_type,
            n_targ=str(self.n_targets),
            descr=self.description,
            n_par=str(self.n_params))

        for para, note, val in zip(self.names_of_params,
                                self.descriptions_of_params,
                                self.parameters):
            message += "{pa} : {no} = {va}\n".format(pa=para,no=note,va=val)

        message += "\ncutoff = {cut}\n".format(cut=str(self.cutoff))

        message += "\ntargeted symbols:"
        if not self.symbols is None:
            for ele in self.symbols:
                message += " {0} ".format(str(ele))
            message += "\n"

        if not self.tags is None:
            message += "\ntargeted tags:"
            for tag in self.tags:
                message += " {0} ".format(str(tag))
            message += "\n"

        if not self.indices is None:
            message += "\ntargeted indices:"
            for indy in self.indices:
                message += " {0} ".format(str(indy))
            message += "\n"
            
        print message

