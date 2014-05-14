#! /usr/bin/env python
"""A module for defining different embedding methods in hybrid calculations"""

from pysic.utility.error import *
from abc import ABCMeta, abstractmethod
from ase import Atom
import numpy as np

#===============================================================================
class EmbeddingMethod(object):
    """An abstract base class for all classes describing an interaction between
    two different subsystems. 
    
    Represents a connection between two groups of ASE Atoms. This virtual base
    class only defines the interface that should be used when creating
    connections. It cannot be instantiated directly.

    You can easily create new embedding schemes by subclassing this class, and
    using the MEHL-class as an example.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def is_ready(self):
        """This method should return whether the connection is ready.
        
        The connection is ready when all modifications to the connected subsystems
        are made and the connection energy can be calculated.
        """
        pass

    @abstractmethod
    def initialize(self):
        """This method should setup the connection scheme. Any alterations to
        the subsystems are made here.
        """
        pass

    @abstractmethod
    def get_connection_energy(self):
        """Returns the energy related to the connection between the two
        subsystems.
        """
        pass


#===============================================================================
class MEHL(EmbeddingMethod):
    """A class derived from EmbeddingMethod that describes an mechanical
    embedding scheme with hydrogen link atoms.
    """
    def __init__(self, primary, secondary, parameters, pbc):
        """
        Args:
            primary: dictionary
               A dictionary with indices of the atoms in the original structure
               as keys, and the atoms themselves as values.
            secondary: dictionary
               A dictionary with indices of the atoms in the original structure
               as keys, and the atoms themselves as values.
            parameters: dictionary
                Contains the parameters for the calculation. 'links' is a list
                of lists contianing the cut bonds as pairs of atoms. The first
                index refers to the primary, quantum mechanical system and the
                second index refers to the secondary system 

        """
        # Check that parameters are OK
        if 'links' not in parameters:
            warn("The list of links not specified in parameters", 5)
        if 'CHL' not in parameters:
            warn("CHL not specified in parameters", 5)

        self.primary_system = primary
        self.secondary_system = secondary
        self.links = parameters["links"]
        self.chl = parameters['CHL']
        self.pbc = pbc
        self.ready = False

    def is_ready(self):
        """docstring for ready"""
        return self.ready

    def initialize(self):
        """docstring for initialize"""
        # Setup the hydrogen link atoms inside the primary subsystem
        hydrogen_counter = -1
        for link in self.links:

            # Extract the position of the boundary atoms  
            q1_index = link[0]
            m1_index = link[1]
            q1 = self.primary_system[q1_index]
            m1 = self.secondary_system[m1_index]
            rq1 = np.array(q1.position)
            rm1 = np.array(m1.position)

            # Calculate position for the hydrogen atom
            distance = rm1-rq1
            rh = rq1+self.chl*distance

            # Create a hydrogen atom at the specified position. These link
            # hydrogens have negative indices
            hydrogen = Atom('H', rh.tolist())
            self.primary_system[hydrogen_counter] = hydrogen
            hydrogen_counter -= 1

        self.ready = True

    def get_connection_energy(self):
        """Calculates the electrostatic energy involved in binding the two
        subsystems together.
        """
        if not self.ready:
            warn("The connection is not ready!", 4)
            return 0
        if self.pbc[0] or self.pbc[1] or self.pbc[2]:
            # Periodic system, use Ewald sums
            
            return 0
        else:
            # Finite system, traditional Coulomb potential energy

            return 0
