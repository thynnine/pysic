#! /usr/bin/env python
"""A module for defining different embedding methods in hybrid calculations"""

from pysic.utility.error import *
from abc import ABCMeta, abstractmethod

#===============================================================================
class EmbeddingMethod(object):
    """An abstract base class for all classes describing interactions between
    subsystems.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def ready(self):
        True

    @abstractmethod
    def initialize(self):
        pass

#===============================================================================
class MechanicalEmbeddingWithHydrogenLinks(EmbeddingMethod):
    """A class derived from EmbeddingMethod that describes an mechanical
    embedding scheme with hydrogen link atoms.
    """
    def __init__(self, primary, secondary, parameters):
        """ 
        Parameters:

        primary: atoms object
            A reference to the primary system in a HybridCalculation.
        secondary: atoms object
            A reference to the secondary system in a HybridCalculation.
        """
        self.primary_system = primary
        self.secondary_system = secondary
        self.links = parameters[0]
        self.CHL = parameters[1]
        self.ready = False
        self.initialize()

    def ready(self):
        """docstring for ready"""
        return self.ready

    def initialize(self):
        """docstring for initialize"""
        self.ready = True

    def get_connection_energy(self):
        """ """
        return 0
