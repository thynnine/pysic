#! /usr/bin/env python
"""A module for defining different embedding methods in hybrid calculations"""

from pysic.utility.error import *
from abc import ABCMeta, abstractmethod
from ase import Atom, Atoms
from ase.visualize import view
from pysic import Pysic, CoulombSummation
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
    def __init__(
            self,
            primary_system,
            secondary_system,
            primary_calculator,
            secondary_calculator):
        """@todo: Docstring for __init__.

        :arg1: @todo
        :returns: @todo

        """
        __metaclass__ = ABCMeta
        self.ready = False
        self.primary_system = primary_system
        self.secondary_system = secondary_system
        self.primary_calculator = primary_calculator
        self.secondary_calculator = secondary_calculator

    @abstractmethod
    def is_ready(self):
        """This method should return whether the connection is ready.
        
        The connection is ready when all modifications to the connected subsystems
        are made and the connection energy can be calculated.
        """
        return self.ready

    @abstractmethod
    def initialize(self):
        """This method should setup the connection scheme. Any alterations to
        the subsystems are made here.

        The calculators for the subsystems are also set here. Notice that this
        should be done after altering the system, because the calculators make a
        copy of the system.
        """
        self.primary_system.set_calculator(self.primary_calculator)
        self.secondary_system.set_calculator(self.secondary_calculator)
        self.ready = True

    @abstractmethod
    def get_connection_energy(self):
        """Returns the energy related to the connection between the two
        subsystems.
        """
        pass

#===============================================================================
class MEHL(EmbeddingMethod):
    """Describes a mechanical embedding scheme with hydrogen link atoms.
    """
    def __init__(
            self,
            primary_system,
            secondary_system,
            primary_calculator,
            secondary_calculator,
            primary_indices,
            secondary_indices,
            parameters):
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
        
        Attributes:
        """
        # Check that parameters are OK
        if 'links' not in parameters:
            warn("The list of links not specified in parameters", 5)
        if 'CHL' not in parameters:
            warn("CHL not specified in parameters", 5)
        
        EmbeddingMethod.__init__(
                self,
                primary_system,
                secondary_system,
                primary_calculator,
                secondary_calculator)
        self.primary_system_original = primary_system.copy()
        self.secondary_system_original = secondary_system.copy()
        self.primary_system_index_map = primary_indices
        self.secondary_system_index_map = secondary_indices
        self.links = parameters["links"]
        self.chl = parameters['CHL']

    def initialize(self):
        """docstring for initialize"""
        # Setup the hydrogen link atoms inside the primary subsystem
        for link in self.links:

            # Extract the position of the boundary atoms  
            q1_index = link[0]
            m1_index = link[1]
            q1 = self.primary_system_index_map[q1_index]
            m1 = self.secondary_system_index_map[m1_index]
            rq1 = np.array(q1.position)
            rm1 = np.array(m1.position)

            # Calculate position for the hydrogen atom
            distance = rm1-rq1
            rh = rq1+self.chl*distance

            # Create a hydrogen atom at the specified position.
            hydrogen = Atom('H', rh.tolist())
            self.primary_system.append(hydrogen)

        # Call base class initialize() only after modifying the system
        EmbeddingMethod.initialize(self)

    def get_connection_energy(self):
        """Calculates the electrostatic energy involved in binding the two
        subsystems together.
        """
        if not self.ready:
            warn("The connection is not ready. Call initialize() first.", 4)
            return 0

        # Ask the charges from the calculators
        primary_charge_list = self.primary_calculator.get_charges()
        secondary_charge_list = self.secondary_calculator.get_charges()
        
        primary_charges = []
        primary_positions = []
        secondary_charges = []
        secondary_positions = []

        for index, charge in enumerate(primary_charge_list):
            if not np.allclose(charge, 0):
                primary_charges.append(charge)
                primary_positions.append(np.array(self.primary_system[index].position))
        for index, charge in enumerate(secondary_charge_list):
            if not np.allclose(charge, 0):
                secondary_charges.append(charge)
                secondary_positions.append(np.array(self.secondary_system[index].position))
        
        # Do calculations only if there is charge on both subsystems
        if (len(primary_charges) is not 0) and (len(secondary_charges) is not 0):

            # Periodic system, use Ewald sums
            pbc = self.primary_system.get_pbc()
            if pbc[0] or pbc[1] or pbc[2]:

                # Setup the pysic calculator. It is used for calculating the
                # electrostatic energies with Ewald sums. Because the pysic
                # implementation can't target certain atoms in the system, we do the
                # calculation in three pieces.
                calc = Pysic()
                ewald = CoulombSummation()
                ewald.set_parameter_value('epsilon',0.00552635)
                ewald.set_parameter_value('k_cutoff',0.7)
                ewald.set_parameter_value('real_cutoff',10.0) #
                ewald.set_parameter_value('sigma',1.4)
                calc.set_coulomb_summation(ewald)

                # Combined system Coulomb energy
                combined = self.primary_system_original + self.secondary_system_original
                calc.set_atoms(combined)
                combined_energy = calc.get_potential_energy()

                # Primary system coulomb energy
                calc.set_atoms(self.primary_system_original)
                primary_energy = calc.get_potential_energy()

                # Secondary system coulomb energy
                calc.set_atoms(self.secondary_system_original)
                secondary_energy = calc.get_potential_energy()

                # Return the binding energy
                connection_energy = combined_energy - primary_energy - secondary_energy
                warn("The connecting electrostatic energy between the subsystems is: "+str(connection_energy)+". Calculated in a periodic system with Ewald summation.", 5)
                return connection_energy

            # Finite system, traditional Coulomb potential energy
            else:
                # Numpy arrays for vectorized calculations
                charge_p = np.array(primary_charges)
                charge_s = np.array(secondary_charges)
                r_p = np.array(primary_positions)
                r_s = np.array(secondary_positions)

                # Calculate coulomb energy between the charges
                connection_energy = 0
                for idx, charge in enumerate(charge_p):
                    r_temp = np.array([r_p[idx,:],]*len(r_s))
                    r = np.linalg.norm(r_s - r_temp, axis=1)
                    connection_energy += np.sum(1/(4*np.pi*0.00552635)*charge*charge_s/r)
                    
                warn("The connecting electrostatic energy between the subsystems is: "+str(connection_energy)+". Calculated in a finite system with Coulomb potentials.", 5)
                return connection_energy
        else:
            warn("There is no electrostatic interaction between the subsystems", 5)
            return 0
