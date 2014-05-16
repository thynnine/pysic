#! /usr/bin/env python
"""A module for defining different embedding methods in hybrid calculations"""

from pysic.utility.error import *
from abc import ABCMeta, abstractmethod
from ase import Atom, Atoms
from ase.visualize import view
from pysic import Pysic, CoulombSummation
import numpy as np
import copy

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
            secondary_system):
        """Docstring for __init__.
        """
        __metaclass__ = ABCMeta
        self.initialized = False
        self.primary_system = primary_system
        self.secondary_system = secondary_system
        self.connection_energy = None

    def get_connection_energy(self):
        """Returns the connection energy.
        """

        self.initialize()
        self.calculate_connection_energy()
        return self.connection_energy

    @abstractmethod
    def initialize(self):
        """This method should setup the connection scheme. Any alterations to
        the subsystems are made here.
        """
        #self.primary_system.modified_atoms.set_calculator(self.primary_system.calculator)
        #self.secondary_system.modified_atoms.set_calculator(self.secondary_system.calculator)
        self.primary_system.calculator.set_atoms(self.primary_system.modified_atoms)
        self.secondary_system.calculator.set_atoms(self.secondary_system.modified_atoms)
        self.initialized = True
        self.primary_system.connection_ready = True
        self.secondary_system.connection_ready = True

    @abstractmethod
    def calculate_connection_energy(self):
        """Calculates the energy related to the connection between the two
        subsystems and stores it as self.connection_energy.
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
        if 'epsilon' not in parameters:
            warn("epsilon not specified in parameters", 5)
        if 'k_cutoff' not in parameters:
            self.k_cutoff = None
        else:
            self.k_cutoff = parameters['k_cutoff']
        if 'real_cutoff' not in parameters:
            self.real_cutoff = None
        else:
            self.real_cutoff = parameters['real_cutoff']
        if 'sigma' not in parameters:
            self.sigma = None
        else:
            self.sigma = parameters['sigma']

        EmbeddingMethod.__init__(
                self,
                primary_system,
                secondary_system)
        self.links = parameters["links"]
        self.chl = parameters['CHL']
        self.epsilon = parameters['epsilon']

    def initialize(self):
        """docstring for initialize"""
            
        # Create a system consisting of the link atom only. The energy
        # between the link atoms is unphysical, so it has to be removed.
        link_system = Atoms()
        link_system.set_pbc(self.primary_system.original_atoms.get_pbc())
        link_system.set_cell(self.primary_system.original_atoms.get_cell())

        # Setup the hydrogen link atoms inside the primary subsystem
        for link in self.links:

            # Extract the position of the boundary atoms  
            q1_index = link[0]
            m1_index = link[1]
            q1 = self.primary_system.index_map[q1_index]
            m1 = self.secondary_system.index_map[m1_index]
            rq1 = np.array(q1.position)
            rm1 = np.array(m1.position)

            # Calculate position for the hydrogen atom
            distance = rm1-rq1
            rh = rq1+self.chl*distance

            # Create a hydrogen atom at the specified position.
            hydrogen = Atom('H', rh.tolist())
            self.primary_system.modified_atoms.append(hydrogen)

            # Add the hydrogen to the link system
            link_system.append(hydrogen)

        # Remove the potential energy in the link system from the energy of the
        # primary system. We do a deep copy of the calculator to ensure that the
        # calculator for the subsystem doesn't break
        link_calculator = copy.deepcopy(self.primary_system.calculator)
        link_system.set_calculator(link_calculator)
        extra_energy = link_system.get_potential_energy()
        self.primary_system.embedding_correction = -extra_energy
            
        # Call base class initialize() only after modifying the system
        EmbeddingMethod.initialize(self)

    def calculate_connection_energy(self):
        """Calculates the electrostatic energy involved in binding the two
        subsystems together.
        """
        if not self.initialized:
            warn("The connection is not ready. Call initialize() first.", 2)
            return None

        # Asks the charges from the calculators
        primary_charge_list = self.primary_system.calculator.get_charges()
        secondary_charge_list = self.secondary_system.calculator.get_charges()
        
        primary_charges = []
        primary_positions = []
        secondary_charges = []
        secondary_positions = []

        for index, charge in enumerate(primary_charge_list):
            if not np.allclose(charge, 0):
                primary_charges.append(charge)
                primary_positions.append(np.array(self.primary_system.original_atoms[index].position))
        for index, charge in enumerate(secondary_charge_list):
            if not np.allclose(charge, 0):
                secondary_charges.append(charge)
                secondary_positions.append(np.array(self.secondary_system.original_atoms[index].position))
        
        # Do calculations only if there is charge on both subsystems
        if (len(primary_charges) is not 0) and (len(secondary_charges) is not 0):

            # Periodic system, use Ewald sums
            pbc = self.primary_system.original_atoms.get_pbc()
            if pbc[0] or pbc[1] or pbc[2]:

                if self.k_cutoff is None:
                    warn("k_cutoff not specified in parameters", 2)
                if self.real_cutoff is None:
                    warn("real_cutoff not specified in parameters", 2)
                if self.sigma is None:
                    warn("sigma not specified in parameters", 2)

                # Setup the pysic calculator. It is used for calculating the
                # electrostatic energies with Ewald sums. Because the pysic
                # implementation can't target certain atoms in the system, we do the
                # calculation in three pieces.
                calc = Pysic()
                ewald = CoulombSummation()
                ewald.set_parameter_value('epsilon', self.epsilon)
                ewald.set_parameter_value('k_cutoff', self.k_cutoff)
                ewald.set_parameter_value('real_cutoff', self.real_cutoff)
                ewald.set_parameter_value('sigma', self.sigma)
                calc.set_coulomb_summation(ewald)

                # Combined system Coulomb energy
                combined = self.primary_system.original_atoms + self.secondary_system.original_atoms
                calc.set_atoms(combined)
                combined_energy = calc.get_potential_energy()

                # Primary system coulomb energy
                calc.set_atoms(self.primary_system.original_atoms)
                primary_energy = calc.get_potential_energy()

                # Secondary system coulomb energy
                calc.set_atoms(self.secondary_system.original_atoms)
                secondary_energy = calc.get_potential_energy()

                # Return the binding energy
                self.connection_energy = combined_energy - primary_energy - secondary_energy

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
                self.connection_energy = connection_energy
        else:
            warn("There is no electrostatic interaction between the subsystems", 5)
            self.connection_energy = 0
