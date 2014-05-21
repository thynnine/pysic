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
        self.primary_system.modified_atoms.set_calculator(self.primary_system.calculator)
        self.secondary_system.modified_atoms.set_calculator(self.secondary_system.calculator)
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

        pbc = primary_system.original_atoms.get_pbc()
        if pbc[0] or pbc[1] or pbc[2]:
            self.k_cutoff = parameters.get('k_cutoff')
            if 'k_cutoff' is None: 
                warn("k_cutoff not specified in parameters", 5)
            self.real_cutoff = parameters.get('real_cutoff')
            if 'real_cutoff' is None: 
                warn("real_cutoff not specified in parameters", 5)
            self.sigma = parameters.get('sigma')
            if 'sigma' is None: 
                warn("sigma not specified in parameters", 5)

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
            q1 = self.primary_system.index_map.get(q1_index)
            m1 = self.secondary_system.index_map.get(m1_index)

            if q1 is None:
                warn("Invalid link: "+str(q1_index)+"-"+str(m1_index)+". The first index does not point to an atom in the primary system.", 2)
            if m1 is None:
                warn("Invalid link: "+str(q1_index)+"-"+str(m1_index)+". The second index does not point to an atom in the secondary system.", 2)

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
        # primary system. We do a copy of the calculator to ensure that the
        # calculator for the subsystem doesn't break
        link_calculator = copy.copy(self.primary_system.calculator)
        link_system.set_calculator(link_calculator)
        warn("Calculating potential energy between link atoms", 5)
        extra_energy = link_system.get_potential_energy()
        self.primary_system.embedding_correction = -extra_energy
            
        # Call base class initialize() only after modifying the system
        EmbeddingMethod.initialize(self)

    def calculate_connection_energy(self):
        """Calculates the electrostatic energy involved in binding the two
        subsystems together.

        If the system is periodic, Ewald summation provided by pysic
        is used. Because the Ewald summation can't target specific
        atoms in the system, the binding energy is calculated by substracting
        the subsystem Coulomb energies from the combined system Coulomb energy.

        If the system is finite...
        """
        if not self.initialized:
            warn("The connection is not ready. Call initialize() first.", 2)
            return None

        primary_charges = []
        primary_positions = []
        secondary_charges = []
        secondary_positions = []

        for atom in self.primary_system.original_atoms:
            if not np.allclose(atom.charge, 0):
                primary_charges.append(atom.charge)
                primary_positions.append(np.array(atom.position))
        for atom in self.secondary_system.original_atoms:
            if not np.allclose(atom.charge, 0):
                secondary_charges.append(atom.charge)
                secondary_positions.append(np.array(atom.position))

        ## Asks the charges from the calculators
        #primary_charge_list = self.primary_system.calculator.get_charges()
        #secondary_charge_list = self.secondary_system.calculator.get_charges()
        
        #for index, charge in enumerate(primary_charge_list):
            #if not np.allclose(charge, 0):
                #primary_charges.append(charge)
                #primary_positions.append(np.array(self.primary_system.original_atoms[index].position))
        #for index, charge in enumerate(secondary_charge_list):
            #if not np.allclose(charge, 0):
                #secondary_charges.append(charge)
                #secondary_positions.append(np.array(self.secondary_system.original_atoms[index].position))
        
        # Do calculations only if there is charge on both subsystems
        if (len(primary_charges) is not 0) and (len(secondary_charges) is not 0):
            
            warn("Calculating electrostatic binding energy between subsystems", 5)

            # Periodic system, use Ewald sums
            pbc = self.primary_system.original_atoms.get_pbc()
            if pbc[0] or pbc[1] or pbc[2]:

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
