#! /usr/bin/env python
"""A module for handling hybrid QM/MM calculations in pysic."""

from pysic.utility.error import *
from pysic import embeddingmethods
from ase import Atoms
from ase.visualize import view # ASEs internal viewer

#===============================================================================
class HybridCalculation():
    """Used to handle hybrid calculations.

    Attributes:

        structure: Atoms
            a copy of the system structure
        subsystem_calculators: dictionary
            apping between set_subsystem name and calculator
        subsystems: dictionary
            mapping between subsystem name and Atoms
        subsystem_index_map: dictionary
            mapping between subsystem name and a dictionary containing mapping 
            between atom index and atom
        embedding_modes: dictionary
            mapping between pair of subsystems and their EmbeddingMethod
        allowed_embedding_modes: tuple
            the allowed embedding modes as strings
    """
    def __init__(self):
        self.structure = None 
        self.subsystem_calculators = {}
        self.subsystem_index_map = {}
        self.subsystems = {}
        self.subsystem_connections = {}
        self.allowed_embedding_modes = ('MEHL')

    def set_system(self, atoms):
        """Set the entire system for hybrid calculations.
        
        Use :meth:`~pysic.hybridcalculation.set_subsystem` for setting up the
        subsystems.
        """
        self.structure = atoms.copy()

    def get_system(self):
        """Return the structure involved in the calculations. Notice that this
        can be very different from the initial structure as it contains also the
        possible link atoms created by the embedding schemes."""
        return self.structure

    def set_subsystem(self, name, atom_indices=None, special_set=None):
        """Add a subsystem to the structure.
            
        Create a named subsystem consisting of a list of indices in an ase
        Atoms object. For the typical subdivision with two subsystems please use
        :meth:'set_primary_system()' and :meth:'set_secondary_system()' 

        Args:
            name: string
                the name for the subsystem.
            atom_indices: list of ints
                a list of atomic indices
        """
        if atom_indices == None:
            if special_set == "remaining":
                atom_indices = self.get_unsubsystemized_atoms()
            else:
                warn("Invalid special set", 5)
                return

        if self.check_subsystem_existence(atom_indices):
            if not self.check_if_overlaps(atom_indices):
                if self.check_subsystem_name_exists(name):
                    warn("Overriding an existing subsystem", 5)

                # Create a copy of the subsystem
                temp_system = Atoms()
                for index in atom_indices:
                    atom = self.structure[index]
                    temp_system.append(atom)
                subsystem = temp_system.copy()
                subsystem.set_pbc(self.structure.get_pbc())
                subsystem.set_cell(self.structure.get_cell())
                self.subsystems[name] = subsystem
    
                # Create a dictionary containing mapping between index in the
                # original system and atom in a subsystem
                index_and_atom = {}
                for index in atom_indices:
                    atom = self.structure[index]
                    index_and_atom[index] = atom
                self.subsystem_index_map[name] = index_and_atom

            else:
                warn("The subsystem overlaps with another system", 5)
                return
        else:
            warn("The subsystem does not exist in the system", 5)
            return

    def set_primary_system(self, atom_indices=None, special_set=None):
        """Set a primary subsystem."""
        self.set_subsystem('primary', atom_indices)

    def set_secondary_system(self, atom_indices=None, special_set=None):
        """Set a secondary subsystem."""
        self.set_subsystem('secondary', atom_indices, special_set)

    def check_subsystem_existence(self, atom_indices):
        """Check that the defined subsystem exists"""
        if self.structure == None:
            warn("The total system is not set", 5)
            return False
        n_atoms = len(self.structure)
        for index in atom_indices:
            if index > n_atoms or index < 0:
                warn("Invalid index in the index list", 5)
                return False
        return True

    def check_subsystem_name_exists(self, subsystem_name):
        """Checks that there is a subsystem with the given name."""
        if self.subsystems == None:
            return False
        return subsystem_name in self.subsystems.keys()

    def get_unsubsystemized_atoms(self):
        """Return a list of indices for the atoms not already in a subsystem."""
        n_atoms = len(self.structure)
        used_indices = []
        unsubsystemized_atoms = []
        for subsystem in self.subsystem_index_map.values():
            for index in subsystem.keys():
             used_indices.append(index)
        for index in range(n_atoms):
            if index not in used_indices:
                unsubsystemized_atoms.append(index)
        return unsubsystemized_atoms

    def get_subsystem_indices(self, subsystem_name):
        """Return a list of atomic indices for the subsystem."""
        if not self.check_subsystem_name_exists(subsystem_name):
            warn("No such subsystem", 5)
            return
        return self.subsystem_index_map[subsystem_name].keys()

    def check_if_overlaps(self, atom_indices):
        """Check that the subsystems don't overlap"""
        if self.subsystems == None:
            return False
        new_indices = []
        new_indices += atom_indices
        for subsystem in self.subsystem_index_map.values():
            new_indices += subsystem.keys()
        if len(new_indices) > len(set(new_indices)):
            return True
        return False

    def set_subsystem_calculator(self, name, calculator):
        """Set a specific calculator to a subsystem."""
        if not self.check_subsystem_name_exists(name):
            warn("Subsystem "+name+" not defined", 5)
            return
        self.subsystem_calculators[name] = calculator

    def set_primary_calculator(self, calculator):
        """Set a calcutor for the primary subsystem."""
        self.set_subsystem_calculator('primary', calculator)

    def set_secondary_calculator(self, calculator):
        """Set a calcutor for the secondary subsystem."""
        self.set_subsystem_calculator('secondary', calculator)

    def set_embedding(self, mode, primary, secondary, parameters):
        """Sets an embedding mode between two subsystems."""
        if mode not in self.allowed_embedding_modes:
            warn("The given mode is not defined.", 5)
            return
        if not self.check_subsystem_name_exists(primary):
            warn("The given subsystem "+primary+" does not exist.", 5)
            return
        if not self.check_subsystem_name_exists(secondary):
            warn("The given subsystem "+secondary+" does not exist.", 5)
            return
        if mode == "MEHL":
            self.subsystem_connections[(primary, secondary)] = \
                    embeddingmethods.MEHL(self.subsystems[primary],
                            self.subsystems[secondary],
                            self.subsystem_index_map[primary],
                            self.subsystem_index_map[secondary],
                            parameters)

    def get_potential_energy(self):
        """Returns the potential energy of the hybrid system.
        
        This method will use a spedific method defined by
        :meth:'~pysic.hybridcalculation.set_embedding_mode' to combine and
        calculate the total potential energy of the combined subsystems.
        """
        total_potential_energy = 0

        # The connections have to be initialized before energies are calculated!
        for connection in self.subsystem_connections.values():
            connection.initialize()
            total_potential_energy += connection.get_connection_energy()

        # The connection initialization may alter the subsystem
        # (e.g. add a hydrogen link atom) so make sure that the connections are
        # initialized before using any calculators
        for name in self.subsystems.keys():
            calculator = self.subsystem_calculators[name]
            subsystem = self.subsystems[name] 

            subsystem.set_calculator(calculator)
            total_potential_energy += subsystem.get_potential_energy()

        return total_potential_energy
    
    def view_subsystems(self):
        """@todo: Docstring for view_subsystems.
        :returns: @todo

        """
        for subsystem in self.subsystems.values():
            view(subsystem)
