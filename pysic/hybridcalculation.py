#! /usr/bin/env python
"""A module for handling hybrid QM/MM calculations in pysic."""

from pysic.utility.error import *
from pysic import embeddingmethods
from ase import Atoms
from ase.visualize import view # ASEs internal viewer

#===============================================================================
class HybridSystem():
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
    def __init__(self, system=None):
        if system is not None:
            self.system = system.copy()
        else:
            self.system = None
        self.subsystem_calculators = {}
        self.subsystem_index_map = {}
        self.subsystems = {}
        self.subsystem_connections = {}
        self.allowed_embedding_modes = ('MEHL')
        self.subsystem_energies = {}

    def set_system(self, system):
        """Set the entire system for hybrid calculations.
        
        Use :meth:`~pysic.hybridcalculation.set_subsystem` for setting up the
        subsystems.
        """
        self.system = system.copy()

    def get_system(self):
        """Return the original structure involved in the calculations."""
        return self.system

    def set_subsystem(self, name, atom_indices=None, special_set=None, calculator=None):
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
                temp_atoms = Atoms()
                for index in atom_indices:
                    atom = self.system[index]
                    temp_atoms.append(atom)
                atoms = temp_atoms.copy()
                atoms.set_pbc(self.system.get_pbc())
                atoms.set_cell(self.system.get_cell())

                # Create a dictionary containing mapping between index in the
                # original system and atom in a subsystem
                index_and_atom = {}
                for index in atom_indices:
                    atom = self.system[index]
                    index_and_atom[index] = atom

                # Create the SubSystem
                subsystem = SubSystem(atoms, calculator, index_and_atom)
                self.subsystems[name] = subsystem

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
        if self.system == None:
            warn("The total system is not set", 5)
            return False
        n_atoms = len(self.system)
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
        n_atoms = len(self.system)
        used_indices = []
        unsubsystemized_atoms = []
        for subsystem in self.subsystems.values():
            for index in subsystem.index_map.keys():
                used_indices.append(index)
        for index in range(n_atoms):
            if index not in used_indices:
                unsubsystemized_atoms.append(index)
        return unsubsystemized_atoms

    def get_subsystem_indices(self, name):
        """Return a list of atomic indices for the subsystem."""
        if not self.check_subsystem_name_exists(name):
            warn("No such subsystem", 5)
            return
        return self.subsystems[name].index_map.keys()

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
        self.subsystems[name].calculator = calculator

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
                    embeddingmethods.MEHL(
                            self.subsystems[primary],
                            self.subsystems[secondary],
                            parameters)

    def get_potential_energy(self):
        """Returns the potential energy of the hybrid system.

        This method will use a spedific method defined by
        :meth:'~pysic.hybridcalculation.set_embedding_mode' to combine and
        calculate the total potential energy of the combined subsystems.
        """
        total_potential_energy = 0

        # The connection energies have to be calculated before the subsystem
        # energies are calculated!
        for connection in self.subsystem_connections.values():
            total_potential_energy += connection.get_connection_energy()

        # The connection initialization may alter the subsystem
        # (e.g. add a hydrogen link atom) so make sure that the connections are
        # initialized before using any calculators
        for name, subsystem in self.subsystems.iteritems():
            warn("Calculating potential energy in subsystem "+name, 5)
            subsystem_energy = subsystem.get_potential_energy()
            total_potential_energy += subsystem_energy
        return total_potential_energy
 
    def view_subsystems(self):
        """@todo: Docstring for view_subsystems.
        :returns: @todo

        """
        for subsystem in self.subsystems.values():
            view(subsystem.modified_atoms)

    def print_potential_energies(self):
        """@todo: Docstring for print_energies.
        :returns: @todo

        """
        for name, subsystem in self.subsystems.iteritems():
            if subsystem.potential_energy is None:
                warn("Potential energy not calculated for subsystem \""+name+"\"", 3)
            else:
                print "Potential energy in subsystem \""+name+"\": "+str(subsystem.potential_energy)
        for pair, connection in self.subsystem_connections.iteritems():
            if connection.connection_energy is None:
                warn("Potential energy not calculated for connection between \""+pair[0]+"\" and \""+pair[1], 3)+"\""
            else:
                print "Connection energy between \""+pair[0]+"\" and \""+pair[1]+"\": "+str(connection.connection_energy)


#===============================================================================
class SubSystem(object):

    """Docstring for SubSystem. """

    def __init__(self, atoms=None, calculator=None, index_map=None):
        """@todo: to be defined1. """

        self.original_atoms = atoms.copy()
        self.modified_atoms = atoms.copy()
        self.calculator = calculator
        self.index_map = index_map
        self.potential_energy = None
        self.connection_ready = False
        self.embedding_correction = 0
        
    def get_potential_energy(self):
        """@todo: Docstring for get_potential_energy.
        :returns: @todo

        """
        if not self.connection_ready:
            warn("""Please setup the embedding method for the subsystem
before calculating it's potential energy""", 2)
            return None
        self.potential_energy = self.modified_atoms.get_potential_energy() + self.embedding_correction
        return self.potential_energy

