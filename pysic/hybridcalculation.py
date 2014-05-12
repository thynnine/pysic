#! /usr/bin/env python
"""A module for handling hybrid QM/MM calculations in pysic."""

from pysic.utility.error import *
from pysic import embeddingmethods

#===============================================================================
class HybridCalculation():
    """Used to handle hybrid calculations. """

    def __init__(self):
        self.structure = None
        self.subsystem_calculators = {}
        self.subsystem_indices = {}
        self.subsystems = {}
        self.embedding_modes = {}
        self.allowed_embedding_modes = ('MEHL')

    def set_system(self, atoms):
        """Set the entire system for hybrid calculations.
        
        Use :meth:`~pysic.hybridcalculation.set_subsystem` for setting up the
        subsystems.
        """
        self.structure = atoms

    def set_subsystem(self, name, atom_indices=None, special_set=None):
        """Add a subsystem to the structure.
            
        Create a named subsystem consisting of a list of indices in an ase
        Atoms object. For the typical subdivision with two subsystems please use
        "set_primary_system()" and "set_secondary_system()" 

        Parameters:

        name: string
            the name for the subsystem.
        atom_indices: list of ints
            a list of atomic indices
        """
        if atom_indices == None:
            if special_set == "remaining":
                atom_indices = self.get_unsubsystemized_atoms()
            else:
                warn("No subsystem specified", 5)
                return
        if self.check_subsystem_existence(atom_indices):
            if not self.check_if_overlaps(atom_indices):
                if name in self.subsystem_indices:
                    warn("Overriding an existing subsystem", 5)
                self.subsystem_indices[name] = atom_indices
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

    def check_subsystem_name(self, subsystem_name):
        """Checks that there is a subsystem with the given name."""
        return subsystem_name in self.subsystem_indices.keys()

    def get_unsubsystemized_atoms(self):
        """Return a list of indices for the atoms not already in a subsystem."""
        n_atoms = len(self.structure)
        used_indices = []
        unsubsystemized_atoms = []
        for subsystem in self.subsystem_indices.values():
            used_indices += subsystem
        for index in range(n_atoms):
            if index not in used_indices:
                unsubsystemized_atoms.append(index)
        return unsubsystemized_atoms

    def get_subsystem_indices(self, subsystem_name):
        """Return a list of atomic indices for the subsystem."""
        if subsystem_name not in self.subsystem_indices:
            warn("No such subsystem", 5)
            return
        return self.subsystem_indices[subsystem_name]

    def check_if_overlaps(self, atom_indices):
        """Check that the subsystems don't overlap"""
        new_indices = []
        new_indices += atom_indices
        for subsystem in self.subsystem_indices.values():
            new_indices += subsystem
        if len(new_indices) > len(set(new_indices)):
            return True
        return False

    def set_subsystem_calculator(self, subsystem_name, calculator):
        """Set a specific calculator to a subsystem."""
        if subsystem_name not in self.subsystem_indices:
            warn("Subsystem "+subsystem_name+" not defined", 5)
            return
        self.subsystem_calculators[subsystem_name] = calculator

    def set_primary_calculator(self, calculator):
        """Set a calcutor for the primary subsystem."""
        self.set_subsystem_calculator('primary', calculator)

    def set_secondary_calculator(self, calculator):
        """Set a calcutor for the secondary subsystem."""
        self.set_subsystem_calculator('secondary', calculator)

    def materialize_subsystems(self):
        """Create a new atoms object for all the specified subsystems."""
        subsystems = {}
        for name, subsystem in self.subsystem_indices.iteritems():
            atoms = self.structure[subsystem]
            subsystems[name] = atoms
        self.subsystems = subsystems

    def set_embedding(self, mode, primary, secondary, parameters):
        """Sets an embedding mode between two subsystems."""
        if mode not in self.allowed_embedding_modes:
            warn("The given mode is not defined.", 5)
            return
        if not self.check_subsystem_name(primary):
            warn("The given subsystem "+primary+" does not exist.", 5)
            return
        if not self.check_subsystem_name(secondary):
            warn("The given subsystem "+secondary+" does not exist.", 5)
            return
        if mode == "MEHL":
            self.embedding_modes[(primary, secondary)] = embeddingmethods.MechanicalEmbeddingWithHydrogenLinks(self.subsystems[primary], self.subsystems[secondary], parameters)

    def get_potential_energy(self):
        """Returns the potential energy of the hybrid system.
        
        This method will use a spedific method defined by
        :meth:'~pysic.hybridcalculation.set_embedding_mode' to combine and
        calculate the total potential energy of the combined subsystems.
        """
        total_potential_energy = 0
        for subsystem_name in self.subsystem_indices.keys():
            calculator = self.subsystem_calculators[subsystem_name]
            atoms = self.subsystems[subsystem_name]
            atoms.set_calculator(calculator)
            total_potential_energy += atoms.get_potential_energy()
        for connection in self.embedding_modes.values():
            connection.initialize()
            total_potential_energy += connection.get_connection_energy()
        return total_potential_energy
