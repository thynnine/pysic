#! /usr/bin/env python
"""A module for handling hybrid QM/MM calculations in pysic."""

from pysic.utility.error import *
from ase import Atoms
from ase.visualize import view # ASEs internal viewer
from pysic.subsystem import *
from pysic.binding import *

#==============================================================================
class HybridCalculator(object):
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

    #--------------------------------------------------------------------------
    # Hybrid methods

    def __init__(self, atoms=None):
        self.subsystem_energies = {}
        self.subsystems = {}
        self.subsystem_bindings = {}
        self.subsystem_info = {}
        self.binding_info = {}
        self.forces = None
        self.potential_energy = None
        self.system_initialized = False
        if atoms is not None:
            self.atoms = atoms.copy()
        else:
            self.atoms = None

    def set_atoms(self, atoms):
        """Set the entire system for hybrid calculations.

        Use :meth:`~pysic.hybridcalculation.set_subsystem` for setting up the
        subsystems.
        """
        # Update the self.atoms and subsystems if necessary
        if not self.identical_atoms(atoms):
            self.update_system(atoms)

        # Initialize the subsystems and bindings if necessary
        if self.system_initialized is False:
            self.initialize_system()

    def get_atoms(self):
        """Return a copy of the full system."""
        return self.atoms.copy()

    def add_subsystem(self, name, indices=None, tag=None, special_set=None, calculator=None):
        """Used to define subsystems
        
        You can define a subsystem with a oneliner, or then call setters on the
        SubSystemInfo object returned by this function.

        Parameters:
            name: string
                The name of the subsystem, used to identify the subsystem.
            indices: tuple or list
                The indices of the atoms to include in the subsystem.
            tag: int
                The tag of the atoms to include in the subsystem.
            special_set: string
                One of the following strings: "remaining" = All the not
                subsystemized atoms.
            calculator: An ASE compatible calculator
                The calculator for the subsystem.

        Returns:
            Returns the SubSystemInfo object for further modification. This
            object later used for initializing the subsystem.
        """
        if name in self.subsystem_info:
            warn("Overriding an existing subsystem", 2)

        subsystem = SubSystemInfo(name, indices, tag, special_set, calculator)
        self.subsystem_info[name] = subsystem
        return subsystem

    def initialize_system(self):
        """ Initializes the subsystems and bindings.
        
        Called once during the lifetime of the calculator. Typically when
        calling set_atoms for the first time, or when calculating any
        quantity for the first time. If the atoms in the simulation need to be
        updated, update_system() is used.
        """
        # Initialize subsystems
        for subsystem_info in self.subsystem_info.itervalues():
            self.initialize_subsystem(subsystem_info)

        # Check that the subsystems cover the entire system
        if len(self.get_unsubsystemized_atoms()) is not 0:
            warn("Subsystems do not cover the entire system", 2)
            return

        # Initialize bindings
        for binding_info in self.binding_info.itervalues():
            self.initialize_binding(binding_info)

        self.system_initialized = True

    def add_binding(
            self,
            primary,
            secondary,
            link_parameters=None,
            coulomb_parameters=None,
            potentials=[]
            ):
        """Used to specify a binding between two subsystems."""
        binding = BindingInfo(primary, secondary, link_parameters, coulomb_parameters, potentials)
        self.binding_info[(primary, secondary)] = binding
        return binding

    def initialize_binding(self, info):
        """@todo: Docstring for initialize_binding."""

        primary = info.primary
        secondary = info.secondary

        # Check subsystem existence
        if not self.check_subsystem_name_exists(primary):
            warn("The given subsystem "+primary+" does not exist.", 2)
            return
        if not self.check_subsystem_name_exists(secondary):
            warn("The given subsystem "+secondary+" does not exist.", 2)
            return
        binding = Binding(
                    self.subsystems[primary],
                    self.subsystems[secondary],
                    info
                    )
        self.subsystem_bindings[(primary, secondary)] = binding

    def initialize_subsystem(self, info):
        """Create the subsystems according to information in
        self.subsystem_info
        """
        name = info.name
        indices = info.indices
        tag = info.tag
        special_set = info.special_set
        calculator = info.calculator

        if calculator is None:
            warn("Calculator is not defined for subsystem "+name, 2)
            return

        # Determine how the atoms are apecified: indices, tag or special set
        if (indices != None) and (tag == None) and (special_set == None):
            real_indices = indices
        elif (indices == None) and (tag != None) and (special_set == None):
            real_indices = []
            for i, t in enumerate(self.atoms.get_tags()):
                if t == tag:
                    real_indices.append(i)
        elif (indices == None) and (tag == None) and (special_set != None):
            if special_set == "remaining":
                real_indices = self.get_unsubsystemized_atoms()
            else:
                warn("Invalid special set", 2)
                return
        else:
            warn("Provide system as indices, tags or a special set", 2)
            return

        # Check subsystem existence and overlap
        if self.check_subsystem_existence(real_indices):
            if not self.check_if_overlaps(real_indices):

                # Create a copy of the subsystem
                temp_atoms = Atoms()
                index_map = {}
                reverse_index_map = {}
                counter = 0
                for index in real_indices:
                    atom = self.atoms[index]
                    temp_atoms.append(atom)
                    index_map[index] = counter
                    reverse_index_map[counter] = index
                    counter += 1
                atoms = temp_atoms.copy()
                atoms.set_pbc(self.atoms.get_pbc())
                atoms.set_cell(self.atoms.get_cell())

                # Create the SubSystem
                subsystem = SubSystem(atoms, calculator, index_map, reverse_index_map)
                self.subsystems[name] = subsystem

            else:
                warn("The subsystem "+name+" overlaps with another system", 4)
                return
        else:
            warn("The subsystem "+name+" contains nonexisting atoms", 2)
            return

    def update_system(self, atoms):
        """Update the subsystem atoms.
        """
        # Replace the old internal atoms
        self.atoms = atoms.copy()
        
        # Update the positions in the subsystems. If link atoms are present,
        # they are also moved. The velocities and momenta are not updated.
        for subsystem in self.subsystems.values():
            for full_index, sub_index in subsystem.index_map.iteritem():
                new_position = atoms[full_index].position
                subsystem.original_atoms[sub_index].set_position(new_position)
                subsystem.modified_atoms[sub_index].set_position(new_position)

        # Update the link atom positions
        for binding in self.subsystem_bindings.values():
            if binding.info.link_parameters is not None:
                binding.update_hydrogen_link_positions()

    def check_subsystem_existence(self, atom_indices):
        """Check that the defined subsystem exists"""
        if self.atoms == None:
            warn("The total system is not set", 5)
            return False
        n_atoms = len(self.atoms)
        for index in atom_indices:
            if index > n_atoms or index < 0:
                warn("Invalid index in the index list", 5)
                return False
        return True

    def check_subsystem_name_exists(self, name):
        """Checks that there is a subsystem with the given name."""
        return name in self.subsystem_info

    def get_unsubsystemized_atoms(self):
        """Return a list of indices for the atoms not already in a subsystem."""
        n_atoms = len(self.atoms)
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
        for subsystem in self.subsystems.values():
            new_indices += subsystem.index_map.keys()
        if len(new_indices) > len(set(new_indices)):
            return True
        return False

    def view_subsystems(self):
        """@todo: Docstring for view_subsystems.
        :returns: @todo

        """
        for subsystem in self.subsystems.values():
            view(subsystem.atoms_for_subsystem)

    def print_energies(self):
        """Print a detailed summary of the different energies in the system.
        This includes the energies in the subsystems, binding energies and
        possible energy corrections.
        """
        print ("\n================================================================================\n"
               "|\n"
               "|\tHYBRIDCALCULATOR ENERGY SUMMARY:\n"
               "|\n"
               "|\tTotal energy: "+str(self.potential_energy)
               )
        for name, subsystem in self.subsystems.iteritems():
            if subsystem.potential_energy is None:
                warn("Potential energy not calculated for subsystem \""+name+"\"", 3)
            else:
                print "|\tPotential energy in subsystem \""+name+"\": "+str(subsystem.potential_energy)
        for pair, binding in self.subsystem_bindings.iteritems():
            if binding.binding_energy is None:
                warn("|\tPotential energy not calculated for connection between \""+pair[0]+"\" and \""+pair[1], 3)+"\""
            else:
                print "|\tBinding energy between \""+pair[0]+"\" and \""+pair[1]+"\": "+str(binding.binding_energy)
        print "|"
        print "================================================================================\n"
    
    def calculate_forces(self):
        """Calculates the forces in the current structure.

        This force includes the forces internal to the subystems and the forces
        that bind the subsystems together.
        """
        # Calculate the binding forces
        forces = np.zeros((len(self.atoms),3))
        binding_forces = []
        for binding in self.subsystem_bindings.values():
            binding_forces += binding.get_binding_forces()
        for pair in binding_forces:
            index = pair[0]
            force = pair[1]
            forces[index,:] += force

        # Calculate the forces internal to the subsystems
        internal_forces = []
        for subsystem in self.subsystems.values():
            internal_forces += subsystem.get_forces()
        for pair in internal_forces:
            index = pair[0]
            force = pair[1]
            forces[index,:] += force

        self.forces = forces

    def calculate_potential_energy(self):
        """Calculates the potential energy in the current structure.
        """
        total_potential_energy = 0

        # The connection initialization may alter the subsystem
        # (e.g. add a hydrogen link atom) so make sure that the connections are
        # initialized before using any calculators
        for name, subsystem in self.subsystems.iteritems():
            warn("Calculating potential energy in subsystem "+name, 5)
            subsystem_energy = subsystem.get_potential_energy()
            total_potential_energy += subsystem_energy

        # Calculate connection energies after the subsystem energies. This way
        # the pseudo_density is already calculated for the primary system 
        for connection in self.subsystem_bindings.values():
            total_potential_energy += connection.get_binding_energy()

        self.potential_energy = total_potential_energy

    def identical_atoms(self, atoms):
        """Compares the given atoms to the stored atoms object.
        """
        # If one of the atoms is None, return false
        if self.atoms is None or atoms is None:
            return False

        # Check if the atoms are identical
        if self.atoms.get_positions() is not atoms.get_positions():
            return True
        if self.atoms.get_atomic_numbers() is not atoms.get_atomic_numbers():
            return True
        if self.atoms.get_pbc() is not atoms.get_pbc():
            return True
        if self.atoms.get_cell() is not atoms.get_pbc():
            return True
        return False

    #--------------------------------------------------------------------------
    # ASE Calculator interface functions

    def calculation_required(self, atoms, quantities=['forces', 'energy', 'stress']):
        """Check if a calculation is required for any of the the given
        properties.

        Check if the quantities in the quantities list have already been
        calculated for the atomic configuration atoms. The quantities can be
        one or more of: 'energy', 'forces', 'stress'.  This method is used to
        check if a quantity is available without further calculations. For this
        reason, calculators should react to unknown/unsupported quantities by
        returning True, indicating that the quantity is not available.

        Two sets of atoms are considered identical if they have the same
        positions, atomic numbers, unit cell and periodic boundary conditions.

        Parameters:
            atoms: `ASE Atoms`_ object
                This structure is compared to the currently stored.
            quantities: list of strings
                list of keywords 'energy', 'forces', 'stress'

        Returns:
            True if the quantities need to be calculated, false otherwise
        """
        # See if internal atoms are present, are atoms given as parameter and
        # whether the given atoms and stored atoms are identical
        if self.atoms is None:
            if atoms is not None:
                return True
            else:
                warn(("No Atoms object given to the calculator. Please provide"
                      "atoms as an argument, or use set_atoms()"), 2)
                return True
        elif atoms is not None:
            if self.identical_atoms(atoms) is False:
                return True

        # Check if the wanted quantities are already calculated
        if type(quantities) is not list:
            quantities = [quantities]
        for quantity in quantities:
            if quantity == 'energy':
                if self.potential_energy is None:
                    return True
            elif quantity == 'forces':
                if self.forces is None:
                    return True
            else:
                return True

        return False

    def get_forces(self, atoms=None):
        """Returns the forces acting on the atoms.

        If the atoms parameter is given, it will be used for updating the
        structure assigned to the calculator prior to calculating the forces.
        Otherwise the structure already associated with the calculator is used.

        The calculator checks if the forces have been calculated already via
        :meth:`~pysic.hybridcalculator.HybridCalculator.calculation_required`.
        If the structure has changed, the forces are calculated using
        :meth:`~pysic.hybridcalculator.HybridCalculator.calculate_forces`

        Parameters:
            atoms: `ASE atoms`_ object
                The structure to calculate the forces on.
        """
        # Can't do calculation without atoms
        if self.atoms is None and atoms is None:
            warn(("No Atoms object given to the calculator."
                  "Please provide atoms as an argument, or use set_atoms()"), 2)

        # See if calculation is required
        if self.calculation_required(atoms, 'forces'):

            # Update the system if necessary
            if atoms is not None:
                if not self.identical_atoms(atoms):
                    self.update_system(atoms)

            # Initialize the subsystems and bindings if necessary
            if self.system_initialized is False:
                self.initialize_system()

            # Calculate the forces and store to self.forces
            self.calculate_forces()

        return np.copy(self.forces)

    def get_potential_energy(self, atoms=None, force_consistent=False):
        """Returns the potential energy of the hybrid system.

        This method will use a spedific method defined by
        :meth:'~pysic.hybridcalculation.set_embedding_mode' to combine and
        calculate the total potential energy of the combined subsystems.
        """
        # Can't do calculation without atoms
        if self.atoms is None and atoms is None:
            warn(("No Atoms object given to the calculator."
                  "Please provide atoms as an argument, or use set_atoms()"), 2)

        # See if calculation is required
        if self.calculation_required(atoms, 'energy'):

            # Update the system if necessary
            if atoms is not None:
                if not self.identical_atoms(atoms):
                    self.update_system(atoms)

            # Initialize the subsystems and bindings if necessary
            if self.system_initialized is False:
                self.initialize_system()

            # Calculate the potential energy and store to self.potential_energy
            self.calculate_potential_energy()

        return np.copy(self.potential_energy)
