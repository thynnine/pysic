#! /usr/bin/env python
"""A module for handling hybrid QM/MM calculations in pysic."""

from pysic.utility.error import *
from ase import Atoms
from ase.visualize import view # ASEs internal viewer
from pysic.subsystem import *
from pysic.binding import *
import colorsys
import ase.data.colors

#===============================================================================
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

    #---------------------------------------------------------------------------
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

    def add_subsystem(self, subsystem):
        """Used to define subsystems
        
        You can define a subsystem with a oneliner, or then call setters on the
        SubSystemInfo object returned by this function.

        Parameters:
            subsystem: SubSystem object
        """
        if subsystem.name in self.subsystem_info:
            warn("Overriding an existing subsystem", 2)

        self.subsystem_info[subsystem.name] = subsystem

    def add_binding(self, binding):
        """Used to add a binding between two subsystems.

        Parameters:
            binding: Binding
                The Binding object containing the information.
        """
        primary = binding.primary
        secondary = binding.secondary

        # Check that the subsystems exist
        if not self.subsystem_defined(primary):
            return
        if not self.subsystem_defined(secondary):
            return

        self.binding_info[(binding.primary, binding.secondary)] = binding

    def initialize_system(self):
        """ Initializes the subsystems and bindings.
        
        Called once during the lifetime of the calculator. Typically when
        calling set_atoms for the first time, or when calculating any
        quantity for the first time. If the atoms in the simulation need to be
        updated, update_system() is used.
        """
        # Can't do calculation without atoms
        if not self.atoms_set():
            return

        # Initialize subsystems
        for subsystem_info in self.subsystem_info.itervalues():
            self.initialize_subsystem(subsystem_info)

        # Check that the subsystems cover the entire system
        if len(self.get_unsubsystemized_atoms()) is not 0:
            warn("Subsystems do not cover the entire system", 2)
            self.system_initialized = True
            return

        # Initialize bindings
        for binding_info in self.binding_info.itervalues():
            self.initialize_binding(binding_info)

        self.system_initialized = True

    def initialize_binding(self, info):
        """@todo: Docstring for initialize_binding."""

        primary = info.primary
        secondary = info.secondary

        # Check subsystem existence
        if not self.subsystem_defined(primary):
            warn("The given subsystem "+primary+" does not exist.", 2)
            return

        if not self.subsystem_defined(secondary):
            warn("The given subsystem "+secondary+" does not exist.", 2)
            return

        binding = BindingInternal(
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
        calculator = info.calculator
        real_indices = self.generate_subsystem_indices(name)

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
        subsystem = SubSystemInternal(atoms, info, index_map, reverse_index_map)
        self.subsystems[name] = subsystem

    def update_system(self, atoms):
        """Update the subsystem atoms.
        """
        # Replace the old internal atoms
        self.atoms = atoms.copy()
        
        # Update the positions in the subsystems. If link atoms are present,
        # they are also moved. The velocities and momenta are not updated.
        for subsystem in self.subsystems.values():
            for full_index, sub_index in subsystem.index_map.iteritems():
                new_position = atoms[full_index].position
                subsystem.atoms_for_binding[sub_index].position = new_position
                subsystem.atoms_for_subsystem[sub_index].position = new_position

        # Update the link atom positions
        for binding in self.subsystem_bindings.values():
            if binding.info.link_parameters is not None:
                binding.update_hydrogen_link_positions()

    def check_subsystem_existence(self, atom_indices, name):
        """Check that the defined subsystem exists"""
        if not self.atoms_set():
            return False

        n_atoms = len(self.atoms)
        for index in atom_indices:
            if index > n_atoms or index < 0:
                warn("The subsystem \"" + name + "\" contains nonexisting atoms", 2)
                return False
        return True

    def subsystem_defined(self, name):
        """Checks that there is a subsystem with the given name."""
        defined = name in self.subsystem_info
        if not defined:
            warn("The subsystem called \"" + name + "\" has not been defined", 2)
        return defined

    def atoms_set(self):
        """Checks that there is a subsystem with the given name."""
        is_set = self.atoms is not None
        if not is_set:
            warn("No Atoms object given to the calculator", 2)
        return is_set

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
        """Return the indices of the atoms in the subsystem in the full system.

        You can ask the indices even if the subsystems have not been
        initialized, but the indices of different subsystems may overlap in
        this case. If the subsystems have been initialized this function will
        only return indices if they are valid.
        """
        # Check that the atoms have been set
        if not self.atoms_set:
            return

        # Check that the subsystem has been defined
        if not self.subsystem_defined(name):
            return

        indices = self.subsystem_info[name].real_indices
        if indices is None:
            return self.generate_subsystem_indices(name)

    def generate_subsystem_indices(self, name):
        """Generates the indices for the given subsystem.
        """
        # Check that the subsystem exists
        if not self.subsystem_defined(name):
            return

        # Check that the atoms have been set
        if not self.atoms_set():
            return

        info = self.subsystem_info[name]
        indices = info.indices
        tag = info.tag
        special_set = info.special_set

        # Determine how the atoms are specified: indices, tag or special set
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
        
        # Check whether the subsystem contains any atoms:
        if len(real_indices) is 0:
            warn("The specified subsystem " + "\""+name+"\""+ " does not contain any atoms.", 2)
            return

        # Check subsystem existence and overlap
        if self.check_subsystem_existence(real_indices, name):
            if not self.check_if_overlaps(real_indices, name):
                return real_indices

    def check_if_overlaps(self, atom_indices, name):
        """Check that the subsystems don't overlap"""
        if self.subsystems == None:
            return False
        new_indices = []
        new_indices += atom_indices
        for subsystem in self.subsystems.values():
            new_indices += subsystem.index_map.keys()
        if len(new_indices) > len(set(new_indices)):
            warn("The subsystem "+name+" overlaps with another system", 4)
            return True
        return False
    
    def calculate_forces(self):
        """Calculates the forces in the current structure.

        This force includes the forces internal to the subsystems and the forces
        that bind the subsystems together.
        """

        forces = np.zeros((len(self.atoms), 3))

        # Calculate the forces internal to the subsystems. These need to be
        # calculated first, so that the pseudo density and grid are available
        # for Bindings
        internal_forces = []
        for subsystem in self.subsystems.values():
            internal_forces += subsystem.get_forces()
        for pair in internal_forces:
            index = pair[0]
            force = pair[1]
            forces[index,:] += force

        # Calculate the binding forces
        binding_forces = []
        for binding in self.subsystem_bindings.values():
            binding_forces += binding.get_binding_forces()
        for pair in binding_forces:
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
            subsystem_energy = subsystem.get_potential_energy_without_corrections()
            total_potential_energy += subsystem_energy

        # Calculate connection energies after the subsystem energies. This way
        # the pseudo_density is already calculated for the primary system 
        for binding in self.subsystem_bindings.values():
            total_potential_energy += binding.get_binding_energy()

            #Calculate the link interaction correction if needed
            if binding.info.link_parameters is not None:
                if binding.info.link_parameters['interaction_correction']:
                    binding.calculate_link_atom_correction()
                    total_potential_energy += binding.primary_system.link_interaction_correction

        self.potential_energy = total_potential_energy

    def identical_atoms(self, atoms):
        """Compares the given atoms to the stored atoms object. The Atoms are
        identical if the positions, atomic numbers, periodic boundary
        conditions and cell are the same.
        """
        # If one of the atoms is None, return false
        if self.atoms is None or atoms is None:
            return False

        # Check if the atoms are identical
        if self.atoms.get_positions() is not atoms.get_positions():
            return False
        if self.atoms.get_atomic_numbers() is not atoms.get_atomic_numbers():
            return False
        if self.atoms.get_pbc() is not atoms.get_pbc():
            return False
        if self.atoms.get_cell() is not atoms.get_cell():
            return False
        return True

    #---------------------------------------------------------------------------
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
            warn(("No Atoms object given to the calculator. "
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
            warn(("No Atoms object given to the calculator. "
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

    #---------------------------------------------------------------------------
    # Miscellanous utility functions
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
               "|  HYBRIDCALCULATOR ENERGY SUMMARY:\n"
               "|\n"
               "|  Total energy: "+str(self.potential_energy)+"\n"
               "|"
               )
        for name, subsystem in self.subsystems.iteritems():
            if subsystem.potential_energy is None:
                warn("Potential energy not calculated for subsystem \""+name+"\"", 3)
            else:
                print "|  Subsystem " + name +":"
                print "|      Potential energy: "+str(subsystem.potential_energy)
                print "|      Link interaction correction energy: "+str(subsystem.link_interaction_correction)
        print "|"
        for pair, binding in self.subsystem_bindings.iteritems():
            if binding.binding_energy is None:
                warn("Potential energy not calculated for connection between \""+pair[0]+"\" and \""+pair[1], 3)+"\""
            else:
                print "|  Binding energy between \""+pair[0]+"\" and \""+pair[1]+"\": "+str(binding.binding_energy)
        print "|"
        print "================================================================================\n"

    def get_colors(self):
        """Returns a color set for AtomEyeViewer with different colours for the
        different subsystems.

        When the subsystems have been defined with add_subsystem() and the
        atoms have been set, this function will return a list of colors for
        each atom. The different subsystems have different colours for easy
        identification. You can provide this list of colors to an AtomEyeViewer
        object for visualization.
        """
        # Initialize the system if necessary
        if not self.system_initialized:
            self.initialize_system()

        # Create the different colours for the subsystems
        n_subsystems = len(self.subsystems)
        hue_space = np.linspace(0.0, 2.0/3.0, n_subsystems)
        rgb_values = []
        for hue in hue_space:
            rgb = colorsys.hls_to_rgb(hue, 0.7, 0.7)
            rgb_values.append(rgb)

        # Create the color list
        colors = np.zeros((len(self.atoms), 3), dtype=float)
        for i_ss, ss in enumerate(self.subsystems.values()):
            ss_color = np.array(rgb_values[i_ss])
            index_map = ss.index_map
            for system_index in index_map:
                number = self.atoms[system_index].number
                atom_color = np.array(ase.data.colors.cpk_colors[number])
                colors[system_index, :] = 0.1*atom_color+0.9*ss_color

        return colors.tolist()



