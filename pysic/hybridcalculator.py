#! /usr/bin/env python
"""A module for handling hybrid QM/MM calculations with Pysic."""

from pysic.utility.error import *
from ase import Atoms
from pysic.subsystem import *
from pysic.interaction import *
from ase.parallel import parprint
from pysic.utility.timer import Timer
import colorsys
import ase.data.colors
import copy


#===============================================================================
class HybridCalculator(object):
    """Used to create and perform hybrid calculations.

    This class is a fully compatible ASE calculator that provides the
    possibility of dividing the Atoms object into subsystems with different ASE
    calculators attached to them.

    You can also define hydrogen link atoms in the interfaces of these
    subsystems or define any Pysic Potentials through which the subsystems
    interact with each other.

    Attributes:
        total_timer: :class:'~pysic.utility.timer.Timer'
            Keeps track of the total time spent in the calculations.
        subsystem_energies: dictionary
            Name to energy.
        subsystems: dictionary
            Name to SubSystemInternal.
        subsystem_interactions: dictionary
            Pair of names to InteractionInternal.
        subsystem_info: dictionary
            Name to SubSystem.
        interaction_info: dictionary
            Pair of names to Interaction.
        forces: numpy array
            -
        stress: None
            -
        potential_energy: float
            -
        system_initialized: bool
            Indicates whether the SubSystemInternals and InteractionInternals
            have been constructed.
        self.atoms: ASE Atoms
            The full system.
    """
    #---------------------------------------------------------------------------
    # Hybrid methods

    def __init__(self, atoms=None):
        """
        Parameters:
            atoms: ASE Atoms
                The full system. You can later on define it with set_atoms, or
                provide it directly in get_potential_energy()/get_forces().
        """
        self.total_timer = Timer(["total"])
        self.total_timer.start("total")
        self.subsystem_energies = {}
        self.subsystems = {}
        self.subsystem_interactions = {}
        self.subsystem_info = {}
        self.interaction_info = {}
        self.forces = None
        self.stress = None
        self.potential_energy = None
        self.system_initialized = False
        if atoms is not None:
            self.atoms = atoms.copy()
        else:
            self.atoms = None

    def set_atoms(self, atoms):
        """Set the full system for hybrid calculations.

        Use :meth:`~pysic.hybridcalculation.add_subsystem` for setting up the
        subsystems.

        Parameters:
            atoms: ASE Atoms
                The full system.
        """
        # Update the self.atoms and subsystems if necessary
        if not self.identical_atoms(atoms):
            self.update_system(atoms)

        # Initialize the subsystems and interactions if necessary
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

    def add_interaction(self, interaction):
        """Used to add a interaction between two subsystems.

        Parameters:
            interaction: Interaction
                The Interaction object containing the information.
        """
        primary = interaction.primary
        secondary = interaction.secondary

        # Check that the subsystems exist
        if not self.subsystem_defined(primary):
            return
        if not self.subsystem_defined(secondary):
            return

        self.interaction_info[(interaction.primary, interaction.secondary)] = interaction

    def initialize_system(self):
        """ Initializes the subsystems and interactions.

        Called once during the lifetime of the calculator. Typically when
        calling set_atoms for the first time, or when calculating any
        quantity for the first time. If the atoms in the simulation need to be
        updated, update_system() is used.
        """
        # Can't do calculation without atoms
        if not self.full_system_set():
            return

        # Initialize subsystems
        for subsystem_info in self.subsystem_info.itervalues():
            self.initialize_subsystem(subsystem_info)

        # Check that the subsystems cover the entire system
        if len(self.get_unsubsystemized_atoms()) is not 0:
            warn("Subsystems do not cover the entire system", 2)
            self.system_initialized = True
            return

        # Initialize interactions
        for interaction_info in self.interaction_info.itervalues():
            self.initialize_interaction(interaction_info)

        self.system_initialized = True

    def initialize_interaction(self, info):
        """Initializes a InteractionInternal from the given Interaction.

        Parameters:
            info: Interaction
        """

        primary = info.primary
        secondary = info.secondary

        # Check subsystem existence
        if not self.subsystem_defined(primary):
            warn("The given subsystem "+primary+" does not exist.", 2)
            return

        if not self.subsystem_defined(secondary):
            warn("The given subsystem "+secondary+" does not exist.", 2)
            return

        interaction = InteractionInternal(
            self.atoms,
            self.subsystems[primary],
            self.subsystems[secondary],
            info)
        self.subsystem_interactions[(primary, secondary)] = interaction

    def initialize_subsystem(self, info):
        """Initializes a SubsystemInternal from the given Subsystem.

        Parameters:
            info: SubSystem
        """
        name = info.name
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
        subsystem = SubSystemInternal(
            atoms,
            info,
            index_map,
            reverse_index_map,
            len(self.atoms))
        self.subsystems[name] = subsystem

    def update_system(self, atoms):
        """Update the subsystem atoms.
        """
        # Replace the old internal atoms
        self.atoms = atoms.copy()

        # Update the positions in the subsystems. If link atoms are present,
        # they are also moved. The velocities and momenta are not updated in
        # the subsystems.
        for subsystem in self.subsystems.values():
            for full_index, sub_index in subsystem.index_map.iteritems():
                new_position = self.atoms[full_index].position
                subsystem.atoms_for_interaction[sub_index].position = copy.copy(new_position)
                subsystem.atoms_for_subsystem[sub_index].position = copy.copy(new_position)

        # Update the link atom positions
        for interaction in self.subsystem_interactions.values():
            interaction.full_system = self.atoms
            if len(interaction.info.links) != 0:
                interaction.update_hydrogen_link_positions()

    def check_subsystem_indices(self, atom_indices, name):
        """Check that the atomic indices of the subsystem are present."""
        if not self.full_system_set():
            return False

        n_atoms = len(self.atoms)
        for index in atom_indices:
            if index > n_atoms or index < 0:
                warn("The subsystem \"" + name + "\" contains nonexisting atoms", 2)
                return False
        return True

    def subsystem_defined(self, name):
        """Checks that there is a subsystem with the given name.
        """
        defined = name in self.subsystem_info
        if not defined:
            warn("The subsystem called \"" + name + "\" has not been defined", 2)
        return defined

    def full_system_set(self):
        """Checks whether the full system has been defined.
        """
        is_set = self.atoms is not None
        if not is_set:
            warn("No Atoms object given to the calculator", 2)
        return is_set

    def get_unsubsystemized_atoms(self):
        """Return a list of indices for the atoms not already in a subsystem.
        """
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
        if not self.full_system_set:
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
        if not self.full_system_set():
            return

        info = self.subsystem_info[name]
        indices = info.indices
        tag = info.tag

        # Determine how the atoms are specified: indices or tag
        if indices == "remaining":
                real_indices = self.get_unsubsystemized_atoms()
        elif (indices is not None) and (tag is None):
            real_indices = indices
        elif (indices is None) and (tag is not None):
            real_indices = []
            for i, t in enumerate(self.atoms.get_tags()):
                if t == tag:
                    real_indices.append(i)

        # Check whether the subsystem contains any atoms:
        if len(real_indices) is 0:
            warn("The specified subsystem " + "\""+name+"\"" + " does not contain any atoms.", 2)
            return

        # Check subsystem indices and overlap
        if self.check_subsystem_indices(real_indices, name):
            if not self.check_subsystem_overlap(real_indices, name):
                return real_indices

    def check_subsystem_overlap(self, atom_indices, name):
        """Check that the subsystem doesn't overlap with another one.
        """
        if self.subsystems is None:
            return False
        new_indices = []
        new_indices += atom_indices
        for subsystem in self.subsystems.values():
            new_indices += subsystem.index_map.keys()
        if len(new_indices) > len(set(new_indices)):
            warn("The subsystem \""+name+"\" overlaps with another system", 2)
            return True
        return False

    def calculate_potential_energy(self):
        """Calculates the potential energy in the current structure.
        """
        total_potential_energy = 0

        # The connection initialization may alter the subsystem
        # (e.g. add a hydrogen link atom) so make sure that the connections are
        # initialized before using any calculators
        for name, subsystem in self.subsystems.iteritems():

            subsystem_energy = subsystem.get_potential_energy()
            total_potential_energy += subsystem_energy

        # Calculate connection energies after the subsystem energies. This way
        # the pseudo_density is already calculated for the primary system
        for interaction in self.subsystem_interactions.values():
            total_potential_energy += interaction.get_interaction_energy()

        self.potential_energy = total_potential_energy

    def calculate_forces(self):
        """Calculates the forces in the current structure.

        This force includes the forces internal to the subsystems and the forces
        that bind the subsystems together.
        """
        forces = np.zeros((len(self.atoms), 3))

        # Calculate the forces internal to the subsystems. These need to be
        # calculated first, so that the pseudo density and grid are available
        # for Interactions
        for subsystem in self.subsystems.values():
            forces += subsystem.get_forces()

        # Calculate the interaction forces
        for interaction in self.subsystem_interactions.values():
            forces += interaction.get_interaction_forces()

        self.forces = forces

    def identical_atoms(self, atoms):
        """Compares the given atoms to the stored atoms object. The Atoms are
        identical if the positions, atomic numbers, periodic boundary
        conditions and cell are the same.
        """
        # If one of the atoms is None, return false
        if self.atoms is None or atoms is None:
            return False

        # Check if the atoms are identical
        if not np.array_equal(self.atoms.get_positions(), atoms.get_positions()):
            return False
        if not np.array_equal(self.atoms.get_atomic_numbers(), atoms.get_atomic_numbers()):
            return False
        if not np.array_equal(self.atoms.get_pbc(), atoms.get_pbc()):
            return False
        if not np.array_equal(self.atoms.get_cell(), atoms.get_cell()):
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
            atoms: ASE Atoms
                This structure is compared to the currently stored.
            quantities: list of strings
                list of keywords 'energy', 'forces', 'stress'

        Returns:
            bool: True if the quantities need to be calculated, false otherwise.
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
            elif quantity == 'stress':
                if self.stress is None:
                    return True
            else:
                return True

        return False

    def get_potential_energy(self, atoms=None, force_consistent=False):
        """Returns the potential energy of the hybrid system.
        """
        # Can't do force_consistent calculations
        if force_consistent:
            warn("Can't do force consistent calculations. Returning the energy extrapolated to zero kelvin.", 2)

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

            # Initialize the subsystems and interactions if necessary
            if self.system_initialized is False:
                self.initialize_system()

            # Calculate the potential energy and store to self.potential_energy
            self.calculate_potential_energy()

        return np.copy(self.potential_energy)

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
            atoms: ASE Atoms object
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

            # Initialize the subsystems and interactions if necessary
            if self.system_initialized is False:
                self.initialize_system()

            # Calculate the forces and store to self.forces
            self.calculate_forces()

        return np.copy(self.forces)

    def get_stress(self, atoms=None, skip_charge_relaxation=False):
        """This function has not been implemented.
        """
        warn("Stress has no been implemented yet.", 2)
        return None

    #---------------------------------------------------------------------------
    # Miscellanous utility functions
    def view_subsystems(self):
        """Views the subsystems with ASE:s built in viewer.
        """
        # Initialize the system if necessary
        if not self.system_initialized:
            self.initialize_system()

        if rank == 0:
            for subsystem in self.subsystems.values():
                view(subsystem.atoms_for_subsystem)

    def get_subsystem_pseudo_density(self, name):
        """Returns the electron pseudo density for the given subsystem.
        """
        if self.subsystem_defined(name):
            return self.subsystems[name].get_pseudo_density()

    def calculate_subsystem_interaction_charges(self, name):
        """Returns the calculated interaction charges for the given subsystem.

        Before the charges can be calulated, one energy calculation has to be
        made so that the electron density is available
        """
        if self.subsystem_defined(name):

            # Initialize the subsystems and interactions if necessary
            if not self.system_initialized:
                self.initialize_system()

            self.subsystems[name].get_potential_energy()
            self.subsystems[name].update_charges()

    def print_interaction_charge_summary(self):
        """Print a summary of the atomic charges that are used in the
        electrostatic interaction between subsystems.
        """
        message = []

        for name, subsystem in self.subsystems.iteritems():
            message.append("Subsystem \"" + name + "\":")

            for i_atom, atom in enumerate(subsystem.atoms_for_interaction):
                symbol = atom.symbol
                charge = atom.charge
                message.append("    Number: " + str(i_atom) + ", Symbol: " + symbol + ", Charge: " + str(charge))

            message.append("")

        str_message = style_message("HYBRIDCALCULATOR SUMMARY OF CHARGES USED IN INTERACTIONS", message)
        parprint(str_message)

    def print_energy_summary(self):
        """Print a detailed summary of the different energies in the system.
        This includes the energies in the subsystems, interaction energies and
        possible energy corrections.
        """
        message = []
        message.append("Total energy: "+str(self.potential_energy))
        message.append("")

        for name, subsystem in self.subsystems.iteritems():

            message.append("Subsystem \"" + name + "\":")

            if subsystem.potential_energy is None:
                ss_energy = "Not calculated"
            else:
                ss_energy = str(subsystem.potential_energy)
            message.append("    Potential energy: " + ss_energy)
            message.append("")

        for pair, interaction in self.subsystem_interactions.iteritems():

            message.append("Interaction between \""+pair[0]+"\" and \""+pair[1] + ":")

            # Total interaction energy
            if interaction.interaction_energy is None:
                b_energy = "Not calculated"
            else:
                b_energy = str(interaction.interaction_energy)
            message.append("    Total interaction energy: "+b_energy)

            # Link atom correction energy
            if interaction.link_atom_correction_energy is None:
                link_atom_correction_energy = "Not calculated"
            else:
                link_atom_correction_energy = str(interaction.link_atom_correction_energy)
            message.append("        Link atom correction energy: "+link_atom_correction_energy)

        str_message = style_message("HYBRIDCALCULATOR ENERGY SUMMARY", message)
        parprint(str_message)

    def print_time_summary(self):
        """Print a detailed summary of the time usage.
        """
        if rank == 0:
            self.total_timer.stop()
            total_time = self.total_timer.get_total_time()
            known_time = 0

            for interaction in self.subsystem_interactions.values():
                known_time += interaction.timer.get_total_time()
            for subsystem in self.subsystems.values():
                known_time += subsystem.timer.get_total_time()

            unknown_time = total_time - known_time
            message = []

            for name, subsystem in self.subsystems.iteritems():

                message.append("Subsystem \"" + name + "\":")
                subsystem_time = subsystem.timer.get_total_time()

                if np.isclose(subsystem_time, 0):
                    message.append("    Time usage: 0 %")
                else:
                    message.append("    Time usage: " + "{0:.1f}".format(subsystem_time/total_time*100.0) + " %")
                    for name, time in subsystem.timer.sections.iteritems():
                        message.append("        " + name + ": " + "{0:.1f}".format(time/total_time*100.0) + " %")

            for pair, interaction in self.subsystem_interactions.iteritems():

                message.append("Interaction between \""+pair[0]+"\" and \""+pair[1] + ":")
                interaction_time = interaction.timer.get_total_time()

                if np.isclose(interaction_time, 0):
                    message.append("    Time usage: 0 %")
                else:
                    message.append("    Time usage: " + "{0:.1f}".format(interaction_time/total_time*100.0) + " %")
                    for name, time in interaction.timer.sections.iteritems():
                        message.append("        " + name + ": " + "{0:.1f}".format(time/total_time*100.0) + " %")

            message.append("Other tasks: " + "{0:.1f}".format(unknown_time/total_time*100.0) + " %")

            str_message = style_message("HYBRIDCALCULATOR TIME SUMMARY", message)
            print(str_message)

    def print_force_summary(self):
        """Print a detailed summary of forces in the system.
        """
        message = []
        message.append("Total forces:")
        if self.forces is None:
            message.append("  Not calculated")
        else:
            for i, force in enumerate(self.forces):
                message.append("  " + str(i) + ": " + str(force))
        message.append("")

        # Forces in subsystems
        for name, subsystem in self.subsystems.iteritems():

            message.append("Subsystem \"" + name + "\":")

            if subsystem.forces is None:
                message.append("  Not calculated")
            else:
                for i, force in enumerate(subsystem.forces):
                    message.append("  " + str(i) + ": " + str(force))
            message.append("")

        # Forces in interactions
        for pair, interaction in self.subsystem_interactions.iteritems():

            message.append("Interaction forces between \""+pair[0]+"\" and \""+pair[1] + ":")

            # Total interaction force
            message.append("  Total:")

            if interaction.interaction_forces is None:
                message.append("    Not calculated")
            else:
                for i, force in enumerate(interaction.interaction_forces):
                    message.append("    " + str(i) + ": " + str(force))
            message.append("")

            # Link atom correction
            message.append("  Link atom correction:")
            if interaction.link_atom_correction_forces is None:
                message.append("      Not calculated")
            else:
                for i, force in enumerate(interaction.link_atom_correction_forces):
                    message.append("    " + str(i) + ": " + str(force))
            message.append("")

        str_message = style_message("HYBRIDCALCULATOR FORCE SUMMARY", message)
        parprint(str_message)

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
            rgb = colorsys.hls_to_rgb(hue, 0.5, 0.7)
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

    def get_subsystem(self, name):
        """Returns a copy of the ASE Atoms object for a certain subsystem.

        The returned subsystem can be used for e.g. visualization or debugging.
        """
        # Check if name defined
        if self.subsystem_defined(name):

            # Initialize the system if necessary
            if not self.system_initialized:
                self.initialize_system()

            return self.subsystems[name].atoms_for_subsystem
