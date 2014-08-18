#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Defines classes used for creating and storing information about energetical
interactions between subsystems in a QM/MM hybrid calculation created with the
HybridCalculator-class.
"""
from pysic.utility.error import *
from ase import Atom, Atoms
from pysic import Pysic, CoulombSummation, Potential, ProductPotential
from pysic.interactions.comb import CombPotential
import numpy as np
from pysic.utility.timer import Timer
from copy import copy


#===============================================================================
class Interaction(object):
    """Used to store information about a interaction between subsystems.

    The end user can create and manipulate these objects when defining
    interactions between subsystems. The interactions between subsystems are
    added to the calculation by calling
    :meth:`~pysic.hybridcalculator.HybridCalculator.add_interaction` which adds
    a Interaction object to the calculation.

    When the HybridCalculator sees fit, the subsystem bindings are materialized
    by converting the stored Interactions to InteractionInternals.

    Attributes:
        primary: string
            Name of the primary subsystem.
        secondary: string
            Name of the secondary subsystem.
        links: list of tuples
            Contains a list of different link types. Each list item is a tuple
            with two items: the first is a list of link pairs, the second one
            is the CHL parameter for these links.
        electrostatic_parameters: dictionary
            Contains all the parameters used for creating a Coulomb potential.
        coulomb_potential_enabled: bool
            -
        comb_potential_enabled: bool
            -
        link_atom_correction_enabled: bool
            -
        potentials: list of :class:`~pysic.interactions.local.Potential`
            -
    """
    def __init__(self,
                 primary,
                 secondary,
                 ):
        """
        Parameters:
            primary: string
                Name of the primary system.
            secondary: string
                Name of the secondary system.
        """
        self.primary = primary
        self.secondary = secondary
        self.links = []
        self.electrostatic_parameters = None
        self.coulomb_potential_enabled = False
        self.comb_potential_enabled = False
        self.link_atom_correction_enabled = True
        self.potentials = []

    def set_potentials(self, potentials):
        """Used to set a list of additional Pysic potentials between the
        subsystems.

        Parameters:
            potentials: list of or single
            :class:`~pysic.interactions.local.Potential`
        """
        if type(potentials) is list:
            self.potentials = potentials
        else:
            self.potentials = [potentials]

    def add_potential(self, potential):
        """Used to add a Pysic potential between the subsystems.

        Parameters:
            potential: :class:`~pysic.interactions.local.Potential`
        """
        self.potentials.append(potential)

    def enable_coulomb_potential(
            self,
            epsilon=0.00552635,
            real_cutoff=None,
            k_cutoff=None,
            sigma=None):
        """Enables the electrical Coulomb interaction between the subsystems.

        Parameters:
            epsilon: float
                The vacumm permittivity. Has a default value from the
                literature, in Atomic units.
            r_cutoff: float
                Real space cutoff radius. Provide if system has periodic
                boundary conditions.
            k_cutoff: float
                Reciprocal space cutoff radius. Provide if system has periodic
                boundary conditions.
            sigma: float
                Ewald summation split parameter. Provide if system has periodic
                boundary conditions.
        """
        self.coulomb_potential_enabled = True
        parameters = {'epsilon': epsilon}
        if sigma is not None:
            parameters['sigma'] = sigma
        if k_cutoff is not None:
            parameters['k_cutoff'] = k_cutoff
        if real_cutoff is not None:
            parameters['real_cutoff'] = real_cutoff
        self.electrostatic_parameters = parameters

    def enable_comb_potential(self):
        """Enable the COMB potential between the subsystems. Valid for Si-Si
        and Si-O interactions.
        """
        self.comb_potential_enabled = True

    def add_hydrogen_links(self, pairs, CHL):
        """Defines hydrogen links in the system.

        Parameters:
            pairs: tuple or list of tuples
                A set of tuples containing two indices, the first index in the
                tuple should point to the subsystem in which the hydrogen link
                atom is actually added to (=primary system).
            CHL: float
                Indicates the position of the hydrogen atom on the line defined
                by the atom pairs. Values should be between 0-1. 0 means that
                the link atom overlaps with the atom in the primary system, 1
                means that the link atom overlaps with the atom in the
                secondary system. The scale is linear.
        """
        if type(pairs) is tuple:
            pairs = [pairs]
        self.links.append((pairs, CHL))

    def set_link_atom_correction(self, value):
        """Sets whether the link interaction correction energy is calculated or
        not.

        By default the correction is enabled, but you have to provide a
        secondary calculator that can calculate the interaction energies
        between the link atoms themselves and between link atoms and the
        primary system.

        The link atom correction is defined as:
        .. math::

            E^\text{link} &= - E^\text{int}_\text{MM}\text{(PS, HL)} - E^\text{tot}_\text{MM}\text{(HL)}

        """
        self.link_atom_correction_enabled = value


#===============================================================================
class InteractionInternal(object):
    """The internal version of the Interaction-class.

    This class is meant only for internal use, and should not be accessed by
    the end-user.

    Attributes:
        info: :class:'~Pysic.interaction.Interaction'
            Contains all the info about the interaction given by the user.
        full_system: ASE Atoms
            -
        primary_subsystem: :class:`~pysic.subsystem.SubSystem`
            -
        secondary_subsystem: :class:`~pysic.subsystem.SubSystem`
            -
        uncorrected_interaction_energy: float
            The interaction energy without the link atom correction.
        uncorrected_interaction_forces: numpy array
            The interaction forces without the link atom correction.
        link_atom_correction_energy: float
            -
        link_atom_correction_forces: numpy array
            -
        interaction_energy: float
            The total interaction energy = uncorrected_interaction_energy +
            link_atom_correction_energy
        interaction_forces: numpy array
            The total interaction forces = uncorrected_interaction_forces +
            link_atom_correction_forces
        has_interaction_potentials: bool
            True if any potentials are defined.
        calculator: :class:'~pysic.calculator.Pysic'
            The pysic calculator used for non-pbc systems.
        pbc_calculator: :class:'~pysic.calculator.Pysic'
            The pysic calculator used for pbc systems.
        timer: :class:'~pysic.utility.timer.Timer'
            Used for tracking time usage.
        has_pbc: bool
            -
        link_atoms: ASE Atoms
            Contains all the hydrogen link atoms. Needed when calculating link
            atom correction.
        n_primary: int
            Number of atoms in primary subsystem.
        n_secondary: int
            Number of atoms in secondary subsystem.
        n_full: int
            Number of atoms in full system.
        n_links: int
            Number of link atoms.
    """
    def __init__(
            self,
            full_system,
            primary_subsystem,
            secondary_subsystem,
            info):
        """
        Parameters:
            full_system: ASE Atoms
            primary_subsystem: :class:`~pysic.subsystem.SubSystem`
            secondary_subsystem: :class:`~pysic.subsystem.SubSystem`
            info: :class:`~pysic.interaction.Interaction`
        """
        self.info = info
        self.full_system = full_system
        self.primary_subsystem = primary_subsystem
        self.secondary_subsystem = secondary_subsystem

        self.uncorrected_interaction_energy = None
        self.uncorrected_interaction_forces = None
        self.link_atom_correction_energy = None
        self.link_atom_correction_forces = None
        self.interaction_energy = None
        self.interaction_forces = None

        # Determine if any potentials have been set
        self.has_interaction_potentials = False
        if self.info.comb_potential_enabled:
            self.has_interaction_potentials = True
        if self.info.coulomb_potential_enabled:
            self.has_interaction_potentials = True
        if len(self.info.potentials) != 0:
            self.has_interaction_potentials = True

        self.calculator = Pysic()
        self.pbc_calculator = Pysic()

        self.timer = Timer([
            "Interaction",
            "Interaction (PBC)",
            "Forces",
            "Forces (PBC)",
            "Link atom correction energy",
            "Link atom correction forces"])

        # The interaction needs to know if PBC:s are on
        pbc = primary_subsystem.atoms_for_interaction.get_pbc()
        if pbc[0] or pbc[1] or pbc[2]:
            self.has_pbc = True
        else:
            self.has_pbc = False

        # Initialize hydrogen links
        self.link_atoms = None
        self.setup_hydrogen_links(info.links)

        # Store the number of atoms in different systems
        self.n_primary = len(primary_subsystem.atoms_for_interaction)
        self.n_secondary = len(secondary_subsystem.atoms_for_interaction)
        self.n_full = self.n_primary + self.n_secondary
        if self.link_atoms is not None:
            self.n_links = len(self.link_atoms)
        else:
            self.n_links = 0

        # Initialize the COMB potential first (set_potentials(COMB) is used,
        # because it isn' the same as add_potential(COMB))
        if info.comb_potential_enabled:
            self.setup_comb_potential()

        # Initialize the coulomb interaction
        if info.electrostatic_parameters is not None:
            self.setup_coulomb_potential()

        # Add the other potentials
        self.setup_potentials()

        # Can't enable link atom correction on system without link atoms
        if len(info.links) == 0:
            self.info.link_atom_correction_enabled = False

    def setup_hydrogen_links(self, links):
        """Setup the hydrogen link atoms to the primary system.

        Parameters:
            link_parameters: list of tuples
                Contains the link atom parameters from the Interaction-object.
                Each tuple in the list is a link atom specification for a
                covalent bond of different type. The first item in the tuple is
                a list of tuples containing atom index pairs. The second item
                in the tuple is the CHL parameter for these links.
        """
        self.link_atoms = Atoms(pbc=copy(self.full_system.get_pbc()),
                                cell=copy(self.full_system.get_cell()))
        for bond in links:

            pairs = bond[0]
            CHL = bond[1]

            for link in pairs:

                # Extract the position of the boundary atoms
                q1_index = link[0]
                m1_index = link[1]
                iq1 = self.primary_subsystem.index_map.get(q1_index)
                im1 = self.secondary_subsystem.index_map.get(m1_index)

                if iq1 is None:
                    error(("Invalid link: "+str(q1_index)+"-"+str(m1_index)+":\n"
                          "The first index does not point to an atom in the primary system."))
                if im1 is None:
                    error(("Invalid link: "+str(q1_index)+"-"+str(m1_index)+":\n"
                          "The second index does not point to an atom in the secondary system."))

                q1 = self.primary_subsystem.atoms_for_subsystem[iq1]
                rq1 = q1.position

                # Calculate the separation between the host atoms from the full
                # system. We need to do this in the full system because the
                # subsystems may not have the same coordinate systems due to cell
                # size minimization.
                frq1 = self.full_system[q1_index].position
                frm1 = self.full_system[m1_index].position
                distance = CHL*(frm1 - frq1)

                # Calculate position for the hydrogen atom in both the primary
                # subsystem and in the system containing only link atoms. The link
                # atom system is used when calculating link atom corrections.
                r_primary = rq1 + distance
                r_link = frq1 + distance

                # Create a hydrogen atom in the primary subsystem and in the
                # full subsystem
                hydrogen_primary = Atom('H', position=r_primary)
                hydrogen_link = Atom('H', position=r_link)
                self.primary_subsystem.atoms_for_subsystem.append(hydrogen_primary)
                self.link_atoms.append(hydrogen_link)

        # Update the cell size after adding link atoms
        if self.primary_subsystem.cell_size_optimization_enabled:
            self.primary_subsystem.optimize_cell()

    def update_hydrogen_link_positions(self):
        """Used to update the position of the hydrogen link atoms specified in
        this interaction

        It is assumed that the position of the host atoms have already been
        updated.
        """
        print "UPDATE"
        counter = 0
        for bond in self.info.links:

            pairs = bond[0]
            CHL = bond[1]

            for i, link in enumerate(pairs):

                # Extract the position of the boundary atoms
                q1_index = link[0]
                m1_index = link[1]
                iq1 = self.primary_subsystem.index_map.get(q1_index)
                q1 = self.primary_subsystem.atoms_for_subsystem[iq1]
                rq1 = q1.position

                # Calculate the separation between the host atoms from the full
                # system. We need to do this in the full system because the
                # subsystems may not have the same coordinate systems due to
                # cell size minimization.
                frq1 = self.full_system[q1_index].position
                frm1 = self.full_system[m1_index].position
                print frq1
                print frm1
                distance = CHL*(frm1 - frq1)

                # Calculate position for the hydrogen atom in both the primary
                # subsystem and in the system containing only link atoms. The
                # link atom system is used when calculating link atom
                # corrections.
                r_primary = rq1 + distance
                r_link = frq1 + distance

                # Update hydrogen atom positions
                j = i + counter
                self.link_atoms[j].position = r_link
                self.primary_subsystem.atoms_for_subsystem[self.n_primary+j].position = r_primary

            # Add the number of links in this bond type to the counter
            counter += len(pairs)

    def setup_coulomb_potential(self):
        """Setups a Coulomb potential between the subsystems.

        Ewald calculation is automatically used for pbc-systems. Non-pbc
        systems use the ProductPotential to reproduce the Coulomb potential.
        """
        parameters = self.info.electrostatic_parameters

        # Check that that should Ewald summation be used and if so, that all the
        # needed parameters are given
        if self.has_pbc:
            needed = ["k_cutoff", "real_cutoff", "sigma"]
            for param in needed:
                if param not in parameters:
                    error(param + " not specified in Interaction.enable_coulomb_potential(). It is needed in order to calculate electrostatic interaction energy with Ewald summation in systems with periodic boundary conditions.")
                    return

            ewald = CoulombSummation()
            ewald.set_parameter_value('epsilon', parameters['epsilon'])
            ewald.set_parameter_value('k_cutoff', parameters['k_cutoff'])
            ewald.set_parameter_value('real_cutoff', parameters['real_cutoff'])
            ewald.set_parameter_value('sigma', parameters['sigma'])
            self.pbc_calculator.set_coulomb_summation(ewald)

        else:
            # Add coulomb force between all charged particles in secondary
            # system, and all atoms in primary system. It is assumed that the
            # combined system will be made so that the primary system comes
            # before the secondary in indexing.
            coulomb_pairs = []
            for i, ai in enumerate(self.primary_subsystem.atoms_for_interaction):
                for j, aj in enumerate(self.secondary_subsystem.atoms_for_interaction):
                    if not np.allclose(aj.charge, 0):
                        coulomb_pairs.append([i, j+self.n_primary])

            # There are no charges in the secondary system, and that can't
            # change unlike the charges in primary system
            if len(coulomb_pairs) is 0:
                warn("There cannot be electrostatic interaction between the "
                     "subsystems, because the secondary system does not have "
                     "any initial charges", 2)
                return

            # The first potential given to the productpotential defines the
            # targets and cutoff
            kc = 1.0/(4.0*np.pi*parameters['epsilon'])
            max_cutoff = np.linalg.norm(np.array(self.primary_subsystem.atoms_for_interaction.get_cell()))
            coul1 = Potential('power', indices=coulomb_pairs, parameters=[1, 1, 1], cutoff=max_cutoff)
            coul2 = Potential('charge_pair', parameters=[kc, 1, 1])
            coulomb_potential = ProductPotential([coul1, coul2])

            self.calculator.add_potential(coulomb_potential)

        self.has_coulomb_interaction = True

    def setup_comb_potential(self):
        """Setups a COMB-potential between the subsystems.
        """
        COMB = CombPotential(excludes=[])
        COMB.set_calculator(self.pbc_calculator, True)
        # Notice that set_potentials is used here instead of add_potential.
        # This means that enable_comb_potential has to be called first in the
        # constructor.
        self.pbc_calculator.set_potentials(COMB)

    def setup_potentials(self):
        """Setups the additional Pysic potentials given in the
        Interaction-object.
        """
        # If pbcs are not on, the targets of the potential are modified and the
        # calculator for finite systems is used
        if not self.has_pbc:

            primary_atoms = self.primary_subsystem.atoms_for_interaction
            secondary_atoms = self.secondary_subsystem.atoms_for_interaction

            # Make the interaction potentials
            for potential in self.info.potentials:
                symbols = potential.get_symbols()
                for pair in symbols:
                    element1 = pair[0]
                    element2 = pair[1]

                    pairs = []

                    for i_a, a in enumerate(primary_atoms):
                        for i_b, b in enumerate(secondary_atoms):
                            if (a.symbol == element1 and b.symbol == element2) or (a.symbol == element2 and b.symbol == element2):
                                pairs.append([i_a, i_b])
                trimmed_potential = copy(potential)
                trimmed_potential.set_symbols(None)
                trimmed_potential.set_indices(pairs)
                self.calculator.add_potential(trimmed_potential)

        # If pbcs are on, the interactions need to be calculated with pbc
        # calculator
        else:
            for potential in self.info.potentials:
                self.pbc_calculator.add_potential(potential)

    def calculate_link_atom_correction_energy(self):
        """Calculates the link atom interaction energy defined as

        .. math::

            E^\\text{link} = -E^\\text{tot}_\\text{MM}(\\text{HL})-E^\\text{int}_\\text{MM}(\\text{PS, HL})

        """
        self.timer.start("Link atom correction energy")

        primary_atoms = self.primary_subsystem.atoms_for_interaction
        link_atoms = self.link_atoms
        primary_and_link_atoms = primary_atoms + link_atoms
        secondary_calculator = self.secondary_subsystem.calculator

        # The Pysic calculators in one simulation all share one CoreMirror
        # object that contains the data about the Atoms. This object should be
        # automatically updated if the number of atoms changes. This is however
        # not happening for some reason, so we temporarily force the updation here. TODO:
        # Find out why this is the case
        secondary_calculator.force_core_initialization = True

        E1 = secondary_calculator.get_potential_energy(link_atoms)
        E2 = secondary_calculator.get_potential_energy(primary_and_link_atoms)
        E3 = secondary_calculator.get_potential_energy(primary_atoms)

        secondary_calculator.force_core_initialization = False

        link_atom_correction_energy = -E1 - (E2 - E1 - E3)

        self.link_atom_correction_energy = link_atom_correction_energy
        self.timer.stop()

        return copy(self.link_atom_correction_energy)

    def calculate_link_atom_correction_forces(self):
        """Calculates the link atom correction forces defined as
        
        .. math::

            F^\\text{link} = -\\nabla(-E^\\text{tot}_\\text{MM}(\\text{HL})-E^\\text{int}_\\text{MM}(\\text{PS, HL}))

        """
        self.timer.start("Link atom correction forces")

        primary_atoms = self.primary_subsystem.atoms_for_interaction
        link_atoms = self.link_atoms
        primary_and_link_atoms = primary_atoms + link_atoms
        secondary_calculator = self.secondary_subsystem.calculator

        # The Pysic calculators in one simulation all share one CoreMirror
        # object that contains the data about the Atoms. This object should be
        # automatically updated if the number of atoms changes. This is however
        # not happening for some reason, so we temporarily force the updation here. TODO:
        # Find out why this is the case
        secondary_calculator.force_core_initialization = True

        # The force arrays from the individual subsystems need to be extended
        # to the size of the combined system
        forces = np.zeros((len(primary_and_link_atoms), 3))
        primary_postfix = np.zeros((len(link_atoms), 3))
        link_prefix = np.zeros((len(primary_atoms), 3))

        # Calculate the force that binds the link atoms and primary atoms. We
        # don't need to calculate the force between the link atoms, although
        # the associated energy had to be calculated.
        F1 = secondary_calculator.get_forces(link_atoms)
        F2 = secondary_calculator.get_forces(primary_and_link_atoms)
        F3 = secondary_calculator.get_forces(primary_atoms)

        secondary_calculator.force_core_initialization = False

        forces += F2
        forces -= np.concatenate((F3, primary_postfix), axis=0)
        forces -= np.concatenate((link_prefix, F1), axis=0)

        # Ignore the forces acting on link atoms, and negate the force
        # direction. The energy can't be ignored, force can
        forces = -forces[0:self.n_primary, :]

        # Store the forces in a suitable numpy array that can be added to the
        # forces of the whole system
        link_atom_correction_forces = np.zeros((self.n_full, 3))
        for sub_index in range(self.n_primary):
            full_index = self.primary_subsystem.reverse_index_map[sub_index]
            force = forces[sub_index, :]
            link_atom_correction_forces[full_index, :] = force

        self.link_atom_correction_forces = link_atom_correction_forces
        self.timer.stop()

        return copy(link_atom_correction_forces)

    def calculate_uncorrected_interaction_energy(self):
        """Calculates the interaction energy of a non-pbc system without the
        link atom correction.
        """
        self.timer.start("Interaction")
        finite_energy = 0
        primary = self.primary_subsystem.atoms_for_interaction
        secondary = self.secondary_subsystem.atoms_for_interaction
        combined = primary + secondary

        finite_energy = self.calculator.get_potential_energy(combined)
        self.uncorrected_interaction_energy = finite_energy
        self.timer.stop()

        return copy(finite_energy)

    def calculate_uncorrected_interaction_energy_pbc(self):
        """Calculates the interaction energy of a pbc system without the link
        atom correction.
        """
        self.timer.start("Interaction (PBC)")
        pbc_energy = 0
        primary = self.primary_subsystem.atoms_for_interaction
        secondary = self.secondary_subsystem.atoms_for_interaction
        combined = primary + secondary

        pbc_energy += self.pbc_calculator.get_potential_energy(combined)
        pbc_energy -= self.pbc_calculator.get_potential_energy(primary)
        pbc_energy -= self.pbc_calculator.get_potential_energy(secondary)

        self.uncorrected_interaction_energy = pbc_energy
        self.timer.stop()

        return copy(pbc_energy)

    def calculate_uncorrected_interaction_forces(self):
        """Calculates the interaction forces of a non-pbc system without the
        link atom correction.
        """
        self.timer.start("Forces")
        primary = self.primary_subsystem.atoms_for_interaction.copy()
        secondary = self.secondary_subsystem.atoms_for_interaction.copy()
        combined_system = primary + secondary

        ## Forces due to finite coulomb potential and other pysic potentials
        forces = np.array(self.calculator.get_forces(combined_system))
        self.uncorrected_interaction_forces = forces
        self.timer.stop()

        return copy(forces)

    def calculate_uncorrected_interaction_forces_pbc(self):
        """Calculates the interaction forces of a pbc system without the link
        atom correction.
        """
        self.timer.start("Forces (PBC)")
        primary = self.primary_subsystem.atoms_for_interaction
        secondary = self.secondary_subsystem.atoms_for_interaction
        combined = primary + secondary

        primary_postfix = np.zeros((self.n_secondary, 3))
        secondary_prefix = np.zeros((self.n_primary, 3))

        # The force arrays from the individual subsystems need to be extended
        # to the size of the combined system
        forces = self.pbc_calculator.get_forces(combined)
        forces -= np.concatenate((self.pbc_calculator.get_forces(primary), primary_postfix), axis=0)
        forces -= np.concatenate((secondary_prefix, self.pbc_calculator.get_forces(secondary)), axis=0)
        self.uncorrected_interaction_forces = forces
        self.timer.stop()

        return copy(forces)

    def update_subsystem_charges(self):
        """Updates the charges in the subsystems involved in the interaction (if
        charge update is enabled in them).
        """
        if self.info.coulomb_potential_enabled:
            self.primary_subsystem.update_charges()
            self.secondary_subsystem.update_charges()

    def get_interaction_energy(self):
        """Returns the total interaction energy which consists of the
        uncorrected energies and possible the link atom correction.

        Returns:
            float: the interaction energy
        """
        # Try to update the charges
        self.update_subsystem_charges()

        interaction_energy = 0

        # Calculate the interaction energies
        if self.has_interaction_potentials:
            if self.has_pbc:
                interaction_energy += self.calculate_uncorrected_interaction_energy_pbc()
            else:
                interaction_energy += self.calculate_uncorrected_interaction_energy()

        # Calculate the link atom correction energy if needed
        if self.info.link_atom_correction_enabled is True:
            interaction_energy += self.calculate_link_atom_correction_energy()

        self.interaction_energy = interaction_energy

        return copy(self.interaction_energy)

    def get_interaction_forces(self):
        """Return a numpy array of total 3D forces for each atom in the whole
        system.

        The forces consists of the uncorrected forces and possibly the link
        atom correction forces. The row index refers to the atom index in the
        original structure.

        Returns:
            numpy array: the forces for each atom in the full system
        """
        # Try to update the charges
        self.update_subsystem_charges()

        interaction_forces = np.zeros((self.n_full, 3))

        # Calculate the uncorrected interaction forces
        if self.has_interaction_potentials:
            if self.has_pbc:
                interaction_forces += self.calculate_uncorrected_interaction_forces_pbc()
            else:
                interaction_forces += self.calculate_uncorrected_interaction_forces()

        # Calculate the link atom correction forces if needed
        if self.info.link_atom_correction_enabled is True:
            interaction_forces += self.calculate_link_atom_correction_forces()

        self.interaction_forces = interaction_forces
        return copy(self.interaction_forces)
