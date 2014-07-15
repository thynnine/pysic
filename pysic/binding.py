#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Defines classes used for creating and storing information about energetical
bindings between subsystems in a HybridCalculator."""

from pysic.utility.error import *
from ase import Atom, Atoms
from pysic import Pysic, CoulombSummation, Potential, ProductPotential
import numpy as np
from pysic.utility.timer import Timer
import copy


#===============================================================================
class Binding(object):

    """Used to store information about a binding between subsystems.

    The end user can create and manipulate these objects when defining bindings
    between subsystems. The bindings between subsystems are added to the
    calculation by calling
    :meth:`~pysic.hybridcalculator.HybridCalculator.add_binding` which adds a
    BindingInfo object to the calculation and returns a reference to this
    object for further manipulation.

    When the HybridCalculator starts the actual calculations, or when atoms are
    set, the subsystem bindigngs are materialized by converting the stored
    BindingInfos to Bindings.
    """

    def __init__(self,
                 primary,
                 secondary,
                 link_atom_correction=False,
                 potentials=[]
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
        self.link_parameters = None
        self.electrostatic_parameters = None
        self.potentials = potentials
        self.has_coulomb_interaction = False
        self.link_atom_correction = link_atom_correction

    def set_hydrogen_links(self, pairs, CHL):
        """Defines hydrogen link positions in the system.

        Parameters:
            pairs: tuple or list of tuples
                A set of tuples containing two indices, the first index in the
                tuple should point to the subsystem in which the hydrogen link
                atom is actually added to.
            CHL: float
                Indicates the position of the hydrogen atom on the line defined
                by the atom pairs. Values should be between 0-1. 1 means that
                the link atom overlaps with the atom in the primary system, 0
                means that the link atom overlaps with the atom in the
                secondary system. The scale is linear.
            interaction_correction: bool
                Is the link self interaction energy corrected. The correction
                may significantly increase computation time because the primary
                calculator has to be used to calculate the correction.
        """
        if type(pairs) is tuple:
            pairs = [pairs]
        self.link_parameters = {
            'pairs': pairs,
            'CHL': CHL
        }

    def set_link_atom_correction(self, value):
        """Sets whether the link interaction correction energy is calculated or
        not.
        """
        self.link_atom_correction = value

    def set_coulomb_interaction(
            self,
            epsilon=0.00552635,
            real_cutoff=None,
            k_cutoff=None,
            sigma=None):
        """Enables the electrical binding between the subsystems.

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
        parameters = {'epsilon': epsilon}
        if sigma is not None:
            parameters['sigma'] = sigma
        if k_cutoff is not None:
            parameters['k_cutoff'] = k_cutoff
        if real_cutoff is not None:
            parameters['real_cutoff'] = real_cutoff
        self.electrostatic_parameters = parameters
        self.has_coulomb_interaction = True

    def set_potentials(self, potentials):
        """Adds the provided Pysic Potential between the subsystems.
        """
        if type(potentials) is list:
            self.potentials = potentials
        else:
            self.potentials = [potentials]


#===============================================================================
class BindingInternal(object):

    """Materialization of a Binding object"""

    def __init__(
            self,
            full_system,
            primary_system,
            secondary_system,
            info,
            record_time_usage):
        """
        :primary_system: @todo
        :secondary_system: @todo
        :returns: @todo

        """
        self.full_system = full_system
        self.info = info
        self.primary_system = primary_system
        self.secondary_system = secondary_system

        self.uncorrected_binding_energy = None
        self.link_atom_correction_energy = None
        self.binding_energy = None

        self.finite_calculator = Pysic()
        self.pbc_calculator = Pysic()

        self.link_atoms = Atoms(pbc=full_system.get_pbc(), cell=full_system.get_cell())
        self.timer = Timer(record_time_usage,
                           {
                               "Interaction": 0,
                               "Interaction (PBC)": 0,
                               "Forces": 0,
                               "Forces (PBC)": 0,
                               "Link atom correction": 0
                           }
                           )

        # The binding needs to know if PBC:s are on
        pbc = primary_system.atoms_for_binding.get_pbc()
        if pbc[0] or pbc[1] or pbc[2]:
            self.has_pbc = True
        else:
            self.has_pbc = False

        # Initialize binding with Binding
        if info.link_parameters is not None:
            self.set_hydrogen_links(info.link_parameters)
        if info.electrostatic_parameters is not None:
            self.set_electrostatic_binding(info.electrostatic_parameters)

        # Add the other potentials
        self.add_binding_potentials()

    def set_hydrogen_links(self, link_parameters):
        """@todo: Docstring for set_hydrogen_links.

        :indices: @todo
        :CHL: @todo
        :returns: @todo

        """
        pairs = link_parameters['pairs']
        CHL = link_parameters['CHL']
        for link in pairs:

            # Extract the position of the boundary atoms
            q1_index = link[0]
            m1_index = link[1]
            iq1 = self.primary_system.index_map.get(q1_index)
            im1 = self.secondary_system.index_map.get(m1_index)

            if iq1 is None:
                warn(("Invalid link: "+str(q1_index)+"-"+str(m1_index)+":\n"
                      "The first index does not point to an atom in the primary system."
                      ), 2)
            if im1 is None:
                warn(("Invalid link: "+str(q1_index)+"-"+str(m1_index)+":\n"
                      "The second index does not point to an atom in the secondary system."), 2)

            q1 = self.primary_system.atoms_for_subsystem[iq1]
            rq1 = np.array(q1.position)

            # Calculate the separation between the host atoms from the full
            # system. We need to do this in the full system because the
            # subsystems may not have the same coordinate systems due to cell
            # size minimization.
            frq1 = np.array(self.full_system[q1_index].position)
            frm1 = np.array(self.full_system[m1_index].position)
            distance = frm1-frq1

            # Calculate position for the hydrogen atom
            rh = rq1+CHL*distance

            # Create a hydrogen atom at the specified position.
            hydrogen = Atom('H', rh.tolist())
            self.primary_system.link_atom_indices.append(len(self.primary_system.atoms_for_subsystem))
            self.primary_system.atoms_for_subsystem.append(hydrogen)
            self.link_atoms.append(hydrogen)

        # Update the cell size after adding link atoms
        if self.primary_system.info.cell_size_optimization:
            self.primary_system.minimize_cell()

    def update_hydrogen_link_positions(self):
        """Used to update the position of the hydrogen link atoms specified in
        this binding,

        It is assumed that the position of the host atoms have already been
        updated.
        """
        pairs = self.info.link_parameters['pairs']
        CHL = self.info.link_parameters['CHL']

        for i, link in enumerate(pairs):

            # Extract the position of the boundary atoms
            q1_index = link[0]
            m1_index = link[1]
            iq1 = self.primary_system.index_map.get(q1_index)
            q1 = self.primary_system.atoms_for_subsystem[iq1]
            rq1 = np.array(q1.position)

            # Calculate the separation between the host atoms from the full
            # system. We need to do this in the full system because the
            # subsystems may not have the same coordinate systems due to cell
            # size minimization.
            frq1 = np.array(self.full_system[q1_index].position)
            frm1 = np.array(self.full_system[m1_index].position)
            distance = frm1-frq1

            # Calculate position for the hydrogen atom
            rh = rq1+CHL*distance

            # Update hydrogen atom position.
            self.link_atoms[i].position = rh.tolist()

    def set_electrostatic_binding(self, parameters):
        """@todo: Docstring for set_electrostatic_binding.

        :epsilon: @todo
        :returns: @todo

        """
        # Check that that should Ewald summation be used and if so, that all the
        # needed parameters are given
        if self.has_pbc:
            needed = ["k_cutoff", "real_cutoff", "sigma"]
            for param in needed:
                if param not in parameters:
                    warn(param + " not specified. It is needed in order to calculate electrostatic binding energy with Ewald summation in systems with periodic boundary conditions", 2)
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
            n_primary_atoms = len(self.primary_system.atoms_for_binding)
            coulomb_pairs = []
            for i, ai in enumerate(self.primary_system.atoms_for_binding):
                for j, aj in enumerate(self.secondary_system.atoms_for_binding):
                    if not np.allclose(aj.charge, 0):
                        coulomb_pairs.append([i, j+n_primary_atoms])

            # There are no charges in the secondary system, and that can't
            # change unlike the charges in primary system
            if len(coulomb_pairs) is 0:
                warn("There cannot be electrostatic binding between the "
                     "subsystems, because the secondary system does not have "
                     "any initial charges", 2)
                return

            # The first potential given to the productpotential defines the targets and cutoff
            kc = 1.0/(4.0*np.pi*parameters['epsilon'])
            max_cutoff = np.linalg.norm(np.array(self.primary_system.atoms_for_binding.get_cell()))
            coul1 = Potential('power', indices=coulomb_pairs, parameters=[1, 1, 1], cutoff=max_cutoff)
            coul2 = Potential('charge_pair', parameters=[kc, 1, 1])
            coulomb_potential = ProductPotential([coul1, coul2])

            self.finite_calculator.add_potential(coulomb_potential)

        self.has_coulomb_interaction = True

    def calculate_link_atom_correction_energy(self):
        """Calculates the link atom interaction energy defined as

        .. math:: E^\text{link} = -E^\text{tot}_\text{MM}(\text{HL})-E^\text{int}_\text{MM}(\text{PS, HL})
        :label: link_atom_correction
        """
        self.timer.start("Link atom correction")

        primary_and_link_atoms = self.primary_system.atoms_for_subsystem
        primary_atoms = self.primary_system.atoms_for_binding
        link_atoms = self.link_atoms
        secondary_calculator = self.secondary_system.calculator

        E1 = secondary_calculator.get_potential_energy(link_atoms)
        E2 = secondary_calculator.get_potential_energy(primary_and_link_atoms) - secondary_calculator.get_potential_energy(primary_atoms) - secondary_calculator.get_potential_energy(link_atoms)
        self.timer.end()
        return -E1 - E2

    def add_binding_potentials(self):
        """@todo: Docstring for add_binding_potential.
        :potential: @todo
        :returns: @todo
        """
        # If pbcs are not on, the interaction potentials can be modified and
        # calculated with the finite calculator
        if not self.has_pbc:

            primary_atoms = self.primary_system.atoms_for_binding
            secondary_atoms = self.secondary_system.atoms_for_binding

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
                trimmed_potential = copy.copy(potential)
                trimmed_potential.set_symbols(None)
                trimmed_potential.set_indices(pairs)
                self.finite_calculator.add_potential(trimmed_potential)

        # If pbcs are on, the interactions need to be calculated with pbc calculator
        else:
            self.pbc_calculator.set_potentials(self.info.potentials)

    def calculate_binding_energy(self):

        self.timer.start("Interaction")

        finite_energy = 0
        primary = self.primary_system.atoms_for_binding
        secondary = self.secondary_system.atoms_for_binding
        combined = primary + secondary

        finite_energy = self.finite_calculator.get_potential_energy(combined)
        self.timer.end()

        return finite_energy

    def calculate_binding_energy_pbc(self):

        self.timer.start("Interaction (PBC)")

        pbc_energy = 0
        primary = self.primary_system.atoms_for_binding
        secondary = self.secondary_system.atoms_for_binding
        combined = primary + secondary

        pbc_energy += self.pbc_calculator.get_potential_energy(combined)
        pbc_energy -= self.pbc_calculator.get_potential_energy(primary)
        pbc_energy -= self.pbc_calculator.get_potential_energy(secondary)
        self.timer.end()

        return pbc_energy

    def calculate_binding_forces(self):

        self.timer.start("Forces")
        # Form the combined system from the original atoms
        combined_system = self.primary_system.atoms_for_binding + self.secondary_system.atoms_for_binding
        forces = np.zeros((len(combined_system), 3))

        # Forces due to finite coulomb potential and other pysic potentials
        forces += np.array(self.finite_calculator.get_forces(combined_system))
        self.timer.end()

        return forces

    def calculate_binding_forces_pbc(self):

        self.timer.start("Forces (PBC)")

        primary = self.primary_system.atoms_for_binding
        secondary = self.secondary_system.atoms_for_binding
        combined = primary + secondary

        forces = np.zeros((len(combined), 3))
        primary_postfix = np.zeros((len(secondary), 3))
        secondary_prefix = np.zeros((len(primary), 3))

        # The force arrays from the individual subsystems need to be extended
        # to the size of the combined system
        forces += self.pbc_calculator.get_forces(combined)
        forces -= np.concatenate((self.pbc_calculator.get_forces(primary), primary_postfix), axis=0)
        forces -= np.concatenate((secondary_prefix, self.pbc_calculator.get_forces(secondary)), axis=0)
        self.timer.end()

        return forces

    def get_binding_energy(self):
        """@todo: Docstring for get_binding_energy.

        :arg1: @todo
        :returns: @todo

        """
        # Update the charges if necessary
        if self.info.has_coulomb_interaction:

            if hasattr(self.primary_system.calculator, "get_pseudo_density"):
                self.primary_system.update_charges()

            if hasattr(self.secondary_system.calculator, "get_pseudo_density"):
                self.secondary_system.update_charges()

        # Calculate the interaction energies
        if self.has_pbc:
            self.uncorrected_binding_energy = self.calculate_binding_energy_pbc()
        else:
            self.uncorrected_binding_energy = self.calculate_binding_energy()

        # Calculate the link atom correction energy if needed
        if self.info.link_atom_correction is True:
            self.link_atom_correction_energy = self.calculate_link_atom_correction_energy()
            self.binding_energy = self.uncorrected_binding_energy + self.link_atom_correction_energy
        else:
            self.link_atom_correction_energy = None
            self.binding_energy = self.uncorrected_binding_energy

        return self.binding_energy

    def get_binding_forces(self):
        """Return a numpy array of 3D forces. The index refers to the atom
        index in the original structure.
        """
        # Update the charges if necessary
        if self.info.has_coulomb_interaction:

            if hasattr(self.primary_system.calculator, "get_pseudo_density"):
                self.primary_system.update_charges()

            if hasattr(self.secondary_system.calculator, "get_pseudo_density"):
                self.secondary_system.update_charges()

        # Calculate the binding forces
        if self.has_pbc:
            forces = self.calculate_binding_forces_pbc()
        else:
            forces = self.calculate_binding_forces()

        # A covalent QM force between the PS and HL is present. This same
        # force, but in opposite direction should also act on the SS atoms that
        # are part of the covalent bonds.

        # The binding between the atoms that are covalently bonded is entirely
        # modeled with the link atoms. Thus we need to eliminate any MM forces
        # between covalently linked atoms. This is done according to the link
        # atom energy correction by removing the MM force that would be between
        # the PS and HL.

        # Create a list of tuples, containing index in the total system, and
        # the force affecting the atom in that index
        force_map = []
        n_primary_atoms = len(self.primary_system.atoms_for_binding)
        n_secondary_atoms = len(self.secondary_system.atoms_for_binding)
        for sub_index in range(n_primary_atoms):
            full_index = self.primary_system.reverse_index_map[sub_index]
            force = forces[sub_index, :]
            index_force_pair = (full_index, force)
            force_map.append(index_force_pair)
        for sub_index in range(n_secondary_atoms):
            full_index = self.secondary_system.reverse_index_map[sub_index]
            force = forces[sub_index+n_primary_atoms, :]
            index_force_pair = (full_index, force)
            force_map.append(index_force_pair)

        return force_map
