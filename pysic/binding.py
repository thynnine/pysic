#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Defines classes used for creating and storing information about energetical
bindings between subsystems in a HybridCalculator."""

from pysic.utility.error import *
from ase import Atom, Atoms
from ase.visualize import view
from pysic import Pysic, CoulombSummation, Potential, ProductPotential
import numpy as np
import copy

#==============================================================================
class BindingInfo(object):

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

    def __init__(
            self,
            primary,
            secondary,
            link_parameters=None,
            coulomb_parameters=None,
            potentials=[]
            ):
        """@todo: to be defined1. """

        self.primary = primary
        self.secondary = secondary

        # Check that the given parameter dictionaries are ok
        if link_parameters is not None:
            if 'pairs' not in link_parameters:
                warn("No pairs defined in link parameters", 5)
            if 'CHL' not in link_parameters: 
                warn("CHL not defined in link parameters", 5)
            if type(link_parameters["pairs"]) is tuple:
                link_parameters["pairs"] = [link_parameters["pairs"]]
            
        if coulomb_parameters is not None:
            if 'epsilon' not in link_parameters:
                warn("epsilon not specified in parameters, "
                     "using default value of 0.0052635 (in units e^2/(eV Å))", 5)
                coulomb_parameters['epsilon'] = 0.0052635
                self.has_electrostatic_binding = True
        else:
            self.has_electrostatic_binding = False

        self.link_parameters = link_parameters
        self.electrostatic_parameters = coulomb_parameters
        self.potentials = potentials

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
        """
        if type(pairs) is tuple:
            pairs = [pairs]
        self.link_parameters = {'pairs': pairs, 'CHL': CHL}

    def set_electrostatic_binding(
            self,
            epsilon=0.0052635,
            sigma=None,
            k_cutoff=None,
            real_cutoff=None
            ):
        """@todo: Docstring for set_electrostatic_bingding
        """
        parameters = {'epsilon': epsilon}
        if sigma is not None:
            parameters['sigma'] = sigma
        if k_cutoff is not None:
            parameters['k_cutoff'] = k_cutoff
        if real_cutoff is not None:
            parameters['real_cutoff'] = real_cutoff
        self.electrostatic_parameters = parameters
        self.has_electrostatic_binding = True

    def add_binding_potential(self, potential):
        """@todo: Docstring for set_potential.
        :returns: @todo

        """
        if type(potential) is list:
            self.potentials += potential
        else:
            self.potentials.append(potential)

#==============================================================================
class Binding(object):

    """Materialization of a BindingInfo object"""

    def __init__(self, primary_system, secondary_system, info):
        """
        :primary_system: @todo
        :secondary_system: @todo
        :returns: @todo

        """
        self.info = info
        self.primary_system = primary_system
        self.secondary_system = secondary_system
        self.binding_energy = None
        self.ewald_calculator = None
        self.potential_calculator = None
        pbc = primary_system.atoms_for_binding.get_pbc()
        self.link_atoms = []
        if pbc[0] or pbc[1] or pbc[2]:
            self.has_pbc = True
        else:
            self.has_pbc = False

        # Initialize with info
        if info.link_parameters is not None:
            self.set_hydrogen_links(info.link_parameters)
        if info.electrostatic_parameters is not None:
            self.set_electrostatic_binding(info.electrostatic_parameters)
        for potential in info.potentials:
            self.add_binding_potential(potential)

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
            q1 = self.primary_system.atoms_for_binding[self.primary_system.index_map[q1_index]]
            m1 = self.secondary_system.atoms_for_binding[self.secondary_system.index_map[m1_index]]

            if q1 is None:
                warn(("Invalid link: "+str(q1_index)+"-"+str(m1_index)+":\n"
                      "The first index does not point to an atom in the primary system."
                      ), 2)
            if m1 is None:
                warn(("Invalid link: "+str(q1_index)+"-"+str(m1_index)+":\n"
                      "The second index does not point to an atom in the secondary system."), 2)

            rq1 = np.array(q1.position)
            rm1 = np.array(m1.position)

            # Calculate position for the hydrogen atom
            distance = rm1-rq1
            rh = rq1+CHL*distance

            # Create a hydrogen atom at the specified position.
            hydrogen = Atom('H', rh.tolist())
            self.primary_system.atoms_for_subsystem.append(hydrogen)
            self.link_atoms.append(hydrogen)
            #self.link_atom_indices.append(len(self.primary_system.atoms_for_subsystem))
            
    def update_hydrogen_link_positions(self):
        """Used to update the position of the hydrogen link atoms specified in
        this binding,

        It is assumed that the position of the host atoms has been already
        updated.
        """
        pairs = self.info.link_parameters['pairs']
        CHL = self.info.link_parameters['CHL']

        for i, link in enumerate(pairs):

            # Extract the position of the boundary atoms  
            q1_index = link[0]
            m1_index = link[1]
            q1 = self.primary_system.atoms_for_binding[self.primary_system.index_map[q1_index]]
            m1 = self.secondary_system.atoms_for_binding[self.secondary_system.index_map[m1_index]]

            rq1 = np.array(q1.position)
            rm1 = np.array(m1.position)

            # Calculate position for the hydrogen atom
            distance = rm1-rq1
            rh = rq1+CHL*distance

            # Update hydrogen atom position.
            self.link_atoms[i].set_position(rh.tolist())

    def set_electrostatic_binding(self, parameters):
        """@todo: Docstring for set_electrostatic_binding.

        :epsilon: @todo
        :returns: @todo

        """
        if self.has_pbc:
            if 'k_cutoff' not in parameters: 
                warn("k_cutoff not specified in parameters", 2)
            if 'real_cutoff' not in parameters: 
                warn("real_cutoff not specified in parameters", 2)
            if 'sigma' not in parameters: 
                warn("sigma not specified in parameters", 2)

            self.ewald_calculator = Pysic()
            ewald = CoulombSummation()
            ewald.set_parameter_value('epsilon', parameters['epsilon'])
            ewald.set_parameter_value('k_cutoff', parameters['k_cutoff'])
            ewald.set_parameter_value('real_cutoff', parameters['real_cutoff'])
            ewald.set_parameter_value('sigma', parameters['sigma'])
            self.ewald_calculator.set_coulomb_summation(ewald)

        else:
            # Add coulomb force between all charged particles in secondary
            # system, and all atoms in primary system. It is assumed that the
            # combined system will be made so that the primary system comes
            # before the secondary in indexing.
            n_primary_atoms = len(self.primary_system.atoms_for_binding)
            coulomb_pairs = []
            for i, ai in enumerate(self.primary_system.atoms_for_binding):
                for j, aj in enumerate(self.secondary_system.atoms_for_binding):
                    if not np.allclose(ai.charge, 0):
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
            
            # Create the potential calculator if necessary
            if self.potential_calculator is None:
                self.potential_calculator = Pysic()

            self.potential_calculator.add_potential(coulomb_potential)

        self.has_electrostatic_binding = True

    def update_charges(self):
        """If DFT calculators are used, updates the atom charges according to
        the pseudodensity.
        """
        pass
        
    def add_binding_potential(self, potential):
        """@todo: Docstring for add_binding_potential.
        :potential: @todo
        :returns: @todo
        """
        # Create the potential calculator if necessary
        if self.potential_calculator is None:
            self.potential_calculator = Pysic()

        potential_pairs = []

    def get_binding_energy(self):
        """@todo: Docstring for get_binding_energy.

        :arg1: @todo
        :returns: @todo

        """
        # Form the combined system from the original atoms
        combined = self.primary_system.atoms_for_binding + self.secondary_system.atoms_for_binding

        # Calculate Ewald energy if needed
        ewald_energy = 0

        if self.ewald_calculator is not None:

            # Update charges
            self.update_charges()

            # Combined system Coulomb energy
            ewald_energy += self.ewald_calculator.get_potential_energy(combined)

            # Primary system coulomb energy
            ewald_energy -= self.ewald_calculator.get_potential_energy(self.primary_system.atoms_for_binding)

            # Secondary system coulomb energy
            ewald_energy -= self.ewald_calculator.get_potential_energy(self.secondary_system.atoms_for_binding)

        # Other potential energies (including Coulomb in system without PBC:s)
        potential_energy = 0
        if self.potential_calculator is not None:
            potential_energy = self.potential_calculator.get_potential_energy(combined)

        self.binding_energy = ewald_energy + potential_energy;
        return self.binding_energy

    def get_binding_forces(self):
        """Return a numpy array of 3D forces. The index refers to the atom
        index in the original structure.
        """
        # Update the charges if necessary
        if self.info.has_electrostatic_binding:
            self.update_charges()

        combined_system = self.primary_system.atoms_for_binding + self.secondary_system.atoms_for_binding
        forces = np.zeros(len(combined_system, 3))

        # Forces due to finite coulomb potential and other pysic potentials
        if self.potential_calculator is not None:
            forces += np.array(self.potential_calculator.get_forces(combined_system))

        # Calculate ewald force if necessary
        if self.ewald_calculator is not None:

            # Combined system Coulomb energy
            forces += self.ewald_calculator.get_forces(combined_system)

            # Primary system coulomb energy
            forces -= self.ewald_calculator.get_forces(self.primary_system.atoms_for_binding)

            # Secondary system coulomb energy
            forces -= self.ewald_calculator.get_forces(self.secondary_system.atoms_for_binding)

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
