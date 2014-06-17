#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Defines classes used for creating and storing information about energetical
bindings between subsystems in a HybridCalculator."""

from pysic.utility.error import *
from ase import Atom, Atoms
from ase.visualize import view
import ase.data
from pysic import Pysic, CoulombSummation, Potential, ProductPotential
import numpy as np
import copy

#===============================================================================
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
                warn("'pairs' not defined in link parameters", 5)
            if 'CHL' not in link_parameters: 
                warn("'CHL' not defined in link parameters", 5)
            if 'interaction_correction' not in link_parameters: 
                warn("'interaction_correction' not defined in link parameters", 5)
            if type(link_parameters["pairs"]) is tuple:
                link_parameters["pairs"] = [link_parameters["pairs"]]

        if coulomb_parameters is not None:
            if 'epsilon' not in link_parameters:
                warn("epsilon not specified in parameters, "
                     "using default value of 0.00552635 (in units e^2/(eV Å))", 5)
                coulomb_parameters['epsilon'] = 0.00552635
                self.has_electrostatic_binding = True
        else:
            self.has_electrostatic_binding = False

        self.link_parameters = link_parameters
        self.electrostatic_parameters = coulomb_parameters
        self.potentials = potentials

    def set_hydrogen_links(self, pairs, CHL, interaction_correction=False):
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
                'CHL': CHL,
                'interaction_correction': interaction_correction
                }

    def set_electrostatic_binding(
            self,
            epsilon=0.00552635,
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

#===============================================================================
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
        self.link_atom_indices = []

         #Used when calculating charges on DFT systems with pseudo density
        #self.primary_density_grid = None
        #self.secondary_density_grid = None
        #self.primary_projected_charges = None
        #self.secondary_projected_charges = None

        # The binding needs to know is PBC:s are on
        if pbc[0] or pbc[1] or pbc[2]:
            self.has_pbc = True
        else:
            self.has_pbc = False

        # Initialize binding with BindingInfo
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
            self.link_atom_indices.append(len(self.primary_system.atoms_for_subsystem))
            self.primary_system.atoms_for_subsystem.append(hydrogen)
            self.link_atoms.append(hydrogen)
            
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
            
            # Create the potential calculator if necessary
            if self.potential_calculator is None:
                self.potential_calculator = Pysic()

            self.potential_calculator.add_potential(coulomb_potential)

        self.has_electrostatic_binding = True

    def calculate_link_atom_correction(self):
        """Calculates the error due to the interaction between link atoms.
        
        If there are multiple link atoms (notice also systems with one link
        attom and pbc), they may start to interact with each other. This
        interaction energy is already included in the secondary system. By
        calling this function this energy is added to the correction term in
        self.primary_system.
        """
        primary_atoms = self.primary_system.atoms_for_binding
        primary_calc = self.primary_system.calculator

        links = Atoms(pbc=primary_atoms.get_pbc(), cell=primary_atoms.get_cell())
        for atom in self.link_atoms:
            links.append(atom)
        view(links)
        view(self.primary_system.atoms_for_subsystem)
        correction = primary_calc.get_potential_energy(links)
        self.primary_system.link_interaction_correction = -correction

    def create_density_grid(self, subsystem):
        """Precalculates a grid 3D grid of 3D points for the charge pseudo
        density calculation.
        """
        calc = subsystem.calculator
        atoms = subsystem.atoms_for_binding

        # We need to make one calculation to get the grid points, at least with
        # GPAW
        calc.get_potential_energy(atoms)
        grid_dim = calc.get_number_of_grid_points()

        dx = grid_dim[0]
        dy = grid_dim[1]
        dz = grid_dim[2]
        cell = atoms.get_cell()
        cx = cell[0]
        cy = cell[1]
        cz = cell[2]

        cux = cx/float((grid_dim[0]-1))
        cuy = cy/float((grid_dim[1]-1))
        cuz = cz/float((grid_dim[2]-1))

        # Create a 3D array of 3D points. Each point is a  xyz-coordinate to a
        # position where the density has been calculated.
        dt = np.dtype('float, float, float')
        grid = np.empty((dx, dy, dz), dtype=dt)
        for x in range(grid_dim[0]):
            for y in range(grid_dim[1]):
                for z in range(grid_dim[2]):
                    r = x*cux+y*cuy+z*cuz
                    grid[x, y, z] = r

        subsystem.density_grid = grid

    def update_charges(self, subsystem):
        """If DFT calculators are used, updates the atomic charges according to
        the electron pseudodensity.
        """
        # The pseudo-density is calculated from the atoms that contain the
        # possible link atoms. This takes into account the effect of the bond
        # on the electronic density.
        atoms = subsystem.atoms_for_subsystem
        calc = subsystem.calculator
        density_grid = np.array(calc.get_pseudo_density(atoms))
        grid = subsystem.density_grid

        initial_charges = subsystem.initial_charges
        atomic_numbers = atoms.get_atomic_numbers()
        total_electron_charge = -np.sum(np.array(atoms.get_atomic_numbers()))

        projected_charges = []

        for atom in atoms:
            atom_charge = 0
            z = atom.number
            # Use the Van der Vaals radius of the atom as a cutoff length. This
            # might be more reasonable than a single cutoff radii for all
            # elements.
            R = ase.data.vdw.vdw_radii[z]
            for x in range(grid.shape[0]):
                for y in range(grid.shape[1]):
                    for z in range(grid.shape[2]):
                        d = density_grid[x, y, z]
                        r_i = np.array(list(grid[x, y, z]))
                        r_atom = np.array(atom.position)
                        if np.linalg.norm(r_i - r_atom) <= R:
                            atom_charge += -d
            projected_charges.append(atom_charge)

        # Normalize the projected charges according to the electronic charge in
        # the whole system
        total_charge = np.sum(np.array(projected_charges))
        projected_charges = np.array(projected_charges)
        projected_charges *= total_electron_charge/total_charge

        # Add the nuclear charges and initial charges
        projected_charges += atomic_numbers
        projected_charges += initial_charges

        # Delete the link atom charges
        projected_charges = np.delete(projected_charges, self.link_atom_indices)

        # Set the calculated charges to the atoms
        subsystem.atoms_for_binding.set_initial_charges(projected_charges.tolist())

        #print projected_charges
        #print initial_charges
        
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
        # Update the charges if necessary
        if self.info.has_electrostatic_binding:
            if hasattr(self.primary_system.calculator, "get_pseudo_density"):
                if self.primary_system.density_grid is None:
                    self.create_density_grid(self.primary_system)
                self.update_charges(self.primary_system)

            if hasattr(self.secondary_system.calculator, "get_pseudo_density"):
                if self.secondary_system.density_grid is None:
                    self.create_density_grid(self.secondary_system)
                self.update_charges(self.secondary_system)

        # Form the combined system from the original atoms
        combined = self.primary_system.atoms_for_binding + self.secondary_system.atoms_for_binding

        # Calculate Ewald energy if needed
        ewald_energy = 0

        if self.ewald_calculator is not None:

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
            if hasattr(self.primary_system.calculator, "get_pseudo_density"):
                if self.primary_system.density_grid is None:
                    self.create_density_grid(self.primary_system)
                self.update_charges(self.primary_system)

            if hasattr(self.secondary_system.calculator, "get_pseudo_density"):
                if self.secondary_system.density_grid is None:
                    self.create_density_grid(self.secondary_system)
                self.update_charges(self.secondary_system)

        # Form the combined system from the original atoms
        combined_system = self.primary_system.atoms_for_binding + self.secondary_system.atoms_for_binding
        forces = np.zeros((len(combined_system), 3))

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
