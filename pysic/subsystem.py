#! /usr/bin/env python
"""Defines classes used for creating and storing information about subsystems
in a HybridCalculator."""

import numpy as np
from pysic.utility.error import *
import ase.data

#===============================================================================
class SubSystem(object):

    """Used to store information about a subsystem.
    
    The end user can create and manipulate these objects when defining
    subsystems. The subsystems are added to the calculation by calling
    :meth:`~pysic.hybridcalculator.HybridCalculator.add_subsystem` which adds a
    SubSystemInfo object to the calculation and returns a reference to this
    object for further manipulation.
    
    When the HybridCalculator starts the actual calculations, or when atoms are
    set, the subsystems are materialized by converting the stored
    SubSystemInfos to SubSystems. 
    """

    def __init__(self,
            name,
            indices=None,
            tag=None,
            special_set=None,
            calculator=None
            ):
        """@todo: to be defined1. """
        self.name = name
        self.calculator = calculator
        self.cell_size_optimization = False
        self.dft_padding = None
        self.set_atoms(indices, tag, special_set)
        self.is_valid()

    def set_calculator(self, calculator):
        """Set the calculator for the subsystem.

        Parameters:
            calculator: ASE Calculator
        """
        self.calculator = calculator

    def set_atoms(self, indices=None, tag=None, special_set=None):
        """Set the atoms that belong to this subsystem. Give only one of the
        specifiers: indices, tag or special_set.

        Parameters:
            indices: list of ints
                A list of atom indices in the full system.
            tag: int
                The tag that is used to identify the subsystem.
            special_set: string
                A special string indicator. One of the following:
                    "remaining": All the the atoms that have not yet been
                    linked to a subsystem.
        """
        if indices is not None and type(indices) is int:
                self.indices = (indices,)
        else:
            self.indices = indices
        self.tag = tag
        self.special_set = special_set

    def is_valid(self):
        """Checks that the given atom specifiers are correctly given. Does not
        yet check that they exist or don't overlap with othe subssytems
        """
        indices = self.indices
        tag = self.tag
        special_set = self.special_set

        # Determine how the atoms are specified: indices, tag or special set
        if (indices != None) and (tag == None) and (special_set == None):
            return True
        elif (indices == None) and (tag != None) and (special_set == None):
            return True
        elif (indices == None) and (tag == None) and (special_set != None):
            if special_set == "remaining":
                return True
            else:
                warn("Invalid special set", 2)
                return False
        else:
            warn("Provide system as indices, tags or a special set", 2)
            return False
        
    def enable_cell_size_optimization(self, padding):
        """Enable cell size optimization.
        
        The optimization is off by default. If there are no periodic boundary
        conditions, the cell size of a DFT subsystem can be minimized. This
        speeds up the calculations in systems where the DFT specific subsystem
        covers only a small portion of the entire system.

        The padding indicates the minimum distance between the cell boundaries
        and the subsystem atoms. If this value is too low, the DFT-calculator
        might not work properly!
        """
        self.cell_size_optimization = True
        self.dft_padding = padding

#===============================================================================
class SubSystemInternal(object):

    """A materialization of a SubSystem object.

    This class is materialised from a SubSystem, and should not be
    accessible to the end user.

    Attributes:
        self.info
        self.atoms_for_binding
        self.atoms_for_subsystem
        ...
    """
    def __init__(self, atoms, info, index_map, reverse_index_map):
        """
        Parameters:
            atoms: ASE Atoms
                The subsystem atoms.
            info: SubSystem object
                Contains all the information about the subsystem
            index_map: dictionary of int to int
                The keys are the atom indices in the full system, values are indices in the subssystem.
            reverse_index_map: dicitonary of int to int
                The keys are the atom indices in the subsystem, values are the keys in the full system.


        """
        self.info = info
        self.calculator = info.calculator
        self.atoms_for_binding = atoms.copy()
        self.atoms_for_subsystem = atoms.copy()
        self.index_map = index_map
        self.reverse_index_map = reverse_index_map
        self.potential_energy = None
        self.forces = None
        self.link_interaction_correction = 0
        self.density_grid = None
        self.link_atom_indices = []
        self.initial_charges = atoms.get_initial_charges()
        self.dft_system = hasattr(self.calculator, "get_pseudo_density")

        # If the cell size minimization flag has been enabled, then try to reduce the
        # cell size
        if info.cell_size_optimization:
            pbc = atoms.get_pbc()
            if pbc[0] or pbc[1] or pbc[2]:
                warn(("Cannot optimize cell size when periodic boundary"
                      "condition have been enabled, disabling optimization."), 3)
                self.info.cell_size_optimization = False
            else:
                self.minimize_cell()

    def minimize_cell(self):
        """@todo: Docstring for minimize_cell.
        :returns: @todo

        """
        cell = self.atoms_for_subsystem.get_cell()
        x_min, y_min, z_min = self.atoms_for_subsystem[0].position
        x_max, y_max, z_max = self.atoms_for_subsystem[0].position

    def update_density_grid(self):
        """Precalculates a grid 3D grid of 3D points for the charge pseudo
        density calculation.
        """
        calc = self.calculator
        atoms = self.atoms_for_subsystem

        # One calculation has to be made before the grid points can be asked
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

        self.density_grid = grid

    def update_charges(self):
        """If DFT calculators are used, updates the atomic charges according to
        the electron pseudodensity.

        The charge for each atom in the system is integrated from the electron
        density inside the van Der Waals radius of the atom in hand. The link
        atoms will affect the distribution of the electron density, but through
        normalization it is ensured that link atoms will not change the total
        charge in the system
        """
        # The pseudo-density is calculated from the atoms that contain the
        # possible link atoms. This takes into account the effect of the bond
        # on the electronic density.
        atoms = self.atoms_for_subsystem
        calc = self.calculator
        density_grid = np.array(calc.get_pseudo_density(atoms))
        grid = self.density_grid

        initial_charges = self.initial_charges
        atomic_numbers = self.atoms_for_binding.get_atomic_numbers()
        total_electron_charge = -np.sum(np.array(self.atoms_for_binding.get_atomic_numbers()))

        projected_charges = []

        for i_atom, atom in enumerate(atoms):
            if i_atom not in self.link_atom_indices:
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
        # the whole system excluding the link atom to preserve charge neutrality
        total_charge = np.sum(np.array(projected_charges))
        projected_charges = np.array(projected_charges)
        projected_charges *= total_electron_charge/total_charge

        # Add the nuclear charges and initial charges
        projected_charges += atomic_numbers
        projected_charges += initial_charges

        # Delete the link atom charges
        projected_charges = np.delete(projected_charges, self.link_atom_indices)

        # Set the calculated charges to the atoms
        self.atoms_for_binding.set_initial_charges(projected_charges.tolist())

    def get_potential_energy_without_corrections(self):
        """@todo: Docstring for get_potential_energy.
        :returns: @todo
        """
        # Update the cell size if minimization is on
        if self.info.cell_size_optimization:
            self.minimize_cell()

        # Ask the energy from the modified atoms (which include possible link
        # atoms), and and possible corrections
        self.potential_energy = self.calculator.get_potential_energy(
            self.atoms_for_subsystem)

        # Update the calculation grid on dft systems, so that it is available
        # for charge updating
        if self.dft_system:
            self.update_density_grid()

        return self.potential_energy

    def get_forces(self):
        """@todo: Docstring for get_potential_energy.
        :returns: @todo
        """
        # Update the cell size if minimization is on
        if self.info.cell_size_optimization:
            self.minimize_cell()

        # Calculate the forces, ignore forces on link atoms
        forces = self.calculator.get_forces(self.atoms_for_subsystem)

        force_map = []
        for full_index, sub_index in self.index_map.iteritems():
            force = forces[sub_index, :]
            index_force_pair = (full_index, force)
            force_map.append(index_force_pair)

        # Update the calculation grid on dft systems, so that it is available
        # for charge updating
        if self.dft_system:
            self.update_density_grid()

        self.forces = force_map
        return self.forces

