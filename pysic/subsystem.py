#! /usr/bin/env python
"""Defines classes used for creating and storing information about subsystems
in a HybridCalculator."""

import numpy as np
from pysic.utility.error import warn
import ase.data
from pysic.utility.timer import Timer
from ase import Atom, Atoms
from ase.visualize import view
from ase.io import write
import copy


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
        self.calculator = copy.deepcopy(calculator)

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
        if (indices is not None) and (tag is None) and (special_set is None):
            return True
        elif (indices is None) and (tag is not None) and (special_set is None):
            return True
        elif (indices is None) and (tag is None) and (special_set is not None):
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
    def __init__(self, atoms, info, index_map, reverse_index_map, record_time_usage):
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
        self.density_grid = None
        self.calculated_charges = None
        self.pseudo_density = None
        self.link_atom_indices = []
        self.timer = Timer(record_time_usage,
                           {
                               "Charge update": 0,
                               "Calculator": 0,
                               "Density grid update": 0,
                               "Cell minimization": 0
                           })

        # The older ASE versions do not support get_initial_charges()
        try:
            charges = np.array(atoms.get_initial_charges())
        except:
            charges = np.array(atoms.get_charges())
        self.initial_charges = charges

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
        """Creates a new cell for the subsystem. This cell is made
        :returns: @todo

        """
        self.timer.start("Cell minimization")
        padding = self.info.dft_padding
        x_min, y_min, z_min = self.atoms_for_subsystem[0].position
        x_max, y_max, z_max = self.atoms_for_subsystem[0].position

        for atom in self.atoms_for_subsystem:
            r = atom.position
            x = r[0]
            y = r[1]
            z = r[2]

            if x > x_max:
                x_max = x
            if y > y_max:
                y_max = y
            if z > z_max:
                z_max = z

            if x < x_min:
                x_min = x
            if y < y_min:
                y_min = y
            if z < z_min:
                z_min = z

        optimized_cell = np.array([2*padding + x_max - x_min, 2*padding + y_max - y_min, 2*padding + z_max - z_min])
        self.atoms_for_subsystem.set_cell(optimized_cell)
        self.atoms_for_subsystem.center()
        self.timer.end()

    def update_density_grid(self):
        """Precalculates a grid 3D grid of 3D points for the charge pseudo
        density calculation.
        """
        self.timer.start("Density grid update")
        calc = self.calculator
        atoms = self.atoms_for_subsystem

        # One calculation has to be made before the grid points can be asked
        grid_dim = calc.get_number_of_grid_points()

        nx = grid_dim[0]
        ny = grid_dim[1]
        nz = grid_dim[2]
        cell = atoms.get_cell()
        cx = cell[0]
        cy = cell[1]
        cz = cell[2]

        cux = cx / float((grid_dim[0] - 1))
        cuy = cy / float((grid_dim[1] - 1))
        cuz = cz / float((grid_dim[2] - 1))

        # Create a 3D array of 3D points. Each point is a  xyz-coordinate to a
        # position where the density has been calculated.
        density_grid = np.zeros((nx, ny, nz, 3))
        for x in range(grid_dim[0]):
            for y in range(grid_dim[1]):
                for z in range(grid_dim[2]):
                    r = x * cux + y * cuy + z * cuz
                    density_grid[x, y, z, :] = r

        self.density_grid = density_grid
        self.timer.end()

    def update_charges(self):
        """If DFT calculators are used, updates the atomic charges according to
        the electron pseudodensity.

        The charge for each atom in the system is integrated from the electron
        density inside the van Der Waals radius of the atom in hand. The link
        atoms will affect the distribution of the electron density, but through
        normalization it is ensured that link atoms will not change the total
        charge in the system
        """
        self.timer.start("Charge update")

        # Turn debugging on or off here
        debugging = False

        atoms_without_links = self.atoms_for_binding
        atoms_with_links = self.atoms_for_subsystem
        calc = self.calculator

        # The electron density is calculated from the system with link atoms.
        # This way the link atoms can modify the charge distribution
        calc.set_atoms(atoms_with_links)
        density = np.array(calc.get_pseudo_density())

        # Write the charge density as .cube file for VMD
        if debugging:
            write('nacl.cube', atoms_with_links, data=density)

        initial_charges = self.initial_charges
        atomic_numbers = np.array(atoms_without_links.get_atomic_numbers())
        total_electron_charge = -np.sum(atomic_numbers)

        grid = self.density_grid

        if debugging:
            debug_list = []

        # The link atoms are at the end of the list
        n_atoms = len(atoms_without_links)
        projected_charges = np.zeros((1, n_atoms))

        for i_atom in range(n_atoms):
            atom = atoms_with_links[i_atom]
            r_atom = atom.position
            z = atom.number

            # Get the van Der Waals radius
            R = ase.data.vdw.vdw_radii[z]

            # Create a 3 x 3 x 3 x 3 array that can be used for vectorized
            # operations with the density grid
            r_atom_array = np.tile(r_atom, (grid.shape[0], grid.shape[1], grid.shape[2], 1))

            diff = grid - r_atom_array
            diff = np.linalg.norm(diff, axis=3)
            indices = np.where(diff <= R)
            densities = density[indices]
            atom_charge = np.sum(densities)
            projected_charges[0, i_atom] = atom_charge

            if debugging:
                debug_list.append((atom, indices, densities))

        #DEBUG: Visualize the grid and contributing grid points as atoms
        if debugging:
            d = Atoms()
            d.set_cell(atoms_with_links.get_cell())

            # Visualize the integration spheres with atoms
            for point in debug_list:
                atom = point[0]
                indices = point[1]
                densities = point[2]
                d.append(atom)
                print "Atom: " + str(atom.symbol) + ", Density sum: " + str(np.sum(densities))
                print "Density points included: " + str(len(densities))
                for i in range(len(indices[0])):
                    x = indices[0][i]
                    y = indices[1][i]
                    z = indices[2][i]
                    a = Atom('H')
                    a.position = grid[x, y, z, :]
                    d.append(a)
            view(d)

        # Normalize the projected charges according to the electronic charge in
        # the whole system excluding the link atom to preserve charge neutrality
        total_charge = np.sum(np.array(projected_charges))
        projected_charges *= total_electron_charge/total_charge

        # Add the nuclear charges and initial charges
        projected_charges += atomic_numbers
        projected_charges += initial_charges

        # Set the calculated charges to the atoms
        self.atoms_for_binding.set_initial_charges(projected_charges[0, :].tolist())
        self.calculated_charges = projected_charges
        self.pseudo_density = density
        self.timer.end()

    def get_potential_energy(self):
        """@todo: Docstring for get_potential_energy.
        :returns: @todo
        """
        # Update the cell size if minimization is on
        if self.info.cell_size_optimization:
            self.minimize_cell()

        # Ask the energy from the modified atoms (which include possible link
        # atoms), and and possible corrections
        self.timer.start("Calculator")
        self.potential_energy = self.calculator.get_potential_energy(
            self.atoms_for_subsystem)
        self.timer.end()

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
        self.timer.start("Calculator")
        forces = self.calculator.get_forces(self.atoms_for_subsystem)
        self.timer.end()

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

    def get_pseudo_density(self):
        """Returns the electron pseudo density if available."""
        if self.pseudo_density is not None:
            return self.pseudo_density
        else:
            if hasattr(self.calculator, "get_pseudo_density"):
                return self.calculator.get_pseudo_density(self.atoms_for_subsystem)
            else:
                warn("The pseudo density for subsystem \"" + self.info.name + "\" is not available.", 2)
