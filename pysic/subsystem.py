#! /usr/bin/env python
"""Defines classes used for creating and storing information about subsystems
in a Pysic QM/MM hybrid calculation created with the HybridCalculator-class.
"""
import numpy as np
from pysic.utility.error import warn, error
from pysic.utility.bader_charges import get_bader_charges
import ase.data
from pysic.utility.timer import Timer
from ase import Atom, Atoms
from ase.visualize import view
from ase.io import write
import copy


#===============================================================================
class SubSystem(object):
    """Used to create and store information about a subsystem.

    The end user can create and manipulate these objects when defining
    subsystems. The subsystems are added to the calculation by calling
    :meth:`~pysic.hybridcalculator.HybridCalculator.add_subsystem` which adds a
    SubSystem-object to the calculation.

    When the HybridCalculator sees fit, the subsystems are materialized by
    converting the stored SubSystems into SubSystemInternals.

    Attributes:
        name: string
            The unique name for this subsystem.
        calculator: ASE Calculator
            The calculator used.
        cell_size_optimization_enabled: bool
            -
        cell_padding: float
            The padding used when optimizing the cell size.
        charge_calculation_enabled: bool
            -
        charge_source: string
            Indicates the electron density that is used in charge calculation.
            Can be "pseudo" or "all-electron".
        division: string
            Indicates the division algorithm used in charge caluclation. Can be
            "Bader" or "van Der Waals".
        gridrefinement: int
            The factor by which the calculation grid is densified in charge
            calculation.
        indices: list of ints
            -
        tag: int
            -
    """

    def __init__(self,
                 name,
                 indices=None,
                 tag=None,
                 calculator=None
                 ):
        """
        Parameters:
            name: string
                A unique identifier for this substring.
            indices: list of integers or string
                A list of atom indices or a special string "remaining", which
                assigns all the yet unassigned atoms to this subsystem.
            tag: int
                The atoms with this tag belong to this subsystem.
            calculator: ASE compatible calculator
                The calculator that is used for this subsystem.
        """
        self.name = name
        self.calculator = calculator
        self.cell_size_optimization_enabled = False
        self.cell_padding = None
        self.charge_calculation_enabled = False
        self.charge_source = None
        self.division = None
        self.gridrefinement = None

        self.set_atoms(indices, tag)

    def is_valid(self, indices, tag):
        """Checks that the given atom specifiers are correctly given. Does not
        yet check that they exist or don't overlap with other subssystems.
        """
        # Determine how the atoms are specified: indices, tag or special set
        if isinstance(indices, str):
            if indices != "remaining":
                error("Use \"remaining\" if you want to assign the yet unassigned atoms to a subsystem")
        elif (indices is None) and (tag is None):
            error("Provide system as indices or tag",)

    def set_calculator(self, calculator):
        """Set the calculator for the subsystem.

        Parameters:
            calculator: ASE compatible calculator
        """
        self.calculator = copy.copy(calculator)

    def set_atoms(self, indices=None, tag=None):
        """Set the atoms that belong to this subsystem. Give only one of the
        specifiers: indices or tag.

        Parameters:
            indices: list of integers or string
                A list of atom indices or a special string "remaining", which
                assigns all the yet unassigned atoms to this subsystem.
            tag: int
                The atoms with this tag belong to this subsystem.
        """
        self.is_valid(indices, tag)
        if indices is not None and type(indices) is int:
                self.indices = (indices,)
        else:
            self.indices = indices
        self.tag = tag

    def enable_cell_optimization(self, padding):
        """Enable cell size optimization.

        A subsystem might spatially reside in only a small portion of the
        entire system. DFT calculators will then waste time doing calculations
        in empty space, where almost none of electron density reaches.

        This optimization minimizes the cell size, so that the atoms in the
        subsystem fit the cell with the given padding. If the padding is too
        small, the DFT-calculator might not work properly!

        The optimization is off by default. It cannot be turned on in systems
        with periodic boundary conditions. The new optimized cell is always
        ortorhombic, regardless of the shape of the original, unoptimized cell.

        Parameters:
            padding: float
                The minimum distance between the subsystem atoms and the cell
                walls.
        """
        self.cell_size_optimization_enabled = True
        self.cell_padding = padding

    def enable_charge_calculation(self, division="Bader", source="all-electron", gridrefinement=4):
        """Enable the dynamic calculation of atom-centered charges with the
        specified algorithm and from the specified electron density. These
        charges are only used for the interaction between other subsystems.

        Parameters:
            division: string
                Indicates the division algorithm that is used. Available options are:
                    
                    - "Bader": Bader algorithm
                    - "van Der Waals": Spheres with van Der Waals radius

            source: string
                Indicates what type of electron density is used. Available
                options are:
                
                    - "pseudo": Use the pseudo electron density provided by all ASE DFT calculators
                    - "all-electron": Use the all-electron density provided by at least GPAW

            gridrefinement: int
                Indicates the subdivision that is used for the all-electron
                density.  Can be other than unity only for Bader algorithm with
                all-electron density.
        """
        divisions = ["Bader", "van Der Waals"]
        charge_sources = ["pseudo", "all-electron"]

        if division not in divisions:
            error("Invalid division algorithm: " + division)

        if source not in charge_sources:
            error("Invalid source for electron density: " + source)

        if gridrefinement != 1:
            if (division != "Bader") or (division == "Bader" and source != "all-electron"):
                warn("The gridrefinement is available only for the Bader algorithm with all-electron density, it is ignored.", 3)

        self.charge_calculation_enabled = True
        self.division = division
        self.charge_source = source
        self.gridrefinement = gridrefinement


#===============================================================================
class SubSystemInternal(object):
    """A materialization of a SubSystem object.

    This class is materialised from a SubSystem, and should not be
    accessible to the end user.

    Attributes:
        name: string
            The unique name for this subsystem.
        calculator: ASE Calculator
            The calculator used.
        cell_size_optimization_enabled: bool
            -
        cell_padding: float
            -
        charge_calculation_enabled: bool
            -
        charge_source: string
            -
        division: string
            -
        gridrefinement: int
            -
        n_atoms: int
            Number of atoms in the full system.
        atoms_for_interaction: ASE Atoms
            The copy of subsystems atoms used in interaction calculations.
        atoms_for_subsystem: ASE Atoms
            The copy of subsystems atoms used in calculating subsystem energies
            etc.
        index_map: dictionary of int to int
            The keys are the atom indices in the full system, values are
            indices in the subssystem.
        reverse_index_map: dicitonary of int to int
            The keys are the atom indices in the subsystem, values are the
            keys in the full system.
        potential_energy: float
            -
        forces: numpy array
            -
        density_grid: numpy array
            Stored if spherical division is used in charge calculation.
        pseudo_density: numpy array
            -
        link_atom_indices: list
            -
        timer: :class:'~pysic.utility.timer.Timer'
            Used to keep track of time usage.
    """
    def __init__(self, atoms, info, index_map, reverse_index_map, n_atoms):
        """
        Parameters:
            atoms: ASE Atoms
                The subsystem atoms.
            info: SubSystem object
                Contains all the information about the subsystem
            index_map: dictionary of int to int
                The keys are the atom indices in the full system, values are
                indices in the subssystem.
            reverse_index_map: dicitonary of int to int
                The keys are the atom indices in the subsystem, values are the
                keys in the full system.
            n_atoms: int
                Number of atoms in the full system.
        """
        # Extract data from info
        self.name = info.name
        self.calculator = copy.copy(info.calculator)
        self.cell_size_optimization_enabled = info.cell_size_optimization_enabled
        self.cell_padding = info.cell_padding
        self.charge_calculation_enabled = info.charge_calculation_enabled
        self.charge_source = info.charge_source
        self.division = info.division
        self.gridrefinement = info.gridrefinement

        self.n_atoms = n_atoms
        self.atoms_for_interaction = atoms.copy()
        self.atoms_for_subsystem = atoms.copy()
        self.index_map = index_map
        self.reverse_index_map = reverse_index_map
        self.potential_energy = None
        self.forces = None
        self.density_grid = None
        self.pseudo_density = None
        self.link_atom_indices = []
        self.timer = Timer([
            "Bader charge calculation",
            "van Der Waals charge calculation",
            "Energy",
            "Forces",
            "Density grid update",
            "Cell minimization"])

        # The older ASE versions do not support get_initial_charges()
        try:
            charges = np.array(atoms.get_initial_charges())
        except:
            charges = np.array(atoms.get_charges())
        self.initial_charges = charges

        ## Can't enable charge calculation on non-DFT calculator
        self.dft_system = hasattr(self.calculator, "get_pseudo_density")
        if self.charge_calculation_enabled is True and not self.dft_system:
            error("Can't enable charge calculation on non-DFT calculator!")

        # If the cell size minimization flag has been enabled, then try to reduce the
        # cell size
        if self.cell_size_optimization_enabled:
            pbc = atoms.get_pbc()
            if pbc[0] or pbc[1] or pbc[2]:
                warn(("Cannot optimize cell size when periodic boundary"
                      "condition have been enabled, disabling optimization."), 2)
                self.cell_size_optimization_enabled = False
            else:
                self.optimize_cell()

    def optimize_cell(self):
        """Tries to optimize the cell of the subsystem so only the atoms in
        this subsystem fit in it with the defined padding to the edges. The new
        cell is always ortorhombic.
        """
        self.timer.start("Cell minimization")
        padding = self.cell_padding
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
        self.timer.stop()

    def update_density_grid(self):
        """Precalculates a grid of 3D points for the charge calculation with
        van Der Waals radius.
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
        self.timer.stop()

    def update_charges(self):
        """Updates the charges in the system. Depending on the value of
        self.division, calls either
        :meth:`~pysic.subsystem.SubSystemInternal.update_charges_bader`, or
        :meth:`~pysic.subsystem.SubSystemInternal.update_charges_van_der_waals`
        """
        if self.charge_calculation_enabled:
            if self.division == "van Der Waals":
                self.update_charges_van_der_waals()
            if self.division == "Bader":
                self.update_charges_bader()

    def update_charges_bader(self):
        """Updates the charges in the atoms used for interaction with the Bader
        algorithm.

        This function uses an external Bader charge calculator from
        http://theory.cm.utexas.edu/henkelman/code/bader/. This tool is
        provided also in pysic/tools. Before using this function the bader
        executable directory has to be added to PATH.
        """
        self.timer.start("Bader charge calculation")

        # The charges are calculated from the system that includes the link atoms
        bader_charges = get_bader_charges(
            self.atoms_for_subsystem,
            self.calculator,
            self.charge_source,
            self.gridrefinement)

        # Set the calculated charges to the interaction atoms. The call for
        # charges was changed between ASE 3.6 and 3.7. Ignore the link atoms
        # from the list of Bader charges
        n_limit = len(self.atoms_for_interaction)

        try:
            self.atoms_for_interaction.set_initial_charges(bader_charges[0:n_limit])
        except:
            self.atoms_for_interaction.set_charges(bader_charges[0:n_limit])

        self.timer.stop()

    def update_charges_van_der_waals(self):
        """Updates the atomic charges by using the electron density within a
        sphere of van Der Waals radius.

        The charge for each atom in the system is integrated from the electron
        density inside the van Der Waals radius of the atom in hand. The link
        atoms will affect the distribution of the electron density.
        """
        self.timer.start("van Der Waals charge calculation")

        # Turn debugging on or off here
        debugging = False

        atoms_with_links = self.atoms_for_subsystem
        calc = self.calculator

        # The electron density is calculated from the system with link atoms.
        # This way the link atoms can modify the charge distribution
        calc.set_atoms(atoms_with_links)

        if self.charge_source == "pseudo":
            try:
                density = np.array(calc.get_pseudo_density())
            except AttributeError:
                error("The DFT calculator on subsystem \"" + self.name + "\" doesn't provide pseudo density.")

        if self.charge_source == "all-electron":
            try:
                density = np.array(calc.get_all_electron_density(gridrefinement=1))
            except AttributeError:
                error("The DFT calculator on subsystem \"" + self.name + "\" doesn't provide all electron density.")

        # Write the charge density as .cube file for VMD
        if debugging:
            write('nacl.cube', atoms_with_links, data=density)

        grid = self.density_grid

        if debugging:
            debug_list = []

        # The link atoms are at the end of the list
        n_atoms = len(atoms_with_links)
        projected_charges = np.zeros((1, n_atoms))

        for i_atom, atom in enumerate(atoms_with_links):
            r_atom = atom.position
            z = atom.number

            # Get the van Der Waals radius
            R = ase.data.vdw.vdw_radii[z]

            # Create a 3 x 3 x 3 x 3 array that can be used for vectorized
            # operations with the density grid
            r_atom_array = np.tile(r_atom, (grid.shape[0], grid.shape[1], grid.shape[2], 1))

            diff = grid - r_atom_array

            # Numpy < 1.8 doesn't recoxnize axis argument on norm. This is a
            # workaround for diff = np.linalg.norm(diff, axis=3)
            diff = np.apply_along_axis(np.linalg.norm, 3, diff)
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
        atomic_numbers = np.array(atoms_with_links.get_atomic_numbers())
        total_electron_charge = -np.sum(atomic_numbers)
        total_charge = np.sum(np.array(projected_charges))
        projected_charges *= total_electron_charge/total_charge

        # Add the nuclear charges and initial charges
        projected_charges += atomic_numbers

        # Set the calculated charges to the atoms.  The call for charges was
        # changed between ASE 3.6 and 3.7
        try:
            self.atoms_for_interaction.set_initial_charges(projected_charges[0, :].tolist())
        except:
            self.atoms_for_interaction.set_charges(projected_charges[0, :].tolist())

        self.pseudo_density = density
        self.timer.stop()

    def get_potential_energy(self):
        """Returns the potential energy contained in this subsystem.
        """
        # Update the cell size if minimization is on
        if self.cell_size_optimization_enabled:
            self.optimize_cell()

        # Ask the energy from the modified atoms (which include possible link
        # atoms)
        self.timer.start("Energy")
        self.potential_energy = self.calculator.get_potential_energy(
            self.atoms_for_subsystem)
        self.timer.stop()

        # Update the calculation grid if charge calculation with van Der Waals
        # division is enabled:
        if self.charge_calculation_enabled:
                if self.division == "van Der Waals":
                    self.update_density_grid()

        return copy.copy(self.potential_energy)

    def get_forces(self):
        """Returns a 3D numpy array that contains forces for this subsystem.

        The returned array contains a row for each atom in the full system, but
        there is only an entry for the atoms in this subsystem. This makes it
        easier to calculate the total forces later on.
        """
        # Update the cell size if minimization is on
        if self.cell_size_optimization_enabled:
            self.optimize_cell()

        # Calculate the forces
        self.timer.start("Forces")
        forces = self.calculator.get_forces(self.atoms_for_subsystem)
        self.timer.stop()

        # Ignore forces on link atoms, link atoms are at the end of the list
        forces = forces[0:len(self.atoms_for_interaction), :]

        # Store the forces in a suitable numpy array that can be added to the forces of
        # the whole system
        full_forces = np.zeros((self.n_atoms, 3))
        for sub_index in range(len(self.atoms_for_interaction)):
            full_index = self.reverse_index_map[sub_index]
            force = forces[sub_index, :]
            full_forces[full_index, :] = force

        # Update the calculation grid if charge calculation with van Der Waals
        # division is enabled:
        if self.charge_calculation_enabled:
                if self.division == "van Der Waals":
                    self.update_density_grid()

        self.forces = full_forces
        return copy.copy(self.forces)

    def get_pseudo_density(self):
        """Returns the electron pseudo density if available.
        """
        if self.pseudo_density is not None:
            return copy.copy(self.pseudo_density)
        else:
            if hasattr(self.calculator, "get_pseudo_density"):
                return copy.copy(self.calculator.get_pseudo_density(self.atoms_for_subsystem))
            else:
                warn("The pseudo density for subsystem \"" + self.name + "\" is not available.", 2)
