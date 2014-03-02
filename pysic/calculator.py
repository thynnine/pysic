#! /usr/bin/env python
"""The main module of Pysic.
    
This module defines the user interface in Pysic for setting up potentials
and calculators.
"""

from pysic.core import *
from pysic.utility.error import *
from pysic.interactions.local import Potential, ProductPotential
from pysic.interactions.compound import CompoundPotential
from pysic.interactions.bondorder import Coordinator, BondOrderParameters
from pysic.interactions.coulomb import CoulombSummation
from pysic.charges.relaxation import ChargeRelaxation

import pysic.pysic_fortran as pf
import pysic.utility.f2py as pu

import numpy as np
import numpy.linalg as npla
import ase.calculators.neighborlist as nbl
from itertools import permutations
import copy
import math

import pysic.utility.debug as d



class FastNeighborList(nbl.NeighborList):
    """ASE has a neighbor list class built in, `ASE NeighborList`_, but its implementation is
        currently inefficient, and building of the list is an :math:`O(n^2)`
        operation. This neighbor list class overrides the 
        :meth:`~pysic.calculator.FastNeighborList.build` method with
        an :math:`O(n)` time routine. The fast routine is based on a
        spatial partitioning algorithm.
        
        The way cutoffs are handled is also somewhat different to the original
        ASE list. In ASE, the distances for two atoms are compared against
        the sum of the individual cutoffs + neighbor list skin. This list, however,
        searches for the neighbors of each atom at a distance of the cutoff of the
        given atom only, plus skin.
        
        
        .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
        .. _ASE NeighborList: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#building-neighbor-lists
        """

    neighbor_marginal = 0.5
    """Default skin width for the neighbor list"""

    
    def __init__(self, cutoffs, skin=None):
        if skin is None:
            skin = FastNeighborList.neighbor_marginal
        nbl.NeighborList.__init__(self,
                              cutoffs=cutoffs, 
                              skin=skin, 
                              sorted=False, 
                              self_interaction=False,
                              bothways=True)    
    
    def build(self,atoms):
        """Builds the neighbor list.
            
            The routine requires that the given atomic structure matches
            the one in the core. This is because the method invokes the
            Fortran core to do the neighbor search.
            The method overrides the similar
            method in the original ASE neighborlist class, which directly operates
            on the given structure, so this method also takes the atomic structure 
            as an argument. However, in order to keep the core modification routines in
            the :class:`~pysic.calculator.Pysic` class, this method does not change the core
            structure. It does raise an error if the structures do not match, though.
            
            The neighbor search is done via the :meth:`generate_neighbor_lists` routine.
            The routine builds the neighbor list in the core, after which the list is
            fed back to the :class:`~pysic.calculator.FastNeighborList` object by looping over all
            atoms and saving the lists of neighbors and offsets.

            Parameters:
            
            atoms: `ASE Atoms`_ object
                the structure for which the neighbors are searched
            """
        
        if not Pysic.core.atoms_ready(atoms):
            raise MissingAtomsError("Neighbor list building: Atoms in the core do not match.")
        if Pysic.core.get_atoms() != atoms:
            raise MissingAtomsError("Neighbor list building: Atoms in the core do not match.")
        
        self.positions = atoms.get_positions()
        self.pbc = atoms.get_pbc()
        self.cell = atoms.get_cell()
        
        
        pf.pysic_interface.generate_neighbor_lists(self.cutoffs)
                
        self.neighbors = [np.empty(0, int) for a in range(len(atoms))]
        self.displacements = [np.empty((0, 3), int) for a in range(len(atoms))]
        
        for i in range(len(atoms)):
            n_nbs = pf.pysic_interface.get_number_of_neighbors_of_atom(i)
            if n_nbs > 0:
                (self.neighbors[i], self.displacements[i]) = pf.pysic_interface.get_neighbor_list_of_atom(i,n_nbs)
                # the offsets are in Fortran array format, so they need to be transposed
                self.displacements[i] = np.transpose(self.displacements[i])
    
        self.nupdates += 1
    
    
    def get_neighbors(self, index, atoms=None, sort=False):
        """Returns arrays containing the indices and offsets of the neighbors of the given atom.
        
        Overrides the method in `ASE NeighborList`_.
        
        Parameters:
        
        index: integer
            the index of the central atom
        atoms: ASE Atoms
            the atoms object containing the absolute coordinates - needed only if sorting is necessary
        sort: boolean
            if True, the list will be sorted according to distance
        """
        n_nbs = pf.pysic_interface.get_number_of_neighbors_of_atom(index)
        nbors = self.neighbors[index][0:n_nbs]
        displ = self.displacements[index][0:n_nbs]
        
        if sort and atoms is not None and n_nbs > 0:
            dists = self.get_neighbor_distances(index, atoms)
            sorted_lists = np.array([np.append(o,np.array(n)) for (d,n,o) in sorted(zip(dists, nbors, displ), key=lambda dis: dis[0])])
            return sorted_lists[:,3], sorted_lists[:,0:3]
        
        return nbors, displ
    
    
    def get_neighbor_separations(self, index, atoms, sort=False):
        """Returns an array of atom-atom separation vectors between the given atom and its neighbors.
        
        Parameters:
        
        index: integer
            the index of the central atom
        atoms: ASE Atoms
            the atoms object containing the absolute coordinates
        sort: boolean
            if True, the list will be sorted according to distance
        """
        
        indices, offsets = self.get_neighbors(index)
        separations = np.array(len(indices)*[[0.0,0.0,0.0]])
        dists = np.array(len(separations)*[0.0])
        running = 0
        for i, offset in zip(indices, offsets):
            sep = atoms.positions[i] + np.dot(offset, atoms.get_cell()) - atoms.positions[index]
            
            separations[running] = sep
            if sort:
                dists[running] = math.sqrt(np.dot(sep,sep))                
            running += 1

        if sort:
            return np.array([s for (d,s) in sorted(zip(dists, separations), key=lambda dis: dis[0])])
        
        return separations

    def get_neighbor_distances(self, index, atoms, sort=False):
        """Returns a list of atom-atom distances between the given atom and its neighbors.
        
        Parameters:
        
        index: integer
            the index of the central atom
        atoms: ASE Atoms
            the atoms object containing the absolute coordinates
        sort: boolean
            if True, the list will be sorted according to distance
        """
        separations = self.get_neighbor_separations(index, atoms)
        dists = np.array(len(separations)*[0.0])
        running = 0
        for sep in separations:
            dists[running] = math.sqrt(np.dot(sep,sep))
            running += 1

        if sort:
            return np.array(sorted(dists))
            
        return dists


class Pysic:
    """A calculator class providing the necessary methods for interfacing with `ASE`_.

    Pysic is a calculator for evaluating energies and forces for given atomic structures
    according to the given :class:`~pysic.interactions.local.Potential` set. Neither the geometry nor the
    potentials have to be specified upon creating the calculator, as they can be specified
    or changed later. They are necessary for actual calculation, of course.

    Simulation geometries must be defined as `ASE Atoms`_. This object contains both the
    atomistic coordinates and supercell parameters.

    Potentials must be defined as a list of :class:`~pysic.interactions.local.Potential` objects. 
    The total potential of the system is then the sum of the individual potentials.
    
    .. _ASE: https://wiki.fysik.dtu.dk/ase/
    .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

    Parameters:

    atoms: `ASE Atoms`_ object
        an Atoms object containing the full simulation geometry
    potentials: list of :class:`~pysic.interactions.local.Potential` objects
        list of potentials for describing interactions
    force_initialization: boolean
        If true, calculations always fully initialize the Fortran core.
        If false, the Pysic tries to evaluate what needs updating by
        consulting the :data:`~pysic.calculator.Pysic.core` instance of :class:`~pysic.core.CoreMirror`.
    """

    core = CoreMirror()
    """An object storing the data passed to the core.

    Whenever a :class:`~pysic.calculator.Pysic` calculator alters the Fortran core,
    it should also modify the :data:`~pysic.calculator.Pysic.core` object so that
    it is always a valid representation of the actual core.
    Then, whenever :class:`~pysic.calculator.Pysic` needs to check if the
    representation in the core is up to date, it only needs to compare
    against :data:`~pysic.calculator.Pysic.core` instead of accessing the
    Fortran core itself.
    """
    def __init__(self,atoms=None,potentials=None,charge_relaxation=None,
                 coulomb=None,full_initialization=False):
        
        self.neighbor_lists_ready = False
        self.saved_cutoffs = None
        
        self.structure = None
        self.neighbor_list = None
        self.potentials = None
        self.charge_relaxation = None
        self.coulomb = None
        self.charges = None
        
        self.set_atoms(atoms)
        self.set_potentials(potentials)
        self.set_charge_relaxation(charge_relaxation)
        self.set_coulomb_summation(coulomb)
        
        self.forces = None
        self.stress = None
        self.energy = None
        self.electronegativities = None

        self.force_core_initialization = full_initialization
    
        self.extra_calculators = []


    def __eq__(self,other):
        try:
            if self.structure != other.structure:
                return False
            if any(self.structure.get_initial_charges() != other.structure.get_initial_charges()):
                return False
            if self.neighbor_list != other.neighbor_list:
                return False
            if self.potentials != other.potentials:
                return False
            if self.extra_calculators != other.extra_calculators:
                return False
        except:
            try:
                if self.structure != other.structure:
                    return False
                if any(self.structure.get_charges() != other.structure.get_charges()):
                    return False
                if self.neighbor_list != other.neighbor_list:
                    return False
                if self.potentials != other.potentials:
                    return False
                if self.extra_calculators != other.extra_calculators:
                    return False
            except:
                return False

        return True

    def __ne__(self,other):
        return not self.__eq__(other)
            

    def __repr__(self):
        return "Pysic(atoms={atoms},potentials={pots},full_initialization={init})".format(atoms=str(self.structure),
                                                                                          pots=str(self.potentials),
                                                                                          init=str(self.force_core_initialization))


    def core_initialization_is_forced(self):
        """Returns true if the core is always fully initialized, false otherwise."""

        return self.force_core_initialization


    def force_core_initialization(self,new_mode):
        """Set the core initialization mode.

        Parameters:

        new_mode: logical
            true if full initialization is required, false if not
        """
        
        self.force_core_initialization = new_mode

    
    def calculation_required(self, atoms=None, 
                             quantities=['forces','energy','stress','electronegativities']):
        """Check if a calculation is required.
        
        When forces or energy are calculated, the calculator saves the
        result in case it is needed several times. This method tells
        if a wanted quantity is not yet calculated for the current
        structure and needs to be calculated explicitly. If a list of
        several quantities is given, the method returns true if any one of
        them needs to be calculated.
        
        Parameters:
        
        atoms: `ASE Atoms`_ object
            ignored at the moment
        quantities: list of strings
            list of keywords 'energy', 'forces', 'stress', 'electronegativities'
        """
        
        do_it = []
        try:
            assert isinstance(quantities, list)
            list_of_quantities = quantities
        except:
            list_of_quantities = [ quantities ]
        
        for mark in list_of_quantities:
            if mark == 'energy':
                do_it.append(self.energy is None)
            elif mark == 'forces':
                do_it.append(self.forces is None)
            elif mark == 'electronegativities':
                do_it.append(self.electronegativities is None)
            elif mark == 'stress':
                do_it.append(self.stress is None)
            else:
                do_it.append(False)
        
        # If the core does not match the Pysic calculator,
        # we may have changed the system or potentials
        # associated with the calculator without telling it.
        # In that case the quantities need to be recalculated.
        # It is of course possible that we have several Pysics
        # changing the core which would lead to unnecessary
        # recalculations.
        if(not Pysic.core.atoms_ready(self.structure)):
            #print "atoms"
            do_it.append(True)
        if(not Pysic.core.charges_ready(self.structure)):
            #print "charges"
            do_it.append(True)
        if(not Pysic.core.cell_ready(self.structure)):
            #print "cell"
            do_it.append(True)
        if(not Pysic.core.potentials_ready(self.potentials)):
            #print "potentials"
            do_it.append(True)
            
        return any(do_it)


    def get_atoms(self):
        """Returns the `ASE Atoms`_ object assigned to the calculator."""
        return self.structure


    def get_neighbor_lists(self):
        """Returns the :class:`~pysic.calculator.FastNeighborList` or `ASE NeighborList`_ 
        object assigned to the calculator.

        The neighbor lists are generated according to the given `ASE Atoms`_ object
        and the :class:`~pysic.interactions.local.Potential` objects of the calculator. Note that the lists
        are created when the core is set or if the method 
        :meth:`~pysic.calculator.Pysic.create_neighbor_lists` is called.
        """
        return self.neighbor_list


    def get_potentials(self):
        """Returns the list of potentials assigned to the calculator."""
        return self.potentials

    def get_calculators(self):
        """Returns the list of other calculators to be used with Pysic."""
        return self.extra_calculators
    
    def get_electronegativities(self, atoms=None):
        """Returns the electronegativities of atoms.
        """
        self.set_atoms(atoms)
        if self.calculation_required(atoms,'electronegativities'):
            self.calculate_electronegativities()
        
        return np.copy(self.electronegativities)
    

    def get_electronegativity_differences(self, atoms=None):
        """Returns the electronegativity differences of atoms from the average of the entire system.
        """
        enegs = self.get_electronegativities(atoms)
        average_eneg = enegs.sum()/len(enegs)
        return enegs - average_eneg

    
    def get_forces(self, atoms=None, skip_charge_relaxation=False):
        """Returns the forces.

        If the atoms parameter is given, it will be used for updating the
        structure assigned to the calculator prior to calculating the forces.
        Otherwise the structure already associated with the calculator is used.

        The calculator checks if the forces have been calculated already
        via :meth:`~pysic.calculator.Pysic.calculation_required`. If the structure
        has changed, the forces are calculated using :meth:`~pysic.calculator.Pysic.calculate_forces`

        Parameters:

        atoms: `ASE atoms`_ object
            the structure for which the forces are determined
        """
        self.set_atoms(atoms)
        if self.calculation_required(atoms,'forces'):
            self.calculate_forces(skip_charge_relaxation=skip_charge_relaxation)

        return np.copy(self.forces)


    def get_potential_energy(self, atoms=None, force_consistent=False,
                             skip_charge_relaxation=False):
        """Returns the potential energy.

        If the atoms parameter is given, it will be used for updating the
        structure assigned to the calculator prior to calculating the energy.
        Otherwise the structure already associated with the calculator is used.

        The calculator checks if the energy has been calculated already
        via :meth:`~pysic.calculator.Pysic.calculation_required`. If the structure
        has changed, the energy is calculated using :meth:`~pysic.calculator.Pysic.calculate_energy`

        Parameters:

        atoms: `ASE atoms`_ object
            the structure for which the energy is determined
        force_consistent: logical
            ignored at the moment
        """
        self.set_atoms(atoms)
        if self.calculation_required(atoms,'energy'):
            self.calculate_energy(skip_charge_relaxation=skip_charge_relaxation)

        return self.energy


    def get_stress(self, atoms=None, skip_charge_relaxation=False):
        """Returns the stress tensor in the format 
        :math:`[\sigma_{xx},\sigma_{yy},\sigma_{zz},\sigma_{yz},\sigma_{xz},\sigma_{xy}]`

        If the atoms parameter is given, it will be used for updating the
        structure assigned to the calculator prior to calculating the stress.
        Otherwise the structure already associated with the calculator is used.

        The calculator checks if the stress has been calculated already
        via :meth:`~pysic.calculator.Pysic.calculation_required`. If the structure
        has changed, the stress is calculated using :meth:`~pysic.calculator.Pysic.calculate_stress`

        Stress (potential part) and force are evaluated in tandem. 
        Therefore, invoking the evaluation of
        one automatically leads to the evaluation of the other. Thus, if you have just
        evaluated the forces, the stress will already be known.
    
        This is because the
        stress tensor is formally defined as
            
        .. math::
        
            \\sigma_{AB} = -\\frac{1}{V} \\sum_i \\left[ m_i (v_i)_A (v_i)_B + (r_i)_A (f_i)_B \\right],
        
            
        where :math:`m`, :math:`v`, :math:`r`, and :math:`f` are mass, velocity,
        position and force of atom :math:`i`, and :math:`A`, :math:`B` denote the
        cartesian coordinates :math:`x,y,z`. 
        (The minus sign is there just to be consistent with the NPT routines in `ASE`_.) 
        However, if periodic boundaries are used,
        the absolute coordinates cannot be used (there would be discontinuities at the
        boundaries of the simulation cell). Instead, the potential energy terms 
        :math:`(r_i)_A (f_i)_B` must be evaluated locally for pair, triplet, and many
        body forces using the relative coordinates of the particles involved in the
        local interactions. These coordinates are only available during the actual force
        evaluation when the local interactions are looped over. Thus, calculating the stress
        requires doing the full force evaluation cycle. On the other hand, calculating the
        stress is not a great effort compared to the force evaluation, so it is convenient
        to evaluate the stress always when the forces are evaluated.
                        
        Parameters:

        atoms: `ASE atoms`_ object
            the structure for which the stress is determined
        """
        self.set_atoms(atoms)
        if self.calculation_required(atoms,'stress'):
            self.calculate_stress(skip_charge_relaxation=skip_charge_relaxation)
        
        # self.stress contains the potential contribution to the stress tensor
        # but we add the kinetic contribution on the fly
        momenta = self.structure.get_momenta()
        masses = self.structure.get_masses()
        velocities = np.divide( momenta, np.array([masses,masses,masses]).transpose() )

        kinetic_stress = np.array([0.0]*6)
        
        # s_xx, s_yy, s_zz, s_yz, s_xz, s_xy
        kinetic_stress[0] = np.dot( momenta[:,0], velocities[:,0] )
        kinetic_stress[1] = np.dot( momenta[:,1], velocities[:,1] )
        kinetic_stress[2] = np.dot( momenta[:,2], velocities[:,2] )
        kinetic_stress[3] = np.dot( momenta[:,1], velocities[:,2] )
        kinetic_stress[4] = np.dot( momenta[:,0], velocities[:,2] )
        kinetic_stress[5] = np.dot( momenta[:,0], velocities[:,1] )
                
        # ASE NPT simulator wants the pressure with an inversed sign
        return np.copy(-( kinetic_stress + self.stress ) / self.structure.get_volume())

    
    def set_atoms(self, atoms=None):
        """Assigns the calculator with the given structure.
            
        This method is always called when any method is given the
        atomic structure as an argument. If the argument is missing
        or None, nothing is done. Otherwise a copy of the given structure
        is saved (according to the instructions in 
        `ASE API <https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#calculator-interface>`_.)
            
        If a structure is already in memory and it is different to the given
        one (as compared with ``__ne__``), it is noted that all quantities
        are unknown for the new system. If the structure is the same as the
        one already known, nothing is done.
        This is because if one wants to
        access the energy of forces of the same system repeatedly, it is unnecessary
        to always calculate them from scratch. Therefore the calculator saves
        the computed values along with a flag stating that the values have been
        computed.
            
        Parameters:

        atoms: `ASE atoms`_ object
            the structure to be calculated
        """
        if atoms == None:
            pass
        else:
            atoms_changed = False
            
            # the call for charges was changed between ASE 3.6 and 3.7
            try:
                atoms_changed = self.structure != atoms or \
                    (self.structure.get_initial_charges() != atoms.get_initial_charges()).any()
            except:
                atoms_changed = self.structure != atoms or \
                    (self.structure.get_charges() != atoms.get_charges()).any()

            if(atoms_changed):
                self.forces = None
                self.energy = None
                self.stress = None
                self.electronegativities = None
                

                # NB: this avoids updating the potential lists every time an atom moves
                try:
                    if((self.structure.get_atomic_numbers() != atoms.get_atomic_numbers()).any()):
                        Pysic.core.potential_lists_ready = False
                        self.neighbor_lists_waiting = False

                    if((self.structure.get_tags() != atoms.get_tags()).any()):
                        Pysic.core.potential_lists_ready = False
                        self.neighbor_lists_waiting = False                

                    if(not Pysic.core.potentials_ready(self.potentials)):
                        Pysic.core.potential_lists_ready = False
                        self.neighbor_lists_waiting = False

                except:
                    Pysic.core.potential_lists_ready = False
                    self.neighbor_lists_waiting = False
            

                self.structure = atoms.copy()


    def set_potentials(self, potentials):
        """Assign a list of potentials to the calculator.
        
        Also a single potential object can be given, instead of a list.
        Note that this method does not conserve any potentials that were
        already known by the calculator. To add potentials to the list
        of known potentials, use :meth:`~pysic.calculator.Pysic.add_potential`.

        Parameters:

        potentials: list of :class:`~pysic.interactions.local.Potential` objects
            a list of potentials to describe interactinos
        """
        self.potentials = []
        if potentials is None:
            pass
        else:
            self.forces = None
            self.energy = None
            self.stress = None
            self.electronegativities = None
            
            new_cutoffs = self.get_individual_cutoffs(1.0)
            self.neighbor_lists_waiting = not self.neighbor_lists_expanded(new_cutoffs)

            try:
                assert isinstance(potentials,list)
                potlist = potentials
            except:
                potlist = [potentials]

            for pot in potlist:
                self.add_potential(pot)
    
    
    def add_calculator(self,calculator):
        """Add a calculator to the list of external calculators.
        
        Parameters:
        
        calculator: an ASE calculator object
            a new calculator to describe interactions
        """
        if self.extra_calculators is None:
            self.extra_calculators = []
            
        if isinstance(calculator,list):
            calcs = calculator
        else:
            calcs = [calculator]
            
        for calc in calcs:
            self.extra_calculators.append(calc)

        self.forces = None
        self.energy = None
        self.stress = None
        self.electronegativities = None



    def add_potential(self, potential):
        """Add a potential to the list of potentials.
        
        Also a list of potentials can be given as an argument.
        In that case, the potentials are added one by one.

        If a :class:`~pysic.interactions.compound.CompoundPotential` is given,
        the addition is done through its 
        :meth:`~pysic.interactions.compound.CompoundPotential.build` method.

        Parameters:

        potential: :class:`.interactions.local.Potential` object
            a new potential to describe interactions
        """

        if self.potentials is None:
            self.potentials = []

        if isinstance(potential,list):
            pots = potential
        else:
            pots = [potential]
    
        for pot in pots:
            if isinstance(pot,CompoundPotential):
                pot.build(self)
            else:
                self.potentials.append(pot)
    
        self.forces = None
        self.energy = None
        self.stress = None
        self.electronegativities = None

        new_cutoffs = self.get_individual_cutoffs(1.0)
        self.neighbor_lists_waiting = not self.neighbor_lists_expanded(new_cutoffs)
    
    
    def remove_calculator(self, calculator):
        """Remove a calculator from the list of calculators.
        
        Parameters:
        
        calculator: an ASE calculator object
            the calculator to be removed
        """
        self.extra_calculators.remove(calculator)

    
    def remove_potential(self, potential):
        """Remove a potential from the list of potentials.
            
        If a :class:`~pysic.interactions.compound.CompoundPotential` is given,
        the removal is done through its 
        :meth:`~pysic.interactions.compound.CompoundPotential.remove` method.
            
        Parameters:
            
        potential: :class:`~pysic.interactions.local.Potential` object
            the potential to be removed
        """

        # cannot just remove a potential due to the cutoffs and other
        # changing factors
        if isinstance(potential, CompoundPotential):
            potential.remove(self)
        else:
            new_pots = self.get_potentials().copy()
            new_pots.remove(potential)
            self.set_potentials(new_pots)

    
    def set_coulomb_summation(self,coulomb):
        """Set the Coulomb summation algorithm for the calculator.
            
            If a Coulomb summation algorithm is set, the Coulomb interactions
            between all charged atoms are evaluated automatically during
            energy and force evaluation. If not, the charges do not directly
            interact.
            
            Parameters:
            
            coulomb: :class:`.interactions.coulomb.CoulombSummation`
                the Coulomb summation algorithm
            """
        self.coulomb = coulomb
        new_cutoffs = self.get_individual_cutoffs(1.0)
        self.neighbor_lists_waiting = not self.neighbor_lists_expanded(new_cutoffs)
    

    def get_coulomb_summation(self):
        """Returns the Coulomb summation algorithm of this calculator.
            """
        return self.coulomb
    
    
    def set_charge_relaxation(self,charge_relaxation):
        """Add a charge relaxation algorithm to the calculator.
            
            If a charge relaxation scheme has been added to the :class:`~pysic.calculator.Pysic`
            calculator, it will be automatically asked to do the charge relaxation 
            before the calculation of energies or forces via 
            :meth:`~pysic.charges.relaxation.ChargeRelaxation.charge_relaxation`.
            
            It is also possible to pass the :class:`~pysic.calculator.Pysic` calculator to the 
            :class:`~pysic.charges.relaxation.ChargeRelaxation` algorithm without creating the opposite
            link using :meth:`~pysic.charges.relaxation.ChargeRelaxation.set_calculator`. 
            In that case, the calculator does not automatically relax the charges, but
            the user can manually trigger the relaxation with 
            :meth:`~pysic.charges.relaxation.ChargeRelaxation.charge_relaxation`.
            
            If you wish to remove automatic charge relaxation, just call this method
            again with None as argument.
            
            Parameters:
            
            charge_relaxation: :class:`~pysic.charges.relaxation.ChargeRelaxation` object
                the charge relaxation algorithm
            """

        try:
            charge_relaxation.set_calculator(self, reciprocal=False)
        except:
            pass
        self.charge_relaxation = charge_relaxation

                
    def get_charge_relaxation(self):
        """Returns the :class:`~pysic.charges.relaxation.ChargeRelaxation` object connected to the calculator.
            """
        return self.charge_relaxation
    
    
    def create_neighbor_lists(self,cutoffs=None,marginal=None):
        """Initializes the neighbor lists.

        In order to do calculations at reasonable speed, the calculator needs 
        a list of neighbors for each atom. For this purpose, either the list
        is built in the Fortran core or the `ASE NeighborList`_
        are used. The ASE lists are very slow, but they will work even for small
        systems where the cutoff is longer than cell size. The custom list
        uses :class:`pysic.calculator.FastNeighborList`, but the build succeeds only
        for cutoffs shorter than cell size. The type of list is determined automatically.
        
        This method initializes these lists according to the given
        cutoffs.

        .. _ASE NeighborList: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#building-neighbor-lists

        Parameters:

        cutoffs: list of doubles
            a list containing the cutoff distance for each atom
        marginal: double
            the skin width of the neighbor list
        """
        
        if marginal is None:
            marginal = FastNeighborList.neighbor_marginal

        
        fastlist = True
        if cutoffs == None:
            cutoffs = self.get_individual_cutoffs(1.0)
        max_cut = np.max(cutoffs)+marginal
        
            
        for i in range(3):
            vec = self.structure.get_cell()[i]
            other_vec1 = self.structure.get_cell()[(i+1)%3]
            other_vec2 = self.structure.get_cell()[(i+2)%3]
            normal = np.cross(other_vec1,other_vec2)
            length = math.fabs(np.dot(vec,normal))/math.sqrt(np.dot(normal,normal))
            if length < max_cut:                
                fastlist = False
                
        if fastlist:
            try:
                self.neighbor_list = FastNeighborList(cutoffs,skin=marginal)
            except:
                fastlist = False

        if not fastlist:
            self.neighbor_list = nbl.NeighborList(cutoffs,skin=marginal,sorted=False,self_interaction=False,bothways=True)

        self.neighbor_lists_waiting = True
        self.set_cutoffs(cutoffs)
    
    
    def get_neighbor_list(self):
        """Returns the neighbor list object.
        """
        return self.neighbor_list
    
            
    def get_individual_cutoffs(self,scaler=1.0):
        """Get a list of maximum cutoffs for all atoms.

        For each atom, the interaction with the longest cutoff is found and
        the associated maximum cutoffs are returned as a list. In case the a list
        of scaled values are required, the scaler can be adjusted. E.g., scaler = 0.5
        will return the cutoffs halved.

        Parameters:

        scaler: double
            a number for scaling all values in the generated list
        """
        if self.structure == None:
            return None
        elif self.potentials == None:
            if self.coulomb == None:
                return self.structure.get_number_of_atoms()*[0.0]
            else:
                return self.structure.get_number_of_atoms()*[self.coulomb.get_realspace_cutoff()]
        else:
            cuts = []            
            # loop over all atoms, with symbol, tags, index containing the corresponding
            # info for a single atom at a time
            for symbol, tags, index in zip(self.structure.get_chemical_symbols(),
                                           self.structure.get_tags(),
                                           range(self.structure.get_number_of_atoms())):
            
                if self.coulomb == None:
                    max_cut = 0.0
                else:
                    max_cut = self.coulomb.get_realspace_cutoff()
                
                for potential in self.potentials:
                    active_potential = False

                    if potential.get_different_symbols().count(symbol) > 0 or potential.get_different_tags().count(tags) > 0 or potential.get_different_indices().count(index) > 0:
                        active_potential = True

                    if active_potential and potential.get_cutoff() > max_cut:
                        max_cut = potential.get_cutoff()

                    try:
                        for bond in potential.get_coordinator().get_bond_order_parameters():
                            active_bond = False
                            if bond.get_different_symbols().count(symbol) > 0:
                                active_bond = True
                                
                            if active_bond:
                                if bond.get_cutoff() > max_cut:
                                    max_cut = bond.get_cutoff()
                    except:
                        pass
            
                cuts.append(max_cut*scaler)
            return cuts


    def calculate_electronegativities(self):
        """Calculates electronegativities.
            
        Calls the Fortran core to calculate forces for the currently assigned structure.
        """
        self.set_core()
        n_atoms = pf.pysic_interface.get_number_of_atoms()
        self.electronegativities = pf.pysic_interface.calculate_electronegativities(n_atoms).transpose()
        
    
    def calculate_forces(self, skip_charge_relaxation=False):
        """Calculates forces (and the potential part of the stress tensor).

        Calls the Fortran core to calculate forces for the currently assigned structure.
            
        If a link exists to a :class:`~pysic.charges.relaxation.ChargeRelaxation`, it is first made to
        relax the atomic charges before the forces are calculated.
        """
        self.set_core()
        if self.charge_relaxation is not None and skip_charge_relaxation == False:
            self.charge_relaxation.charge_relaxation()
        n_atoms = pf.pysic_interface.get_number_of_atoms()
        self.forces, self.stress = pf.pysic_interface.calculate_forces(n_atoms)#.transpose()
        self.forces = self.forces.transpose()

        if not self.extra_calculators is None:
            if len(self.extra_calculators) > 0:
                # calculators should copy the system themselves, so the
                # original system could be passed here. Maybe it is best to
                # give a copy though, just in case? Will take more memory though.
                system_copy = copy.deepcopy(self.structure)
                for calc in self.extra_calculators:
                    self.forces = self.forces + calc.get_forces(system_copy)
                    self.stress = self.stress + calc.get_stress(system_copy)
        
        

    def calculate_energy(self, skip_charge_relaxation=False):
        """Calculates the potential energy.

        Calls the Fortran core to calculate the potential energy for the currently assigned structure.
 
        If a link exists to a :class:`~pysic.charges.relaxation.ChargeRelaxation`, it is first made to
        relax the atomic charges before the forces are calculated.
        """
        self.set_core()

        if self.charge_relaxation is not None and skip_charge_relaxation == False:
            self.charge_relaxation.charge_relaxation()
        #n_atoms = pf.pysic_interface.get_number_of_atoms()
        self.energy = pf.pysic_interface.calculate_energy()

        if not self.extra_calculators is None:
            if len(self.extra_calculators) > 0:
                system_copy = copy.deepcopy(self.structure)
                for calc in self.extra_calculators:
                    self.energy = self.energy + calc.get_potential_energy(system_copy)


    def calculate_stress(self, skip_charge_relaxation=False):
        """Calculates the potential part of the stress tensor (and forces).

        Calls the Fortran core to calculate the stress tensor for the currently assigned structure.
        """

        if self.charge_relaxation is not None and skip_charge_relaxation == False:            self.charge_relaxation.charge_relaxation()
        
        self.set_core()
        n_atoms = pf.pysic_interface.get_number_of_atoms()
        self.forces, self.stress = pf.pysic_interface.calculate_forces(n_atoms)
        self.forces = self.forces.transpose()

        if not self.extra_calculators is None:
            if len(self.extra_calculators) > 0:
                system_copy = copy.deepcopy(self.system)
                for calc in self.extra_calculators:
                    self.forces = self.forces + calc.get_forces(system_copy)
                    self.stress = self.stress + calc.get_stress(system_copy)


    def set_core(self):
        """Sets up the Fortran core for calculation.

        If the core is not initialized, if the number of atoms has changed, or
        if full initialization is forced, the core is initialized from scratch.
        Otherwise, only the atomic coordinates and momenta are updated.
        Potentials, neighbor lists etc. are also updated if they have been edited.
        """        
        
        do_full_init = False
        if self.force_core_initialization:
            do_full_init = True
        elif not Pysic.core.mpi_ready:
            do_full_init = True
        elif Pysic.core.get_atoms() == None:
            do_full_init = True
        elif self.structure.get_number_of_atoms() != Pysic.core.structure.get_number_of_atoms():
            do_full_init = True
        elif self.structure.get_number_of_atoms() != pf.pysic_interface.get_number_of_atoms():
            do_full_init = True
            
        if do_full_init:
            self.initialize_fortran_core()
        else:
            if not Pysic.core.cell_ready(self.structure):
                self.update_core_supercell()
            
            if not Pysic.core.atoms_ready(self.structure):
                self.update_core_coordinates()
            
            if not Pysic.core.charges_ready(self.structure):
                self.update_core_charges()
                    
            if not Pysic.core.potentials_ready(self.potentials):
                self.update_core_potentials()

            if self.coulomb != None:
                if not Pysic.core.coulomb_summation_ready(self.coulomb):
                    self.update_core_coulomb()
            
            if not Pysic.core.potential_lists_ready:
                self.update_core_potential_lists()

            if not self.neighbor_lists_waiting:
                self.create_neighbor_lists(self.get_individual_cutoffs(1.0))

            if not Pysic.core.neighbor_lists_ready(self.neighbor_list):
                self.update_core_neighbor_lists()
                

    def update_core_potential_lists(self):
        """Initializes the potential lists.

        Since one often runs :class:`~pysic.calculator.Pysic` with a set of potentials,
        the core pre-analyzes which potentials affect each atom and saves a list
        of such potentials for every particle. This method asks the core to
        generate these lists.
        """
        if not Pysic.core.atoms_ready(self.structure):
            raise MissingAtomsError("Creating potential lists before updating atoms in core.")
        pf.pysic_interface.create_potential_list()
        pf.pysic_interface.create_bond_order_factor_list()
        Pysic.core.potential_lists_ready = True


    def update_core_potentials(self):
        """Generates potentials for the Fortran core."""
        
        Pysic.core.potential_lists_ready = False
        if self.potentials == None:
            pf.pysic_interface.allocate_potentials(0)
            pf.pysic_interface.allocate_bond_order_factors(0)
            
            n_atoms = pf.pysic_interface.get_number_of_atoms()
            pf.pysic_interface.allocate_bond_order_storage(n_atoms,0,0)
            warn("There are no potentials associated with the Pysic calculator!",2)
            return

        if len(self.potentials) == 0:
            pf.pysic_interface.allocate_potentials(0)
            pf.pysic_interface.allocate_bond_order_factors(0)
            n_atoms = pf.pysic_interface.get_number_of_atoms()
            pf.pysic_interface.allocate_bond_order_storage(n_atoms,0,0)
            warn("There are no potentials associated with the Pysic calculator!",2)
            return
        
        n_pots = 0
        coord_list = []
        pot_index = 0
        # count the number of separate potentials
        for pot in self.potentials:

            # grab the coordinators associated with the potentials
            # but check that each coordinator is included only once
            # to prevent multiple evaluations of the same bond order factors
            coord = pot.get_coordinator()
            if(not coord is None):
                repeat_coord = False
                for crd in coord_list:
                    if crd[0] == coord:
                        repeat_coord = True
                if not repeat_coord:
                    coord_list.append([coord,pot_index])
                pot_index += 1
            
                # check the collection of bond order parameters for this potential
                bond_level = 0
                scaling_elements = []
                noted_scaling = False

                for bond in coord.get_bond_order_parameters():

                    # !!!: check for overriding scaling
                    if len(scaling_elements) == 0:
                        if bond.includes_scaling():
                            for elems in bond.get_symbols():
                                scaling_elements.append(elems[0])
                    else:
                        if bond.includes_scaling() and not noted_scaling:

                            # gather all targets of the new scaler
                            new_elems = []
                            for elems in bond.get_symbols():
                                new_elems.append(elems[0])

                            override = False

                            # if the target already exists, we are overriding scaling
                            for ne in new_elems:
                                for elems in scaling_elements:
                                    if ne == elems:
                                        override = True

                            # add the elements to the list of scaling targets
                            for ne in new_elems:
                                scaling_elements.append(ne)
                                
                            if override:
                                noted_scaling = True
                                warn("You are overriding bond order scaling in \n"+str(coord),3)
                                    
                    # !!!: check for redundancy (elements not in the potential)
                    for elems in bond.get_symbols():
                        for symbols in pot.get_symbols():
                            found_symbol = False
                            for symb in symbols:
                                if elems[0] == symb:
                                    found_symbol = True
                            if not found_symbol:
                                warn("You are applying a bond order factor on "+elems[0]+\
                                    " to a potential which does not target this element.\n\n"+\
                                    str(pot)+"\n\n"+str(bond),2)

                    # check for different levels of bond factors
                    if bond_level > 0:
                        if bond_level != bond.get_level():
                            warn("You are mixing per-atom and per-bond bond order factors for potential \n"+\
                                str(pot),2)
                    else:
                        bond_level = bond.get_level()
                
                    # check that pairwise factors are only applied on pair potentials
                    if bond.get_level() > 1:
                        if pot.get_number_of_targets() != bond.get_level():
                            warn("You are applying a level "+str(bond.get_level())+\
                                " bond order factor on a "+str(pot.get_number_of_targets())+\
                                "-body potential: it will be zero!",1)
            
            try:
                alltargets = pot.get_symbols()
                for targets in alltargets:
                    perms = permutations(targets)
                    different = set(perms)
                    n_pots += len(different)
            except:
                if not pot.get_symbols() is None:
                    raise InvalidPotentialError("Invalid potential symbols: "+str(pot.get_symbols()))
            try:
                alltargets = pot.get_tags()
                for targets in alltargets:
                    perms = permutations(targets)
                    different = set(perms)
                    n_pots += len(different)
            except:
                if not pot.get_tags() is None:
                    raise InvalidPotentialError("Invalid potential tags: "+str(pot.get_tags()))
            try:
                alltargets = pot.get_indices()
                for targets in alltargets:
                    perms = permutations(targets)
                    different = set(perms)
                    n_pots += len(different)
            except:
                if not pot.get_indices() is None:
                    raise InvalidPotentialError("Invalid potential indices: "+str(pot.get_indices()))
                
        pf.pysic_interface.allocate_potentials(n_pots)

        elemental_potentials = []
        is_multiplier = []
        master_potentials = []
        for pots in self.potentials:
        
            # warn for missing cutoffs
            if pots.get_number_of_targets() > 1 and pots.get_cutoff() < 0.01:
                warn("Potential with zero cutoff present: \n\n"+str(pots),1)
            if pots.get_symbols() is None and \
                pots.get_tags() is None and \
                pots.get_indices() is None:
                warn("Potential with no targets present: \n\n"+str(pots),1)
        
            # reverse the lists so that the master potentials come last
            for rpot in reversed(pots.get_potentials()):
                elemental_potentials.append(rpot)
            for rmul in reversed(pots.is_multiplier()):
                is_multiplier.append(rmul)
            master_potentials.extend([pots.get_potentials()[0]]*\
                                     len(pots.get_potentials()))
                    
        pot_index = 0
        for pot, mul, mpot in zip(elemental_potentials,is_multiplier, master_potentials):
    
            multiplier_added = False
            group_index = -1
            if not mpot.get_coordinator() is None:
                # pick the potential index of the first coordinator that matches
                # that of the master potential - this is done to prevent
                # creation of multiple equal bond factors that would then be
                # evaluated repeatedly
                for crd in coord_list:
                    # only check the bond parameters since the group indices are
                    # manipulated here so the coordinators don't exactly match
                    if mpot.get_coordinator().get_bond_order_parameters() == crd[0].get_bond_order_parameters() and group_index < 0:
                        group_index = crd[1]
                #group_index = pot_index # the old method assigns a new coordinator to each potential
                mpot.get_coordinator().set_group_index(group_index)
            if not mul:
                pot_index += 1

            n_targ = mpot.get_number_of_targets()
            no_symbs = np.array( n_targ*[pu.str2ints('xx',2)] ).transpose()
            no_tags = np.array( n_targ*[-9] )
            no_inds = np.array( n_targ*[-9] )

            try:
                        
                if mul:
                    alltargets = [mpot.get_symbols()[0]]
                else:
                    alltargets = mpot.get_symbols()
                for targets in alltargets:
                    int_orig_symbs = []
                    for orig_symbs in targets:
                        int_orig_symbs.append( pu.str2ints(orig_symbs,2) )

                    if mul:
                        different = [targets]
                    else:
                        perms = permutations(targets)
                        different = set(perms)

                    for symbs in different:
                        int_symbs = []
                        for label in symbs:
                            int_symbs.append( pu.str2ints(label,2) )

                        if not mul or not multiplier_added:
                            success = pf.pysic_interface.add_potential(pot.get_potential_type(),
                                                         np.array( pot.get_parameter_values() ),
                                                         mpot.get_cutoff(),
                                                         mpot.get_soft_cutoff(),
                                                         np.array( int_symbs ).transpose(),
                                                         no_tags,
                                                         no_inds,
                                                         np.array( int_orig_symbs ).transpose(),
                                                         no_tags,
                                                         no_inds,
                                                         group_index,
                                                         mul )
                            multiplier_added = True
                        else:
                            success = True

                        if not success:
                            raise InvalidPotentialError("")
            except:
                if not mpot.get_symbols() is None:
                    raise InvalidPotentialError("Failed to create a potential in the core: "+str(mpot))
            try:
                if mul:
                    alltargets = [mpot.get_tags()[0]]
                else:
                    alltargets = mpot.get_tags()
                for targets in alltargets:
                    orig_tags = targets

                    if mul:
                        different = [targets]
                    else:
                        perms = permutations(targets)
                        different = set(perms)

                    for tags in different:
                                                
                        if not mul or not multiplier_added:
                        
                            success = pf.pysic_interface.add_potential(pot.get_potential_type(),
                                                         np.array( pot.get_parameter_values() ),
                                                         mpot.get_cutoff(),
                                                         mpot.get_soft_cutoff(),
                                                         no_symbs,
                                                         np.array( tags ),
                                                         no_inds,
                                                         no_symbs,
                                                         np.array(orig_tags),
                                                         no_inds,
                                                         group_index,
                                                         mul)
                        
                            multiplier_added = True
                        else:
                            success = True
                        
                        if not success:
                            raise InvalidPotentialError("")
            except:
                if not mpot.get_tags() is None:
                    raise InvalidPotentialError("Failed to create a potential in the core: "+str(mpot))
            try:
                if mul:
                    alltargets = [mpot.get_indices()[0]]
                else:
                    alltargets = mpot.get_indices()                
                for targets in alltargets:
                    orig_inds = targets
                        
                    if mul:
                        different = [targets]
                    else:
                        perms = permutations(targets)
                        different = set(perms)

                    for inds in different:
                                                
                        if not mul or not multiplier_added:
                        
                            success = pf.pysic_interface.add_potential(pot.get_potential_type(),
                                                         np.array( pot.get_parameter_values() ),
                                                         mpot.get_cutoff(),
                                                         mpot.get_soft_cutoff(),
                                                         no_symbs,
                                                         no_tags,
                                                         np.array( inds ),
                                                         no_symbs,
                                                         no_tags,
                                                         np.array(orig_inds),
                                                         group_index,
                                                         mul )
                            multiplier_added = True
                        else:
                            success = True


                        if not success:
                            raise InvalidPotentialError("")
            except:
                if not mpot.get_indices() is None:
                    raise InvalidPotentialError("Failed to create a potential in the core: "+str(mpot))
                        
                        
            if not mul:
                pf.pysic_interface.clear_potential_multipliers()

        n_bonds = 0
        permutate = False
        for coord in coord_list:
            try:
                allbonds = coord[0].get_bond_order_parameters()
                for bond in allbonds:
                
                    # warn for missing cutoffs
                    if (bond.get_number_of_targets() > 1 and bond.get_cutoff() < 0.01):
                        warn("Bond order factor with zero cutoff present: \n\n"+str(bond),1)
                
                    alltargets = bond.get_symbols()
                    for targets in alltargets:
                    
                        if(permutate):
                            # permutate bond factor symbols
                            perms = permutations(targets)
                            different = set(perms)
                            n_bonds += len(different)

                        else:
                            # do not permutate the bond factor symbols
                            n_bonds += 1
            except:
                raise InvalidParametersError("Invalid bond order parameter symbols: "+str(bond.get_symbols()))

        pf.pysic_interface.allocate_bond_order_factors(n_bonds)

        for coord in coord_list:
            try:
                allbonds = coord[0].get_bond_order_parameters()
                for bond in allbonds:
                    alltargets = bond.get_symbols()
                    for targets in alltargets:

                        int_orig_symbs = []
                        for orig_symbs in targets:
                            int_orig_symbs.append( pu.str2ints(orig_symbs,2) )
                    
                        if(permutate):
                            # permutate bond factor symbols
                            perms = permutations(targets)
                            different = set(perms)
                        else:
                            # do not permutate the bond factor symbols
                            different = [targets]

                        for symbs in different:
                            int_symbs = []
                            for label in symbs:
                                int_symbs.append( pu.str2ints(label,2) )

                            success = pf.pysic_interface.add_bond_order_factor(bond.get_bond_order_type(),
                                                                   np.array( bond.get_parameters_as_list() ),
                                                                   np.array( bond.get_number_of_parameters() ),
                                                                   bond.get_cutoff(),
                                                                   bond.get_soft_cutoff(),
                                                                   np.array( int_symbs ).transpose(),
                                                                   np.array( int_orig_symbs ).transpose(),
                                                                   coord[1])
                            if not success:
                                raise InvalidParametersError("")


            except:
                raise InvalidParametersError("Failed to create a bond order factor in the core: "+str(bond))

        n_atoms = pf.pysic_interface.get_number_of_atoms()
        pf.pysic_interface.allocate_bond_order_storage(n_atoms,
                                                       pot_index,
                                                       len(coord_list))

        Pysic.core.set_potentials(self.potentials)


        self.neighbor_lists_waiting = False


            
    def update_core_coulomb(self):
        """Updates the Coulomb summation parameters in the Fortran core.
            """
        
        if self.coulomb != None:
            if self.coulomb.method == CoulombSummation.summation_modes[0]: # ewald summation
                rcut = self.coulomb.parameters['real_cutoff']
                kcut = self.coulomb.parameters['k_cutoff']
                sigma = self.coulomb.parameters['sigma']
                epsilon = self.coulomb.parameters['epsilon']
                
                scales = self.coulomb.get_scaling_factors()
                
                # calculate the truncation limits for the k-space sum
                reci_cell = np.multiply(2.0*math.pi,self.structure.get_reciprocal_cell())
                volume = np.dot( reci_cell[0], np.cross( reci_cell[1], reci_cell[2] ) )
                k1 = int( kcut * np.linalg.norm( np.cross( reci_cell[1], reci_cell[2] ) ) / volume + 0.5 )
                k2 = int( kcut * np.linalg.norm( np.cross( reci_cell[0], reci_cell[2] ) ) / volume + 0.5 )
                k3 = int( kcut * np.linalg.norm( np.cross( reci_cell[0], reci_cell[1] ) ) / volume + 0.5 )

                if scales == None:
                    scales = [1.0]*self.structure.get_number_of_atoms()
                elif(len(scales) != self.structure.get_number_of_atoms()):
                    raise InvalidParametersError("Length of the scaling factor vector does not match the number of atoms.")
                
                pf.pysic_interface.set_ewald_parameters(rcut,
                                                        kcut,
                                                        np.array([k1,k2,k3]),
                                                        sigma,
                                                        epsilon,
                                                        scales)

                Pysic.core.set_coulomb(self.coulomb)
        
    
    def update_core_coordinates(self):
        """Updates the positions and momenta of atoms in the Fortran core.

        The core must be initialized and the number of atoms must match.
        Upon the update, it is automatically checked if the neighbor lists
        should be updated as well.
        """
        
        if self.structure.get_number_of_atoms() != pf.pysic_interface.get_number_of_atoms():
            raise LockedCoreError("The number of atoms does not match.")
        
        positions = np.array( self.structure.get_positions() ).transpose()
        momenta = np.array( self.structure.get_momenta() ).transpose()

        self.forces = None
        self.energy = None
        self.stress = None
        self.electronegativities = None

        pf.pysic_interface.update_atom_coordinates(positions,momenta)

        Pysic.core.set_atomic_positions(self.structure)
        Pysic.core.set_atomic_momenta(self.structure)
                
        if not self.neighbor_lists_waiting:
            self.create_neighbor_lists(self.get_individual_cutoffs(1.0))
        
        self.update_core_neighbor_lists()


    def get_charges(self, system=None):
        """Update for ASE 3.7"""
    
        if self.charges is not None:
            return self.charges
        else:
            if system is not None:
                try:
                    return system.get_initial_charges()
                except:
                    return system.get_charges()
            else:
                try:
                    return self.structure.get_initial_charges()
                except:
                    return self.structure.get_charges()

    def update_core_charges(self):
        """Updates atomic charges in the core."""
        
        try:
            charges = np.array( self.structure.get_initial_charges() )
        except:
            charges = np.array( self.structure.get_charges() )

        self.forces = None
        self.energy = None
        self.stress = None
        self.electronegativities = None
        
        total_charge = np.sum(charges)
        if total_charge > 0.001 or total_charge < -0.001:
            warn("There is non-zero total charge in the system: {tc}".format(tc=total_charge),5)
        
        pf.pysic_interface.update_atom_charges(charges)
        self.charges = charges
        
        Pysic.core.set_charges(charges)
            
            
    def update_core_supercell(self):
        """Updates the supercell in the Fortran core."""
        vectors = np.array( self.structure.get_cell() ).transpose()
        inverse = np.linalg.inv(np.array( self.structure.get_cell() )).transpose()
        periodicity = np.array( self.structure.get_pbc() )
        
        pf.pysic_interface.create_cell(vectors,inverse,periodicity)
        
        Pysic.core.set_cell(self.structure)
        Pysic.core.set_neighbor_lists(None)
            

    def update_core_neighbor_lists(self):
        """Updates the neighbor lists in the Fortran core.

         If uninitialized, the lists are created first via :meth:`~pysic.calculator.Pysic.create_neighbor_lists`.
         """

        if not Pysic.core.atoms_ready(self.structure):
            raise MissingAtomsError("Creating neighbor lists before updating atoms in the core.")
        cutoffs = self.get_individual_cutoffs(1.0)

        if not self.neighbor_lists_waiting:
            self.create_neighbor_lists(cutoffs)
            self.set_cutoffs(cutoffs)
            self.neighbor_lists_waiting = True

        self.neighbor_list.update(self.structure)
        if isinstance(self.neighbor_list,FastNeighborList):
            # if we used the fast list, the core is already updated
            pass
        else:
            # if we have used the ASE list, it must be passed on to the core
            for index in range(self.structure.get_number_of_atoms()):
                [nbors,offs] = self.neighbor_list.get_neighbors(index)                
                pf.pysic_interface.create_neighbor_list(index+1,np.array(nbors),np.array(offs).transpose())

        Pysic.core.set_neighbor_lists(self.neighbor_list)


    def initialize_fortran_core(self):
        """Fully initializes the Fortran core, creating the atoms, supercell, potentials, and neighbor lists."""
        
        
        masses = np.array( self.structure.get_masses() )
        
        # the call for charges was changed between ASE 3.6 and 3.7
        try:
            self.charges = np.array( self.structure.get_initial_charges() )
        except:
            self.charges = np.array( self.structure.get_charges() )
        
        charges = self.charges
        positions = np.array( self.structure.get_positions() ).transpose()
        momenta = np.array( self.structure.get_momenta() ).transpose()
        tags = np.array( self.structure.get_tags() )
        elements = self.structure.get_chemical_symbols()

        for index in range(len(elements)):
            elements[index] = pu.str2ints(elements[index],2)

        elements = np.array( elements ).transpose()

        #self.create_neighbor_lists(self.get_individual_cutoffs(1.0))
        #self.neighbor_lists_waiting = True

        total_charge = np.sum(self.charges)
        if total_charge > 0.001 or total_charge < -0.001:
            warn("There is non-zero total charge in the system: {tc}".format(tc=total_charge),5)


        pf.pysic_interface.create_atoms(masses,charges,positions,momenta,tags,elements)
        Pysic.core.set_atoms(self.structure)

        pf.pysic_interface.distribute_mpi(self.structure.get_number_of_atoms())
        Pysic.core.mpi_ready = True
                
        self.update_core_supercell()
        self.update_core_potentials()
        self.neighbor_lists_waiting = False
        self.update_core_neighbor_lists()
        self.update_core_potential_lists()
        self.update_core_coulomb()

    def get_numerical_energy_gradient(self, atom_index, shift=0.0001, atoms=None):
        """Numerically calculates the negative gradient of energy with respect to moving a single particle.

        This is for debugging the forces."""

        if(atoms == None):
            system = self.structure
            orig_system = self.structure.copy()
        else:
            system = atoms.copy()
            orig_system = atoms.copy()
            self.set_atoms(system)
        
        self.energy == None
        energy_xp = self.get_potential_energy()
        system[atom_index].x += shift
        energy_xp = self.get_potential_energy()
        system[atom_index].x -= 2.0*shift
        energy_xm = self.get_potential_energy()
        system[atom_index].x += shift

        system[atom_index].y += shift
        energy_yp = self.get_potential_energy()
        system[atom_index].y -= 2.0*shift
        energy_ym = self.get_potential_energy()
        system[atom_index].y += shift

        system[atom_index].z += shift
        energy_zp = self.get_potential_energy()
        system[atom_index].z -= 2.0*shift
        energy_zm = self.get_potential_energy()
        system[atom_index].z += shift

        self.energy == None
        self.get_potential_energy(orig_system)
                
        return [ -(energy_xp-energy_xm)/(2.0*shift),
                 -(energy_yp-energy_ym)/(2.0*shift),
                 -(energy_zp-energy_zm)/(2.0*shift) ]


            
    def set_cutoffs(self, cutoffs):
        """Copy and save the list of individual cutoff radii.
            
            Parameters:
            
            cutoffs: list of doubles
            new cutoffs
            """
        self.saved_cutoffs = copy.deepcopy(cutoffs)
            
            
    def neighbor_lists_expanded(self, cutoffs):
        """Check if the cutoffs have been expanded.
                    
        If the cutoffs have been made longer than before,
        the neighbor lists have to be recalculated.
        This method checks the individual cutoffs of all atoms
        to check if the cutoffs have changed.
        
        Parameters:
        
        cutoffs: list of doubles
            new cutoffs
        """
        if self.saved_cutoffs == None:
            return True
        if cutoffs == None:
            return True
                                
        if len(self.saved_cutoffs) != len(cutoffs):
            return True
        for old_cut, new_cut in zip(self.saved_cutoffs, cutoffs):
            if old_cut < new_cut:
                return True
                
        return False

            


    def get_numerical_bond_order_gradient(self, coordinator, atom_index, moved_index, shift=0.001, atoms=None):
        """Numerically calculates the gradient of a bond order factor with respect to moving a single particle.

        This is for debugging the bond orders."""

        if(atoms == None):
            system = self.structure.copy()
            orig_system = self.structure.copy()
        else:
            system = atoms.copy()
            orig_system = atoms.copy()

        self.energy == None
        crd = coordinator
        system[moved_index].x += shift
        self.set_atoms(system)
        self.set_core()
        crd.calculate_bond_order_factors()
        bond_xp = crd.get_bond_order_factors()[atom_index]
        system[moved_index].x -= 2.0*shift
        self.set_atoms(system)
        self.set_core()
        crd.calculate_bond_order_factors()
        bond_xm = crd.get_bond_order_factors()[atom_index]
        system[moved_index].x += shift        

        system[moved_index].y += shift
        self.set_atoms(system)
        self.set_core()
        crd.calculate_bond_order_factors()
        bond_yp = crd.get_bond_order_factors()[atom_index]
        system[moved_index].y -= 2.0*shift
        self.set_atoms(system)
        self.set_core()
        crd.calculate_bond_order_factors()
        bond_ym = crd.get_bond_order_factors()[atom_index]
        system[moved_index].y += shift

        system[moved_index].z += shift
        self.set_atoms(system)
        self.set_core()
        crd.calculate_bond_order_factors()
        bond_zp = crd.get_bond_order_factors()[atom_index]
        system[moved_index].z -= 2.0*shift
        self.set_atoms(system)
        self.set_core()
        crd.calculate_bond_order_factors()
        bond_zm = crd.get_bond_order_factors()[atom_index]
        system[moved_index].z += shift

        self.energy == None
        self.set_atoms(orig_system)
        self.set_core()

        
        return [ (bond_xp-bond_xm)/(2.0*shift),
                 (bond_yp-bond_ym)/(2.0*shift),
                 (bond_zp-bond_zm)/(2.0*shift) ]



    
    def get_numerical_electronegativity(self, atom_index, shift=0.001, atoms=None):
        """Numerically calculates the derivative of energy with respect to charging a single particle.
            
            This is for debugging the electronegativities."""
        
        if(atoms == None):
            system = self.structure.copy()
            orig_system = self.structure.copy()
        else:
            system = atoms.copy()
            orig_system = self.structure.copy()
        
        try:
            charges = system.get_charges()
        except:
            charges = system.get_initial_charges()
        
        self.energy == None
        self.set_atoms(system)
        self.set_core()
        charges[atom_index] += 1.0*shift

        # the call for charges was changed between ASE 3.6 and 3.7
        try:
            system.set_charges(charges)
        except:
            system.set_initial_charges(charges)

        energy_p = self.get_potential_energy(system)
        charges[atom_index] -= 2.0*shift
        
        # the call for charges was changed between ASE 3.6 and 3.7
        try:
            system.set_charges(charges)
        except:
            system.set_initial_charges(charges)

        energy_m = self.get_potential_energy(system)
        charges[atom_index] += 1.0*shift
        
        # the call for charges was changed between ASE 3.6 and 3.7
        try:
            system.set_charges(charges)
        except:
            system.set_initial_charges(charges)
                
        self.energy == None
        self.set_atoms(orig_system)
        self.set_core()
        
        return (energy_m-energy_p)/(2.0*shift)






