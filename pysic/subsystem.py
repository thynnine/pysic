#! /usr/bin/env python
"""Defines classes used for creating and storing information about subsystems
in a HybridCalculator."""

#==============================================================================
class SubSystemInfo(object):

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

    def __init__(self, name, indices=None, tag=None, special_set=None, calculator=None):
        """@todo: to be defined1. """
        self.name = name
        self.calculator = calculator
        self.set_atoms(indices, tag, special_set)

    def set_calculator(self, calculator):
        """@todo: Docstring for set_calculator.

        :calculator: @todo
        :returns: @todo

        """
        self.calculator = calculator

    def set_atoms(self, indices=None, tag=None, special_set=None):
        """@todo: Docstring for set_calculator.

        :calculator: @todo
        :returns: @todo

        """
        if indices is not None and type(indices) is int:
                self.indices = (indices,)
        else:
            self.indices = indices
        self.tag = tag
        self.special_set = special_set

#==============================================================================
class SubSystem(object):

    """Represents a subsystem used in HybridCalculator.

    This class is materialised from a SubSystemInfo, and should not be
    accessible to the end user.
    """
    def __init__(self, atoms=None, calculator=None, index_map=None, reverse_index_map=None):
        """@todo: to be defined1. """

        self.atoms_for_binding = atoms.copy()
        self.atoms_for_subsystem = atoms.copy()
        self.calculator = calculator
        self.index_map = index_map
        self.reverse_index_map = reverse_index_map
        self.potential_energy = None
        self.forces = None
        self.link_interaction_correction = 0
        self.density_grid = None
        self.initial_charges = atoms.get_initial_charges()

    def get_potential_energy_without_corrections(self):
        """@todo: Docstring for get_potential_energy.
        :returns: @todo
        """
        # Ask the energy from the modified atoms (which include possible link
        # atoms), and and possible corrections
        self.potential_energy = self.calculator.get_potential_energy(
            self.atoms_for_subsystem)
        return self.potential_energy

    def get_forces(self):
        """@todo: Docstring for get_potential_energy.
        :returns: @todo
        """
        # Calculate the forces, ignore forces on link atoms
        forces = self.calculator.get_forces(self.atoms_for_subsystem)

        force_map = []
        for full_index, sub_index in self.index_map.iteritems():
            force = forces[sub_index, :]
            index_force_pair = (full_index, force)
            force_map.append(index_force_pair)

        self.forces = force_map
        return force_map




