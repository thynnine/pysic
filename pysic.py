#! /usr/bin/env python
"""Pysic calculator and atomistic potential."""
import pysic_fortran as pf
import pysic_utility as pu
import numpy as np
import ase.calculators.neighborlist as nbl
from itertools import permutations
import random as rnd

#
# PySiC : Pythonic Simulation Code
#
# Teemu Hynninen 2011
#
# pysic.py is the main interface module of Pysic.
# The module contains methods allowing the user
# to interface with the rest of the simulation tools.
#
# The module defines potentials and a calculator class
# ready to be used with the ASE simulation package.
# 


class InvalidPotentialError(Exception):
    """An error raised when an invalid potential is about to be created or used.

    Parameters:

    message: string
        information describing why the error occurred
    potential: :class:`~pysic.Potential`
        the errorneous potential
    """
    def __init__(self,message='',potential=None):
        self.message = message
        self.potential = potential

    def __str__(self):
        if(self.potential == None):
            return self.message
        else:
            return self.messsage + " - the Potential: " + repr(potential)


class MissingAtomsError(Exception):
    """An error raised when the core is being updated with per atom information before updating the atoms.

    Parameters:

    message: string
        information describing why the error occurred
    """
    def __init__(self,message=''):
        self.message = message

    def __str__(self):
        return self.message

    
class MissingNeighborsError(Exception):
    """An error raised when a calculation is initiated without initializing the neighbor lists.

    In principle :class:`~pysic.Pysic` should always take care of handling the neighbors automatically.
    This error is an indication that there is loophole in the built-in preparations.

    Parameters:

    message: string
        information describing why the error occurred
    """
    def __init__(self,message=''):
        self.message = message

    def __str__(self):
        return self.message

class LockedCoreError(Exception):
    """An error raised when a :class:`~pysic.Pysic` tries to access the core which is locked
    by another calculator.

    Parameters:

    message: string
        information describing why the error occurred
    """
    def __init__(self,message=''):
        self.message = message

    def __str__(self):
        return self.message

# Below are general methods for accessing information about the built-in potentials.
# Note that the methods directly inquire the Fortran core, and thus all changes made
# to the core should automatically be taken into account, as long as the proper
# arrays are updated in Potentials.f90

# automatically call the initialization routine in Potentials.f90
pf.pysic_interface.start_potentials()

# automatically initialize mpi - fortran side decides if we really are in mpi
pf.pysic_interface.start_mpi()

# automatically initialize the fortran rng
# it needs to be possible to force a seed by the user
pf.pysic_interface.start_rng(5301)#rnd.randint(1, 999999))


def finish_mpi():
    """Terminates the MPI framework.

    If the Fortran core is compiled in MPI mode, :mod:`~pysic` will automatically
    initialize MPI upon being imported. As the environment won't know
    when the user is done with his simulation, terminating the MPI is
    left as a manual operation. If one terminates Python running :mod:`~pysic`
    without finishing the MPI first, the MPI framework in the Fortran core will
    print an error.
    """
    pf.pysic_interface.finish_mpi()


def sync_mpi():
    """Calls MPI barrier from the Fortran core MPI framework."""

    pf.pysic_interface.sync_mpi()

def get_number_of_cpus():
    """Gets the number of cpus from the Fortran MPI.
    """
    return pf.pysic_interface.get_number_of_cpus()


def get_cpu_id():
    """Gets the cpu ID from the Fortran MPI.
    """
    return pf.pysic_interface.get_cpu_id()


#def distribute_mpi_load(pysic_calculator):
#    """Distributes atoms among processors in the Fortran MPI core.
#
#    Parameters:
#
#    pysic_calculator: :class:`pysic.Pysic~` object
#        the calculator object whose structure is used for distributing the MPI load.
#    """
#    pf.pysic_interface.distribute_mpi(pysic_calculator.get_atoms().get_number_of_atoms())


def list_valid_potentials():
    """A list of names of potentials currently known by the core.

    The method retrieves from the core a list of the names of different potentials
    currently implemented. Since the fortran core is directly accessed, any
    updates made in the core source code should get noticed automatically.
    """
    n_pots = pf.pysic_interface.number_of_potentials()
    pot_codes = pf.pysic_interface.list_valid_potentials(n_pots).transpose()
    pot_names = []
    for code in pot_codes:
        pot_names.append(pu.ints2str(code))

    return pot_names

# Checks if the given keyword defines a potential
def is_valid_potential(potential_name):
    """Tells if the given string is the name of a potential.

    Parameters:

    potential_name: string
        the name of the potential
    """
    return pf.pysic_interface.is_potential(potential_name)

def number_of_targets(potential_name):
    """Tells how many targets a potential acts on, i.e., is it pair or many-body.

    Parameters:

    potential_name: string
        the name of the potential
    """
    return pf.pysic_interface.number_of_targets_of_potential(potential_name)

def number_of_parameters(potential_name):
    """Tells how many parameters a potential incorporates.

    Parameters:

    potential_name: string
        the name of the potential
    """
    return pf.pysic_interface.number_of_parameters_of_potential(potential_name)

def names_of_parameters(potential_name):
    """Lists the names of the parameters of a potential.

    Parameters:

    potential_name: string
        the name of the potential
    """
    param_codes = pf.pysic_interface.names_of_parameters_of_potential(potential_name).transpose()
    param_names = []
    n_params = number_of_parameters(potential_name)
    index = 0
    for code in param_codes:
        if index < n_params:
            param_names.append(pu.ints2str(code).strip())
            index += 1

    return param_names

def index_of_parameter(potential_name, parameter_name):
    """Tells the index of a parameter of a potential in the list of parameters the potential uses.

    Parameters:

    potential_name: string
        the name of the potential
    parameter_name: string
        the name of the parameter
    """
    param_names = names_of_parameters(potential_name)
    index = 0
    for name in param_names:
        if name == parameter_name:
            return index
        index += 1


def descriptions_of_parameters(potential_name):
    """Returns a list of strings containing physical names of the parameters of a potential,
    e.g., 'spring constant' or 'decay length'.
    
    Parameters:

    potential_name: string
        the name of the potential
    """
    param_codes = pf.pysic_interface.descriptions_of_parameters_of_potential(potential_name).transpose()
    param_notes = []
    n_params = number_of_parameters(potential_name)
    index = 0
    for code in param_codes:
        if index < n_params:
            param_notes.append(pu.ints2str(code).strip())
            index += 1

    return param_notes

def description_of_potential(potential_name, parameter_values=None, cutoff=None,
                             elements=None, tags=None, indices=None):
    """Prints a brief description of a potential. If optional arguments are provided,
    they are incorporated in the description. That is, by default the method describes
    the general features of a potential, but it can also be used for describing a
    particular potential with set parameters.

    Parameters:

    potential_name: string
        the name of the potential
    """
    param_names = names_of_parameters(potential_name)
    param_notes = descriptions_of_parameters(potential_name)
    n_targets = number_of_targets(potential_name)
    n_params = number_of_parameters(potential_name)
    description = pf.pysic_interface.description_of_potential(potential_name)

    message = """
potential '{pot}'
    
{n_targ}-body interaction

{descr}
parameters ({n_par}):
""".format(pot=potential_name,n_targ=n_targets,descr=description,n_par=str(n_params))

    if parameter_values == None:
        for para, note in zip(param_names, param_notes):
            message += "{pa} : {no}\n".format(pa=para,no=note)
    else:
        for para, note, val in zip(param_names, param_notes, parameter_values):
            message += "{pa} : {no} = {va}\n".format(pa=para,no=note,va=val)

    if cutoff != None:
        message += "\ncutoff = {cut}\n".format(cut=str(cutoff))

    if elements != None:
        message += "\ntargeted symbols:"
        for ele in elements:            
            message += " {0} ".format(str(ele))
        message += "\n"

    if tags != None:
        message += "\ntargeted tags:"
        for tag in tags:
            message += " {0} ".format(str(tag))
        message += "\n"

    if indices != None:
        message += "\ntargeted indices:"
        for indy in indices:
            message += " {0} ".format(str(indy))
        message += "\n"
            
    print message



class Potential:
    """Class for representing a potential.

    Several types of potentials can be defined by specifying the type of the potential
    as a keyword. The potentials contain a host of parameters and information on
    what types of particles they act on. To view a list of available potentials, use the method
    :meth:`~pysic.list_valid_potentials`.

    A potential may be a pair or many-body potential: here, the bodies a potential
    acts on are called targets. Thus specifying the number of targets of a potential also
    determines if the potential is a many-body potential.
    
    A potential may be defined for atom types or specifically for certain atoms. These are
    specified by the symbols, tags, and indices. Each of these should be either 'None' or
    a list of lists of n values where n is the number of targets. For example, if::

        indices = [[0, 1], [1, 2], [2, 3]]

    the potential will be applied between atoms 0 and 1, 1 and 2, and 2 and 3.

    Parameters:

    symbols: list of string
        the chemical symbols (elements) on which the potential acts
    tags: integer
        atoms with specific tags on which the potential acts
    indices: list of integers
        atoms with specific indices on which the potential acts
    potential_type: string
        a keyword specifying the type of the potential
    parameters: list of doubles
        a list of parameters for characterizing the potential; their meaning depends on the type of potential
    cutoff: double
        the maximum atomic separation at which the potential is applied
    """

    def __init__(self,potential_type,symbols=None,tags=None,indices=None,parameters=None,cutoff=0.0,cutoff_margin=0.0):
        if(is_valid_potential(potential_type)):
            self.symbols = None
            self.tags = None
            self.indices = None
            self.potential_type = potential_type
            self.cutoff = cutoff
            self.cutoff_margin = 0.0
            self.set_cutoff_margin(cutoff_margin)
            self.n_targets = number_of_targets(potential_type)
            self.names_of_params = names_of_parameters(potential_type)
            self.set_symbols(symbols)
            self.set_tags(tags)
            self.set_indices(indices)
            if(parameters == None):
                self.parameters = len(self.names_of_params)*[0.0]
            else:
                self.parameters = parameters
                
            if(len(self.names_of_params) != len(self.parameters)):
                raise InvalidPotentialError(
                    'The potential "{pot}" requires {num} parameters: {par}'.format(
                    pot=potential_type,
                    num=str(number_of_parameters(potential_type)),
                    par=str(parameters(potential_type))
                    ) )
        else:
            raise InvalidPotentialError('There is no potential called "{pot}".'.format(pot=potential_type))


    def __eq__(self,other):
        try:
            if other == None:
                return False
            if self.symbols != other.symbols:
                return False
            if self.tags != other.tags:
                return False
            if self.indices != other.indices:
                return False
            if self.potential_type != other.potential_type:
                return False
            if self.parameters != other.parameters:
                return False
            if self.cutoff != other.cutoff:
                return False            
        except:
            return False
        
        return True


    def __repr__(self):
        return "Potential({name},symbols={symbs},tags={tags},indices={inds},parameters={params},cutoff={cut},cutoff_margin={marg})".format(name=self.potential_type,
                                                                                                                      symbs=str(self.symbols),
                                                                                                                      tags=str(self.tags),
                                                                                                                      inds=str(self.indices),
                                                                                                                      params=str(self.parameters),
                                                                                                                      cut=str(self.cutoff),
                                                                                                                      marg=str(self.cutoff_margin))
            

    def get_symbols(self):
        """Return a list of the chemical symbols (elements) on which the potential
        acts on."""
        return self.symbols

    def get_tags(self):
        """Return the tags on which the potential acts on."""
        return self.tags

    def get_indices(self):
        """Return a list of indices on which the potential acts on. """
        return self.indices

    def get_different_symbols(self):
        """Returns a list containing each symbol the potential affects once."""
        all_symbols = []
        if self.symbols == None:
            return all_symbols

        for symblist in self.symbols:
            for symb in symblist:
                if all_symbols.count(symb) == 0:
                    all_symbols.append(symb)
        return all_symbols

    def get_different_tags(self):
        """Returns a list containing each tag the potential affects once."""
        all_tags = []
        if self.tags == None:
            return all_tags
            
        for taglist in self.tags:
            for tag in taglist:
                if all_tags.count(tag) == 0:
                    all_tags.append(tag)
        return all_tags

    def get_different_indices(self):
        """Returns a list containing each index the potential affects once."""
        all_indices = []
        if self.indices == None:
            return all_indices

        for indicelist in self.indices:
            for index in indicelist:
                if all_indices.count(index) == 0:
                    all_indices.append(index)
        return all_indices
    
    def get_parameter_values(self):
        """Returns a list containing the current parameter values of the potential."""
        return self.parameters

    def get_parameter_names(self):
        """Returns a list of the names of the parameters of the potential."""
        return self.names_of_params

    def get_parameter_value(self,param_name):
        """Returns the value of the given parameter.

        Parameters:

        param_name: string
            name of the parameter
        """
        return self.parameters[index_of_parameter(self.potential_type,param_name)]

    def get_potential_type(self):
        """Returns the keyword specifying the type of the potential."""
        return self.potential_type

    def get_cutoff(self):
        """Returns the cutoff."""
        return self.cutoff

    def get_cutoff_margin(self):
        """Returns the margin for a smooth cutoff."""
        return self.cutoff_margin

    def get_soft_cutoff(self):
        """Returns the lower limit for a smooth cutoff."""
        return self.cutoff-self.cutoff_margin

    def get_number_of_targets(self):
        """Returns the number of targets."""
        return self.n_targets

    def accepts_target_list(self,targets):
        """Tests whether a list is suitable as a list of targets, i.e., symbols, tags, or indices
        and returns True or False accordingly.

        A list of targets should be of the format::

            targets = [[a, b], [c, d]]

        where the length of the sublists must equal the number of targets.

        It is not tested that the values contained in the list are valid.

        Parameters:

        targets: list of strings or integers
            a list whose format is checked
        """
        n_targets = self.get_number_of_targets()
        try:
            for a_set in targets:
                assert isinstance(a_set,list)
                if len(a_set) != n_targets:
                    return False
            return True
        except:
            return False

    def set_symbols(self,symbols):
        """Sets the list of symbols to equal the given list.

        Parameters:

        symbols: list of strings
            list of symbols on which the potential acts
        """
        if symbols == None:
            return
        if(self.accepts_target_list(symbols)):
            self.symbols = symbols
        elif(self.accepts_target_list([symbols])):
            self.symbols = [symbols]
        else:
            raise InvalidPotentialError("Invalid number of targets.")

    def set_tags(self,tags):
        """Sets the list of tags to equal the given list.

        Parameters:

        tags: list of integers
            list of tags on which the potential acts
        """
        if tags == None:
            return
        if(self.accepts_target_list(tags)):
            self.tags = tags
        elif(self.accepts_target_list([tags])):
            self.tags = [tags]
        else:
            raise InvalidPotentialError("Invalid number of targets.")

    def set_indices(self,indices):
        """Sets the list of indices to equal the given list.

        Parameters:

        indices: list of integers
            list of integers on which the potential acts
        """
        if indices == None:
            return
        if(self.accepts_target_list(indices)):
            self.indices = indices
        elif(self.accepts_target_list([indices])):
            self.indices = [indices]
        else:
            raise InvalidPotentialError("Invalid number of targets.")

    def add_symbols(self,symbols):
        """Adds the given symbols to the list of symbols.

        Parameters:

        symbols: list of strings
            list of additional symbols on which the potential acts
        """
        if(self.accepts_target_list(symbols)):
            if self.symbols == None:
                self.symbols = []

            for stuff in symbols:
                self.symbols.append(stuff)
        elif(self.accepts_target_list([symbols])):
            if self.symbols == None:
                self.symbols = []

            self.symbols.append(symbols)
        else:
            raise InvalidPotentialError("Invalid number of targets.")
            

    def add_tags(self,tags):
        """Adds the given tags to the list of tags.

        Parameters:

        tags: list of integers
            list of additional tags on which the potential acts
        """
        if(self.accepts_target_list(tags)):
            if self.tags == None:
                self.tags = []

            for stuff in tags:
                self.tags.append(stuff)
        elif(self.accepts_target_list([tags])):
            if self.tags == None:
                self.tags = []

            self.tags.append(tags)
        else:
            raise InvalidPotentialError("Invalid number of targets.")
        

    def add_indices(self,indices):
        """Adds the given indices to the list of indices.

        Parameters:

        indices: list of integers
            list of additional indices on which the potential acts
        """
        if(self.accepts_target_list(indices)):
            if self.indices == None:
                self.indices = []

            for stuff in indices:
                self.indices.append(stuff)
        elif(self.accepts_target_list([indices])):
            if self.indices == None:
                self.indices = []

            self.indices.append(indices)
        else:
            raise InvalidPotentialError("Invalid number of targets.")

    def set_parameter_values(self,values):
        """Sets the numeric values of all parameters.

        Parameters:

        values: list of doubles
            list of values to be assigned to parameters
        """
        if len(values) == number_of_parameters(self.potential_type):
            self.parameters = values
        else:
            raise InvalidPotentialError("The potential '{pot}' takes {n_par} parameters, not {n_in}.".
                                        format(pot=self.potential_type,
                                               n_par=number_of_parameters(self.potential_type),
                                               n_in=len(values)))

    def set_parameter_value(self,parameter_name,value):
        """Sets a given parameter to the desired value.

        Parameters:

        parameter_name: string
            name of the parameter
        value: double
            the new value of the parameter
        """
        index = index_of_parameter(self.potential_type,parameter_name)
        if index > -1:
            self.parameters[index] = value
        else:
            raise InvalidPotentialError("The potential '{pot}' does not have a parameter '{para}'".
                                        format(pot=self.potential_type,para=parameter_name))

    def set_cutoff(self,cutoff):
        """Sets the cutoff to a given value.

        This method affects the hard cutoff.
        For a detailed explanation on how to define a soft cutoff, see :meth:`~pysic.Potential.set_cutoff_margin`.
        
        Parameters:

        cutoff: double
            new cutoff for the potential
        """
        self.cutoff = cutoff


    def set_cutoff_margin(self,margin):
        """Sets the margin for smooth cutoff to a given value.

        Many potentials decay towards zero in infinity, but in a numeric simulation
        they are cut at a finite range as specified by the cutoff radius. If the potential
        is not exactly zero at this range, a discontinuity will be introduced.
        It is possible to avoid this by including a smoothening factor in the potential
        to force a decay to zero in a finite interval.

        This method defines the decay interval :math:`r_\mathrm{hard}-r_\mathrm{soft}`.
        Note that if the soft cutoff value is made smaller than 0 or larger than the
        hard cutoff value an :class:`~pysic.InvalidPotentialError` is raised.

        Parameters:

        margin: double
            The new cutoff margin
        """
        if margin < 0.0:
            raise InvalidPotentialError("A negative cutoff margin of {marg} was given.".format(marg=str(margin)))
        if margin > self.cutoff:
            raise InvalidPotentialError("A cutoff margin ({marg}) longer than the hard cutoff ({hard}) was given.".format(marg=str(margin), hard=str(self.cutoff)))
        
        self.cutoff_margin = margin


    def set_soft_cutoff(self,cutoff):
        """Sets the soft cutoff to a given value.

        For a detailed explanation on the meaning of a soft cutoff, see
        :meth:`~pysic.Potential.set_cutoff_margin`.
        Note that actually the cutoff margin is recorded, so changing the
        hard cutoff (see :meth:`~pysic.Potential.set_cutoff`) will also affect the
        soft cutoff.

        Parameters:

        cutoff: double
            The new soft cutoff
        """
        self.set_cutoff_margin(self.cutoff-cutoff)


    def describe(self):
        """Prints a short description of the potential using the method :meth:`~pysic.describe_potential`."""
        description_of_potential(self.potential_type,
                                 self.parameters,
                                 self.cutoff,
                                 self.symbols,
                                 self.tags,
                                 self.indices)

class Pysic:
    """A calculator class providing the necessary methods for interfacing with `ASE`_.

    Pysic is a calculator for evaluating energies and forces for given atomic structures
    according to the given :class:`~pysic.Potential` set. Neither the geometry nor the
    potentials have to be specified upon creating the calculator, as they can be specified
    or changed later. They are necessary for actual calculation, of course.

    Simulation geometries must be defined as `ASE Atoms`_. This object contains both the
    atomistic coordinates and supercell parameters.

    Potentials must be defined as a list of :class:`~pysic.Potential` objects. The full potential
    is a sum of the individual potentials.
    
    .. _ASE: https://wiki.fysik.dtu.dk/ase/
    .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

    Parameters:

    atoms: `ASE Atoms`_ object
        an Atoms object containing the full simulation geometry
    potentials: list of :class:`~pysic.Potential` objects
        list of potentials for describing interactions
    force_initialization: boolean
        If true, calculations always fully initialize the fortran core.
        If false, only the coordinates of atoms and cell vectors are updated
        in the core upon calculation, if possible. Especially the charges of
        atoms will be remembered by the core.
    """

    id_number = 1
    
    def __init__(self,atoms=None,potentials=None,full_initialization=False):
        self.id = Pysic.id_number
        Pysic.id_number += 1
        
        self.neighbor_lists_ready = False
        self.neighbor_lists_waiting = False
        self.potential_lists_ready = False
        self.potentials_ready = False
        self.atoms_ready = False
        self.cell_ready = False
        
        self.structure = None
        self.neighbor_list = None
        self.potentials = None
        self.set_atoms(atoms)
        self.set_potentials(potentials)
        self.forces = 0.0
        self.stress = 0.0
        self.energy = 0.0
        self.forces_calculated = False
        self.energy_calculated = False
        self.stress_calculated = False
        self.force_core_initialization = full_initialization


    def __del__(self):
        self.release_fortran_core()
            

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

    
    def calculation_required(self, atoms=None, quantities=['forces','energy','stress']):
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
            list of keywords 'energy', 'forces', 'stress'
        """
        
        do_it = []
        try:
            assert isinstance(quantities, list)
            list_of_quantities = quantities
        except:
            list_of_quantities = [ quantities ]
        
        for mark in list_of_quantities:
            if mark == 'energy':
                do_it.append(not self.energy_calculated)
            elif mark == 'forces':
                do_it.append(not self.forces_calculated)
            elif mark == 'stress':
                do_it.append(not self.stress_calculated)
            else:
                do_it.append(False)
        return any(do_it)


    def get_atoms(self):
        """Returns the `ASE Atoms`_ object assigned to the calculator."""
        return self.structure


    def get_neighbor_lists(self):
        """Returns the `ASE NeighborList`_ object assigned to the calculator.

        The neighbor lists are generated according to the given `ASE Atoms`_ object
        and the :class:`~pysic.Potential` objects of the calculator. Note that the lists
        are created when the core is set or if the method :meth:`~pysic.Pysic.create_neighbor_lists`
        is called.
        """
        return self.neighbor_list


    def get_potentials(self):
        """Returns the list of potentials assigned to the calculator."""
        return self.potentials

    
    def get_forces(self, atoms=None):
        """Returns the forces.

        If the atoms parameter is given, it will be used for updating the
        structure assigned to the calculator prior to calculating the forces.
        Otherwise the structure already associated with the calculator is used.

        The calculator checks if the forces have been calculated already
        via :meth:`~pysic.Pysic.calculation_required`. If the structure
        has changed, the forces are calculated using :meth:`~pysic.Pysic.calculate_forces`

        Parameters:

        atoms: `ASE atoms`_ object
            the structure for which the forces are determined
        """
        self.set_atoms(atoms)
        if self.calculation_required(atoms,'forces'):
            self.calculate_forces()

        return self.forces


    def get_potential_energy(self, atoms=None, force_consistent=False):
        """Returns the potential energy.

        If the atoms parameter is given, it will be used for updating the
        structure assigned to the calculator prior to calculating the energy.
        Otherwise the structure already associated with the calculator is used.

        The calculator checks if the energy has been calculated already
        via :meth:`~pysic.Pysic.calculation_required`. If the structure
        has changed, the energy is calculated using :meth:`~pysic.Pysic.calculate_energy`

        Parameters:

        atoms: `ASE atoms`_ object
            the structure for which the energy is determined
        force_consistent: logical
            ignored at the moment
        """
        self.set_atoms(atoms)
        if self.calculation_required(atoms,'energy'):
            self.calculate_energy()            
        
        return self.energy


    def get_stress(self, atoms=None):
        """Returns the stress.

        If the atoms parameter is given, it will be used for updating the
        structure assigned to the calculator prior to calculating the stress.
        Otherwise the structure already associated with the calculator is used.

        The calculator checks if the stress has been calculated already
        via :meth:`~pysic.Pysic.calculation_required`. If the structure
        has changed, the stress is calculated using :meth:`~pysic.Pysic.calculate_stress`

        Parameters:

        atoms: `ASE atoms`_ object
            the structure for which the stress is determined
        """
        self.set_atoms(atoms)
        if self.calculation_required(atoms,'stress'):
            self.calculate_stress()
        
        return self.stress

    
    def set_atoms(self, atoms=None):
        """Assigns the calculator with the given structure .

        Parameters:

        atoms: `ASE atoms`_ object
            the structure to be calculated
        """
        if atoms == None:
            pass
        else:
            if(self.structure != atoms):
                self.forces_calculated = False
                self.energy_calculated = False
                self.stress_calculated = False
                self.neighbor_lists_ready = False

                # NB: this avoids updating the potential lists every time an atom moves, but
                # at the same time it is prone to miss an update if the user changes the
                # structure so that the number of atoms remains the same...
                try:
                    if(self.structure.get_number_of_atoms() != atoms.get_number_of_atoms()):
                        self.potential_lists_ready = False
                except:
                    self.potential_lists_ready = False
                
                self.atoms_ready = False
                try:
                    if((self.structure.get_cell() != atoms.get_cell()).any()):
                        self.cell_ready = False
                except:
                    self.cell_ready = False
                
                self.structure = atoms.copy()


    def set_potentials(self, potentials):
        """Assign a list of potentials to the calculator.

        Parameters:

        potentials: list of :class:`~pysic.Potential` objects
            a list of potentials to describe interactinos
        """
        if potentials == None:
            pass
        else:
            self.forces_calculated = False
            self.energy_calculated = False
            self.stress_calculated = False
            self.potential_lists_ready = False
            self.potentials_ready = False
            self.neighbor_lists_ready = False
            self.neighbor_lists_waiting = False

            try:
                assert isinstance(potentials,list)
                self.potentials = potentials
            except:
                self.potentials = [potentials]
    
    
    def add_potential(self, potential):
        """Add a potential to the list of potentials.

        Parameters:

        potential: :class:`~pysic.Potential` object
            a new potential to describe interactions
        """

        if self.potentials == None:
            self.potentials = []

        self.potentials.append(potential)
        self.forces_calculated = False
        self.energy_calculated = False
        self.stress_calculated = False
        self.potential_lists_ready = False
        self.potentials_ready = False
        self.neighbor_lists_ready = False
        self.neighbor_lists_waiting = False
        
    
    def create_neighbor_lists(self,cutoffs=None,marginal=pu.neighbor_marginal):
        """Initializes the neighbor lists.

        In order to do calculations at reasonable speed, the calculator needs 
        a list of neighbors for each atom. For this purpose, the `ASE NeighborList`_
        are used. This method initializes these lists according to the given
        cutoffs.

        .. _ASE NeighborList: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#building-neighbor-lists

        Parameters:

        cutoffs: list of doubles
            a list containing the cutoff distance for each atom
        marginal: double
            the skin width of the neighbor list
        """
        if cutoffs == None:
            cutoffs = self.get_individual_cutoffs(0.5)
        self.neighbor_list = nbl.NeighborList(cutoffs,skin=marginal,sorted=False,self_interaction=False,bothways=True)
        self.neighbor_lists_waiting = True
        self.neighbor_lists_ready = False


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
            return self.structure.get_number_of_atoms()*[0.0]
        else:
            cuts = []
            for symbol, tags, index in zip(self.structure.get_chemical_symbols(),
                                           self.structure.get_tags(),
                                           range(self.structure.get_number_of_atoms())):
                max_cut = 0.0
                
                for potential in self.potentials:
                    active_potential = False
                    
                    if potential.get_different_symbols().count(symbol) > 0 or potential.get_different_tags().count(tags) > 0 or potential.get_different_indices().count(index) > 0:
                        active_potential = True
                    
                    if active_potential and potential.get_cutoff() > max_cut:
                        max_cut = potential.get_cutoff()
                cuts.append(max_cut*scaler)
            return cuts


    def calculate_forces(self):
        """Calculates forces.

        Calls the Fortran core to calculate forces for the currently assigned structure.
        """
        self.set_core()
        self.forces = pf.pysic_interface.calculate_forces(self.structure.get_number_of_atoms()).transpose()
        self.forces_calculated = True
        

    def calculate_energy(self):
        """Calculates the potential energy.

        Calls the Fortran core to calculate the potential energy for the currently assigned structure.
        """
        self.set_core()
        self.energy = pf.pysic_interface.calculate_energy()
        self.energy_calculated = True


    def calculate_stress(self):
        """Calculates the stress tensor.

        Calls the Fortran core to calculate the stress tensor for the currently assigned structure.
        """
        self.stress = 0.0
        self.stress_calculated = True
    

    def get_core_readiness(self):
        """Returns a tuple of booleans used for notifying necessary core updates.

        This is mainly a debug utility. If any of the values are False,
        the corresponding part of the core will be updated when calculation
        is invoked the next time.

        The returned booleans are in order:
        atoms_ready,
        cell_ready,
        potentials_ready,
        potential_lists_ready,
        neighbor_lists_waiting,
        neighbor_lists_ready
        """
        return (self.atoms_ready,
                self.cell_ready,
                self.potentials_ready,
                self.potential_lists_ready,
                self.neighbor_lists_waiting,
                self.neighbor_lists_ready)

    def set_core(self):
        """Sets up the Fortran core for calculation.

        If the core is not initialized, if the number of atoms has changed, or
        if full initialization is forces, the core is initialized from scratch.
        Otherwise, only the atomic coordinates and momenta are updated.
        Potentials, neighbor lists etc. are also updated if they have been edited.
        """
        if not self.core_is_available():
            raise LockedCoreError("core is not available")

        do_full_init = False
        if self.force_core_initialization:
            do_full_init = True
        elif not self.owns_the_core():
            do_full_init = True
        elif self.structure.get_number_of_atoms() != pf.pysic_interface.get_number_of_atoms():
            do_full_init = True
            
        if do_full_init:
            self.initialize_fortran_core()
        else:
            if not self.atoms_ready:
                self.update_core_coordinates()

            if not self.cell_ready:
                self.update_core_supercell()
                    
            if not self.potentials_ready:
                self.update_core_potentials()

            if not self.potential_lists_ready:
                self.update_core_potential_lists()

            if not self.neighbor_lists_waiting:
                self.create_neighbor_lists(self.get_individual_cutoffs(0.5))

            if not self.neighbor_lists_ready:
                self.update_core_neighbor_lists()



    def update_core_potential_lists(self):
        """Initializes the potential lists.

        Since one often runs :class:`~pysic.Pysic` with a set of potentials,
        the core pre-analyzes which potentials affect each atom and saves a list
        of such potentials for every particle. This method asks the core to
        generate these lists.
        """
        if not self.atoms_ready:
            raise MissingAtomsError("Creating potential lists before updating atoms in core.")
        pf.pysic_interface.create_potential_list()
        self.potential_lists_ready = True
        

    def update_core_potentials(self):
        """Generates potentials for the Fortran core."""
        n_pots = 0
        # count the number of separate potentials
        for pot in self.potentials:
            
            try:
                alltargets = pot.get_symbols()
                for targets in alltargets:
                    perms = permutations(targets)
                    different = set(perms)
                    n_pots += len(different)
            except:
                pass
            try:
                alltargets = pot.get_tags()
                for targets in alltargets:
                    perms = permutations(targets)
                    different = set(perms)
                    n_pots += len(different)
            except:
                pass
            try:
                alltargets = pot.get_indices()
                for targets in alltargets:
                    perms = permutations(targets)
                    different = set(perms)
                    n_pots += len(different)
            except:
                pass

        pf.pysic_interface.allocate_potentials(n_pots)
        for pot in self.potentials:
            n_targ = pot.get_number_of_targets()
            no_symbs = np.array( n_targ*[pu.str2ints('xx',2)] ).transpose()
            no_tags = np.array( n_targ*[-9] )
            no_inds = np.array( n_targ*[-9] )

            try:
                alltargets = pot.get_symbols()
                for targets in alltargets:
                    int_orig_symbs = []
                    for orig_symbs in targets:
                        for label in orig_symbs:
                            int_orig_symbs.append( pu.str2ints(label,2) )

                    perms = permutations(targets)
                    different = set(perms)
                    for symbs in different:
                        int_symbs = []
                        for label in symbs:
                            int_symbs.append( pu.str2ints(label,2) )

                        pf.pysic_interface.add_potential(pot.get_potential_type(),
                                                         np.array( pot.get_parameter_values() ),
                                                         pot.get_cutoff(),
                                                         pot.get_soft_cutoff(),
                                                         np.array( int_symbs ).transpose(),
                                                         no_tags,
                                                         no_inds,
                                                         np.array( int_orig_symbs ).transpose(),
                                                         no_tags,
                                                         no_inds )
            except:
                pass
            try:
                alltargets = pot.get_tags()
                for targets in alltargets:
                    orig_tags = targets
                    perms = permutations(targets)
                    different = set(perms)

                    for tags in different:
                        pf.pysic_interface.add_potential(pot.get_potential_type(),
                                                         np.array( pot.get_parameter_values() ),
                                                         pot.get_cutoff(),
                                                         pot.get_soft_cutoff(),
                                                         no_symbs,
                                                         np.array( tags ),
                                                         no_inds,
                                                         no_symbs,
                                                         np.array(orig_tags),
                                                         no_inds )
            except:
                pass
            try:
                alltargets = pot.get_indices()                
                for targets in alltargets:
                    orig_inds = targets
                    perms = permutations(targets)
                    different = set(perms)

                    for inds in different:
                        pf.pysic_interface.add_potential(pot.get_potential_type(),
                                                         np.array( pot.get_parameter_values() ),
                                                         pot.get_cutoff(),
                                                         pot.get_soft_cutoff(),
                                                         no_symbs,
                                                         no_tags,
                                                         np.array( inds ),
                                                         no_symbs,
                                                         no_tags,
                                                         np.array(orig_inds) )
            except:
                pass

        self.potentials_ready = True
        
    
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

        self.forces_calculated = False
        self.energy_calculated = False
        self.stress_calculated = False

        pf.pysic_interface.update_atom_coordinates(positions,momenta)
        self.atoms_ready = True
        if not self.neighbor_lists_waiting:
            self.create_neighbor_lists(self.get_individual_cutoffs(0.5))
        
        self.update_core_neighbor_lists()
        
        
    def update_core_supercell(self):
        """Updates the supercell in the Fortran core."""
        vectors = np.array( self.structure.get_cell() ).transpose()
        inverse = np.linalg.inv(np.array( self.structure.get_cell() )).transpose()
        periodicity = np.array( self.structure.get_pbc() )
        
        pf.pysic_interface.create_cell(vectors,inverse,periodicity)
        
        self.cell_ready = True
        self.neighbor_lists_ready = False


    def update_core_neighbor_lists(self):
        """Updates the neighbor lists in the Fortran core.

         If uninitialized, the lists are created first via :meth:`~pysic.Pysic.create_neighbor_lists`.
         """
        
        if not self.atoms_ready:
            raise MissingAtomsError("Creating neighbor lists before updating atoms in the core.")
        if not self.neighbor_lists_waiting:            
            self.create_neighbor_lists(self.get_individual_cutoffs(0.5))
            #raise MissingNeighborsError("The neighbor lists have not been pre-initialized.")

        self.neighbor_list.update(self.structure)
        for index in range(self.structure.get_number_of_atoms()):
            [nbors,offs] = self.neighbor_list.get_neighbors(index)                
            pf.pysic_interface.create_neighbor_list(index+1,np.array(nbors),np.array(offs).transpose())
        self.neighbor_lists_ready = True
        

    def initialize_fortran_core(self):
        """Fully initializes the Fortran core, creating the atoms, supercell, potentials, and neighbor lists."""
        
        masses = np.array( self.structure.get_masses() )
        charges = np.array( self.structure.get_charges() )
        positions = np.array( self.structure.get_positions() ).transpose()
        momenta = np.array( self.structure.get_momenta() ).transpose()
        tags = np.array( self.structure.get_tags() )
        elements = self.structure.get_chemical_symbols()

        for index in range(len(elements)):
            elements[index] = pu.str2ints(elements[index],2)

        elements = np.array( elements ).transpose()

        self.create_neighbor_lists(self.get_individual_cutoffs(0.5))
        #self.neighbor_list.update(self.structure) #now done in update_core_neighbor_lists

        pf.pysic_interface.set_owner_id(-1)
        pf.pysic_interface.create_atoms(masses,charges,positions,momenta,tags,elements)
        self.atoms_ready = True
        
        self.update_core_supercell()
        self.update_core_potentials()
        self.update_core_neighbor_lists()
        self.update_core_potential_lists()
        pf.pysic_interface.distribute_mpi(self.structure.get_number_of_atoms())
        pf.pysic_interface.set_owner_id(self.id)


    def release_fortran_core(self):
        """Releases the Fortran core, deallocating memory and allowing another calculator to act."""
        if self.owns_the_core():
            pf.pysic_interface.release()        

        self.neighbor_lists_ready = False
        self.potential_lists_ready = False
        self.potentials_ready = False
        self.atoms_ready = False
        self.cell_ready = False


    def view_fortran_core(self):
        """Print some information on the data allocated in the Fortran core."""
        pf.pysic_interface.examine_atoms()
        pf.pysic_interface.examine_cell()
        pf.pysic_interface.examine_potentials()

    def core_is_available(self):
        """True if this calculator can use the Fortran core."""
        return pf.pysic_interface.can_be_accessed(self.id)

    def owns_the_core(self):
        """True if this calculator owns the core.

        One can create several instances of :class:`~pysic.Pysic` within Python, but there is only one
        Fortran core for the actual calculations. Therefore, it is recommended that
        only one :class:`~pysic.Pysic` is created. In case there are several calculators, some basic
        checks are done to assure that the calculators do not use the core simultaneously.
        This method simply checks if the core is currently owned by this particular
        :class:`~pysic.Pysic`."""
        return (self.id == pf.pysic_interface.get_owner_id())


    def get_numerical_energy_gradient(self, atom_index, shift=0.001, atoms=None):
        """Numerically calculates the gradient of energy with respect to moving a single particle.

        This is for debugging the forces."""

        if(atoms == None):
            system = self.structure.copy()
            orig_system = self.structure.copy()
        else:
            system = atoms.copy()
            orig_system = atoms.copy()

        self.energy_calculated = False
        self.atoms_ready = False
        system[atom_index].x += shift
        energy_xp = self.get_potential_energy(system)
        system[atom_index].x -= 2.0*shift
        energy_xm = self.get_potential_energy(system)
        system[atom_index].x += shift

        system[atom_index].y += shift
        energy_yp = self.get_potential_energy(system)
        system[atom_index].y -= 2.0*shift
        energy_ym = self.get_potential_energy(system)
        system[atom_index].y += shift

        system[atom_index].z += shift
        energy_zp = self.get_potential_energy(system)
        system[atom_index].z -= 2.0*shift
        energy_zm = self.get_potential_energy(system)

        self.energy_calculated = False
        self.atoms_ready = False
        self.get_potential_energy(orig_system)
        
        return [ -(energy_xp-energy_xm)/(2.0*shift),
                 -(energy_yp-energy_ym)/(2.0*shift),
                 -(energy_zp-energy_zm)/(2.0*shift) ]
