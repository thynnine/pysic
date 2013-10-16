#! /usr/bin/env python

import pysic.pysic_fortran as pf
import pysic.utility.f2py as pu
from pysic.charges.relaxation import ChargeRelaxation
from pysic.interactions.coulomb import CoulombSummation
import numpy as np
import copy
import atexit
from pysic.utility.mpi import *

# !!!: automatic initialization functions

# Below are general methods for accessing information about the built-in potentials.
# Note that the methods directly inquire the Fortran core, and thus all changes made
# to the core should automatically be taken into account, as long as the proper
# arrays are updated in Potentials.f90

# automatically call the initialization routine in Potentials.f90
pf.pysic_interface.start_potentials()
pf.pysic_interface.start_bond_order_factors()

# automatically initialize mpi - fortran side decides if we really are in mpi
pf.pysic_interface.start_mpi()

# automatically initialize the fortran rng
# it needs to be possible to force a seed by the user - to be implemented
pf.pysic_interface.start_rng(5301)


# !!!: automatic termination functions

def termination():
    """Function for deallocating memory in the core and finalizing the MPI.
    
    The function is automatically called upon program termination by 'atexit'.
    """
    pf.pysic_interface.release()
    pf.pysic_interface.finish_mpi()

atexit.register(termination)


def list_potentials():
    """Same as :meth:`~pysic.core.list_valid_potentials`
    """
    return list_valid_potentials()

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


def list_bond_order_factors():
    """Same as :meth:`~pysic.core.list_valid_bond_order_factors`
    """
    return list_valid_bond_order_factors()

def list_valid_bond_order_factors():
    """A list of names of bond order factors currently known by the core.

    The method retrieves from the core a list of the names of different bond factors
    currently implemented. Since the fortran core is directly accessed, any
    updates made in the core source code should get noticed automatically.
    """
    n_bonds = pf.pysic_interface.number_of_bond_order_factors()
    bond_codes = pf.pysic_interface.list_valid_bond_order_factors(n_bonds).transpose()
    bond_names = []
    for code in bond_codes:
        bond_names.append(pu.ints2str(code))

    return bond_names


def is_potential(potential_name):
    """Same as :meth:`~pysic.core.is_valid_potential`
    
    Parameters:
    
    potential_name: string
        the name of the potential
    """
    return pf.pysic_interface.is_potential(potential_name)

def is_valid_potential(potential_name):
    """Tells if the given string is the name of a potential.

    Parameters:

    potential_name: string
        the name of the potential
    """
    return pf.pysic_interface.is_potential(potential_name)


def is_bond_order_factor(bond_order_name):
    """Same as :meth:`~pysic.core.is_valid_bond_order_factor`

    Parameters:

    bond_order_name: string
        the name of the bond order factor
    """
    return pf.pysic_interface.is_bond_order_factor(bond_order_name)

def is_valid_bond_order_factor(bond_order_name):
    """Tells if the given string is the name of a bond order factor.

    Parameters:

    bond_order_name: string
        the name of the bond order factor
    """
    return pf.pysic_interface.is_bond_order_factor(bond_order_name)


def is_charge_relaxation(relaxation_name):
    """Same as :meth:`~pysic.core.is_valid_charge_relaxation`
        
        Parameters:
        
        relaxation_name: string
            the name of the relaxation mode
        """
    for mode in ChargeRelaxation.relaxation_modes:
        if(relaxation_name == mode):
            return True
    return False


def is_valid_charge_relaxation(relaxation_name):
    """Tells if the given string is the name of a charge relaxation mode.
        
        Parameters:
        
        relaxation_name: string
            the name of the relaxation mode
        """
    for mode in ChargeRelaxation.relaxation_modes:
        if(relaxation_name == mode):
            return True
    return False


def is_coulomb_summation(summation_name):
    """Same as :meth:`~pysic.core.is_valid_coulomb_summation`
        
        Parameters:
        
        summation_name: string
            the name of the summation mode
        """
    for mode in CoulombSummation.summation_modes:
        if(summation_name == mode):
            return True
    return False



def is_valid_coulomb_summation(summation_name):
    """Tells if the given string is the name of a coulomb summation mode.
        
        Parameters:
        
        summation_name: string
            the name of the summation mode
        """
    for mode in CoulombSummation.summation_modes:
        if(summation_name == mode):
            return True
    return False
    



def number_of_targets(potential_name):
    """Tells how many targets a potential or bond order factor acts on, i.e., is it pair or many-body.

    Parameters:

    potential_name: string
        the name of the potential or bond order factor
    """

    if(is_potential(potential_name)):
        return pf.pysic_interface.number_of_targets_of_potential(potential_name)
    elif(is_bond_order_factor(potential_name)):
        return pf.pysic_interface.number_of_targets_of_bond_order_factor(potential_name)
    else:
        return 0


def level_of_factor(bond_order_name):
    """Tells the level of a bond order factor, i.e., is it a per-atom or per-pair factor.
    
    Parameters:

    bond_order_name: string
        the name of the bond order factor
    """
    if(is_bond_order_factor(bond_order_name)):
        return pf.pysic_interface.level_of_bond_order_factor(bond_order_name)
    else:
        return 0


def number_of_parameters(potential_name,as_list=False):
    """Tells how many parameters a potential, bond order factor, charge relaxation mode or coulomb summation mode incorporates.

    A potential has a simple list of parameters and thus the function returns
    by default a single number. A bond order factor can incorporate parameters for
    different number of targets (some for single elements, others for pairs, etc.), and
    so a list of numbers is returned, representing the number of single, pair etc.
    parameters. If the parameter 'as_list' is given and is True, the result is
    a list containing one number also for a potential.

    Parameters:

    potential_name: string
        the name of the potential or bond order factor
    as_list: logical
        should the result always be a list
    """
    
    if(is_potential(potential_name)):
        n_params = pf.pysic_interface.number_of_parameters_of_potential(potential_name)
        if as_list:
            return [n_params]
        else:
            return n_params
    elif(is_bond_order_factor(potential_name)):
        n_targets = pf.pysic_interface.number_of_targets_of_bond_order_factor(potential_name)
        n_params = [ [""] ]*2
        for i in range(2):
            n_params[i] = pf.pysic_interface.number_of_parameters_of_bond_order_factor(potential_name,i+1)
        return n_params
    elif(is_charge_relaxation(potential_name)):
        n_params = len(ChargeRelaxation.relaxation_parameters[potential_name])
        if as_list:
            return [n_params]
        else:
            return n_params
    elif(is_coulomb_summation(potential_name)):
        n_params = len(CoulombSummation.summation_parameters[potential_name])
        if as_list:
            return [n_params]
        else:
            return n_params
    else:
        return 0

def names_of_parameters(potential_name):
    """Lists the names of the parameters of a potential, bond order factor, charge relaxation mode or coulomb summation mode.

    For a potential, a simple list of names is returned. For a bond order factor, the
    parameters are categorised according to the number of targets they apply to
    (single element, pair, etc.). So, for a bond order factor, a list of lists is
    returned, where the first list contains the single element parameters, the second
    list the pair parameters etc.

    Parameters:

    potential_name: string
        the name of the potential or bond order factor
    """
    
    if(is_potential(potential_name)):
        param_codes = pf.pysic_interface.names_of_parameters_of_potential(potential_name).transpose()

        param_names = []
        n_params = number_of_parameters(potential_name)
        index = 0
        for code in param_codes:
            if index < n_params:
                param_names.append(pu.ints2str(code).strip())
                index += 1
                
        return param_names
    
    elif(is_bond_order_factor(potential_name)):

        n_targets = pf.pysic_interface.number_of_targets_of_bond_order_factor(potential_name)
        param_codes = [ [0] ]*2
        param_names = [ [] ]*2
        
        n_params = number_of_parameters(potential_name)
        
        for i in range(2):
            param_codes[i] = pf.pysic_interface.names_of_parameters_of_bond_order_factor(potential_name,i+1).transpose()
            param_names[i] = [""]*n_params[i]


        target_index = 0
        for target_codes in param_codes:
            index = 0
            for code in target_codes:
                if index < n_params[target_index]:
                    param_names[target_index][index] = pu.ints2str(code).strip()
                    index += 1
            target_index += 1
        return param_names

    elif(is_charge_relaxation(potential_name)):
        return ChargeRelaxation.relaxation_parameters[potential_name]

    elif(is_coulomb_summation(potential_name)):
        return CoulombSummation.summation_parameters[potential_name]
            
    else:
        return []

def index_of_parameter(potential_name, parameter_name):
    """Tells the index of a parameter of a potential or bond order factor in the list of parameters the potential uses.

    For a potential, the index of the specified parameter is given.
    For a bond order factor, a list of two integers is given.
    These give the number of targets (single element, pair etc.) the parameter
    is associated with and the list index.
    
    Note especially that an index is returned, and these start counting from 0.
    So for a bond order factor, a parameter for pairs (2 targets) will return 1
    as the index for number of targets.
    
    Parameters:
    
    potential_name: string
        the name of the potential or bond order factor
    parameter_name: string
        the name of the parameter
    """
    
    if(is_potential(potential_name) or is_charge_relaxation(potential_name)):
        param_names = names_of_parameters(potential_name)
        index = 0
        for name in param_names:
            if name == parameter_name:
                return index
            index += 1
    elif(is_bond_order_factor(potential_name)):
        param_names = names_of_parameters(potential_name)
        target_index = 0
        for listed in param_names:
            index = 0
            for name in listed:
                if name == parameter_name:
                    return [target_index,index]
                index += 1
            target_index += 1

    else:
        return None



def descriptions_of_parameters(potential_name):
    """Returns a list of strings containing physical names of the parameters of a potential, bond order factor, or charge relaxation mode,
    e.g., 'spring constant' or 'decay length'.
    
    For a potential, a simple list of descriptions is returned. For a bond order factor, the
    parameters are categorised according to the number of targets they apply to
    (single element, pair, etc.). So, for a bond order factor, a list of lists is
    returned, where the first list contains the single element parameters, the second
    list the pair parameters etc.
    
    Parameters:

    potential_name: string
        the name of the potential or bond order factor
    """
    
    if(is_potential(potential_name)):
        param_codes = pf.pysic_interface.descriptions_of_parameters_of_potential(potential_name).transpose()
        param_notes = []
        n_params = number_of_parameters(potential_name)
        index = 0
        for code in param_codes:
            if index < n_params:
                param_notes.append(pu.ints2str(code).strip())
                index += 1
    
    elif(is_bond_order_factor(potential_name)):
        n_targets = pf.pysic_interface.number_of_targets_of_bond_order_factor(potential_name)
        param_codes = [ [0] ]*n_targets
        param_names = [ [] ]*n_targets
        n_params = number_of_parameters(potential_name)

        for i in range(n_targets):
            param_codes[i] = pf.pysic_interface.descriptions_of_parameters_of_bond_order_factor(potential_name,i+1).transpose()
            param_names[i] = [""]*n_params[i]

        target_index = 0
        for target_codes in param_codes:
            index = 0
            for code in target_codes:
                if index < n_params[target_index]:
                    param_names[target_index][index] = pu.ints2str(code).strip()
                    index += 1
            target_index += 1
        return param_names
            
    
    elif(is_charge_relaxation(potential_name)):
        return ChargeRelaxation.relaxation_parameter_descriptions[potential_name]
    
    elif(is_coulomb_summation(potential_name)):
        return CoulombSummation.summation_parameter_descriptions[potential_name]

    else:
        return []

    return param_notes

#
# ToDo: implement descriptions of bond order factor and charge relaxation types
#

def description_of_potential(potential_name, parameter_values=None, cutoff=None,
                             elements=None, tags=None, indices=None):
    """Prints a brief description of a potential.

    If optional arguments are provided,
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


class CoreMirror:
    """A class representing the status of the core.

    Whenever data is being passed over to the core for calculation,
    it should also be saved in the CoreMirror. This makes the CoreMirror
    reflect the current status of the core. Then, when something needs to be
    calculated, the :class:`~pysic.calculator.Pysic` calculator can simply check that
    it contains the same system as the CoreMirror to ensure that the
    core operates on the correct data.

    All data given to CoreMirror is saved as deep copies, i.e., not as
    the objects themselves but objects with exactly the same contents.
    This way if the original objects are modified, the ones in CoreMirror
    are not. This is the proper way to work, since the Fortran core
    obviously does not change without pushing the changes in the Python
    side to the core first.

        
    Since exactly one CoreMirror should exist during the simulation, 
    deletion of the instance (which should happen at program termination automatically)
    will automatically trigger release of memory in the Fortran core
    as well as termination of the MPI framework.
        
    .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
    .. _ASE NeighborList: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#building-neighbor-lists

    """

    def __init__(self):
        self.structure = None
        self.potentials = None
        self.neighbor_lists = None
        self.coulomb = None
        self.cutoffs = None
        self.potential_lists_ready = False
        self.bond_order_factor_lists_ready = False
        self.mpi_ready = False        

    def __repr__(self):
        return "CoreMirror()"


    def get_atoms(self):
        """Returns the `ASE Atoms`_ structure stored in the CoreMirror.
        """
        return self.structure


    def view_fortran(self):
        """Print some information on the data allocated in the Fortran core.

        This is mainly a debugging utility for directly seeing the status of the core.
        It can be accessed through::

           >>> pysic.Pysic.core.view_fortran()

        The result is a bunch of data dumped to stdout. The function does
        not return anything.
        """
        pf.pysic_interface.examine_atoms()
        pf.pysic_interface.examine_cell()
        pf.pysic_interface.examine_potentials()
        pf.pysic_interface.examine_bond_order_factors()


    def set_atoms(self, atoms):
        """Copies and stores the entire `ASE Atoms`_ instance.

        Parameters:

        atoms: `ASE Atoms`_ object
            atomic structure to be saved"""
        self.structure = copy.deepcopy(atoms)
        del self.structure.constraints
        self.potential_lists_ready = False

    def set_charges(self, charges):
        """Copies and stores the charges of atoms in the `ASE Atoms`_ instance.
            
            Parameters:
            
            atoms: `ASE Atoms`_ object
                atomic structure containing the positions to be saved.
            """
        # the call for charges was changed between ASE 3.6 and 3.7
        try:
            self.structure.set_charges(charges)
        except:
            self.structure.set_initial_charges(charges)


    def set_atomic_positions(self, atoms):
        """Copies and stores the positions of atoms in the `ASE Atoms`_ instance.

        Parameters:

        atoms: `ASE Atoms`_ object
            atomic structure containing the positions to be saved.
        """
        self.structure.set_positions(atoms.get_positions())
        
    def set_atomic_momenta(self, atoms):
        """Copies and stores the momenta of atoms in the `ASE Atoms`_ instance.

        Parameters:

        atoms: `ASE Atoms`_ object
            atomic structure containing the momenta to be saved.
        """
        self.structure.set_momenta(atoms.get_momenta())
    
    
    def set_cell(self, atoms):
        """Copies and stores the supercell in the `ASE Atoms`_ instance.

        Parameters:

        atoms: `ASE Atoms`_ object
            atomic structure containing the supercell to be saved.
        """
        self.structure.set_cell(atoms.get_cell())
        self.structure.set_pbc(atoms.get_pbc())
        
    def set_potentials(self, potentials):
        """Copies and stores :class:`~pysic.interactions.local.Potential` potentials.

        The :class:`~pysic.interactions.local.Potential` instances are copied as a whole,
        so any possible :class:`~pysic.interactions.bondorder.Coordinator` and
        :class:`~pysic.interactions.bondorder.BondOrderParameters` objects are also stored.

        Parameters:

        atoms: list of :class:`~pysic.interactions.local.Potential` objects
            Potentials to be saved.
        """
        self.potentials = copy.deepcopy(potentials)
        self.potential_lists_ready = False

    def set_neighbor_lists(self, lists):
        """Copies and stores the neighbor lists.

        Parameters:

        atoms: `ASE NeighborList`_ object
            Neighbor lists to be saved.
        """
        self.neighbor_lists = copy.deepcopy(lists)

        
    def set_coulomb(self,coulomb):
        """Copies and stores the Coulomb summation algorithm.
        
        Parameters:
            
        coulomb: :class:`~pysic.interactions.coulomb.CoulombSummation`
            Coulomb summation algorithm to be saved            
            """
        self.coulomb = copy.deepcopy(coulomb)
    
    
    def atoms_ready(self, atoms):
        """Checks if the positions and momenta of the given atoms match those in the core.

        True is returned if the structures match, False otherwise.

        Parameters:

        atoms: `ASE Atoms`_ object
            The atoms to be compared.
        """
        if self.structure is None:
            return False
        if(len(self.structure) != len(atoms)):
            return False
        if((self.structure.get_atomic_numbers() != atoms.get_atomic_numbers()).any()):
            return False
        if((self.structure.get_positions() != atoms.get_positions()).any()):
            return False
        if((self.structure.get_momenta() != atoms.get_momenta()).any()):
            return False
    
        return True
                        
    def charges_ready(self, atoms):
        """Checks if the charges of the given atoms match those in the core.
            
            True is returned if the charges match, False otherwise.
            
            Parameters:
            
            atoms: `ASE Atoms`_ object
                The atoms to be compared.
            """
        
        if self.structure == None:
            return False
        if len(self.structure) != len(atoms):
            return False
            
        # the call for charges was changed between ASE 3.6 and 3.7
        try:
            if ((self.structure.get_initial_charges() != atoms.get_initial_charges()).any()):
                return False
        except:
            if ((self.structure.get_charges() != atoms.get_charges()).any()):
                return False

        return True
            
    def cell_ready(self,atoms):
        """Checks if the given supercell matches that in the core.

        True is returned if the structures match, False otherwise.

        Parameters:

        atoms: `ASE Atoms`_ object
            The cell to be compared.
        """
        if self.structure is None:
            return False
        if((self.structure.get_cell() != atoms.get_cell()).any()):
            return False
        if((self.structure.get_pbc() != atoms.get_pbc()).any()):
            return False
        return True

    def potentials_ready(self, pots):
        """Checks if the given potentials match those in the core.

        True is returned if the potentials match, False otherwise.

        Parameters:

        atoms: list of :class:`~pysic.interactions.local.Potential` objects
            The potentials to be compared.
        """
        if self.potentials == None:
            if pots == None:
                return True
            else:
                return False
        return (self.potentials == pots)

    def neighbor_lists_ready(self, lists):
        """Checks if the given neighbor lists match those in the core.

        True is returned if the structures match, False otherwise.

        Parameters:

        atoms: `ASE NeighborList`_ object
            The neighbor lists to be compared.
        """
        if self.neighbor_lists == None:
            return False
        return self.neighbor_lists == lists

            
            

            
    def coulomb_summation_ready(self,coulomb):
        """Checks if the given Coulomb summation matches that in the core.
            
            True is returned if the summation algorithms match, False otherwise.
            
            Parameters: :class:`~pysic.interactions.coulomb.CoulombSummation`
                the summation algorithm to be compared
            """
        if self.coulomb == None:
            return False
        return self.coulomb == coulomb


