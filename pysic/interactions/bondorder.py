#! /usr/bin/env python

from pysic.core import *
from pysic.utility.error import InvalidParametersError, InvalidCoordinatorError, warn
import pysic.pysic_fortran as pf

class BondOrderParameters:
    """Class for representing a collection of parameters for bond order calculations.

    Calculating bond order factors using Tersoff-like methods defined in
    :class:`~pysic.interactions.bondorder.Coordinator` requires several parameters per element and
    element pair. To facilitate the handling of all these parameters, they are
    wrapped in a BondOrderParameters object.

    The object can be created empty and filled later with the parameters. Alternatively,
    a list of parameters can be given upon initialization in which case it is passed
    to the :meth:`~pysic.interactions.bondorder.BondOrderParameter.set_parameters` method.

    Parameters:

    bond_order_type: string
        a keyword specifying the type of the bond order factor
    soft_cut: double
        The soft cutoff for calculating partial coordination.
        Any atom closer than this is considered a full neighbor.
    hard_cut: double
        The hard cutoff for calculating partial coordination.
        Any atom closer than this is considered (at least) a partial neighbor
        and will give a fractional contribution to the total coordination.
        Any atom farther than this will not contribute to the neighbor count.
    parameters: list of doubles
        a list of parameters to be contained in the parameter object
    symbols: list of strings
        a list of elements on which the factor is applied
    """


    def __init__(self,bond_order_type,cutoff=0.0,cutoff_margin=0.0,parameters=None,symbols=None):
        self.params = []

        if(not is_valid_bond_order_factor(bond_order_type)):
            raise InvalidParametersError('There is no bond order factor called "{bof}"'.format(bof=bond_order_type))

        self.bond_order_type = bond_order_type
        self.cutoff = cutoff
        self.cutoff_margin = 0.0
        self.set_cutoff_margin(cutoff_margin)
        self.n_targets = number_of_targets(bond_order_type)
        self.names_of_params = names_of_parameters(bond_order_type)
        self.n_params = number_of_parameters(bond_order_type)
        self.level = level_of_factor(bond_order_type)

        if parameters == None:
            self.parameters = 2*[[]]
            for index in range(len(self.parameters)): 
                self.parameters[index] = self.n_params[index]*[0.0]
        else:  
            self.set_parameters(parameters)

        self.symbols = None
        
        self.set_symbols(symbols)

        
    def __repr__(self):
        return "BondOrderParameters( bond_order_type='"+self.bond_order_type+\
               "',cutoff="+str(self.cutoff)+",cutoff_margin="+str(self.cutoff_margin)+\
               ",parameters="+str(self.parameters)+",symbols="+str(self.symbols)+" )"


    def __eq__(self,other):
        try:
            if self.bond_order_type != other.bond_order_type:
                return False
            if self.cutoff != other.cutoff:
                return False
            if self.cutoff_margin != other.cutoff_margin:
                return False
            if self.parameters != other.parameters:
                return False
            if self.symbols != other.symbols:
                return False
        except:
            return False

        return True

    def __ne__(self,other):
        return not self.__eq__(other)
        

            
    def accepts_parameters(self,params):
        """Test if the given parameters array has the correct dimensions.

        A bond order parameter can contain separate parameters for single,
        pair etc. elements and each class can have a different number of
        parameters. This method checks if the given list has the correct
        dimensions.

        Parameters:

        params: list of doubles
            list of parameters
        """
        if(len(params) != len(self.n_params)):
            return False
        for i in range(len(self.n_params)):
            if(len(params[i]) != self.n_params[i]):
                return False
        return True


    def get_number_of_parameters(self):
        """Returns the number of parameters the bond order parameter object contains.
        """
        return self.n_params


    def includes_scaling(self):
        """Returns True iff there are scaling paramters.
        """
        return self.n_params[0] > 0

    def get_level(self):
        """Returns the level of the factor, i.e., is it a per-atom or par-pair factor.
        """
        return self.level

    def get_number_of_targets(self):
        """Returns the (maximum) number of targets the bond order factor affects.
        """
        return self.n_targets

    def get_parameters_as_list(self):
        """Returns the parameters of the bond order factor as a single list.

        The generated list first contains the single element parameters, then
        pair parameters, etc.
        """
        allparams = []
        for target_params in self.parameters:
            for param in target_params:
                allparams.append(param)
        return allparams

    def set_symbols(self,symbols):
        """Sets the list of symbols to equal the given list.

        Parameters:

        symbols: list of strings
            list of element symbols on which the bond order factor acts
        """
        if symbols == None:
            return
        if(self.accepts_target_list(symbols)):
            self.symbols = symbols
        elif(self.accepts_target_list([symbols])):
            self.symbols = [symbols]
        else:
            raise InvalidParametersError("Invalid number of targets.")

        
    def add_symbols(self,symbols):
        """Adds the given symbols to the list of symbols.

        Parameters:

        symbols: list of strings
            list of additional symbols on which the bond order factor acts
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
            raise InvalidParametersError("Invalid number of targets.")


    def get_symbols(self):
        """Returns the symbols the bond parameters affect.
        """

        return self.symbols

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
        indices = index_of_parameter(self.bond_order_type,param_name)
        return self.parameters[indices[0]][indices[1]]


    def get_bond_order_type(self):
        """Returns the keyword specifying the type of the bond order factor."""
        return self.bond_order_type

    def get_cutoff(self):
        """Returns the cutoff."""
        return self.cutoff

    def get_cutoff_margin(self):
        """Returns the margin for a smooth cutoff."""
        return self.cutoff_margin

    def get_soft_cutoff(self):
        """Returns the lower limit for a smooth cutoff."""
        return self.cutoff-self.cutoff_margin


    def set_parameter_values(self,values):
        """Sets the numeric values of all parameters.

        Parameters:
        
        params: list of doubles
            list of values to be assigned to parameters
        """
        self.set_parameters(values)

    def set_parameter_value(self,parameter_name,value):
        """Sets a given parameter to the desired value.

        Parameters:

        parameter_name: string
            name of the parameter
        value: double
            the new value of the parameter
        """
        indices = index_of_parameter(self.bond_order_type,parameter_name)
        if indices != None:
            self.parameters[indices[0]][indices[1]] = value
        else:
            raise InvalidParametersError("The bond order factor '{pot}' does not have a parameter '{para}'".
                                        format(pot=self.bond_order_type,para=parameter_name))

    def set_cutoff(self,cutoff):
        """Sets the cutoff to a given value.

        This method affects the hard cutoff.
        
        Parameters:

        cutoff: double
            new cutoff for the bond order factor
        """
        self.cutoff = cutoff


    def set_cutoff_margin(self,margin):
        """Sets the margin for smooth cutoff to a given value.

        This method defines the decay interval :math:`r_\mathrm{hard}-r_\mathrm{soft}`.
        Note that if the soft cutoff value is made smaller than 0 or larger than the
        hard cutoff value an :class:`~pysic.utility.error.InvalidParametersError` is raised.

        Parameters:

        margin: double
            The new cutoff margin
        """
        if margin < 0.0:
            raise InvalidParametersError("A negative cutoff margin of {marg} was given.".format(marg=str(margin)))
        if margin > self.cutoff:
            raise InvalidParametersError("A cutoff margin ({marg}) longer than the hard cutoff ({hard}) was given.".format(marg=str(margin), hard=str(self.cutoff)))
        
        self.cutoff_margin = margin


    def set_soft_cutoff(self,cutoff):
        """Sets the soft cutoff to a given value.

        Note that actually the cutoff margin is recorded, so changing the
        hard cutoff (see :meth:`~pysic.interactions.bondorder.BondOrderParameters.set_cutoff`) will also affect the
        soft cutoff.

        Parameters:

        cutoff: double
            The new soft cutoff
        """
        self.set_cutoff_margin(self.cutoff-cutoff)




    def accepts_target_list(self,targets):
        """Tests whether a list is suitable as a list of targets, i.e., element symbols
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


    def set_parameters(self,params):
        """Sets the numeric values of all parameters.

        Equivalent to :meth:`~pysic.interactions.bondorder.BondOrderParameters.set_parameter_values`.

        Parameters:
        
        params: list of doubles
            list of values to be assigned to parameters
        """
        if self.accepts_parameters(params):
            self.parameters = params
        else:
            new_params = [[],[]]
            if len(params) == self.n_params[0]+self.n_params[1]:
                new_params[0] = params[0:self.n_params[0]]
                new_params[1] = params[self.n_params[0]:self.n_params[1]]
            
            if self.accepts_parameters(new_params):
                warn("Using parameters \n"+str(new_params)+\
                    "\ninstead of \n"+str(params),3)
                self.parameters = new_params
            else:
                raise InvalidParametersError('The bond order factor "{bof}" requires {num} parameters.'.format(bof=self.bond_order_type,num=str(self.n_params)))


        
        
class Coordinator:
    """Class for representing a calculator for atomic coordination numbers and bond order factors.

    Pysic can utilise 'Tersoff-like' potentials which are locally scaled according to
    bond order factors, related to the number of neighbors of each atom.
    The coordination calculator keeps track of updating
    the bond order factors and holds the parameters for calculating the values.

    When calculating forces also the derivatives of the coordination numbers are needed.
    Coordination numbers may be used repeatedly when calculating energies and forces even within
    one evaluation of the forces and therefore they are stored by the calculator. Derivatives
    are not stored since storing them could
    potentially require an N x N matrix, where N is the number of particles.

    The calculation of coordination is an operation on the geometry, not the complete physical
    system including the interactions, and so one can define coordination calculators as
    standalone objects as well. They always operate on the geometry currently allocated in the core.

    Parameters:

    bond_order_parameters: list of :class:`~pysic.interactions.bondorder.BondOrderParameters` objects
        Parameters for calculating bond order factors.
    """

    def __init__(self,bond_order_parameters=None):
        self.bond_order_params = []
        self.bond_order_factors = None
        if(bond_order_parameters != None):
            self.set_bond_order_parameters(bond_order_parameters)
        self.group_index = -1


    def __eq__(self,other):
        try:
            if self.group_index != other.group_index:
                return False
            if self.bond_order_params != other.bond_order_params:
                return False
        except:
            return False

        return True

    def __ne__(self,other):
        return not self.__eq__(other)

    def __repr__(self):
        return "Coordinator(bond_order_parameters={bonds})".format(bonds=str(self.bond_order_params))



    def get_bond_order_parameters(self):
        """Returns the bond order parameters of this Coordinator.
        """
        return self.bond_order_params


    def set_bond_order_parameters(self,params):
        """Assigns new bond order parameters to this Coordinator.

        Parameters:

        params: :class:`~pysic.interactions.bondorder.BondOrderParameters`
            new bond order parameters
        """

        if isinstance(params,list):
            self.bond_order_params = params
        else:
            self.bond_order_params = [params]

    def add_bond_order_parameters(self,params):
        """Adds the given parameters to this Coordinator.

        Parameters:
            
        params: :class:`~pysic.interactions.bondorder.BondOrderParameters`
            new bond order parameters
        """
        if isinstance(params,list):
            for param in params:
                self.bond_order_params.append(param)
        else:
            self.bond_order_params.append([params])
        

    def set_group_index(self,index):
        """Sets an index for the coordinator object.

        In the fortran core, bond order parameters are calculated by
        bond order parameters. Since a coordinator contains many, they are grouped
        to a coordinator via a grouping index when the core is initialized.
        This method is meant to be used for telling the Coordinator of this
        index. That allows the bond orders can be calculated by calling the
        Coordinator itself, since the index tells which bond parameters in the
        core are needed.

        Parameters:

        index: integer
            an index for grouping bond order parameters in the core
        """
        self.group_index = index

    def get_group_index(self):
        """Returns the group index of the Coordinator.
        """
        return self.group_index

    def calculate_bond_order_factors(self):
        """Recalculates the bond order factors for all atoms and stores them.

        Similarly to coordination numbers (:meth:`~pysic.interactions.bondorder.Coordinator.calculate_coordination`),
        this method only calculates the factors and stores them but does not return them.
        """
        n_atoms = pf.pysic_interface.get_number_of_atoms()
        self.bond_order_factors = pf.pysic_interface.calculate_bond_order_factors(n_atoms,self.group_index)


    def get_bond_order_factors(self):
        """Returns an array containing the bond order factors of all atoms.

        This method does not calculate the bond order factors but returns the
        precalculated array.
        """
        return self.bond_order_factors


    
    def get_bond_order_gradients(self,atom_index):
        """Returns an array containing the gradients of bond order factors of all atoms with respect to moving one atom.

        Parameters:

        atom_index: integer
            the index of the atom the position of which is being differentiated with
        """
        n_atoms = pf.pysic_interface.get_number_of_atoms()
        # add 1 to atom index because fortran indexing starts from 1, not 0
        return pf.pysic_interface.calculate_bond_order_gradients(n_atoms,self.group_index,atom_index+1).transpose()



    def get_bond_order_gradients_of_factor(self,atom_index):
        """Returns an array containing the gradients of the bond order factor of one atom with respect to moving any atom.

        Parameters:

        atom_index: integer
            the index of the atom the position of which is being differentiated with
        """
        n_atoms = pf.pysic_interface.get_number_of_atoms()
        # add 1 to atom index because fortran indexing starts from 1, not 0
        return pf.pysic_interface.calculate_bond_order_gradients_of_factor(n_atoms,self.group_index,atom_index+1).transpose()


