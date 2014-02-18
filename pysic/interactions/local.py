#! /usr/bin/env python

from pysic.core import *
from pysic.utility.error import InvalidPotentialError
import copy

class ProductPotential:
    """Class representing an interaction obtained by multiplying several :class:`~pysic.interactions.local.Potential` objects.
        
        Parameters:
        
        potentials: a list of :class:`~pysic.interactions.local.Potential` objects
            the potentials
    """

    def __init__(self,potentials):
        self.potentials = []
        self.n_targets = None
        for pot in potentials:
            self.add_potential(pot)
        

    def __eq__(self,other):
        try:
            if self.n_targets != other.n_targets:
                return False
            if self.potentials != other.potentials:
                return False
        except:
            return False
        return True


    def __ne__(self,other):
        return not self.__eq__(other)
    
    
    def set_potentials(self,potentials):
        """Sets the list of potentials for the product.
            
            Parameters:
            
            potentials: a list of :class:`~pysic.interactions.local.Potential` objects
                the potentials
            """

        self.potentials = []
        self.n_targets = None
        for pot in potentials:
            self.add_potential(pot)
        
    
    def add_potential(self,potential):
        """Adds a potential in the product.
        
        If another ProductPotential is given as an argument, the
        individual potentials in the product are added one by one to this
        product.
        
        Also a list of potentials can be given. Then all the listed potentials
        are added to the product.
            
        Parameters:
            
        potential: a :class:`~pysic.interactions.local.Potential` object
            the potential to be added
            """
        
        if isinstance(potential, list):
            potlist = potential
        else:
            potlist = [potential]
            
        for potti in potlist:
            if self.n_targets is None:
                self.n_targets = potti.get_number_of_targets()
        
            elif self.n_targets != potti.get_number_of_targets():
                raise InvalidPotentialError("You may only multiply potentials with the same number of targets.")

            # in case a product potential is given, add the individual potentials
            for pot in potti.get_potentials():
                self.potentials.append(copy.deepcopy(pot))


    def get_symbols(self):
        """Return a list of the chemical symbols (elements) on which the potential
            acts on."""
        return self.potentials[0].get_symbols()

    def get_tags(self):
        """Return the tags on which the potential acts on."""
        return self.potentials[0].get_tags()

    def get_indices(self):
        """Return a list of indices on which the potential acts on. """
        return self.potentials[0].get_indices()

    def get_different_symbols(self):
        """Returns a list containing each symbol the potential affects once."""
        return self.potentials[0].get_different_symbols()

    def get_different_tags(self):
        """Returns a list containing each tag the potential affects once."""
        return self.potentials[0].get_different_tags()

    def get_different_indices(self):
        """Returns a list containing each index the potential affects once."""
        return self.potentials[0].get_different_indices()

    def get_cutoff(self):
        """Returns the cutoff."""
        return self.potentials[0].get_cutoff()

    def get_cutoff_margin(self):
        """Returns the margin for a smooth cutoff."""
        return self.potentials[0].get_cutoff_margin()

    def get_soft_cutoff(self):
        """Returns the lower limit for a smooth cutoff."""
        return self.potentials[0].get_soft_cutoff()

    def get_number_of_targets(self):
        """Returns the number of targets."""
        return self.n_targets

    def get_coordinator(self):
        """Returns the Coordinator.
            """
        return self.potentials[0].get_coordinator()

    def set_coordinator(self, coordinator):
        """Sets a new Coordinator.
            
            Parameters:
            
            coordinator: a :class:`~pysic.interactions.bondorder.Coordinator` object
            the Coordinator
            """
        self.potentials[0].set_coordinator(coordinator)

    def set_symbols(self, symbols):
        """Sets the list of symbols to equal the given list.
            
            Parameters:
            
            symbols: list of strings
            list of element symbols on which the potential acts
            """
        self.potentials[0].set_symbols(symbols)

    def set_tags(self, tags):
        """Sets the list of tags to equal the given list.
            
            Parameters:
            
            tags: list of integers
            list of tags on which the potential acts
            """
        self.potentials[0].set_tags(tags)

    def set_indices(self, indices):
        """Sets the list of indices to equal the given list.
            
            Parameters:
            
            indices: list of integers
            list of integers on which the potential acts
            """
        self.potentials[0].set_indices(indices)

    def add_symbols(self, symbols):
        """Adds the given symbols to the list of symbols.
            
            Parameters:
            
            symbols: list of strings
            list of additional symbols on which the potential acts
            """
        self.potentials[0].add_symbols(symbols)

    def add_tags(self, tags):
        """Adds the given tags to the list of tags.
            
            Parameters:
            
            tags: list of integers
            list of additional tags on which the potential acts
            """
        self.potentials[0].add_tags(tags)
    
    def add_indices(self, indices):
        """Adds the given indices to the list of indices.
            
            Parameters:
            
            indices: list of integers
            list of additional indices on which the potential acts
            """
        self.potentials[0].add_indices(indices)

    def set_cutoff(self, cutoff):
        """Sets the cutoff to a given value.
            
            This method affects the hard cutoff.
            For a detailed explanation on how to define a soft cutoff, see :meth:`~pysic.interactions.local.Potential.set_cutoff_margin`.
            
            Parameters:
            
            cutoff: double
            new cutoff for the potential
            """
        self.potentials[0].set_cutoff(cutoff)
    
    def set_cutoff_margin(self, cutoff_margin):
        """Sets the margin for smooth cutoff to a given value.
            
            Many potentials decay towards zero in infinity, but in a numeric simulation
            they are cut at a finite range as specified by the cutoff radius. If the potential
            is not exactly zero at this range, a discontinuity will be introduced.
            It is possible to avoid this by including a smoothening factor in the potential
            to force a decay to zero in a finite interval.
            
            This method defines the decay interval :math:`r_\mathrm{hard}-r_\mathrm{soft}`.
            Note that if the soft cutoff value is made smaller than 0 or larger than the
            hard cutoff value an :class:`~pysic.utility.error.InvalidPotentialError` is raised.
            
            Parameters:
            
            margin: double
            The new cutoff margin
            """
        self.potentials[0].set_cutoff_margin(cutoff_margin)

    def set_soft_cutoff(self, cutoff):
        """Sets the soft cutoff to a given value.
            
            For a detailed explanation on the meaning of a soft cutoff, see
            :meth:`~pysic.interactions.local.Potential.set_cutoff_margin`.
            Note that actually the cutoff margin is recorded, so changing the
            hard cutoff (see :meth:`~pysic.interactions.local.Potential.set_cutoff`) will also affect the
            soft cutoff.
            
            Parameters:
            
            cutoff: double
            The new soft cutoff
            """
        self.potentials[0].set_soft_cutoff(cutoff)

    def get_potentials(self):
        """Returns the potentials stored in the :class:`~pysic.interactions.local.ProductPotential`.
            """
        return self.potentials

    def is_multiplier(self):
        """Returns a list of logical values specifying if the potentials are multipliers for a product.
            
            The leading potential is considered not to be a multiplier, while the rest are multipliers.
            Therefore the returned list is [False, True, ..., True], withlength equal to the length of
            potentials in the product.
            """
        return [False] + [True] * (len(self.potentials)-1)

    def accepts_target_list(self,targets):
        """Tests whether a list is suitable as a list of targets, i.e., symbols, tags, or indices and returns True or False accordingly.
            
        A list of targets should be of the format::
            
          targets = [[a, b], [c, d]]
            
        where the length of the sublists must equal the number of targets.
            
        It is not tested that the values contained in the list are valid.
            
        Parameters:
            
        targets: list of strings or integers
            a list whose format is checked
            """
        return self.potentials[0].accepts_target_list(targets)


class Potential(object):
    """Class for representing a potential.

    Several types of potentials can be defined by specifying the type of the potential
    as a keyword. The potentials contain a host of parameters and information on
    what types of particles they act on. To view a list of available potentials, use the method
    :meth:`~pysic.core.list_valid_potentials`.

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
    cutoff_margin: double
        the margin in which the potential is smoothly truncated to zero
    coordinator: :class:`~pysic.interactions.bondorder.Coordinator` object
        the coordinator defining a bond order factor for scaling the potential
    """

    def __init__(self,potential_type,symbols=None,tags=None,indices=None,parameters=None,cutoff=0.0,cutoff_margin=0.0,coordinator=None):
        if(is_valid_potential(potential_type)):
            self.symbols = None
            self.tags = None
            self.indices = None
            self.potential_type = potential_type
            self.cutoff = cutoff
            self.cutoff_margin = 0.0
            self.coordinator = None
            self.set_cutoff_margin(cutoff_margin)
            self.n_targets = number_of_targets(potential_type)
            self.names_of_params = names_of_parameters(potential_type)
            self.set_symbols(symbols)
            self.set_tags(tags)
            self.set_indices(indices)
            if(parameters is None):
                self.parameters = len(self.names_of_params)*[0.0]
            else:
                self.parameters = parameters
                
            if(len(self.names_of_params) != len(self.parameters)):
                raise InvalidPotentialError(
                    'The potential "{pot}" requires {num} parameters: {par}'.format(
                    pot=potential_type,
                    num=str(self.get_number_of_parameters()),
                    par=str(names_of_parameters(potential_type))
                    ) )
            self.set_coordinator(coordinator)
        else:
            raise InvalidPotentialError('There is no potential called "{pot}".'.format(pot=potential_type))


    def __eq__(self,other):
        try:
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
            if self.coordinator != other.coordinator:
                return False               
        except:
            return False
        
        return True

    def __ne__(self,other):
        return not self.__eq__(other)

    def __repr__(self):
        return ("Potential('{name}',symbols={symbs},"+ \
                "tags={tags},indices={inds},parameters={params}"+ \
                ",cutoff={cut},cutoff_margin={marg},"+ \
                "coordinator={coord})").format(name=self.potential_type,
                                               symbs=str(self.symbols),
                                               tags=str(self.tags),
                                               inds=str(self.indices),
                                               params=str(self.parameters),
                                               cut=str(self.cutoff),
                                               marg=str(self.cutoff_margin),
                                               coord=str(self.coordinator))
    
    
    def get_number_of_parameters(self):
        """Return the number of parameters the potential expects.
            """
        return number_of_parameters(self.potential_type)
    
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

    def get_coordinator(self):
        """Returns the Coordinator.
        """
        return self.coordinator

    def set_coordinator(self,coordinator):
        """Sets a new Coordinator.
            
        Parameters:
            
        coordinator: a :class:`~pysic.interactions.bondorder.Coordinator` object
            the Coordinator
        """
        self.coordinator = coordinator
        

    def set_symbols(self,symbols):
        """Sets the list of symbols to equal the given list.

        Parameters:

        symbols: list of strings
            list of element symbols on which the potential acts
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

    def set_parameters(self, values):
        """Sets the numeric values of all parameters.

        Equivalent to :meth:`~pysic.interactions.local.Potential.set_parameter_values`.

        Parameters:

        values: list of doubles
            list of values to be assigned to parameters
        """
        self.set_parameter_values(values)

    
    def set_parameter_values(self,values):
        """Sets the numeric values of all parameters.

        Parameters:

        values: list of doubles
            list of values to be assigned to parameters
        """
        if len(values) == self.get_number_of_parameters():
            self.parameters = values
        else:
            raise InvalidPotentialError("The potential '{pot}' takes {n_par} parameters, not {n_in}.".
                                        format(pot=self.potential_type,
                                               n_par=self.get_number_of_parameters(),
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
        For a detailed explanation on how to define a soft cutoff, see :meth:`~pysic.interactions.local.Potential.set_cutoff_margin`.
        
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
        hard cutoff value an :class:`~pysic.utility.error.InvalidPotentialError` is raised.

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
        :meth:`~pysic.interactions.local.Potential.set_cutoff_margin`.
        Note that actually the cutoff margin is recorded, so changing the
        hard cutoff (see :meth:`~pysic.interactions.local.Potential.set_cutoff`) will also affect the
        soft cutoff.

        Parameters:

        cutoff: double
            The new soft cutoff
        """
        self.set_cutoff_margin(self.cutoff-cutoff)


    def get_potentials(self):
        """Returns a list containing the potential itself.

            This is a method ensuring cross compatibility with the 
            :class:`~pysic.interactions.local.ProductPotential` class.
            """
        return [self]
    
    def is_multiplier(self):
        """Returns a list containing False.

            This is a method ensuring cross compatibility with the 
            :class:`~pysic.interactions.local.ProductPotential` class.
            """
        return [False]
    
    
    def describe(self):
        """Prints a short description of the potential using the method :meth:`~pysic.core.describe_potential`."""
        description_of_potential(self.potential_type,
                                 self.parameters,
                                 self.cutoff,
                                 self.symbols,
                                 self.tags,
                                 self.indices)

