#! /usr/bin/env python

from pysic.utility.error import InvalidSummationError


class CoulombSummation:
    """Class for representing a collection of parameters for evaluating Coulomb potentials.

    Summing :math:`1/r` potentials in periodic systems requires more 
    advanced techniques than just direct summation of pair interactions.
    The starndard method for evaluating these kinds of potentials is through
    Ewald summation, where the long range part of the potential is evaluated
    in reciprocal space.
    
    Instances of this class are used for wrapping the parameters controlling
    the summations. Passing such an instance to the :class:`~pysic.calculator.Pysic`
    calculator activates the evaluation of Coulomb interactions.
    
    Currently, only Ewald summation is available as a calculation method.
        
    Parameters:
        
    method: string
        keyword specifying the method of summation
    parameters: list of doubles
        numeric values of summation parameters
    scaler: list of doubles
        numeric values for scaling the atomic charges in summation
    """

    summation_modes = [ 'ewald' ]
    """Names of the summation methods. These are keywords used for setting up the summation algorithms."""
    summation_parameters = { summation_modes[0] : ['real_cutoff',
                                                   'k_cutoff',
                                                   'sigma',
                                                   'epsilon'] }
    """Names of the parameters of the summation algorithm."""
    summation_parameter_descriptions = { summation_modes[0] : ['real space cutoff radius',
                                                   'reciprocal space cutoff radius',
                                                   'ewald summation split parameter',
                                                   'vacuum permittivity'] }
    """Short descriptions of the parameters of the summation algorithm."""
    
    def __init__(self,method='ewald',parameters=None,scaler=None):
        self.set_summation(method)
        self.set_parameters(parameters)
        self.set_scaling_factors(scaler)
    
    
    def __eq__(self,other):
        try:
            if self.method != other.method:
                return False
            if self.parameters != other.parameters:
                return False
            if self.scaler != other.scaler:
                return False
        except:
            return False
                    
        return True
    
    
    def __ne__(self,other):
        return not self == other


    def __repr__(self):
        return "CoulombSummation(method='{m}',parameters={p},scaler={s})".format(m=self.method,
                                                                                 p=str(self.parameters),
                                                                                 s=str(self.scaler))


    def set_summation(self,method):
        """Sets the summation method.
            
            The method also creates a dictionary of parameters initialized to 0.0
            by invoking :meth:`~pysic.interactions.coulomb.CoulombSummation.initialize_parameters`.
            
            Parameters:
            
            method: string
                a keyword specifying the mode of summation
            """
        if CoulombSummation.summation_modes.count(method) == 1:
            self.method = method
            self.initialize_parameters()
        else:
            raise InvalidSummationError("no such summation mode "+method)


    def initialize_parameters(self):
        """Creates a dictionary of parameters and initializes all values to 0.0.
        """        
        self.parameters = {}
        for param in CoulombSummation.summation_parameters[self.method]:
            self.parameters[param] = 0.0

                

    def set_parameters(self, parameters):
        """Sets the numeric values for all parameters.
        
        Equivalent to :meth:`~pysic.interactions.coulomb..CoulombSummation.set_parameter_values`
        
        Parameters:
        
        parameters: list of doubles
            list of values to be assigned to parameters
        """
        self.set_parameter_values(parameters)
    
    
    def set_parameter_values(self, parameters):
        """Sets the numeric values for all parameters.                    
            
            Parameters:
            
            parameters: list of doubles
                list of values to be assigned to parameters
            """
        if parameters == None:
            return
        if self.parameters == None:
            self.initialize_parameters()
        if len(parameters) != len(self.parameters):
            raise InvalidSummationError("The summation mode "+self.method+" requires "+
                                         len(self.parameters)+" parameters.")
        for key, value in zip(CoulombSummation.summation_parameters[self.method], parameters):
            self.parameters[key] = value
        
    def set_parameter_value(self, parameter_name, value):
        """Sets a given parameter to the desired value.
        
        Parameters:
        
        parameter_name: string
            name of the parameter
        value: double
            the new value of the parameter
        """
        self.parameters[parameter_name] = value


    def set_scaling_factors(self,scaler):
        """Set the list of scaling parameters for atomic charges.

            Parameters:
            
            scaler: list of doubles
                the list of scaling factors
            """
        self.scaler = scaler
                

    def get_scaling_factors(self):
        """Returns the list of scaling parameters for atomic charges.
            """
        return self.scaler
                
                
    def get_summation(self):
        """Returns the mode of summation."""
        return self.method

    def get_parameters(self):
        """Returns a list containing the numeric values of the parameters.
        """
        return self.parameters

    def get_realspace_cutoff(self):
        """Returns the real space cutoff.
            """
        if self.method == CoulombSummation.summation_modes[0]: #ewald
            return self.parameters['real_cutoff']
        return 0.0


def estimate_ewald_parameters(real_cutoff=10.0, accuracy='normal'):
    """Returns a tuple containing a good initial guess for Ewald parameters in the order `real_cutoff, k_cutoff, sigma, epsilon`.
    
    The returned values are (real_cutoff, k_cutoff, sigma, epsilon).
    Epsilon is always 0.00552635 :math:`\\varepsilon_0` in units of :math:`\\frac{e^2}{eV A}`
    Real cutoff can be given by the user and should be short enough to make the
    real space summation efficient (note that this affects neighborlisting and thus
    also other real space summations). Sigma and k_cutoff are determined by simple scaling rules
    where :math:`\\sigma \\sim r_\\mathrm{cut}` and :math:`k_\\mathrm{cut} \\sim \\sigma^{-1}`.
    The scaling constants are determined by the accuracy setting.
    
    Note that the given parameters are not analysed in any way - they are only a first guess.
    You should always test the parameters for accuracy and speed before production simulations.
    
    Parameters:
    
    real_cutoff: double
        the real space cutoff to be used - it should be shorter than the size of the simulation box
    accuracy: string
        either 'low', 'normal', 'high', 'real' or 'reciprocal'
    """

    epsilon = 0.00552635
    if(accuracy == 'low'):
        sigma = real_cutoff/2.0
        k_cutoff = 2.0/sigma
    elif(accuracy == 'normal'):
        sigma = real_cutoff/3.0
        k_cutoff = 3.0/sigma
    elif(accuracy == 'high'):
        sigma = real_cutoff/4.0
        k_cutoff = 4.0/sigma
    elif(accuracy == 'real'):
        sigma = real_cutoff/3.0
        k_cutoff = 2.0/sigma
    elif(accuracy == 'reciprocal'):
        sigma = real_cutoff/2.0
        k_cutoff = 3.0/sigma

    return (real_cutoff, k_cutoff, sigma, epsilon)
