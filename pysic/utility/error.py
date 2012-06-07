#! /usr/bin/env python


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
            return self.messsage + " \n the Potential: " + str(self.potential)



class InvalidCoordinatorError(Exception):
    """An error raised when an invalid coordination calculator is about to be created or used.

    Parameters:
    
    message: string
        information describing why the error occurred
    coordinator: :class:`~pysic.Coordinator`
        the errorneous coordinator
    """
    def __init__(self,message='',coordinator=None):
        self.message = message
        self.coordinator = coordinator

    def __str__(self):
        if(self.coordinator == None):
            return self.message
        else:
            return self.messsage + " \n  the Coordinator: " + str(self.coordinator)




class InvalidParametersError(Exception):
    """An error raised when an invalid set of parameters is about to be created.

    Parameters:

    message: string
        information describing why the error occurred
    params: :class:`~pysic.BondOrderParameters`
        the errorneous parameters
    """
    def __init__(self,message='',params=None):
        self.message = message
        self.params = params

    def __str__(self):
        if(self.params == None):
            return self.message
        else:
            return self.messsage + " \n  the Parameters: " + str(self.params)



class InvalidSummationError(Exception):
    """An error raised when an invalid coulomb summation is about to be created.

    Parameters:

    message: string
        information describing why the error occurred
    params: :class:`~pysic.CoulombSummation`
        the errorneous summation
    """
    def __init__(self,message='',summer=None):
        self.message = message
        self.summer = summer

    def __str__(self):
        if(self.summer == None):
            return self.message
        else:
            return self.messsage + " \n  the Summation: " + str(self.summer)


class InvalidRelaxationError(Exception):
    """An error raised when an invalid charge relaxation is about to be created.
        
        Parameters:
        
        message: string
            information describing why the error occurred
        params: :class:`~pysic.ChargeRelaxation`
            the errorneous parameters
        """
    def __init__(self,message='',relaxation=None):
        self.message = message
        self.relaxation = relaxation
    
    def __str__(self):
        if(self.params == None):
            return self.message
        else:
            return self.messsage + " \n  the Relaxation: " + str(self.relaxation)




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

