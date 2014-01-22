#! /usr/bin/env python

import numpy as np
from scipy.optimize import fmin_slsqp
from pysic.utility.error import InvalidRelaxationError

class ChargeRelaxation:
    """A class for handling charge dynamics and relaxation.
        
        Pysic does not implement molecular dynamics or geometric
        optimization since they are handled by ASE.
        Conceptually, the structural dynamics of the system are
        properties of the atomic geometry and so it makes sense that
        they are handled by ASE, which defines the atomic structure
        in the first place, in the `ASE Atoms`_ class.
        
        On the other hand, charge dynamics are related to the electronic
        structure of the system. Since ASE is meant to use methods such as 
        density functional theory (DFT) in the calculators is employs, all
        electronic properties are left at the responsibilty of the
        calculator. This makes sense since in DFT the electron density is 
        needed for calculations of forces and energies.
        
        Pysic is not a DFT calculator and there is no electron density
        but the atomic charges can be made to develop dynamically. The
        ChargeRelaxation class handles these dynamics.
        
        .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
        
        Parameters:
        
        relaxation: string
            a keyword specifying the mode of charge relaxation
        calculator: :class:`~pysic.calculator.Pysic` object
            a Pysic calculator 
        parameters: list of doubles
            numeric values for parameters        
        atoms: `ASE Atoms`_ object
            The system whose charges are to be relaxed. 
            Note! The relaxation is always done using the atoms copy in 
            :class:`~pysic.calculator.Pysic`, but if the original structure needs to be
            updated as well, the relaxation algorithm must have access to it.
    """

    relaxation_modes = [ 'dynamic', 'potentiostat', 'optimize' ]
    """Names of the charge relaxation algorithms available. 
        
        These are keywords needed when creating the 
        :class:`~pysic.charges.relaxation.ChargeRelaxation` objects as type specifiers."""
    
    relaxation_parameters = { relaxation_modes[0] : ['n_steps', 'timestep', 'inertia', 'friction', 'tolerance'],
        relaxation_modes[1] : ['n_steps', 'timestep', 'inertia', 'potential'],
        relaxation_modes[2] : ['n_steps', 'tolerance'] }
    """Names of the parameters of the charge relaxation algorithms."""
    
    relaxation_parameter_descriptions = { relaxation_modes[0] : ['number of time steps of charge dynamics between molecular dynamics', 
            'time step of charge dynamics', 
            'fictional charge mass', 
            'friction coefficient for damped charge dynamics',
            'convergence tolerance'],
        relaxation_modes[1] : ['number of time steps of charge dynamics between molecular dynamics',
            'time step of charge dynamics',
            'fictional charge mass', 
            'external potential'],
        relaxation_modes[2] : ['maximum number of optimization steps',
            'convergence tolerance'] }
    """Short descriptions of the relaxation parameters."""
    
    def __init__(self, relaxation='dynamic', calculator=None, parameters=None, atoms=None):
        self.set_relaxation(relaxation)
        self.set_calculator(calculator)
        self.set_parameters(parameters)
        self.set_atoms(atoms)
        self.piped_relaxation = None
        self.charge_rates = None
        self.observers = []

    
    def __eq__(self,other):
        
        try:
            if self.relaxation != other.relaxation:
                return False
            if self.calculator != other.calculator:
                return False
            if self.parameters != other.parameters:
                return False
            if self.atoms != other.atoms:
                return False
            if self.piped_relaxation is None:
                if not other.piped_relaxation is None:
                    return False
            elif self.piped_relaxation != other.piped_relaxation:
                return False
        except:
            return False

        return True


    def __ne__(self,other):
        return not self == other

            
    def __repr__(self):
        return "ChargeRelaxation('{m}',calculator={c},parameters={p},atoms={a})".format(
            m=self.relaxation, 
            c=str(self.calculator), 
            p=str(self.parameters),
            a=str(self.atoms) )


    def set_relaxation_pipe(self, relaxation, pass_info=True):
        """Add another :class:`~pysic.charges.relaxation.ChargeRelaxation` to take place immediately after this relaxation.
        
            Sometimes it is most efficient to start the relaxation with a robust but 
            possibly slow algorithm and switch to a more advanced one once closer to
            equilibrium. On the other hand, when relaxation is needed between MD steps,
            this is not possible manually. Therefore, one can link several ChargeRelaxation
            objects together to automatically invoke a sequential relaxation procedure.
            
            You should only give one ChargeRelaxation for piping. If you want to
            pipe more than two relaxation sequences together, you should always link the
            additional relaxation to the currently final ChargeRelaxation object to
            create a nested pipeline::
            
              rel1 = ChargeRelaxation()
              rel2 = ChargeRelaxation()
              rel3 = ChargeRelaxation()
              
              rel1.set_relaxation_pipe(rel2)
              rel2.set_relaxation_pipe(rel3)
              
              rel1.charge_relaxation() # will do rel1 -> rel2 -> rel3
              
            Parameters:
            
            relaxation: :class:`~pysic.charges.relaxation.ChargeRelaxation` object
                the relaxation to be piped after this one
            pass_info: boolean
                If True, the calculator and atoms of this calculator will be automatically
                copied over to the piped calculator.
        """

        # Check that this relaxation does not appear in the pipeline to be attached.
        # Otherwise we would create an infinite loop.
        if relaxation is self:
            raise InvalidRelaxationError("Trying to create a looped relaxation pipe.")
        if not relaxation.piped_relaxation is None:
            pipeline = relaxation.get_relaxation_pipe()
            for rel in pipeline:
                if rel is self:
                    raise InvalidRelaxationError("Trying to create a looped relaxation pipe.")

        
        self.piped_relaxation = relaxation

        if pass_info:
            pipeline = self.get_relaxation_pipe()
            
            for rel in pipeline:
                rel.set_calculator(self.calculator)
                rel.set_atoms(self.atoms)

        if self.calculator != relaxation.calculator:
            warn("The piped ChargeRelaxation has a different calculator.",2)
        if not self.atoms is None:
            if relaxation.atoms is None:
                warn("No Atoms object linked to the piped ChargeRelaxation.",2)


    def clear_relaxation_pipe(self):
        """Equivalent to set_relaxation_pipe(None).
        """
        self.set_relaxation_pipe(None)
        

    def set_calculator(self, calculator, reciprocal=False):
        """Assigns a :class:`~pysic.calculator.Pysic` calculator.
            
            The calculator is necessary for calculation of electronegativities.
            It is also possible to automatically assign the charge relaxation
            method to the calculator by setting ``reciprocal = True``. 
            
            Note though
            that it does make a difference whether the calculator knows the
            charge relaxation or not: If the :class:`~pysic.calculator.Pysic` has a connection
            to the :class:`~pysic.charges.relaxation.ChargeRelaxation`, every time force or energy calculations
            are requested the charges are first relaxed by automatically invoking
            :meth:`~pysic.charges.relaxation.ChargeRelaxation.charge_relaxation`. If there is no link,
            it is up to the user to start the relaxation.
            
            Parameters:
            
            calculator: :class:`~pysic.calculator.Pysic` object
                a Pysic calculator
            reciprocal: logical
                if True, also the :class:`~pysic.charges.relaxation.ChargeRelaxation` is passed to the
                :class:`~pysic.calculator.Pysic` through :meth:`~pysic.calculator.Pysic.set_charge_relaxation`.
            """
        
        # remove possible return link from the calculator
        try:
            if self.calculator != calculator:
                self.calculator.set_charge_relaxation(None)
        except:
            pass
        self.calculator = calculator
        if reciprocal:
            calculator.set_charge_relaxation(self)


    
    def set_relaxation(self, relaxation):
        """Sets the relaxation method.
            
            The method also creates a dictionary of parameters initialized to 0.0
            by invoking :meth:`~pysic.charges.relaxation.ChargeRelaxation.initialize_parameters`.
            
            Parameters:
            
            relaxation: string
                a keyword specifying the mode of charge relaxation
            """
        if ChargeRelaxation.relaxation_modes.count(relaxation) == 1:
            self.relaxation = relaxation
            self.initialize_parameters()
        else:
            raise InvalidRelaxationError("no such relaxation mode "+relaxation)


    def initialize_parameters(self):
        """Creates a dictionary of parameters and initializes all values to 0.0.
        """        
        self.parameters = {}
        for param in ChargeRelaxation.relaxation_parameters[self.relaxation]:
            self.parameters[param] = 0.0


    def set_parameters(self, parameters):
        """Sets the numeric values for all parameters.
            
            Equivalent to :meth:`~pysic.charges.relaxation.ChargeRelaxation.set_parameter_values`
            
            Parameters:
            
            parameters: list of doubles
                list of values to be assigned to parameters
            """
        self.set_parameter_values(parameters)


    def set_parameter_values(self, parameters):
        """Sets the numeric values for all parameters.
            
            Use get_parameters() to see the correct order of values.
            
            Parameters:
                    
            parameters: list of doubles
                list of values to be assigned to parameters
            """
        if parameters == None:
            return
        if self.parameters == None:
            self.initialize_parameters()
        if len(parameters) != len(self.parameters):
            raise InvalidRelaxationError("The relaxation mode "+self.relaxation+" requires "+
                                         len(self.parameters)+" parameters.")
        for key, value in zip(self.parameters.keys(), parameters):
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

    
    def get_calculator(self):
        """Returns the :class:`~pysic.calculator.Pysic` calculator assigned to this :class:`~pysic.charges.relaxation.ChargeRelaxation`.
            """
        return self.calculator
    
    def get_relaxation(self):
        """Returns the keyword specifying the mode of relaxation.
            """
        return self.relaxation

    def get_relaxation_pipe(self, recursive=True):
        """Return the piped :class:`~pysic.charges.relaxation.ChargeRelaxation` objects.
        
        By default, a list is given containing all the ChargeRelaxation objects
        in the pipeline (or an empty list).
        
        If called with the argument False, only the ChargeRelaxation next-in-line
        will be returned (or None).
        
        Parameters:
        
        recursive: boolean
            If True, will return a list containing all the piped objects. Otherwise, only the next-in-line is returned.
        """
        if recursive:
            pipe = []
            if not self.piped_relaxation is None:
                pipe.append(self.piped_relaxation)
                while not pipe[-1].piped_relaxation is None:
                    pipe.append(pipe[-1].piped_relaxation)
            return pipe
        else:
            return self.piped_relaxation
    
    def get_parameters(self):
        """Returns a list containing the numeric values of the parameters.
            """
        return self.parameters
        
    
    def set_atoms(self, atoms, pass_to_calculator=False):
        """Lets the relaxation algorithm know the atomic structure to be updated.
            
            The relaxation algorithm always works with the structure stored in the
            :class:`~pysic.calculator.Pysic` calculator it knows. If ``pass_to_calculator = True``,
            this method also updates the structure known by the calculator. However,
            this is not the main purpose of letting the :class:`~pysic.charges.relaxation.ChargeRelaxation`
            know the structure -
            it is not even necessary that the structure known by the relaxation algorithm
            is the same as that known by the calculator.
            
            The structure given to the algorithm is the structure whose charges it 
            automatically updates after relaxing the charges in
            :meth:`~pysic.charges.relaxation.ChargeRelaxation.charge_relaxation`. In other words, if no
            structure is given, the relaxation will update the charges in the strucure
            known by :class:`~pysic.calculator.Pysic`, but this is always just a copy and so the
            original structure is left untouched.
            
            
            Parameters:
            
            atoms: `ASE Atoms`_ object
                The system whose charges are to be relaxed. 
                Note! The relaxation is always done using the atoms copy in 
                :class:`~pysic.calculator.Pysic`, but if the original structure needs to be
                updated as well, the relaxation algorithm must have access to it.
            pass_to_calculator: logical
                if True, the atoms are also set for the calculator via 
                :meth:`~pysic.Pysic.set_atoms`
            """
        self.atoms = atoms
        if pass_to_calculator:
            self.calculator.set_atoms(atoms)
    
    def get_atoms(self):
        """Returns the atoms object known by the algorithm.
            
        This is the `ASE Atoms`_ which will be automatically updated when
        charge relaxation is invoked.
        """
        return self.atoms
    
    
    def _energy_function(self, charges, orig_charge):

        try:
            self.calculator.get_atoms().set_initial_charges(charges)
        except:
            self.calculator.get_atoms().set_charges(charges)

        return self.calculator.get_potential_energy(skip_charge_relaxation=True)
    
    def _nabla_function(self, charges, orig_charge):
    
        try:
            self.calculator.get_atoms().set_initial_charges(charges)
        except:
            self.calculator.get_atoms().set_charges(charges)

        return -self.calculator.get_electronegativities()
    
    def _charge_constraint(self, charges, orig_charge):
        return np.sum(charges) - orig_charge

    def charge_relaxation(self):
        """Performs the charge relaxation.

            The relaxation is always performed on the system associated with
            the :class:`~pysic.calculator.Pysic` calculator joint with this :class:`~pysic.charges.relaxation.ChargeRelaxation`.
            The calculated equilibrium charges are returned as a numeric array.
            
            If an `ASE Atoms`_ structure is known by the
            :class:`~pysic.charges.relaxation.ChargeRelaxation` 
            (given through :meth:`~pysic.charges.relaxation.ChargeRelaxation.set_atoms`), the charges of
            the structure are updated according to the calculation result.
            If the structure is not known, the charges are updated in the
            structure stored in the :class:`~pysic.calculator.Pysic` calculator, but not in any other
            object. Since :class:`~pysic.calculator.Pysic` only stores a copy of the structure it 
            is given, the original `ASE Atoms`_ object will not be updated.
            """
        
        atoms = self.calculator.get_atoms()
        try:
            charges = atoms.get_initial_charges() # starting charges
        except:
            charges = atoms.get_charges() # starting charges
        
        if self.relaxation == ChargeRelaxation.relaxation_modes[0]: # damped dynamic relaxation
        
            dt = self.parameters['timestep']
            dt2 = dt*dt
            mq = self.parameters['inertia']
            inv_mq = 1.0/mq
            n_steps = self.parameters['n_steps']
            tolerance = self.parameters['tolerance']
            friction = self.parameters['friction']
            future_friction = 1.0 / (1.0 + 0.5 * friction * dt)

            error = 2.0 * tolerance
            step = 0
            
            if self.charge_rates is None:
                self.charge_rates = np.array( len(atoms)*[0.0] ) # starting "velocities" \dot{q}
            charge_forces = self.calculator.get_electronegativity_differences(atoms)
            
            while (error > tolerance and step < n_steps):
                charges += self.charge_rates * dt + \
                           (charge_forces - friction * self.charge_rates) * inv_mq * dt2
                self.charge_rates = (1.0 - 0.5 * friction * dt) * self.charge_rates + \
                            charge_forces * 0.5 * inv_mq * dt
                try:
                    atoms.set_initial_charges(charges)
                except:
                    atoms.set_charges(charges)

                charge_forces = self.calculator.get_electronegativity_differences(atoms)
                self.charge_rates = future_friction * \
                            ( self.charge_rates + charge_forces * 0.5 * inv_mq * dt )
                step += 1
                error = charge_forces.max()
                self.call_observers(step)

    
        elif self.relaxation == ChargeRelaxation.relaxation_modes[1]: # potentiostat
        
            dt = self.parameters['timestep']
            dt2 = dt*dt
            mq = self.parameters['inertia']
            inv_mq = 1.0/mq
            n_steps = self.parameters['n_steps']
            step = 0
            ext_potential = self.parameters['potential']

            if self.charge_rates is None:
                self.charge_rates = np.array( len(atoms)*[0.0] ) # starting "velocities" \dot{q}
            charge_forces = self.calculator.get_electronegativities(atoms)-ext_potential
            
            while (step < n_steps):
                charges += self.charge_rates * dt + charge_forces * inv_mq * dt2
                self.charge_rates = self.charge_rates + charge_forces * 0.5 * inv_mq * dt

                try:
                    atoms.set_initial_charges(charges)
                except:
                    atoms.set_charges(charges)

                charge_forces = self.calculator.get_electronegativities(atoms)-ext_potential
            
                self.charge_rates = self.charge_rates + charge_forces * 0.5 * inv_mq * dt

                step += 1
                self.call_observers(step)
        
        elif self.relaxation == ChargeRelaxation.relaxation_modes[2]: #
                
            charge0 = charges
            orig_charge = np.sum(charges),
            
            charges = fmin_slsqp(self._energy_function, charge0, fprime=self._nabla_function, eqcons=[self._charge_constraint], acc=self.parameters['tolerance'], iter=self.parameters['n_steps'], args=orig_charge)
        
        else:
            pass

        if self.atoms != None:
            #print "ch", charges
            #print "atoms", self.atoms
            try:
                self.atoms.set_initial_charges(charges)
            except:
                self.atoms.set_charges(charges)


        if not self.piped_relaxation is None:
            charges = self.piped_relaxation.charge_relaxation()


        return charges



    def add_observer(self, function, interval=1, *args, **kwargs):
        """Attach a callback function.
    
        Call the given function with the given arguments at constant
        intervals, after a given number of relaxation steps.
        
        Parameters:
        
        function: Python function
            the function to be called
        interval: integer
            the number of steps between callback
        args: string
            arguments for the callback function
        kwargs: string
            keyword arguments for the callback function
        """
    
        if not hasattr(function, '__call__'):
            function = function.write
        self.observers.append((function, interval, args, kwargs))


    def clear_observers(self):
        """Removes all observers.
        """
        self.observers = []
    
    def call_observers(self,step):
        """Calls all the attached callback functions.
        
        Parameters:
        
        step: integer
            the number of dynamics steps taken during the relaxation
        """
        for function, interval, args, kwargs in self.observers:
            try:
                if step % interval == 0:
                    function(*args, **kwargs)
            except:
                function(*args, **kwargs)