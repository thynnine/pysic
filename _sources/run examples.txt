.. file:run examples

.. _run examples:



Examples
--------

Next, we go through some basic examples on how to use Pysic, starting from launching Python, finishing with relatively advanced scripting techniques.


.. file:python setup

.. _python setup:



First steps with Python
_______________________

Once you have installed Python - chances are it is already preinstalled in your system - you can launch the Python interpreter with the command ``python`` on the command line::
 
  > python
  Python 2.7.2 (#1, Feb 28 2010, 00:02:06)
  Type "help", "copyright", "credits" or "license" for more information.
  >>>

This will start up Python in interactive mode and you can start writing commands. For instance::

  >>> a = 1+1
  >>> a
  2
  >>> b = 'hello'
  >>> c = b + ' world!'
  >>> c
  'hello world!'

Another way to run Python is by feeding a source file to the interpreter. Say you have the file ``first_step.py`` containing::

  left_first = False
  if left_first:
     print "Left, right, left!"
  else:
     print "Right, left, right!"

This can be run with either::

  > python first_step.py
  Right, left, right!

or simply::

  > ./first_step.py
  Right, left, right!

assuming you have execution permissions for the script.
You can also pass command line arguments to a Python script when launching one.

That's it, basically! 
In the following examples it is shown how Pysic is set up within Python
and how simulations can be run.
However, as Python is a powerful language, you will certainly get the most
out of Pysic by learning some of the basic features of Python.
Taking a look at the official `Python tutorial`_ is a good place to start.

.. _Python tutorial: http://docs.python.org/tutorial/index.html

.. file:importing pysic

.. _importing pysic:



Importing Pysic to Python
_________________________

Modules are accessed in Python by importing them with the ``import`` command. To get access to Pysic, simply write::

  >>> import pysic

Then, you can access all the functionality in the module :mod:`pysic`::

  >>> pysic_calculator = pysic.Pysic()
  >>> pysic.get_names_of_parameters('LJ')
  ['epsilon', 'sigma']

The Fortran interface module is imported with :mod:`pysic` as ``pysic.pf.pysic_interface``, but it is strongly recommended you do not touch the Fortran core directly.

Altogether, when the module :mod:`pysic` is imported, it imports the following non-standard modules::

  >>> import pysic.pysic_fortran as pf
  >>> import numpy as np
  >>> import ase.calculators.neighborlist as nbl

and defines the functions and classes as documented in the syntax description of :mod:`pysic`. In addition, the :func:`start_potentials` and :func:`start_bond_order_factors` routines of the Fortran core are automatically invoked in order to initialize the descriptors in the core (see :ref:`potentials`). In an MPI environment, the MPI initialization routines :func:`start_mpi` are also called from the Fortran core. Finally, the random number generator is initialized in the Fortran core by :func:`start_rng`.


.. file:minimal setup

.. _minimal setup:




Minimal example of running Pysic
________________________________

Here is an example of setting up a Pysic calculator for `ASE`_::

    >>> from ase import Atoms
    >>> import pysic
    >>> system = Atoms('He2', [[0.0, 0.0, 0.0], [0.0, 0.0, 3.0]])  
    >>> calc = pysic.Pysic()  
    >>> system.set_calculator(calc)  
    >>> physics = pysic.Potential('LJ', cutoff = 10.0)  
    >>> physics.set_symbols(['He', 'He'])  
    >>> physics.set_parameter_value('epsilon', 0.1) 
    >>> physics.set_parameter_value('sigma', 2.5)
    >>> calc.add_potential(physics) 


The example above creates a system of two helium atoms interacting via a Lennard-Jones
potential 

.. math::

   V(r) = \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^{6} \right]
   
   \varepsilon = 0.1

   \sigma = 2.5

In the code above, ``system`` is an `ASE Atoms`_ object containing the structure of the system to be calculated -
two He atoms in this case. The object ``calc`` is an instance of :class:`~pysic.calculator.Pysic`, the `ASE calculator`_ class
defined by Pysic. The interactions governing the system are defined by the ``physics`` object, which is an instance of the
:class:`~pysic.interactions.local.Potential` class of Pysic.

Now, the potential energy of the system and the forces acting on the atoms can be calculated with::

  >>> system.get_potential_energy()
  -0.022274132189576905
  >>> system.get_forces()
  array([[ 0.        ,  0.        ,  0.02211693],
         [ 0.        ,  0.        , -0.02211693]])

These commands are addressed to the ``system`` object, but under the hood ``system`` asks ``calc``, i.e., Pysic,
to do the actual calculations. In order to evaluate the requested quantities,
Pysic needs the parameters contained in ``physics``.

A more compact way to create the calculator would be::

    >>> physics = pysic.Potential('LJ', 
    ...                           cutoff = 10.0, 
    ...                           symbols = ['He','He'],
    ...                           parameters = [0.1, 2.5])
    >>> calc = pysic.Pysic(system, physics)

Setting up more complicated interactions works similarly, as is shown in later examples.

.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/
.. _ASE calculator: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html

.. file:scripting - generate systems

.. _scripting - generate systems:



.. file:running MD

.. _running MD:



Running molecular dynamics
__________________________

This example shows how to set up a dynamic simulation with `ASE`_ (also see the `ASE MD tutorial`_).

First set up a MgO crystal::

  from ase.lattice.compounds import Rocksalt

  # Set up a crystal
  size = 4    
  lattice_constant = 4.212
  system = Rocksalt(size=(size,size,size), symbol=("Mg", "O"), 
                    latticeconstant=lattice_constant, pbc=True)
  # Set charges for the atoms
  for i in range(len(system)):
      if system[i].get_symbol() == "Mg":
          system[i].set_charge(2)
      elif system[i].get_symbol() == "O":
          system[i].set_charge(-2)

Create a Pysic calculator. We set up pairwise interactions with the :ref:`Buckingham potential` and a Coulomb interaction::

  import pysic

  # Set up a calculator
  calc = pysic.Pysic()

  # Mg-O pair potential
  pot_mgo = pysic.Potential('Buckingham', symbols=['Mg','O'], cutoff=8.0, cutoff_margin=1.0)
  pot_mgo.set_parameter_value('A', 1284.38)
  pot_mgo.set_parameter_value('C', 0.0)
  pot_mgo.set_parameter_value('sigma', 0.2997)
  calc.add_potential(pot_mgo)

  # O-O pair potential
  pot_oo = pysic.Potential('Buckingham', symbols=['O','O'], cutoff=10.0, cutoff_margin=1.0)
  pot_oo.set_parameter_value('A', 9574.96)
  pot_oo.set_parameter_value('C', 288474.00)
  pot_oo.set_parameter_value('sigma', 0.2192)
  calc.add_potential(pot_oo)

  # Coulomb interaction
  ewald = pysic.CoulombSummation()
  ewald.set_parameter_value('epsilon',0.00552635) # epsilon_0 in e**2/(eV A)
  ewald.set_parameter_value('k_cutoff',0.7)     # one should always test the Ewald parameters for convergence
  ewald.set_parameter_value('real_cutoff',10.0) # and optimal speed - long cutoffs can substantially increase
  ewald.set_parameter_value('sigma',1.4)        # needed cpu time while short cutoffs give incorrect results
  calc.set_coulomb_summation(ewald)
  
  system.set_calculator(calc)


Now, let us set up the dynamic simulation.

We initialize the velocity distribution::

  from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
  from ase import units

  # Set the momenta corresponding to T=300K
  MaxwellBoltzmannDistribution(system, 300*units.kB)

And define the dynamics as a ``VelocityVerlet`` object::

  from ase.md.verlet import VelocityVerlet

  # We want to run MD with constant energy using the VelocityVerlet algorithm.
  dyn = VelocityVerlet(system, 5*units.fs)  # 5 fs time step.

In order to record the `ASE trajectory`_ and print information during simulation,
we also define observers and attach them to the dynamics. 
In a similar fashion we can collect whatever information we need::

  # Function to print the potential, kinetic and total energy.
  def print_energy(a=system):
      epot = a.get_potential_energy() / len(a)
      ekin = a.get_kinetic_energy() / len(a)
      print ("Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  Etot = %.3feV" %
             (epot, ekin, ekin/(1.5*units.kB), epot+ekin))
  dyn.attach(print_energy, interval=100) # print after every 100 steps

  # Save a trajectory
  from ase.io.trajectory import PickleTrajectory

  traj = PickleTrajectory("example.traj", "w", system) # write a trajectory to the file 'example.traj'
  dyn.attach(traj.write, interval=10) # save every 10 step

Finally, we run the dynamics::

  dyn.run(1000) # run for 1000 steps
  traj.close()  # close the trajectory recording


.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _ASE MD tutorial: https://wiki.fysik.dtu.dk/ase/tutorials/md/md.html
.. _ASE trajectory: https://wiki.fysik.dtu.dk/ase/ase/trajectory.html


.. file:scripting - collect data

.. _scripting - collect data:



.. file:potential targeting

.. _potential targeting:



.. file:using bond order potentials

.. _using bond order potentials:



.. file:running MPI

.. _running MPI:



.. file:plotting utilities

.. _plotting utilities:



.. file:hybrid calculation

.. _hybrid calculation:



A hybrid calculation
__________________________


The following python script demonstrates the Python interface for creating a simple hybrid QM/MM simulation with Pysic. The potentials and calculator parameters used in this example are not physically justified and serve only to illustrate the usage::

	#! /usr/bin/env python
	from ase.structure import molecule
	from pysic import *
	from gpaw import GPAW

	# Create the atomic system as an ASE Atoms object
	ethane = molecule(’C2H6’, cell=(6, 6, 6))
	ethane.center()

	# Create the hybrid calculator
	hybrid_calc = HybridCalculator()

	# Setup the primary system
	gpaw_calc = GPAW(h=0.25, txt=None)
	PS = SubSystem("primary", indices=(0, 2, 3, 4), calculator=gpaw_calc)
	hybrid_calc.add_subsystem(PS)

	# Setup the secondary system
	1pysic_calc = Pysic()
	pysic_calc.add_potential(Potential(’LJ’, symbols=[[’H’, ’C’]], parameters=[0.05, 2.5], cutoff=5))
	SS = SubSystem("secondary", indices="remaining", calculator=pysic_calc)
	hybrid_calc.add_subsystem(SS)
	
	# Setup the interaction between the subsystems
	interaction = Interaction("primary", "secondary")
	interaction.set_potentials(Potential(’LJ’, symbols=[[’C’, ’C’]], parameters=[0.05, 3], cutoff=5))
	interaction.add_hydrogen_links((0, 1), CHL=0.5)
	hybrid_calc.add_interaction(interaction)

	# Do some calculations
	ethane.set_calculator(hybrid_calc)
	ethane.get_potential_energy()
	ethane.get_forces()
	
	# Output results
	hybrid_calc.print_energy_summary()
	hybrid_calc.print_force_summary()
	hybrid_calc.print_time_summary()

The code is now split into sections that are further elaborated:

* First we define the entire structure that we want to analyze. It is a regular ASE Atoms object::

	ethane = molecule(’C2H6’, cell=(6, 6, 6))
	ethane.center()

* A HybridCalculator object is then created. It is the core class in Pysic’s hybrid calculations. The end-user may treat is like any other ASE compatible calculator::

	hybrid_calc = HybridCalculator()

* A SubSystem object is created for each desired subsystem. A subsystem is defined by giving it a unique name, the indices of the atoms in the subsystem and the calculator which is to be used for it. The newly created SubSystem object is then passed to the HybridCalculator::

	# Setup the primary system
	gpaw_calc = GPAW(h=0.25, txt=None)
	PS = SubSystem("primary", indices=(0, 2, 3, 4), calculator=gpaw_calc)
	hybrid_calc.add_subsystem(PS)

	# Setup the secondary system
	pysic_calc = Pysic()
	pysic_calc.add_potential(Potential(’LJ’, symbols=[[’H’, ’C’]], parameters=[0.05, 2.5], cutoff=5))
	SS = SubSystem("secondary", indices="remaining", calculator=pysic_calc)
	hybrid_calc.add_subsystem(SS)

* Interaction potentials and hydrogen link atoms between any two subsystems are set up by creating an Interaction object. The names of the interacting subsystems are given in the constructor, the first being the primary system, and the latter being the secondary system. Potentials and link atoms are created through simple function calls. Finally the Interaction object is passed to the HybridCalculator::

	interaction = Interaction("primary", "secondary")
	interaction.set_potentials(Potential(’LJ’, symbols=[[’C’, ’C’]], parameters=[0.05, 3], cutoff=5))
	interaction.add_hydrogen_links((0, 1), CHL=0.5)
	hybrid_calc.add_interaction(interaction)

* When the subsystems and interactions are successfully configured, the HybridCalculator can be used like any regular ASE calculator to calculate the potential energy and forces in the system::

	ethane.set_calculator(hybrid_calc)
	ethane.get_potential_energy()
	ethane.get_forces()

* The HybridCalculator provides useful tools for inspecting the energies, calculation times and forces related to the different subsystems and interactions::

	hybrid_calc.print_energy_summary()
	hybrid_calc.print_force_summary()
	hybrid_calc.print_time_summary()



.. file:screencasts

.. _screencasts:




Screencasts
-----------

We have also prepared screencasts demonstrating the use of Pysic.



.. file:Pysic basics

.. _Pysic basics:



Pysic basics
____________

This is a demonstration of a very simple simulation setup while running Pysic in an interactive Python session.

.. raw:: html

  <iframe width="853" height="480" src="http://www.youtube-nocookie.com/embed/f3VKiclBvQY?rel=0" frameborder="0" allowfullscreen></iframe>

`YouTube link <http://www.youtube.com/watch?v=f3VKiclBvQY>`_

.. file:Pysic in a script

.. _Pysic in a script:



Pysic in a script
_________________

This is a demonstration of running Pysic with a script, including an MPI parallel run.

.. raw:: html

  <iframe width="853" height="480" src="http://www.youtube-nocookie.com/embed/cz9sL5Zh1xg?rel=0" frameborder="0" allowfullscreen></iframe>

`YouTube link (high definition) <http://www.youtube.com/watch?v=cz9sL5Zh1xg>`_
