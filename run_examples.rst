.. file:run_examples

.. _examples:

Examples
--------

Next, we go through some basic examples on how to use Pysic, starting from launching Python, finishing with relatively advanced scripting techniques.


.. file:python setup

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

.. file:importing_pysic

Importing Pysic to Python
_________________________

Modules are accessed in Python by importing them with the ``import`` command. To get access to Pysic, simply write::

  >>> import pysic

Then, you can access all the functionality in the module :mod:`pysic`::

  >>> pysic_calculator = pysic.Pysic()
  >>> pysic.get_names_of_parameters('LJ')
  ['epsilon', 'sigma']

The Fortran interface module is imported with :mod:`pysic` as ``pysic.pf.pysic_interface``, but it is strongly recommended you do not touch the Fortran core directly.

The utility module :mod:`pysic_utility` is also imported as ``pysic.pu``, but it is also convenient to import it separately::

  >>> import pysic_utility
  >>> pysic_utility.plot_energy_on_plane(0,system,[[1,0,0],[0,1,0]],[10,10])

A shorthand can be introduced through the ``as`` keyword::

  >>> import pysic_utility as pu
  >>> pu.plot_energy_on_plane(0,system,[[1,0,0],[0,1,0]],[10,10])


Altogether, when the module :mod:`pysic` is imported, it imports the following non-standard modules::

  >>> import pysic_fortran as pf
  >>> import pysic_utility as pu
  >>> import numpy as np
  >>> import ase.calculators.neighborlist as nbl

and defines the functions and classes as documented in the syntax description of :mod:`pysic`. In addition, the :func:`start_potentials` and :func:`start_bond_order_factors` routines of the Fortran core are automatically invoked in order to initialize the descriptors in the core (see :ref:`potentials`). In an MPI environment, the MPI initialization routines :func:`start_mpi` are also called from the Fortran core. Finally, the random number generator is initialized in the Fortran core by :func:`start_rng`.


.. file:minimal setup


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
two He atoms in this case. The object ``calc`` is an instance of :class:`~pysic.Pysic`, the `ASE calculator`_ class
defined by Pysic. The interactions governing the system are defined by the ``physics`` object, which is an instance of the
:class:`~pysic.Potential` class of Pysic.

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

.. file:scripting: generate systems

.. file:running MD

.. file:scripting: collect data

.. file:Potential targets

.. file:bond order potentials

.. file:MPI

.. file:plotting utilities

.. file:syntax
