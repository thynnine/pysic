.. file:pysic class

.. _pysic class:



.. file:pysic class - description

.. _pysic class - description:



.. module:: pysic.calculator

========================
Pysic class
========================

This class provides an interface both to the `ASE`_ environment and the
Fortran core of Pysic through the :mod:`~pysic_fortran` module.
In `ASE`_ terms, Pysic is a `calculator`_, i.e., an object that calculates
forces and energies of given atomistic structures provided as an
`ASE Atoms`_ object.

.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html  
.. _calculator: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html

Pysic is not a fixed potential calculator, where the interactions are determined
solely based on the elements present in the system. Instead, Pysic allows and requires 
the user to specify the interactions governing the system. This is done by providing
the calculator with one or several :class:`~pysic.interactions.local.Potential` objects.


.. file:ASE calculators

.. _ASE calculators:



ASE calculators
-------------------

`ASE calculators <https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html>`_ are a class representing the physics of the atomistic system. Or put in another way, a calculator is the mathematical machine containing the algorithms and routines for calculating, for instance, the total potential energy of the atomistic system represented as an `ASE atoms`_ object (hence the name). As ASE itself contains only some very simple calculators, in practice the calculators are often just Python interfaces to other atomistic simulators. Pysic too is such a calculator.

.. _ASE atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

Pysic implements classical potentials, where the interactions between atomic pairs, triplets etc. are represented via simple mathematical functions or a :ref:`tabulated potential`. These can be combined freely using the :ref:`Potential class <potential class>` and its modifiers and extensions such as the :ref:`Coordinator <coordinator class>`, :ref:`Bond order parameter <bondorderparameters class>`, :ref:`Product potential <productpotential class>`, :ref:`Compound Potential <compoundpotential class>`, and :ref:`Coulomb summation <coulombsummation class>` classes.

Wrapping calculators
------------------------

Although Pysic can be used as a calculator on its own, it is also possible to use it as a wrapper for other ASE calculators.
By doing so, one may combine different simulators and even different descriptions in one calculator, and all the wrapped calculators will be executed when the Pysic calculator needs to calculate energies or forces. This is done simply by adding the external ASE calculators to Pysic using the :meth:`~pysic.calculator.Pysic.add_calculator` method.

As a simple example, one can use Pysic to add simple components like springs in an ab initio calculation.
This example uses the `Gpaw <https://wiki.fysik.dtu.dk/gpaw/>`_ code.::

    import pysic
    from ase import Atoms  
    from gpaw import GPAW  

    # create the system
    system = Atoms( ... )

    # create the pysic calculator
    calc = pysic.Pysic()
    
    # create a classical spring potential between atoms '0' and '1'
    spr = pysic.Potential('spring', indices=[0,1], parameters=[10.0, 3.0], cutoff=10.0, cutoff_margin=1.0)
    calc.add_potential(spr)

    # create an ab initio calculator
    dft = GPAW(h=0.18, nbands=1, xc='PBE', txt='wrapped.out')

    # wrap the calculators
    calc.add_calculator(dft)

    # calculate the energy
    system.set_calculator(calc)
    print system.get_potential_energy()

It is also possible to wrap several instances of Pysic together, but usually there should be no need to do this. 
It is more efficient to just add all the potentials in one Pysic instance. However, as an example, here is how it works::

  import pysic
  import ase
  
  # create some test system
  system = ase.Atoms('COH',[[0.0,0.0,0.0], [2.0,0.0,0.0], [2.0,2.0,2.0]])

  # create two calculators
  calc1 = pysic.Pysic()
  calc2 = pysic.Pysic()

  # potential 1
  pot1 = pysic.Potential('LJ',symbols=['C','O'],parameters=[1,1],cutoff=10.0)
  # potential 2
  pot2 = pysic.Potential('LJ',symbols=['O','H'],parameters=[1,1],cutoff=10.0)

  # add potentials to respective calculators
  calc1.add_potential(pot1)
  calc2.add_potential(pot2)
  system.set_calculator(calc1)

  # write energies and forces for the calculators both separately and combined
  print "Only potential 1"
  print system.get_potential_energy()
  print system.get_forces()
  print "\n"

  print "Only potential 2"
  print calc2.get_potential_energy(system)
  print calc2.get_forces(system)
  print "\n"

  print "Both potentials by calculator wrapping"
  calc1.add_calculator(calc2)
  print system.get_potential_energy()
  print system.get_forces()
  print "\n"

  calc1.remove_calculator(calc2)
  calc1.add_potential(pot2)
  print "Both potentials in one calculator"
  print system.get_potential_energy()
  print system.get_forces()
  print "\n"

This will output::

  Only potential 1
  -0.015380859375
  [[ 0.04541016  0.          0.        ]
   [-0.04541016  0.          0.        ]
   [ 0.          0.          0.        ]]


  Only potential 2
  -0.00194931030273
  [[ 0.          0.          0.        ]
   [ 0.          0.00291824  0.00291824]
   [ 0.         -0.00291824 -0.00291824]]


  Both potentials by calculator wrapping
  -0.0173301696777
  [[ 0.04541016  0.          0.        ]
   [-0.04541016  0.00291824  0.00291824]
   [ 0.         -0.00291824 -0.00291824]]


  Both potentials in one calculator
  -0.0173301696777
  [[ 0.04541016  0.          0.        ]
   [-0.04541016  0.00291824  0.00291824]
   [ 0.         -0.00291824 -0.00291824]]

Notice how in the third printout the calculators `calc1` and `calc2` are wrapped together, and the result is their sum. In the fourth case, `calc2` is discarded and instead `pot2` is added directly to `calc1`. The result is the same, but the latter option is better performance wise, since wrapping the calculators requires both to be separately managed during calculation.


.. file:energy evaluation

.. _energy evaluation:



.. file:force evaluation

.. _force evaluation:



.. file:electronegativity evaluation

.. _electronegativity evaluation:



.. file:stress evaluation

.. _stress evaluation:



.. file:role of other classes

.. _role of other classes:



.. file:pysic class - autogenerated

.. _pysic class - autogenerated:






List of methods
---------------
  
Below is a list of methods in :class:`~pysic.calculator.Pysic`, grouped according to
the type of functionality.
  
Structure handling
__________________
  
- :meth:`~pysic.calculator.Pysic.create_neighbor_lists` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.get_atoms`
- :meth:`~pysic.calculator.Pysic.get_neighbor_lists`
- :meth:`~pysic.calculator.Pysic.neighbor_lists_expanded` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.set_atoms`
  
Potential handling
__________________
  
- :meth:`~pysic.calculator.Pysic.add_potential`
- :meth:`~pysic.calculator.Pysic.get_individual_cutoffs`
- :meth:`~pysic.calculator.Pysic.get_potentials`
- :meth:`~pysic.calculator.Pysic.remove_potential`
- :meth:`~pysic.calculator.Pysic.set_potentials`

Coulomb handling
________________

- :meth:`~pysic.calculator.Pysic.get_coulomb_summation`
- :meth:`~pysic.calculator.Pysic.set_coulomb_summation`

Charge relaxation handling
__________________________

- :meth:`~pysic.calculator.Pysic.get_charge_relaxation`
- :meth:`~pysic.calculator.Pysic.set_charge_relaxation`  

Calculator handling
_______________________

- :meth:`~pysic.calculator.Pysic.add_calculator`
- :meth:`~pysic.calculator.Pysic.get_calculators`
- :meth:`~pysic.calculator.Pysic.remove_calculator`


Calculations
____________
  
- :meth:`~pysic.calculator.Pysic.calculate_electronegativities` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.calculate_energy` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.calculate_forces` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.calculate_stress` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.calculation_required` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.get_electronegativities`
- :meth:`~pysic.calculator.Pysic.get_electronegativity_differences`
- :meth:`~pysic.calculator.Pysic.get_forces`
- :meth:`~pysic.calculator.Pysic.get_numerical_bond_order_gradient` (for testing)
- :meth:`~pysic.calculator.Pysic.get_numerical_energy_gradient` (for testing)
- :meth:`~pysic.calculator.Pysic.get_numerical_electronegativity` (for testing)
- :meth:`~pysic.calculator.Pysic.get_potential_energy`
- :meth:`~pysic.calculator.Pysic.get_stress`

  
Core
____
  
- :data:`~pysic.calculator.Pysic.core`
- :meth:`~pysic.calculator.Pysic.core_initialization_is_forced`
- :meth:`~pysic.calculator.Pysic.force_core_initialization`
- :meth:`~pysic.calculator.Pysic.initialize_fortran_core`
- :meth:`~pysic.calculator.Pysic.set_core`
- :meth:`~pysic.calculator.Pysic.update_core_charges` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.update_core_coordinates` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.update_core_coulomb` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.update_core_neighbor_lists` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.update_core_potential_lists` (meant for internal use)
- :meth:`~pysic.calculator.Pysic.update_core_supercell` (meant for internal use)




Full documentation of the Pysic class
-------------------------------------

.. currentmodule:: pysic.calculator
.. autoclass:: Pysic
   :members:
   :undoc-members:

