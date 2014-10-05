
=============
Documentation
=============

Hybrid QM/MM calculations in Pysic
==================================

Pysic provides a framework for creating and running hybrid QM/MM simulations. Within this framework it is possible to calculate the potential energy and forces in an atomistic configuration which has multiple subsystems and interactions between them. Any external calculator with an ASE interface can be assigned to a subsystem or the classical potentials provided by Pysic can be used. The QM/MM implementation in Pysic uses the mechanical embedding scheme with hydrogen link atoms. It is possible to enable any Pysic-supported interaction potentials between the subsystems, with special functions provided for the easy use of Coulomb and COMB interactions.

The following figure demonstrates the simplified workflow comparison of a regular ASE calculation vs. Pysic's hybrid QM/MM calculation.

.. image:: workflow.png
	:scale: 70 %
	:align: center

Mechanical Embedding
--------------------

The electrostatic interaction between subsystems is treated with the mechanical embedding scheme. The Coulomb interaction is thus calculated in the MM level by assigning point charges for the atoms in the  primary and secondary system. The charges in the secondary system are static and can be obtained from experimental data or theoretical calculations. Rather than also parametrizing the charges in the primary system, in Pysic they can be dynamically calculated by using the electron density provided by the DFT calculator. The positive charge of the nuclei can be accurately modeled as an atom-centered point charge. The effect of the negative electron charge density is approximated also as a single atom-centered point charge. The value of this point charge is calculated by integrating the electron charge density inside a certain volume. Primarily Pysic uses the Bader partitioning scheme to determine these volumes, but a spherical partitioning scheme is also available.

Hydrogen links in Pysic
-----------------------

The covalent bonds between subsystems are treated with the link atom approach, i.e. hydrogen atoms are used to cap the bonds in the primary QM system. The hydrogen link atoms are placed on the line connecting the two atoms that form the bond. The exact position of the link atoms in this line can be controlled with the CHL parameter, given when defining the links in :meth:`~pysic.interaction.add_hydrogen_links`. The link atoms automatically keep their correct position when running dynamics or geometry optimization.

Example
-------

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

HybridCalculator class
======================

The HybridCalculator provides a regular ASE calculator interface, but internally combines results from multiple different calculators and interactions between them. Currently you can use it for potential energy and force calculations (stress calculations not yet implemented). Like with any ASE calculator, you can also run molecular dynamics and local geometry optimizations with it, but global energy optimization schemes that require a stress implementation will naturally not work due to missing stress implementation.

The following methods form the basis for creating and performing hybrid simulations with HybridCalculator:

	- :meth:`~pysic.hybridcalculator.HybridCalculator.add_subsystem`
	- :meth:`~pysic.hybridcalculator.HybridCalculator.add_interaction`
	- :meth:`~pysic.hybridcalculator.HybridCalculator.get_potential_energy`
	- :meth:`~pysic.hybridcalculator.HybridCalculator.get_forces`
	- :meth:`~pysic.hybridcalculator.HybridCalculator.print_energy_summary`
	- :meth:`~pysic.hybridcalculator.HybridCalculator.print_force_summary`
	- :meth:`~pysic.hybridcalculator.HybridCalculator.print_time_summary`
	- :meth:`~pysic.hybridcalculator.HybridCalculator.print_interaction_charge_summary`
	- :meth:`~pysic.hybridcalculator.HybridCalculator.view_subsystems`

Full documentation of the HybridCalculator class
------------------------------------------------

.. currentmodule:: pysic.hybridcalculator
.. autoclass:: HybridCalculator
   :members:
   :undoc-members:

Subsystem class
===============

The SubSystem class works as an interface for defining subsystems in hybrid simulations. Through it you can define the atoms that belong to the subsystem and the calculator that is used for it. SubSystem objects also provides special options for QM subsystems: you can setup dynamical charge calculation with :meth:`~pysic.subsystem.SubSystem.enable_charge_calculation` and cell optimization with :meth:`~pysic.subsystem.SubSystem.enable_cell_optimization`.

You can define the atoms in the subsystem as a list of indices, with a tag or with a special string: "remaining", which means all the atoms that are not yet assigned to a subsystem.

Full documentation of the Subsystem class
-----------------------------------------

.. currentmodule:: pysic.subsystem
.. autoclass:: SubSystem
   :members:
   :undoc-members:

SubSystemInternal class
=======================

This class is meant for internal use only. 

Full documentation of the SubSystemInternal class
-------------------------------------------------

.. currentmodule:: pysic.subsystem
.. autoclass:: SubSystemInternal
   :members:
   :undoc-members:

Interaction class
=================

The Interaction class works as an interface for defining interactions and link atoms between subsystems in hybrid calulations. With the methods :meth:`~pysic.interaction.Interaction.add_potential`, :meth:`~pysic.interaction.Interaction.set_potentials`, :meth:`~pysic.interaction.Interaction.enable_coulomb_potential` and :meth:`~pysic.interaction.Interaction.enable_comb_potential` you can define the Pysic Potentials that are present between the subsystems. With :meth:`~pysic.interaction.Interaction.add_hydrogen_links` you can setup the hydrogen link atoms.

Full documentation of the Interaction class
-------------------------------------------

.. currentmodule:: pysic.interaction
.. autoclass:: Interaction
   :members:
   :undoc-members:

InteractionInternal class
=========================

This class is meant for internal use only. 

Full documentation of the InteractionInternal class
---------------------------------------------------

.. currentmodule:: pysic.interaction
.. autoclass:: InteractionInternal
   :members:
   :undoc-members:

