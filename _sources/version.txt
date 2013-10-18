.. file:version

.. _version:



Version history
===============

This is a list of the main updates in the different versions of Pysic.
Descriptions of all incremental changes can be found in the `commit log <https://github.com/thynnine/pysic/commits/master>`_ 
(see also the `develop version <https://github.com/thynnine/pysic/commits/develop>`_).

The version number is available in Pysic as the variable::

  pysic.version


Version 0.4.8
----------------

- Fixed compatibility with syntax changes in ASE 3.7.0.
- Added piping of :class:`~pysic.charge.relaxation.ChargeRelaxation` objects.
- Added a potentiostat as a :class:`~pysic.charge.relaxation.ChargeRelaxation` type.
- Added utility functions :meth:`~pysic.calculator.FastNeighborList.get_neighbor_distances`, :func:`~pysic.utility.convenience.expand_symbols_string`


Version 0.4.7
---------------

- Added wrapping of several :ref:`ase calculators`


Version 0.4.6
--------------

- Code optimization: updates in the memory allocation and computation routines in the Fortran core for better performance.
- Added a warning system (:class:`pysic.utility.error.Warning`)
- Added support for pairwise bond order factors (in addition to per-atom factors), see :ref:`atomic and pairwise factors`.
- Numerous bug fixes, see the github commit log for details.


Version 0.4.5
--------------

- Added the :class:`~pysic.interactions.compound.CompoundPotential` class for representing complicated potentials.
- Added the :class:`~pysic.interactions.suttonchen.SuttonChenPotential` class describing the Sutton-Chen potential.
- Added the :ref:`shifted power potential`
- Interface enhancements. E.g., it is possible to pass lists of potentials to methods like :meth:`pysic.calculator.Pysic.add_potential`, or pass :class:`~pysic.interactions.local.ProductPotential` objects to other :class:`~pysic.interactions.local.ProductPotential` objects.
- Efficiency fixes. E.g., when bond order factors are created in the core, it is checked that duplicates are not created.

Version 0.4.4
-------------

- Added the ability to calculate products of local potentials (see :class:`pysic.interactions.local.ProductPotential`)
- Added the :ref:`charged-pair potential`
- Separated the old charged exponential potential to :ref:`exponential potential` and :ref:`charge exponential potential`
- Changed the :ref:`bond bending potential` to allow more general expressions.
- Added the :ref:`absolute charged-pair potential`


Version 0.4.3
-------------

- Major restructuring of the Python source code
- Provided a Makefile for compiling
- Added calculation of the stress tensor with the method :meth:`pysic.calculator.Pysic.get_stress`
- Added the :ref:`tabulated potential`
- Added the :ref:`tabulated scaling function`
- Added the :ref:`tabulated bond order factor`
- Bug fix: Fixed an issue with core initialization where changing the size of the supercell would lead to a conflict in neighbor list updating (the list update was tried before the cell update but failed due to the cell having been changed).
- Bug fix: Fixed an issue with the parallel neighbor list building algorithm which did not properly broadcast the calculated lists to all cpus.

Version 0.4.2
-------------

- Restructured the interaction evaluation loops in the Fortran core (:ref:`potentials`)
- Added support for 4-body potentials
- Added the :ref:`dihedral angle potential`
- Added the :ref:`Buckingham potential`
- Added the :ref:`power decay potential`
- Added the :ref:`power decay bond order factor`
- Added the :ref:`square root scaling function`
- Bug fix: fixed a memory issue in Ewald summation :class:`~pysic.interactions.coulomb.CoulombSummation`
- Bug fix: fixed an issue with periodic boundaries in :class:`~pysic.calculator.FastNeighborList`
- Bug fix: fixed an issue with special parameter values in Tersoff bond order factor evaluation
- Bug fix: fixed an issue where the cutoff of a bond order factor could overwrite a longer cutoff a potential
- Bug fix: fixed an indexing error in evaluation of 3-body interaction which gave to incorrect forces
- Bug fix: fixed and indexing error in neighbor offsets in :class:`~pysic.calculator.FastNeighborList`

Version 0.4.1
-------------

- Implemented an order :math:`\mathcal{O}(n)` neighbor finding algorithm in Fortran (see :class:`pysic.calculator.FastNeighborList`)



Version 0.4
-----------

- Implemented the Ewald summation of :math:`\frac{1}{r}` potentials (see :class:`pysic.interactions.coulomb.CoulombSummation`)
- The framework allows for the addition of other summation methods later on, but for now only standard Ewald is available


Version 0.3
-----------

- Implemented framework for charge relaxation (see :class:`pysic.charges.relaxation.ChargeRelaxation`)
- Implemented the :ref:`damped dynamics` charge relaxation algorithm.
- Implemented the :ref:`charge exponential potential` potential.


Version 0.2
-----------

- Implemented bond order factors (see :class:`pysic.interactions.bondorder.Coordinator` and :class:`pysic.interactions.bondorder.BondOrderParameters`) for scaling of potential energy according to local bond structure.
- Implemented a more robust method for tracking the status of the Fortran core (see :class:`pysic.core.CoreMirror`). This makes it less likely that wrong results are produced due to the changes in the user interface not propagating to the core.


Version 0.1
-----------

- First publicly available version
- Python interface

  * :mod:`pysic`
  * :class:`pysic.calculator.Pysic`
  * :class:`pysic.interactions.local.Potential`
  * ``pysic_utility``

- Framework for handling pair- and three-body potentials
- ASE compatibility

  * :meth:`pysic.calculator.Pysic.get_forces`
  * :meth:`pysic.calculator.Pysic.get_potential_energy`

