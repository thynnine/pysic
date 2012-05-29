.. file:version

.. _version:



Version history
===============

This is a list of the main updates in the different versions of Pysic.


Version 0.4.3
-------------

- Restructured the interaction evaluation loops in the Fortran core (:ref:`potentials`)
- Added support for 4-body potentials
- Added the :ref:`dihedral angle potential` potential
- Bug fix: fixed an issue where the cutoff of a bond order factor could overwrite a longer cutoff a potential
- Bug fix: fixed an indexing error in evaluation of 3-body interaction which gave to incorrect forces


Version 0.4.2
-------------

- Some code restructuring in the Fortran core (:ref:`potentials`)
- Added the :ref:`Buckingham potential` potential
- Bug fix: fixed an issue with periodic boundaries in :class:`~pysic.FastNeighborList`
- Bug fix: fixed an issue with special parameter values in Tersoff bond order factor evaluation


Version 0.4.1
-------------

- Implemented an order :math:`\mathcal{O}(n)` neighbor finding algorithm in Fortran (see :class:`pysic.FastNeighborList`)



Version 0.4
-----------

- Implemented the Ewald summation of :math:`\frac{1}{r}` potentials (see :class:`pysic.CoulombSummation`)
- The framework allows for the addition of other summation methods later on, but for now only standard Ewald is available


Version 0.3
-----------

- Implemented framework for charge relaxation (see :class:`pysic.ChargeRelaxation`)
- Implemented the :ref:`damped dynamics` charge relaxation algorithm.
- Implemented the :ref:`charge exponential potential` potential.


Version 0.2
-----------

- Implemented bond order factors (see :class:`pysic.Coordinator` and :class:`pysic.BondOrderParameters`) for scaling of potential energy according to local bond structure.
- Implemented a more robust method for tracking the status of the Fortran core (see :class:`pysic.CoreMirror`). This makes it less likely that wrong results are produced due to the changes in the user interface not propagating to the core.


Version 0.1
-----------

- First publicly available version
- Python interface

  * :mod:`pysic`
  * :class:`pysic.Pysic`
  * :class:`pysic.Potential`
  * :mod:`pysic_utility`

- Framework for handling pair- and three-body potentials
- ASE compatibility

  * :meth:`pysic.Pysic.get_forces`
  * :meth:`pysic.Pysic.get_potential_energy`

