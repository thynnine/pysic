.. file:version

Version history
===============

This is a list of the main updates in the different versions of Pysic.


Version 0.3
-----------

- Implemented framework for charge relaxation (see :class:`pysic.ChargeRelaxation`)
- Implemented the :ref:`damped dynamic charge relaxation` charge relaxation algorithm.
- Implemented the :ref:`charge dependent exponential potential` potential.


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

