.. file:syntax_forewords

.. _syntax:

Syntax of Pysic
===============

Pysic can be used with basic functionality with just a few commands. To give access to all the functionality in Pysic, the full API of Pysic is documented in the following section. All the classes and methods including their arguments are explained in detail. Also the types of potentials and bond order factors available are documented, including their mathematical descriptions and the keywords needed for access.

In addition to the Python documentation, also the variables, types, and routines in the Fortran core are listed with comments and explanations. Besides the interface module :ref:`pysic_interface`, this part of the code is not accessible through Python. Thus, you need not know what the core contains. However, if you plan to modify the core, studying also the Fortran documentation is useful.

.. only:: html

 Modules
 -------

 - :mod:`pysic`
 - :mod:`pysic_utility`
 - :mod:`pysic_fortran`

 Classes
 -------

 - :class:`pysic.Pysic`
 - :class:`pysic.Potential`
 - :class:`pysic.Coordinator`
 - :class:`pysic.BondOrderParameters`
 - :class:`pysic.CoreMirror`

