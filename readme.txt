
Welcome to Pysic
----------------

Pysic is a calculator incorporating various empirical pair and many-body potentials in an object-based Python environment and user interface while implementing an efficient numeric core written in Fortran. The immediate aim of the Pysic project is to implement advanced variable charge potentials.

The code is developed as part of the EU FP7 Mordred project.

Pysic is open source, available at http://github.com/thynnine/pysic/
The official documentation is found at http://thynnine.github.com/pysic/


Contents of the Pysic repository
--------------------------------

When you download Pysic, you should receive a package containing three directories [doc], [fortran], and [pysic]. The directory [doc] contains the source reStructuredText used for generating the documentation (with Sphinx). For readability, check the online version. The directories [fortran] and [pysic] contain the Fortran and Python source code, respectively. In addition to these folders, there should be a readme.txt file containing installation instructions and a Makefile for compiling with make.


Requirements
------------

To compile and run Pysic you need:

- Python 2.x
- Numpy (which includes f2py)
- ASE

Pysic is designed to interface with the Atomic Simulation Environment package 
(ASE, developed by CAMd, DTU), and although in theory ASE is not a prerequisite
for running Pysic, in practice you need it to do anything interesting.


Quick compilation
-----------------

To compile Pysic, go to the directory where you have the Makefile and folders *fortran* and 
*pysic*, and enter the command ``make <target>`` on command line, where ``<target>`` is for instance 
``serial`` or ``parallel``. A list of available targets are shown with the command ``make help``.

The default compilers listed in the Makefile are ``gfortran`` and ``mpif90``. 
You may need to edit the names of the compilers to match those available on your system.

By default, make creates the pysic Python module as a folder named *pysic* (or *pysic_debug* etc., 
depending on the type of compilation) in the folder *build*. To run Pysic, copy the folder *pysic* 
to whereever you wish to run your simulations and import it in Python with ``import pysic``.


Compilation under the hood
--------------------------

The Python part of Pysic does not have to be compiled separately, but the Fortran core must be.
In order to link the Fortran and Python sides together, you need to compile the code with the
*f2py* tool, which is part of the *Numpy* package.

The Fortran part is compiled with the commands::

  f2py -m pysic_fortran -h pysic_fortran.pyf PyInterface.F90
  f2py -c --fcompiler=gfortran --f90exec=mpif90 --f90flags="-D MPI" \
  pysic_fortran.pyf Mersenne.F90 MPI.F90 Quaternions.F90 Utility.F90 \
  Geometry.F90 Potentials.F90 Core.F90 PyInterface.F90

The first run creates the interface file ``pysic_fortran.pyf`` from the Fortran source file 
``PyInterface.F90``. The second run compiles the rest of the code in a native Fortran mode and
creates the Python accessible binary file ``pysic_fortran.so``. 
From this file, the subroutines of ``PyInterface.F90`` are callable as Python functions.

The created ``pysic_fortran.so`` appears in Python as a module ``pysic_fortran.pysic_interface``.
The package structure of pysic expects to find it at ``pysic.pysic_fortran.pysic_interface``, i.e.,
the ``.so`` file should be in the folder ``pysic`` (the root folder of the Python package).
The Makefile should take care of putting it in the right place, but if you decide to compile the
Fortran part manually, keep in mind that the file needs to be moved there.
