.. file:setup

.. _setup:



.. file:obtaining

.. _obtaining:



Getting Pysic
=============

Pysic is currently in development and not yet properly tested.
Nonetheless, the source code is available through `github`_ if you wish to
try it out. The code is provided with no warranty or support.

.. _github: https://github.com/thynnine/pysic/
.. _Gitorious: https://gitorious.org/pysic/pysic


Contents of the Pysic repository
--------------------------------

When you download Pysic, you should receive a package containing three directories *doc*, *fortran*, and *pysic*. The directory *doc* contains the reStructuredText (``.rst``) source files used for generating the documentation (with `Sphinx <http://sphinx.pocoo.org/>`_). For readability, check the `online version <http://thynnine.github.com/pysic/>`_. The directories *fortran* and *pysic* contain the Fortran and Python source code, respectively. In addition to these folders, there should be a ``readme.txt`` file containing installation instructions and a ``Makefile`` for compiling with ``make``.



.. file:compiling

.. _compiling:



Compiling Pysic
===============

This section gives brief instructions on how to compile Pysic.

Requirements
------------

To run Pysic, you will need `Python 2.x`_ and the `numpy`_ and `scipy`_ modules.
For dynamic simulation, the `ASE`_ package is required. Some plotting tools
also require the `matplotlib`_ package, but Pysic will work without it.

To compile the
Fortran core, one needs a Fortran 90 compiler and `f2py`_ (the latter is
part of `numpy`_.) For compiling an MPI compatible version, an MPI-Fortran
compiler is needed.




Quick compilation
-----------------

To compile Pysic, go to the directory where you have the ``Makefile`` and folders *fortran* and 
*pysic*, and enter the command ``make <target>`` on command line, where ``<target>`` is for instance 
``serial`` or ``parallel``. A list of available targets are shown with the command ``make help``.

The default compilers listed in the ``Makefile`` are ``gfortran`` and ``mpif90``. 
You may need to edit the names of the compilers to match those available on your system.

By default, ``make`` creates the pysic Python module as a folder named *pysic* (or *pysic_debug* etc., 
depending on the type of compilation) in the folder *build*. To run Pysic, copy the folder *pysic* 
to wherever you wish to run your simulations and import it in Python with ``import pysic``.


Compilation under the hood
--------------------------

The Python part of Pysic does not have to be compiled separately, but the Fortran core must be.
In order to link the Fortran and Python sides together, you need to compile the code with the
`f2py`_ tool, which is part of the `numpy`_ package.

The Fortran part is compiled with the commands::

  f2py -m pysic_fortran -h pysic_fortran.pyf PyInterface.F90
  f2py -c --fcompiler=gfortran --f90exec=mpif90 --f90flags="-D MPI" \
  pysic_fortran.pyf Mersenne.F90 MPI.F90 Quaternions.F90 Utility.F90 \
  Geometry.F90 Potentials.F90 Core.F90 PyInterface.F90

The first run creates the interface file ``pysic_fortran.pyf`` from the Fortran source file 
``PyInterface.F90``. The second run compiles the rest of the code in a native Fortran mode and
creates the Python accessible binary file ``pysic_fortran.so``. From this file, the subroutines of ``PyInterface.F90`` are callable as Python functions.

Above, the ``gfortran`` and ``mpif90`` compilers are used, but change the names to whichever compilers you wish to call. The flag ``--f90flags="-D MPI"``
invokes the preprocessor which picks the MPI part of the source code. To compile a serial version,
leave this flag out. 

The created ``pysic_fortran.so`` appears in Python as a module ``pysic_fortran.pysic_interface``.
The package structure of pysic expects to find it at ``pysic.pysic_fortran.pysic_interface``, i.e.,
the ``.so`` file should be in the folder *pysic* (the root folder of the Python package).
The ``Makefile`` should take care of putting it in the right place, but if you decide to compile the
Fortran part manually, keep in mind that the file needs to be moved there.

For further information on compiling with `f2py`_, consult the `f2py manual`_.

 .. _f2py: http://www.scipy.org/F2py
 .. _f2py manual: http://cens.ioc.ee/projects/f2py2e/usersguide/
 .. _ASE: https://wiki.fysik.dtu.dk/ase/
 .. _numpy: http://numpy.scipy.org/
 .. _scipy: http://www.scipy.org/
 .. _matplotlib: http://matplotlib.sourceforge.net/
 .. _Python 2.x: http://www.python.org



.. file:resources

.. _resources:



External resources
==================

Below is a list of tools one may find useful or even necessary when using Pysic:

- `Python`_: The Python language
- `Python documentation`_: The official documentation for Python
- `Python tutorial`_: The official tutorials for Python
- `F2py`_: Fortran to Python interface generator
- `F2py manual`_: The (old) official F2py manual
- `github`_: Version control and repository storing Pysic
- `ASE`_: Atomistic Simulation Environment (ASE)
- `ASE tutorials`_: The official ASE tutorials
- `NumPy`_: Numerical Python
- `SciPy`_: Scientific Python
- `Matplotlib`_: Matplotlib Python plotting library
- `iPython`_: iPython enhanced Python interpreter

.. _Python: http://www.python.org/
.. _Python documentation: http://docs.python.org/
.. _Python tutorial: http://docs.python.org/tutorial/index.html
.. _github: https://github.com/thynnine/pysic/
.. _Gitorious: https://gitorious.org/pysic/pysic
.. _F2py: http://www.scipy.org/F2py
.. _F2py manual: http://cens.ioc.ee/projects/f2py2e/usersguide/
.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _ASE tutorials: https://wiki.fysik.dtu.dk/ase/tutorials/tutorials.html
.. _NumPy: http://numpy.scipy.org/
.. _SciPy: http://www.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/
.. _iPython: http://ipython.org/

