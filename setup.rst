.. file:setup

.. file:obtaining

Getting Pysic
=============

Pysic is currently in development and not yet properly tested.
Nonetheless, the source code is available through `Gitorious`_ if you wish to
play with it. However, the code is provided with no warranty or support.

.. _Gitorious: https://gitorious.org/pysic/pysic

.. file:compiling

Compiling Pysic
===============

No installers of makefiles are currently provided with Pysic. Once Pysic is in
a production capable stage, some kind of an installation tool will be produced.

To run Pysic, you will need `Python 2.7`_ and the `numpy`_ and `scipy`_ modules.
For dynamic simulation, the `ASE`_ package is required. Some plotting tools
also require the `matplotlib`_ package, but Pysic will work without it.

To compile the
Fortran core, one needs a Fortran 90 compiler and `f2py`_ (the latter is
part of `numpy`_.) For compiling an MPI compatible version, an MPI-Fortran
compiler is needed.

When compiling Pysic, one must wrap the interface module :ref:`pysic_interface`
with `f2py`_ and compile all other Fortran source directly to ``.mod`` modules.
This can be achieved with::

   > f2py -m pysic_fortran -h pysic_fortran.pyf PyInterface.F90
   > f2py -c pysic_fortran.pyf Mersenne.F90 MPI.F90 Quaternions.F90 \ 
   > Utility.F90 Geometry.F90 Potentials.F90 Core.F90 PyInterface.F90

For MPI, use the conditional compiling flag ``-D MPI``.
This tells the Fortran compiler to include the MPI portions of the code::

   > f2py -m pysic_fortran -h pysic_fortran.pyf PyInterface.F90
   > f2py -c --fcompiler=gfortran --f90exec=mpif90 --f90flags="-D MPI" \
   > pysic_fortran.pyf Mersenne.F90 MPI.F90 Quaternions.F90 Utility.F90 \
   > Geometry.F90 Potentials.F90 Core.F90 PyInterface.F90

Above, the ``gfortran`` and ``mpif90`` compilers are used, but change the names to whichever compilers you wish to call.
For further information on compiling with `f2py`_, consult the `f2py manual`_.

 .. _f2py: http://www.scipy.org/F2py
 .. _f2py manual: http://cens.ioc.ee/projects/f2py2e/usersguide/
 .. _ASE: https://wiki.fysik.dtu.dk/ase/
 .. _numpy: http://numpy.scipy.org/
 .. _scipy: http://www.scipy.org/
 .. _matplotlib: http://matplotlib.sourceforge.net/
 .. _Python 2.7: http://www.python.org



.. file:resources

External resources
==================

Below is a list of tools one may find useful or even necessary when using Pysic:

- `Python`_: The Python language
- `Python documentation`_: The official documentation for Python
- `Python tutorial`_: The official tutorials for Python
- `F2py`_: Fortran to Python interface generator
- `F2py manual`_: The (old) official F2py manual
- `Gitorious`_: Version control and repository storing Pysic
- `ASE`_: Atomistic Simulation Environment (ASE)
- `ASE tutorials`_: The official ASE tutorials
- `NumPy`_: Numerical Python
- `SciPy`_: Scientific Python
- `Matplotlib`_: Matplotlib Python plotting library
- `iPython`_: iPython enhanced Python interpreter

.. _Python: http://www.python.org/
.. _Python documentation: http://docs.python.org/
.. _Python tutorial: http://docs.python.org/tutorial/index.html
.. _Gitorious: https://gitorious.org/pysic/pysic
.. _F2py: http://www.scipy.org/F2py
.. _F2py manual: http://cens.ioc.ee/projects/f2py2e/usersguide/
.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _ASE tutorials: https://wiki.fysik.dtu.dk/ase/tutorials/tutorials.html
.. _NumPy: http://numpy.scipy.org/
.. _SciPy: http://www.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/
.. _iPython: http://ipython.org/

