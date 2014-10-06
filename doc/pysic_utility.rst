.. file:pysic_utility

.. _pysic_utility:



.. module:: pysic.utility

=====================
Pysic Utility modules
=====================

Modules for providing utility tools and parameters.

.. file:plot utility

.. _plot utility:



Plot
----

The plot utility defines a group of tools for exploring and plotting energy and force landscapes. It uses the `matplotlib library <http://matplotlib.org>`_.

.. automodule:: pysic.utility.plot
   :members:

.. file:outliers utility

.. _outliers utility:



Outliers
----------

This module contains tools for statistical analysis of atomic structures. It can be used for calculating the bond length and angle distributions in a given structure according to bond rules given by the user. Log-likelihood distributions are generated based on these distributions, which can be further grouped or used to create false color representations of the structure.

This is an example of how to create the distributions::

 import pysic
 from pysic.utility.outliers import *
 from ase.io import read

 def do_everything(inputfile, outputfile, radii, periodic_directions, cell_lengths):
    """Stitch everything together.
    """

    system = read(inputfile)
    system.set_pbc(periodic_directions)
    system.set_cell(cell_lengths)
    structure = Structure(system)
    
    structure.add_bond(['Si','Si'], radii['Si']*2*0.6)
    structure.add_bond(['Si','O'], (radii['Si']+radii['O'])*0.6)
    structure.add_bond(['O','O'], radii['O']*2*0.6)
    structure.add_bond(['Si','H'], (radii['Si']+radii['H'])*0.6)
    structure.add_bond(['H','O'], (radii['H']+radii['O'])*0.6)
    structure.create_neighbor_lists()
    angles = structure.get_all_angles()
    distances = structure.get_all_distances()
    max_n = np.shape(system)[0]
    
    for a in angles:                                         
        print "a", a.center_index, a.type1, a.type2, a.type3, a.value
    
    for d in distances:                                         
        print "d", d.primary_index, d.type1, d.type2, d.value
        
    angle_distribs, dist_distribs = get_distributions(angles, distances, radii)
    a_logls, d_logls = get_log_likelihoods(angles, distances, angle_distribs,
                                           dist_distribs, max_n)
    write_to_file(outputfile, cell_lengths, system.get_chemical_symbols(), 
                    system.get_positions(), a_logls, d_logls)
    

 inputfile           = 'sio2.xyz'
 outputfile          = 'sio2_analysis.txt'
 periodic_directions = [True, True, True]
 cell_lengths = [14.835, 14.835, 14.835]

 radii = {'O' : 1.52,
          'Si': 2.10}


 # All set? Run.
 do_everything(inputfile, outputfile, radii, periodic_directions, cell_lengths)

The atomic system is first read from an xyz file using ASE read method.
It is then translated to a Structure class, which the outliers tools use,
and bonding rules are created based on element types. Once bonding has been defined,
bond lengths, bond angles, and their distributions are calculated with the
outliers tools. The results are printed in a file in this example, but it would
of course be possible to access them directly in the script for further analysis.


Structure class
_________________________________

.. autoclass:: pysic.utility.outliers.Structure
   :members:


Distance class
_________________________________

.. autoclass:: pysic.utility.outliers.Distance
   :members:


Angle class
_________________________________

.. autoclass:: pysic.utility.outliers.Angle
   :members:


Outliers module
_________________________________


.. automodule:: pysic.utility.outliers
   :members:

.. file:mpi utility

.. _mpi utility:



MPI
----

The mpi utility defines a group of tools for parallel computations.

.. automodule:: pysic.utility.mpi
   :members:

.. file:archive utility

.. _archive utility:



Archive
---------

The archive utility defines a group of functions for archiving and retrieving simulation data in the hdf5 format. This requires `hdf5 <http://www.hdfgroup.org/HDF5/>`_  and the `h5py Python library <http://www.h5py.org/>`_.

.. automodule:: pysic.utility.archive
   :members:
   :undoc-members:

.. file:convenience utility

.. _convenience utility:



Convenience
-----------

The convenience utility defines functions to ease the handling of complicated data structures.

.. automodule:: pysic.utility.convenience
   :members:

.. file:geometry utility

.. _geometry utility:



Geometry
--------

The geometry utility defines tools for geometric operations.

.. automodule:: pysic.utility.geometry
   :members:

.. file:error utility

.. _error utility:



Error 
----------------------------------

The error utility defines a group of intrinsic errors to describe situations where
one tries to use or set up the calculator with errorneous or insufficient
information.

The module also defines the :class:`~pysic.utility.error.Warning` class and related routines
for displaying warnings for the user upon suspicious but non-critical behavior.

.. automodule:: pysic.utility.error
   :members:





.. file:debug utility

.. _debug utility:



Debug
----------------------------------

The debug utility defines debugging tools.

.. automodule:: pysic.utility.debug
   :members:





.. file:f2py utility

.. _f2py utility:



F2py
----------------------------------

The f2py utility defines tools for the Python-Fortran interfacing.

.. automodule:: pysic.utility.f2py
   :members:




