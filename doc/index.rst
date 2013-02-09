.. file:index

.. _index:




Pysic (Pythonic simulation code) 
==============================================

|logo|

`Teemu Hynninen <http://http://physics.aalto.fi/personnel/?id=209>`_ (2011-2013)

`Tampere University of Technology <http://www.tut.fi>`_ 
`Aalto University, Helsinki <http://www.aalto.fi>`_

Contact: `@Pysic_code <https://twitter.com/#!/Pysic_code>`_, teemu.hynninen aalto.fi, `@thynnine <https://twitter.com/#!/thynnine>`_

Pysic is a calculator incorporating various empirical pair and many-body
potentials in an object-based Python environment and user interface while
implementing an efficient numeric core written in Fortran. The immediate aim
of the Pysic project is to implement advanced variable charge potentials. 

Pysic is designed to interface with the `ASE`_ simulation environment.

.. _ASE: https://wiki.fysik.dtu.dk/ase/

The code is developed as part of the `Mordred`_ project.

.. _Mordred: http://webhotel2.tut.fi/fys/mordred/

Pysic is an open source code with emphasis on making the program simple to learn and control as well as the source code readable, extendable and conforming to good programming standards. The source code is freely available at the `github`_ repository. **Note the relocation since June 2012.**

.. _github: https://github.com/thynnine/pysic/
.. _Gitorious: https://gitorious.org/pysic/

Pysic is in development. Simulations can be run with it, but many parts 
of the program are not yet fully tested and so bugs must be expected.
Similarly, also this documentation is being constantly updated.

.. only:: html

    This document is also available as a `pdf`_.

    .. _pdf: pysic.pdf



Physical background
--------------------

.. toctree::
   :maxdepth: 3
   
   physics


Getting Pysic
-------------

.. toctree::
   :maxdepth: 3
   
   setup
   


Performing simulations with Pysic
----------------------------------

.. toctree::
   :maxdepth: 3
   
   running
   run examples
   

Structure and syntax in Pysic
-----------------------------

.. toctree::
   :maxdepth: 3
   
   syntax
   pysic
   pysic_utility
   pysic_fortran


Development of Pysic
--------------------

.. toctree::
   :maxdepth: 3

   version
   roadmap
   issues



.. |logo| image:: ../../../Graphics/mordred/pysic/pysic_logo.png
       :alt: 
       :height: 320

