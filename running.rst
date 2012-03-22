.. file:running

.. file:run_forewords

Running Pysic
-------------

Once you have Python and Pysic working, it's time to learn how to use them. Next we will go through the basic concepts and ideas behind running simulations with Pysic. A collection of examples will demonstrate how to set up simulations in practice.


.. file:python_general

Why Python?
-----------

`Python`_ is a programming language with a simple human readable syntax yet powerful features and
an extensive library. Python is an interpreted language meaning it does not need to be compiled.
One can simply feed the source code to an interpreter which reads and interprets it during run.
The Python interpreter can also perform calculations, read and write files, render graphics,
etc. making Python a powerful tool for scripting. Python is also an object based language
enabling sophisticated :ref:`object oriented programming <objects>`.

In computational physics codes, the most common method of performing the calculations is to have
the program read input files to extract the necessary parameters, perform the simulation based on
the given input, and print output based on the results. The generated data is then analysed using 
separate tools. In some cases, more common in commercial programs, the user can control the
program through a graphical user interface.

In Pysic, another approach is chosen. Pysic is not a "black box" simulator that turns input data to
output data. Instead, Pysic is a *Python module*. In essence, it is a library of tools one can use
within Python to perform complicated calculations. Instead of writing an input file - often
a complicated and error prone collection of numbers - the user needs
to build the simulation in Python. Pysic can then be used in evaluating the energies,
forces and other physical quantities of the given system. Python syntax is in general simple and
simplicity and user friendlyness has also been a key goal in the design of Pysic syntax.

Since Python is a programming language, having Pysic be a part of the language makes it straightforward
to write scripts that generate the physical system to be studied and also extract the needed information
from the simulations.
Instead of generating enormous amounts of data which would then have to be fed to some other analysis
program, one can precisely control what kind of data should be produced in the simulations and even analyse
the results simultaneously as the simulation is run. Even controlling the simulation based on the
produced data is not only possible but easy.

The downside of Pysic being a Python module instead of an independent program is that one has to know the
basics of how to run Python. `Python documentation`_ is the best resource for getting started, but some
simple first step instructions and :ref:`examples` are also given in this document.

.. only:: html 
        
      |!| See the official `Python documentation`_ and `tutorials`_.


.. _Python: http://www.python.org
.. _Python documentation: http://docs.python.org/
.. _tutorials: http://docs.python.org/tutorial/index.html
.. |!| image:: ../../Graphics/mordred/pysic/exclamation.png
       :alt: < ! >
       :height: 32

.. file:ase_general

The Atomic Simulation Environment
---------------------------------

The Atomistic Simulation Environment (`ASE`_ [#]_) is a simulation tool developed originally at `CAMd`_, DTU (Technical University of Denmark). The package defines a full atomistic simulation environment for Python, including utilities such as a `graphical user interface`_. Like Pysic, ASE is a collection of Python modules. These modules allow the user to construct atomistic systems, run `molecular dynamics`_ or `geometry optimization`_, do calculations and analysis, etc. 

.. only:: html 
        
      |!| See the official `ASE tutorials`_.

ASE is easily extendable, and in fact, the common way to use ASE is join it with an external calculator which produces the needed physical quantities based on the structures defined by ASE. Pysic is such an extension. In the terms used in ASE, Pysic is a `calculator`_. The main role of calculators in ASE is to calculate forces and energies, and this is also what Pysic does. In addition, Pysic incorporates variable charge potentials with dynamic charge equilibration routines. Since ASE does not contain such routines, they are also handled by Pysic, making it more than just an extension of ASE. ASE does contain efficient dynamics and optimization routines including `constrained algorithms`_ and `nudged elastic band`_ as well as a choice of `thermostats`_, and so Pysic does not implement any such functionality.

The central object in ASE, also used by Pysic, is the `Atoms`_ class. This class defines the complete system geometry including atomic species, coordinates, momenta, supercell, and boundary conditions. Pysic interprets instances of this class as the simulation geometry.

.. [#] ASE: Comput. Sci. Eng., Vol. 4, 56-66, 2002; https://wiki.fysik.dtu.dk/ase/

.. _CAMd: http://www.camd.dtu.dk/English/Software.aspx
.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html 
.. _graphical user interface: https://wiki.fysik.dtu.dk/ase/ase/gui/gui.html#module-gui
.. _molecular dynamics: https://wiki.fysik.dtu.dk/ase/ase/md.html#module-md
.. _geometry optimization: https://wiki.fysik.dtu.dk/ase/ase/optimize.html#module-optimize
.. _constrained algorithms: https://wiki.fysik.dtu.dk/ase/ase/constraints.html#module-constraints
.. _nudged elastic band: https://wiki.fysik.dtu.dk/ase/ase/neb.html#module-neb
.. _thermostats: https://wiki.fysik.dtu.dk/ase/ase/md.html#constant-nvt-simulations-the-canonical-ensemble
.. _calculator: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html
.. _ASE tutorials: https://wiki.fysik.dtu.dk/ase/tutorials/tutorials.html


.. file:object_oriented

.. _objects:

Thinking in objects
-------------------

If you are already familiar with Python or languages such as Java or C++, you probably know what object oriented programming is. If you are more of a Fortran77 type, maybe not. Since Pysic relies heavily on the object paradigm, a few words should be said about it.

Let's say we want to simulate a bunch of atoms. To do this, we need to know where they are (coordinates), what element they are, what are their charges, momenta, etc. One way to store the information would be to assign indices to the atoms and create arrays containing the data. One for coordinates, another for momenta, third for charges etc. However, this type of bookkeeping approach gets more and more complicated the more strucuted data one needs to store. E.g., if we have a neighbor list as an array of atomic indices and wish to find the coordinates of a neighbor of a given atom, we first need to find the correct entry in the list of neighbors to find the index of the neighbor, and then find the entry corresponding to this index in the list of coordinates. In code, it could look like this::

  neighbor_atom_index = neighbor_lists[atom_index, neighbor_index]
  neighbor_atom_coordinates = coordinates[neighbor_atom_index, 1:3]

For complicated data hierarchies one may need to do this kind of index juggling for several rounds, which leads to code that is difficult to read and very susceptible to bugs. Furthermore, if one would want to edit the lists of atoms by, say, removing an atom from the system, all the arrays containing data related to atoms would have to be checked in case they contain the to be deleted particle and updated accordingly.

In the object oriented approach, one defines a data structure (called a *class* in Python) capable of storing various types of information in a single instance. So one can define an 'atom' datatype which contains the coordinates (as real numbers), momenta, etc. in one neat package. One can also define a 'neighbor list' datatype which contains a list of 'atom' datatypes. And the 'atom' datatype can contain a 'neighbor list' [#]_. Now, the problem of finding the coordinates of a neighbor is solved in a more intuitive way by asking the atom who the neigbor is and the neighbor its coordinates. This might look something like this::

  neighbor_atom_coordinates = atom.get_neighbor(neighbor_index).get_coordinates()

The above example also demonstrates *methods* - object specific functions allowing one to essentially give orders to objects. Objects and methods make it easy to write code that is simple to read and understand, since we humans intuitively see the world as objects, not as arrays of data. Another great benefit of the object based model is that when an object is modified, the changes automatically propagate everywhere where that object is referred.

The classes and their methods defined in Pysic are documented in detail in :ref:`syntax`, and their basic use is shown in the collection of provided :ref:`examples`. The central class in Pysic is :class:`~pysic.Pysic`, which is an energy and force calculator for `ASE`_. The interactions according to which the energies are calculated are constructed through the class :class:`~pysic.Potential`. Utilizing these classes is necessary to run meaningful calculations, though also other classes are defined for special purposes.

.. [#] Although, this is not exactly how ASE handles neighbor lists...

