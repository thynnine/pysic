.. file:syntax

.. _syntax:



.. file:syntax forewords

.. _syntax forewords:




Structure and syntax of Pysic
=============================

Pysic can be used with basic functionality with just a few commands. To give access to all the functionality in Pysic, the full API of Pysic is documented in the following section. All the classes and methods including their arguments are explained in detail. Also the types of potentials and bond order factors available are documented, including their mathematical descriptions and the keywords needed for access.

In addition to the Python documentation, also the variables, types, and routines in the Fortran core are listed with comments and explanations. Besides the interface module :ref:`pysic_interface`, this part of the code is not accessible through Python. Thus, you need not know what the core contains. However, if you plan to modify the core, studying also the Fortran documentation is useful.

The graph below show the main dependancies between the Python and Fortran classes and modules. (This is not a full UML diagram, just a schematic presentation.)

|diagram|

.. |diagram| image:: ../../../Graphics/mordred/pysic/pysic_diagram.pdf
             :alt: class and module graph
             :height: 400

.. comment out graphviz

   .. graphviz::

   digraph pysic {
     fontname=Helvetica;
     size=8;
	node [fontname=Helvetica fontsize=10];
	subgraph cluster_5{
	subgraph cluster_0 {
	        style=filled;
		color=lightgray;
		node [shape=rectangle, fillcolor=white, style=filled];
		CoreMirror -> Pysic;
		ChargeRelaxation -> Pysic;
		CoulombSummation -> Pysic;
		BondOrderParameters -> Coordinator -> Potential -> Pysic;
		FastNeighborList -> Pysic;
		label = "pysic";
	}
	subgraph cluster_1 {
	        style=filled;
		color=lightgray;
		node [shape=rectangle, fillcolor=white, style=filled];
	        Atom -> Atoms -> Pysic;
		Atoms -> NeighborList;
		label = "ase";
		FastNeighborList -> NeighborList [arrowhead=empty];
	}
	subgraph cluster_2 {
	        color=white;
		node [shape=rectangle, fillcolor=white, style=filled];
		pysic_utility -> Pysic;
		pysic_interface -> Pysic;
	}
	label = "Python";
	}
	subgraph cluster_3 {
	        style=filled;
		color=lightgray;
		node [shape=rectangle, fillcolor=white, style=filled];
		PyInterface -> pysic_interface;
		label = "Fortran 90";
	}

   }


   .. graphviz::

   digraph pysic {
     fontname=Helvetica;
     size=8;
	node [fontname=Helvetica fontsize=10];
	subgraph cluster_5{
	subgraph cluster_2 {
	        color=white;
		node [shape=rectangle, fillcolor=white, style=filled];
		pysic_interface;
	}
	label = "Python";
	}
	subgraph cluster_3 {
	        style=filled;
		color=lightgray;
		node [shape=rectangle, fillcolor=white, style=filled];
		PyInterface -> pysic_interface;
		Mersenne -> PyInterface;
		MPI -> PyInterface;
		Utility -> PyInterface;
		Core -> PyInterface;
		Mersenne -> MPI -> Core;
		MPI -> Potentials -> Core;
		Utility -> Geometry -> Core;
		Utility -> Potentials;
		Quaternions -> Potentials;
		Quaternions -> Geometry -> Potentials;
		label = "Fortran 90";
	}

   }

.. only:: html

 Modules
 -------

 - :mod:`pysic`
 - :mod:`pysic_utility`
 - :mod:`pysic_fortran`

 Classes
 -------

 - :class:`pysic.calculator.Pysic`
 - :class:`pysic.Potential`
 - :class:`pysic.CoulombSummation`
 - :class:`pysic.Coordinator`
 - :class:`pysic.BondOrderParameters`
 - :class:`pysic.ChargeRelaxation`
 - :class:`pysic.FastNeighborList`
 - :class:`pysic.CoreMirror`

