.. file:fastneighborlist class

.. _fastneighborlist class:



.. file:fastneighborlist class - description

.. _fastneighborlist class - description:



=========================
FastNeighborList class
=========================

This class extends the `ASE NeighborList`_ class to provide a more efficient neighbor finding tool.
The neighbor finding routine searches the neighborhoods of all atoms and for each atom records which other atoms are closer than a given cutoff distance.



.. file:neighbor lists

.. _neighbor lists:




The benefit of neighbor lists
------------------------------

Atomistic pair and many-body potentials typically depend on the local atomic structure and especially the relative coordinates of the atoms. However, finding the separation vector and distance between coordinates in periodic 3D space is computationally fairly costly operation and the number of atom-atom pairs in the system grows as :math:`\mathcal{O}(n^2)`. Therefore the evaluation of local potentials can be made efficient by storing lists of nearby atoms for all particles to narrow down the scope of search for interacting neighbors.

Typically one chooses a cutoff distance :math:`r_\mathrm{cut}` beyond which the atoms do not see each other. Then, the neighbor lists should always contain all the atoms within this cutoff radius :math:`r_{ij} \le r_\mathrm{cut}`. In dynamic simulations where the atoms move, the typical scheme is to list atoms within a slightly longer radius, :math:`r_\mathrm{cut} + r_\mathrm{skin}` because then the lists need not be updated until an atom has moved by more than :math:`r_\mathrm{skin}`.

.. file:fast search

.. _fast search:





Faster neighbor search
---------------------------

There is a built in neighbor searching tool in ASE, `ASE NeighborList`_. It is, however, a pure Python implementation using a brute-force :math:`\mathcal{O}(n^2)` algorithm making it slow - even prohibitively slow - for large systems especially when periodic boundary conditions are used.

To overcome this performance bottleneck, Pysic implements the :class:`~pysic.calculator.FastNeighborList` class. This class inherits other properties from the built-in ASE class except for the :meth:`~pysic.calculator.FastNeighborList.build` method, which is replaced by a faster algorithm. The fast neighbor search is implemented in Fortran and parallelized with MPI. The algorithm is based on a spatial divisioning, i.e.

- the simulation volume is divided in subvolumes
- for each atom the subvolume where it is contained is found
- for each atom, the neighbors are searched for only in the adjacent subvolumes

For a fixed cutoff, the neighborhood searched for each atom is constant and thus this is an :math:`\mathcal{O}(n)` algorithm. [#]_ The method is also faster the shorter the cutoffs are. For short cutoffs (~ 5 Å), a 10000 atom periodic system is expected to be handled 100 or even 1000 fold faster with :class:`~pysic.calculator.FastNeighborList` than with the ASE method.

.. [#] For very large systems the number of subdivisions is limited to conserve memory so the :math:`\mathcal{O}(n)` scaling is eventually lost. Say we divide the volume in a hundred subvolumes along each axis; we end up with a million subvolumes which is a lot!

.. file:limitations

.. _limitations:




Limitations in the implementation
---------------------------------

Since the fast algorithm is implemented in Fortran, it operates on the structure allocated in the Fortran core. Therefore, even though the :meth:`~pysic.calculator.FastNeighborList.build` method takes an `ASE Atoms`_ object as an argument, it does not analyze the given structure. It does check against :class:`~pysic.CoreMirror` to see if the given structure matches the one in the core and raises an error if not, but accessing the core has to still be done through :class:`~pysic.calculator.Pysic`. When :class:`~pysic.calculator.Pysic` is run normally, this is automatically taken care of. As the implementation is MPI parallelized, it is also necessary that the MPI environment has been set up - especially the distribution of load (i.e. atoms) between processors must be done before the lists can be built.

Another more profound limitation in the current implementation of the algorithm is the fact that it limits the neighbor finding to neighboring subvolumes. Since the subvolumes are not allowed to be larger than the actual simulation volume, the cutoffs cannot be longer than the shortest perpendicular separation between facets of the subvolume. For rectangular cells, this is just the minimum of the lengths of the vectors spanning the cell, :math:`\mathbf{v}_{i,j,k}`. For inclined cell shapes, the perpendicular distance between cell facets, :math:`d`, is

.. math::
  :nowrap:

  \begin{eqnarray}
  d_i & = & \frac{|\mathbf{v}_i \cdot \mathbf{n_i}|}{|\mathbf{n}_{i}|}\\
  \mathbf{n}_i & = & \mathbf{v}_j \times \mathbf{v}_k 
  \end{eqnarray}

where :math:`\mathbf{n}_i` are the normal vectors of the plane spanned by the vectors :math:`\mathbf{v}_{j,k}`. If one wishes to find neighbors in a radius containing the simulation volume several times, the original `ASE NeighborList`_ should be used instead. :class:`~pysic.calculator.Pysic` does this choice automatically when building the neighbor lists. One should usually avoid such long cutoffs in the first place, but if your system is very small that may not be possible.


.. _ASE NeighborList: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#building-neighbor-lists
.. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html



.. file:listing neighbors

.. _listing neighbors:



Accessing neighbor lists
-----------------------------

Although the neighbor lists are automatically used by the calculator, it is also possible to use the list for manually analysing local atomic neighborhoods. This can be done with the methods

- :meth:`~pysic.calculator.FastNeighborList.get_neighbors`
- :meth:`~pysic.calculator.FastNeighborList.get_neighbor_separations`
- :meth:`~pysic.calculator.FastNeighborList.get_neighbor_distances`

which list, for a specified atom, the indices and periodic boundary offsets of neighbors (cf. `get_neighbors <https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#ase.calculators.neighborlist.NeighborList.get_neighbors>`_ of the ASE neighbor list), separation vectors, and separation distances, respectively. The arrays can also be requested pre-sorted according to the atom-atom distance by providing the keyword ``sort=True``.

Since the lists are created in the Fortran core according to specified interaction cutoffs, it is unfortunately not possible to inquire the neighbors within arbitrary radii without touching the potentials. Normally the lists contain the atoms which interact. In order to just analyse structures, a dummy calculator needs to be created::

 system = ase.Atoms(...)
 system.set_cell(5*np.identity(3))
 system.set_pbc([True,True,True])

 # create a dummy calculator
 dummy = pysic.Pysic()
 # Find neighbors at distance 3.0 + 0.5 (0.5 is the default marginal).
 # Note that the list finds all neighbors within the maximum interaction
 # radius of the particular atom.
 pot = pysic.Potential('LJ', cutoff=3.0, symbols=...)
 dummy.add_potential(pot)
 system.set_calculator(dummy)
 # It is important to manually initialize the core since no actual calculations are not carried out.
 dummy.set_core()

 # get the list and access its contents
 nbl = dummy.get_neighbor_list()
 neighbors, offsets = nbl.get_neighbors(0, system, True)
 separations = nbl.get_neighbor_separations(0, system, True)
 distances = nbl.get_neighbor_distances(0, system, True)
 




.. file:fastneighborlist class - autogenerated

.. _fastneighborlist class - autogenerated:




Methods inherited from ASE NeighborList
---------------------------------------

- `update <https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#ase.calculators.neighborlist.NeighborList.update>`_

List of methods
---------------

- :meth:`~pysic.calculator.FastNeighborList.build`
- :meth:`~pysic.calculator.FastNeighborList.get_neighbors`
- :meth:`~pysic.calculator.FastNeighborList.get_neighbor_separations`
- :meth:`~pysic.calculator.FastNeighborList.get_neighbor_distances`


Full documentation of the FastNeighborList class
------------------------------------------------

.. currentmodule:: pysic.calculator
.. autoclass:: FastNeighborList
   :members:
   :undoc-members:

