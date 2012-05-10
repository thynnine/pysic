
.. _geometry:
        
===============================================
geometry (Geometry.f90)
===============================================



A module for handling the geometric structure of the system.


.. only:: html


    Modules used by geometry
    ------------------------
    - :ref:`quaternions`
    - :ref:`utility`

    List of global variables in geometry
    ------------------------------------
    - :data:`label_length`

    List of custom types in geometry
    --------------------------------
    - :data:`atom`
    - :data:`neighbor_list`
    - :data:`subcell`
    - :data:`supercell`

    List of subroutines in geometry
    -------------------------------
        
    - :func:`absolute_coordinates`
    - :func:`assign_bond_order_factor_indices`
    - :func:`assign_neighbor_list`
    - :func:`assign_potential_indices`
    - :func:`divide_cell`
    - :func:`expand_subcell_atom_capacity`
    - :func:`find_subcell_for_atom`
    - :func:`generate_atoms`
    - :func:`generate_supercell`
    - :func:`get_optimal_splitting`
    - :func:`relative_coordinates`
    - :func:`separation_vector`
    - :func:`update_atomic_charges`
    - :func:`update_atomic_positions`
    - :func:`wrapped_coordinates`

    List of functions in geometry
    -----------------------------
        
    - :func:`pick`


Full documentation of global variables in geometry
--------------------------------------------------
        
        
  .. data:: label_length

    integer    *scalar*  *parameter*  

    *initial value* = 2
    
    the number of characters available for denoting chemical symbols
    

Full documentation of custom types in geometry
----------------------------------------------
        
        
  .. data:: atom

    Defines an atomic particle.
    

    Contained data:

    neighbor_list: type(neighbor_list)    *scalar*
        the list of neighbors for the atom
    index: integer    *scalar*
        index of the atom
    n_pots: integer    *scalar*
        number of potentials that may affect the atom
    tags: integer    *scalar*
        integer tag
    potential_indices: integer  *pointer*  *size(:)*
        the indices of the potentials for which this atom is a valid target at first position (see :func:`potential_affects_atom`)
    potentials_listed: logical    *scalar*
        logical tag for checking if the potentials affecting the atom have been listed in potential_indices
    bond_indices: integer  *pointer*  *size(:)*
        the indices of the bond order factors for which this atom is a valid target at first position (see :func:`bond_order_factor_affects_atom`)
    element: character(len=label_length)    *scalar*
        the chemical symbol of the atom
    charge: double precision    *scalar*
        charge of the atom
    subcell_indices: integer    *size(3)*
        indices of the subcell containing the atom, used for fast neighbor searching (see :data:`subcell`)
    mass: double precision    *scalar*
        mass of th atom
    n_bonds: integer    *scalar*
        number of bond order factors that may affect the atom
    bond_order_factors_listed: logical    *scalar*
        logical tag for checking if the bond order factors affecting the atom have been listed in bond_indices
    position: double precision    *size(3)*
        coordinates of the atom
    momentum: double precision    *size(3)*
        momentum of the atom
  .. data:: neighbor_list

    Defines a list of neighbors for a single atom.
    The list contains the indices of the neighboring atoms
    as well as the periodic boundary condition (PBC) offsets.
    
    The offsets are integer
    triplets showing how many times must the supercell vectors
    be added to the position of the neighbor to find the
    neighboring image in a periodic system.
    For example, let the supercell be::
    
     [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]],
    
    i.e., a unit cube, with periodic boundaries.
    Now, if we have particles with coordinates::
    
     a = [1.5, 0.5, 0.5]
     b = [0.4, 1.6, 3.3]
    
    the closest separation vector :math:`\mathbf{r}_b-\mathbf{r}_a` between the particles is::
    
      [-.1, .1, -.2]
    
    obtained if we add the vector of periodicity::
    
      [1.0, -1.0, -3.0]
    
    to the coordinates of particle b. The offset vector
    (for particle b, when listing neighbors of a) is then::
    
      [1, -1, -3]
    
    Note that if the system is small, one atom can in
    principle appear several times in the neighbor list with
    different offsets.

    Contained data:

    neighbors: integer  *pointer*  *size(:)*
        indices of the neighboring atoms
    max_length: integer    *scalar*
        The allocated length of the neighbor lists. To avoid deallocating and reallocating memory, extra space is reserved for the neighbors in case the number of neighbors increases during simulation (due to atoms moving).
    pbc_offsets: integer  *pointer*  *size(:, :)*
        offsets for periodic boundaries for each neighbor
    n_neighbors: integer    *scalar*
        the number of neighbors in the lists
  .. data:: subcell

    Subvolume, which is a part of the supercell containing the simulation.
    
    The subcells are used in partitioning of the simulation space in subvolumes.
    This divisioning of the simulation cell is needed for quickly finding the
    neighbors of atoms (see also :class:`pysic.FastNeighborList`).
    The fast neighbor search is based on dividing the system, locating the subcell
    in which each atom is located, and then searching for neighbors for each atom
    by only checking the adjacent subcells. For small subvolumes (short cutoffs)
    this method is much faster than a brute force algorithm that checks all atom
    pairs. It also scales :math:`\mathcal{O}(n)`.
    

    Contained data:

    neighbors: integer    *size(3, -1:1, -1:1, -1:1)*
        indices of the 3 x 3 x 3 neighboring subcells (note that the neighboring subcell 0,0,0 is the cell itself)
    vector_lengths: double precision    *size(3)*
        lengths of the vectors spanning the subcell
    offsets: integer    *size(3, -1:1, -1:1, -1:1)*
        integer offsets of the neighboring subcells - if a neighboring subcell is beyond a periodic border, the offset records the fact
    max_atoms: integer    *scalar*
        the maximum number of atoms the cell can contain in the currently allocated memory space
    vectors: double precision    *size(3, 3)*
        the vectors spanning the subcell
    atoms: integer  *pointer*  *size(:)*
        indices of the atoms in this subcell
    n_atoms: integer    *scalar*
        the number of atoms contained by the subcell
    indices: integer    *size(3)*
        integer coordinates of the subcell in the subcell divisioning of the supercell
    include: logical    *size(-1:1, -1:1, -1:1)*
        A logical array noting if the neighboring subcells should be included in the neighbor search. Usually all neighbors are included, but in a non-periodic system, there is only a limited number of cells and once the system border is reached, this tag will be set to ``.false.`` to notify that there is no neighbor to be found.
  .. data:: supercell

    Supercell containing the simulation.
    
    The supercell is spanned by three vectors :math:`\mathbf{v}_1,\mathbf{v}_2,\mathbf{v}_3` stored as a
    :math:`3 \times 3` matrix in format
    
    .. math::
    
      \mathbf{M} = \left[
      \begin{array}{ccc}
      v_{1,x} & v_{1,y} & v_{1,z} \\
      v_{2,x} & v_{2,y} & v_{2,z} \\
      v_{3,x} & v_{3,y} & v_{3,z}
      \end{array}
      \right].
    
    Also the inverse cell matrix is kept for transformations between the absolute and fractional coordinates.
    

    Contained data:

    vector_lengths: double precision    *size(3)*
        the lengths of the cell spanning vectors (stored to avoid calculating the vector norms over and over)
    max_subcell_atom_count: integer    *scalar*
        the maximum number of atoms any of the subcells has
    n_splits: integer    *size(3)*
        the number of subcells there are in the subdivisioning of the cell, in the directions of the spanning vectors
    inverse_cell: double precision    *size(3, 3)*
        the inverse of the cell matrix :math:`\mathbf{M}^{-1}`
    subcells: type(subcell)  *pointer*  *size(:, :, :)*
        an array of :data:`subcell` subvolumes which partition the supercell
    vectors: double precision    *size(3, 3)*
        vectors spanning the supercell containing the system as a matrix :math:`\mathbf{M}`
    volume: double precision    *scalar*
        volume of the cell
    periodic: logical    *size(3)*
        logical switch determining if periodic boundary conditions are applied in the directions of the three cell spanning vectors
    reciprocal_cell: double precision    *size(3, 3)*
        the reciprocal cell as a matrix, :math:`\mathbf{M}_R = 2 \pi( \mathbf{M}^{-1} )^T`. That is, if :math:`\mathbf{b}_i` are the reciprocal lattice vectors and :math:`\mathbf{a}_j` the real space lattice vectors, then :math:`\mathbf{b}_i \mathbf{a}_j = 2 \pi \delta_{ij}`.

Full documentation of subroutines in geometry
---------------------------------------------
        
        
            
  .. function:: absolute_coordinates(relative, cell, position)

    Transforms from fractional to absolute coordinates.
    
    Absolute coordinates are the coordinates in the normal
    :math:`xyz` base,
    
    .. math::
    
       \mathbf{r} = x\mathbf{i} + y\mathbf{j} + z\mathbf{k}.
    
    Fractional coordiantes are the coordiantes in the base
    spanned by the vectors defining the supercell,
    :math:`\mathbf{v}_1`, :math:`\mathbf{v}_2`, :math:`\mathbf{v}_3`,
    
    .. math::
    
       \mathbf{r} = \tilde{x}\mathbf{v}_1 + \tilde{y}\mathbf{v}_2 + \tilde{z}\mathbf{v}_3.
    
    Notably, for positions inside the supercell, the fractional
    coordinates fall between 0 and 1.
    
    Transformation between the two bases is given by the cell
    matrix
    
    .. math::
    
       \left[
       \begin{array}{c}
       x \\
       y \\
       z
       \end{array} \right] = \mathbf{M}
       \left[
       \begin{array}{c}
       \tilde{x} \\
       \tilde{y} \\
       \tilde{z}
       \end{array} \right]
    

    Parameters:

    relative: double precision  *intent(in)*    *size(3)*  
        the fractional coordinates
    cell: type(supercell)  *intent(in)*    *scalar*  
        the supercell
    **position**: double precision  **intent(out)**    *size(3)*  
        the absolute coordinates
            
  .. function:: assign_bond_order_factor_indices(n_bonds, atom_in, indices)

    Save the indices of bond order factors affecting an atom.
    
    In bond order factor evaluation, it is important to loop
    over bond parameters quickly. As the evaluation of factors
    goes over atoms, atom pairs etc., it is useful to first
    filter the parameters by the first atom participating
    in the factor. Therefore, the atoms can be given
    a list of bond order parameters for which they are a suitable target
    as a 'first participant' (in a triplet A-B-C, A is the
    first participant).
    

    Parameters:

    n_bonds: integer  *intent(in)*    *scalar*  
        number of bond order factors
    **atom_in**: type(atom)  **intent(inout)**    *scalar*  
        the atom for which the bond order factors are assigned
    indices: integer  *intent(in)*    *size(n_bonds)*  
        the indices of the bond order factors
            
  .. function:: assign_neighbor_list(n_nbs, nbor_list, neighbors, offsets)

    Creates a neighbor list for one atom.
    
    The neighbor list will contain an array of the indices
    of the neighboring atoms as well as periodicity offsets,
    as explained in :data:`neighbor_list`
    
    The routine takes the neighbor_list object to be created
    as an argument. If the list is empty, it is initialized.
    If the list already contains information, the list is emptied and
    refilled. If the previous list has room to contain the new list
    (as in, it has enough allocated memory), no memory reallocation
    is done (since it will be slow if done repeatedly). Only if the
    new list is too long to fit in the reserved memory, the pointers
    are deallocated and reallocated.
    

    Parameters:

    n_nbs: integer  *intent(in)*    *scalar*  
        number of neighbors
    **nbor_list**: type(neighbor_list)  **intent(inout)**    *scalar*  
        The list of neighbors to be created.
    neighbors: integer  *intent(in)*    *size(n_nbs)*  
        array containing the indices of the neighboring atoms
    offsets: integer  *intent(in)*    *size(3, n_nbs)*  
        periodicity offsets
            
  .. function:: assign_potential_indices(n_pots, atom_in, indices)

    Save the indices of potentials affecting an atom.
    
    In force and energy evaluation, it is important to loop
    over potentials quickly. As the evaluation of energies
    goes over atoms, atom pairs etc., it is useful to first
    filter the potentials by the first atom participating
    in the interaction. Therefore, the atoms can be given
    a list of potentials for which they are a suitable target
    as a 'first participant' (in a triplet A-B-C, A is the
    first participant).
    

    Parameters:

    n_pots: integer  *intent(in)*    *scalar*  
        number of potentials
    **atom_in**: type(atom)  **intent(inout)**    *scalar*  
        the atom for which the potentials are assigned
    indices: integer  *intent(in)*    *size(n_pots)*  
        the indices of the potentials
            
  .. function:: divide_cell(cell, splits)

    Split the cell in subcells according to the given number of divisions.
    
    The argument 'splits' should be a list of three integers determining how many
    times the cell is split. For instance, if splits = [3,3,5], the cell is divided in
    3*3*5 = 45 subcells: 3 cells along the first two cell vectors and 5 along the third.
    
    The Cell itself is not changed, but an array 'subcells' is created, containing
    the subcells which are Cell instances themselves. These cells will contain additional
    data arrays 'neighbors' and 'offsets'. These are 3-dimensional arrays with each dimension
    running from -1 to 1. The neighbors array contains references to the neighboring subcell
    Cell instances.
    The offsets contain coordinate offsets with respect to the periodic boundaries. In other words,
    if a subcell is at the border of the original Cell, it will have neighbors at the other side
    of the cell due to periodic boundary conditions. But from the point of view of the subcell,
    the neighboring cell is not on the other side of the master cell, but a periodic image of that
    cell. Therefore, any coordinates in the the subcell to which the neighbors array refers to must
    in fact be shifted by a vector of the master cell. The offsets list contains the multipliers
    for the cell vectors to make these shifts.
    
    Example in 2D for simplicity: ``split = [3,4]`` creates subcells::
    
     (0,3) (1,3) (2,3)
     (0,2) (1,2) (2,2)
     (0,1) (1,1) (2,1)
     (0,0) (1,0) (2,0)
    
    subcell (0,3) will have the neighbors::
     (2,0) (0,0) (1,0)
     (2,3) (0,3) (1,3)
     (2,2) (0,2) (1,2)
    
    and offsets::
     [-1,1] [0,1] [0,1]
     [-1,0] [0,0] [0,0]
     [-1,0] [0,0] [0,0]
    
    Note that the central 'neighbor' is the cell itself.
    
    If a boundary is not periodic, extra subcells with indices 0 and split+1
    are created to pad the simulation cell. These will contain the atoms that
    are outside the simulation cell.

    Parameters:

    **cell**: type(supercell)  **intent(inout)**    *scalar*  
        
    splits: integer  *intent(in)*    *size(3)*  
        
            
  .. function:: expand_subcell_atom_capacity(atoms_list, old_size, new_size)


    Parameters:

    atoms_list: integer  *intent()*  *pointer*  *size(:)*  
        
    old_size: integer  *intent(in)*    *scalar*  
        
    new_size: integer  *intent(in)*    *scalar*  
        
            
  .. function:: find_subcell_for_atom(cell, at)


    Parameters:

    **cell**: type(supercell)  **intent(inout)**    *scalar*  
        
    **at**: type(atom)  **intent(inout)**    *scalar*  
        
            
  .. function:: generate_atoms(n_atoms, masses, charges, positions, momenta, tags, elements, atoms)

    Creates atoms to construct the system to be simulated.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    masses: double precision  *intent(in)*    *size(n_atoms)*  
        array of masses for the atoms
    charges: double precision  *intent(in)*    *size(n_atoms)*  
        array of charges for the atoms
    positions: double precision  *intent(in)*    *size(3, n_atoms)*  
        array of coordinates for the atoms
    momenta: double precision  *intent(in)*    *size(3, n_atoms)*  
        array of momenta for the atoms
    tags: integer  *intent(in)*    *size(n_atoms)*  
        array of integer tags for the atoms
    elements: character(len=label_length)  *intent(in)*    *size(n_atoms)*  
        array of chemical symbols for the atoms
    atoms: type(atom)  *intent()*  *pointer*  *size(:)*  
        array of the atom objects created
            
  .. function:: generate_supercell(vectors, inverse, periodicity, cell)

    Creates the supercell containing the simulation geometry.
    
    The supercell is spanned by three vectors :math:`\mathbf{v}_1,\mathbf{v}_2,\mathbf{v}_3` stored as a
    :math:`3 \times 3` matrix in format
    
    .. math::
    
      \mathbf{M} = \left[
      \begin{array}{ccc}
      v_{1,x} & v_{1,y} & v_{1,z} \\
      v_{2,x} & v_{2,y} & v_{2,z} \\
      v_{3,x} & v_{3,y} & v_{3,z}
      \end{array}
      \right].
    
    Also the inverse cell matrix :math:`\mathbf{M}^{-1}` must be given
    for transformations between the absolute and fractional coordinates.
    However, it is not checked that the given matrix and inverse truly
    fulfill :math:`\mathbf{M}^{-1}\mathbf{M} = \mathbf{I}` - it is the
    responsibility of the caller to give the true inverse.
    
    Also the periodicity of the system in the directions of the
    cell vectors need to be given.
    

    Parameters:

    vectors: double precision  *intent(in)*    *size(3, 3)*  
        the cell spanning matrix :math:`\mathbf{M}`
    inverse: double precision  *intent(in)*    *size(3, 3)*  
        the inverse cell :math:`\mathbf{M}`
    periodicity: logical  *intent(in)*    *size(3)*  
        logical switch, true if the boundaries are periodic
    **cell**: type(supercell)  **intent(out)**    *scalar*  
        the created cell object
            
  .. function:: get_optimal_splitting(cell, max_cut, splits)


    Parameters:

    cell: type(supercell)  *intent(in)*    *scalar*  
        
    max_cut: double precision  *intent(in)*    *scalar*  
        
    **splits**: integer  **intent(out)**    *size(3)*  
        
            
  .. function:: relative_coordinates(position, cell, relative)

    Transforms from absolute to fractional coordinates.
    
    Absolute coordinates are the coordinates in the normal
    :math:`xyz` base,
    
    .. math::
    
       \mathbf{r} = x\mathbf{i} + y\mathbf{j} + z\mathbf{k}.
    
    Fractional coordiantes are the coordiantes in the base
    spanned by the vectors defining the supercell,
    :math:`\mathbf{v}_1`, :math:`\mathbf{v}_2`, :math:`\mathbf{v}_3`,
    
    .. math::
    
       \mathbf{r} = \tilde{x}\mathbf{v}_1 + \tilde{y}\mathbf{v}_2 + \tilde{z}\mathbf{v}_3.
    
    Notably, for positions inside the supercell, the fractional
    coordinates fall between 0 and 1.
    
    Transformation between the two bases is given by the inverse cell
    matrix
    
    .. math::
    
       \left[
       \begin{array}{c}
       \tilde{x} \\
       \tilde{y} \\
       \tilde{z}
       \end{array} \right] = \mathbf{M}^{-1}
       \left[
       \begin{array}{c}
       x \\
       y \\
       z
       \end{array} \right]
    

    Parameters:

    position: double precision  *intent(in)*    *size(3)*  
        the absolute coordinates
    cell: type(supercell)  *intent(in)*    *scalar*  
        the supercell
    **relative**: double precision  **intent(out)**    *size(3)*  
        the fractional coordinates
            
  .. function:: separation_vector(r1, r2, offset, cell, separation)

    Calculates the minimum separation vector between two atoms, :math:`\mathbf{r}_2-\mathbf{r}_1`, including possible periodicity.
    

    Parameters:

    r1: double precision  *intent(in)*    *size(3)*  
        coordiantes of atom 1, :math:`\mathbf{r}_1`
    r2: double precision  *intent(in)*    *size(3)*  
        coordinates of atom 1, :math:`\mathbf{r}_2`
    offset: integer  *intent(in)*    *size(3)*  
        periodicity offset (see :data:`neighbor_list`)
    cell: type(supercell)  *intent(in)*    *scalar*  
        supercell spanning the system
    **separation**: double precision  **intent(out)**    *size(3)*  
        the calculated separation vector, :math:`\mathbf{r}_2-\mathbf{r}_1`
            
  .. function:: update_atomic_charges(n_atoms, charges, atoms)

    Updates the charges of the given atoms.
    Other properties are not altered.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    charges: double precision  *intent(in)*    *size(n_atoms)*  
        new charges for the atoms
    atoms: type(atom)  *intent()*  *pointer*  *size(:)*  
        the atoms to be edited
            
  .. function:: update_atomic_positions(n_atoms, positions, momenta, atoms)

    Updates the positions and momenta of the given atoms.
    Other properties are not altered.
    
    This is meant to be used
    during dynamic simulations or geometry optimization
    where the atoms are only moved around, not changed in other ways.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    positions: double precision  *intent(in)*    *size(3, n_atoms)*  
        new coordinates for the atoms
    momenta: double precision  *intent(in)*    *size(3, n_atoms)*  
        new momenta for the atoms
    atoms: type(atom)  *intent()*  *pointer*  *size(:)*  
        the atoms to be edited
            
  .. function:: wrapped_coordinates(position, cell, wrapped)

    Wraps a general coordinate inside the supercell if the system is periodic.
    
    In a periodic system, every particle has periodic images at intervals
    defined by the cell vectors :math:`\mathbf{v}_1,\mathbf{v}_2,\mathbf{v}_3`.
    That is, for a particle at :math:`\mathbf{r}`, there are periodic
    images at
    
    .. math::
    
       \mathbf{R} = \mathbf{r} + a_1 \mathbf{v}_1 + a_2 \mathbf{v}_2 + a_3 \mathbf{v}_3
    
    for all :math:`a_1, a_2, a_3 \in \mathbf{Z}`.
    These are equivalent positions in the sense that if a particle is
    situated at any of one of them, the set of images is the same.
    Exactly one of the images is inside the cell - this routine gives
    the coordinates of that particular image.
    
    If the system is periodic in only some directions, the wrapping is
    done only along those directions.
    

    Parameters:

    position: double precision  *intent(in)*    *size(3)*  
        the absolute coordinates
    cell: type(supercell)  *intent(in)*    *scalar*  
        the supercell
    **wrapped**: double precision  **intent(out)**    *size(3)*  
        the wrapped absolute coordinates

Full documentation of functions in geometry
---------------------------------------------
        
        
            
  .. function:: pick(index1, index2, offset)

    A utility function for sorting the atoms.
    
    The function return ``true`` if ``index1 < index2`` and ``false`` otherwise.
    If ``index1 == index2``, the comparison is made through the separation vector.
    The vector is examined element at a time, and if a positive number is found,
    ``true`` is returned, if a negative one, ``false``. For values of zero, the next
    element is examined.
    
    The purpose for this function is to sort the atoms to prevent double counting when summing
    over pairs. In principle, a sum over pairs :math:`(i,j)` can be done with
    :math:`\frac{1}{2} \sum_{i \ne j}`, but this leads to evaluation of all elements twice
    (both :math:`(i,j)` and :math:`(j,i)` are considered separately).
    It is more efficient to evaluate :math:`\sum_{i < j}`, where only one of :math:`(i,j)` and :math:`(j,i)`
    fullfill the condition.
    
    A special case arises if interactions are so long ranged that an atom can see its own periodic
    images. Then, one will need to sum terms for atom pairs where both atoms have the same index
    :math:`\sum_\mathrm{images} \sum_{i,j}` if they are in different periodic copies of the actual
    simulation cell. In order to still pick only one of the pairs :math:`(i,i')` and :math:`(i',i)`,
    we compare the offset vectors. If atom :math:`i'` is in the neighboring cell of :math:`i` in the
    first cell vector direction, it has an offset of :math:`[1,0,0]` and vice versa :math:`i` has
    an offset of :math:`[-1,0,0]` from :math:`i'`. Instead of the index, the sorting :math:`i' < i`
    is then done by comparing these offset vectors, element by element.
    

    Parameters:

    index1: integer  *intent(in)*    *scalar*  
        index of first atom
    index2: integer  *intent(in)*    *scalar*  
        index of second atom
    offset: integer  *intent(in)*    *size(3)*  
        pbc offset vector from atom1 to atom2