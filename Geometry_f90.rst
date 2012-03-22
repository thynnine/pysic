
.. _geometry:
        
===============================================
geometry (Geometry.f90)
===============================================



A module for handling the geometric structure of the system.


.. only:: html


    Modules used by geometry
    ------------------------
    - :ref:`quaternions`

    List of global variables in geometry
    ------------------------------------
    - :data:`label_length`

    List of custom types in geometry
    --------------------------------
    - :data:`atom`
    - :data:`neighbor_list`
    - :data:`supercell`

    List of subroutines in geometry
    -------------------------------
        
    - :func:`absolute_coordinates`
    - :func:`assign_bond_order_factor_indices`
    - :func:`assign_neighbor_list`
    - :func:`assign_potential_indices`
    - :func:`generate_atoms`
    - :func:`generate_supercell`
    - :func:`relative_coordinates`
    - :func:`separation_vector`
    - :func:`update_atomic_positions`
    - :func:`wrapped_coordinates`


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
    element: character(len=label_length)    *scalar*
        the chemical symbol of the atom
    charge: double precision    *scalar*
        charge of the atom
    bond_indices: integer  *pointer*  *size(:)*
        the indices of the bond order factors for which this atom is a valid target at first position (see :func:`bond_order_factor_affects_atom`)
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

    periodic: logical    *size(3)*
        
    vector_lengths: double precision    *size(3)*
        
    vectors: double precision    *size(3, 3)*
        vectors spanning the supercell containing the system as a matrix :math:`\mathbf{M}`
    inverse_cell: double precision    *size(3, 3)*
        the inverse of the cell matrix :math:`\mathbf{M}^{-1}`

Full documentation of subroutines in geometry
---------------------------------------------
        
        
            
  .. function:: absolute_coordinates(relative, cell, position)


    Parameters:

    relative: double precision  *intent(in)*    *size(3)*  
        
    cell: type(supercell)  *intent(in)*    *scalar*  
        
    **position**: double precision  **intent(out)**    *size(3)*  
        
            
  .. function:: assign_bond_order_factor_indices(n_bonds, atom_in, indices)


    Parameters:

    n_bonds: integer  *intent(in)*    *scalar*  
        
    **atom_in**: type(atom)  **intent(inout)**    *scalar*  
        
    indices: integer  *intent(in)*    *size(n_bonds)*  
        
            
  .. function:: assign_neighbor_list(n_nbs, nbor_list, neighbors, offsets)


    Parameters:

    n_nbs: integer  *intent(in)*    *scalar*  
        
    **nbor_list**: type(neighbor_list)  **intent(inout)**    *scalar*  
        
    neighbors: integer  *intent(in)*    *size(n_nbs)*  
        
    offsets: integer  *intent(in)*    *size(3, n_nbs)*  
        
            
  .. function:: assign_potential_indices(n_pots, atom_in, indices)


    Parameters:

    n_pots: integer  *intent(in)*    *scalar*  
        
    **atom_in**: type(atom)  **intent(inout)**    *scalar*  
        
    indices: integer  *intent(in)*    *size(n_pots)*  
        
            
  .. function:: generate_atoms(n_atoms, masses, charges, positions, momenta, tags, elements, atoms)

    Creates atoms to construct the system to be simulated

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        
    masses: double precision  *intent(in)*    *size(n_atoms)*  
        
    charges: double precision  *intent(in)*    *size(n_atoms)*  
        
    positions: double precision  *intent(in)*    *size(3, n_atoms)*  
        
    momenta: double precision  *intent(in)*    *size(3, n_atoms)*  
        
    tags: integer  *intent(in)*    *size(n_atoms)*  
        
    elements: character(len=label_length)  *intent(in)*    *size(n_atoms)*  
        
    atoms: type(atom)  *intent()*  *pointer*  *size(:)*  
        
            
  .. function:: generate_supercell(vectors, inverse, periodicity, cell)

    Creates the supercell

    Parameters:

    vectors: double precision  *intent(in)*    *size(3, 3)*  
        
    inverse: double precision  *intent(in)*    *size(3, 3)*  
        
    periodicity: logical  *intent(in)*    *size(3)*  
        
    **cell**: type(supercell)  **intent(out)**    *scalar*  
        
            
  .. function:: relative_coordinates(position, cell, relative)


    Parameters:

    position: double precision  *intent(in)*    *size(3)*  
        
    cell: type(supercell)  *intent(in)*    *scalar*  
        
    **relative**: double precision  **intent(out)**    *size(3)*  
        
            
  .. function:: separation_vector(r1, r2, offset, cell, separation)

    Calculates the minimum separation vector between two atoms, r2-r1, including possible periodicity

    Parameters:

    r1: double precision  *intent(in)*    *size(3)*  
        
    r2: double precision  *intent(in)*    *size(3)*  
        
    offset: integer  *intent(in)*    *size(3)*  
        
    cell: type(supercell)  *intent(in)*    *scalar*  
        
    **separation**: double precision  **intent(out)**    *size(3)*  
        
            
  .. function:: update_atomic_positions(n_atoms, positions, momenta, atoms)

    Updates the positions and momenta of atoms

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        
    positions: double precision  *intent(in)*    *size(3, n_atoms)*  
        
    momenta: double precision  *intent(in)*    *size(3, n_atoms)*  
        
    atoms: type(atom)  *intent()*  *pointer*  *size(:)*  
        
            
  .. function:: wrapped_coordinates(position, cell, wrapped)

    Wraps a general coordinate inside the supercell if the system is periodic

    Parameters:

    position: double precision  *intent(in)*    *size(3)*  
        
    cell: type(supercell)  *intent(in)*    *scalar*  
        
    **wrapped**: double precision  **intent(out)**    *size(3)*  
        