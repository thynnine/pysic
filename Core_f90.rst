
.. _pysic_core:
        
=============================================
pysic_core (Core.f90)
=============================================



Core, true to its name, is the heart of the Fortran core
of Pysic. It contains the data structures defining the simulation
geometry and interactions and defines the central routines for
calculating the total energy of and the forces acting on the system.

Many of the routines in :ref:`pysic_interface` which `f2py`_ interfaces
to Python are simply calling routines here.


.. _f2py: http://www.scipy.org/F2py


.. only:: html


    Modules used by pysic_core
    --------------------------
    - :ref:`geometry`
    - :ref:`mpi`
    - :ref:`potentials`

    List of global variables in pysic_core
    --------------------------------------
    - :data:`atoms`
    - :data:`atoms_created`
    - :data:`bond_factors`
    - :data:`bond_factors_allocated`
    - :data:`bond_storage_allocated`
    - :data:`cell`
    - :data:`evaluate_ewald`
    - :data:`ewald_allocated`
    - :data:`ewald_cutoff`
    - :data:`ewald_epsilon`
    - :data:`ewald_k_cutoffs`
    - :data:`ewald_scaler`
    - :data:`ewald_sigma`
    - :data:`group_index_save_slot`
    - :data:`interactions`
    - :data:`n_bond_factors`
    - :data:`n_interactions`
    - :data:`n_saved_bond_order_factors`
    - :data:`potentials_allocated`
    - :data:`saved_bond_order_factors`
    - :data:`saved_bond_order_gradients`
    - :data:`saved_bond_order_sums`
    - :data:`use_saved_bond_order_factors`
    - :data:`use_saved_bond_order_gradients`

    List of subroutines in pysic_core
    ---------------------------------
        
    - :func:`core_add_bond_order_factor`
    - :func:`core_add_potential`
    - :func:`core_allocate_bond_order_factors`
    - :func:`core_allocate_bond_order_storage`
    - :func:`core_allocate_potentials`
    - :func:`core_assign_bond_order_factor_indices`
    - :func:`core_assign_potential_indices`
    - :func:`core_calculate_bond_order_factors`
    - :func:`core_calculate_bond_order_gradients`
    - :func:`core_calculate_bond_order_gradients_of_factor`
    - :func:`core_calculate_electronegativities`
    - :func:`core_calculate_energy`
    - :func:`core_calculate_forces`
    - :func:`core_clear_atoms`
    - :func:`core_clear_bond_order_factors`
    - :func:`core_clear_bond_order_storage`
    - :func:`core_clear_potentials`
    - :func:`core_create_cell`
    - :func:`core_create_neighbor_list`
    - :func:`core_empty_bond_order_gradient_storage`
    - :func:`core_empty_bond_order_storage`
    - :func:`core_fill_bond_order_storage`
    - :func:`core_generate_atoms`
    - :func:`core_get_bond_order_factor_of_atom`
    - :func:`core_get_bond_order_factors`
    - :func:`core_get_bond_order_gradients`
    - :func:`core_get_bond_order_sums`
    - :func:`core_get_cell_vectors`
    - :func:`core_get_ewald_energy`
    - :func:`core_get_number_of_atoms`
    - :func:`core_post_process_bond_order_factors`
    - :func:`core_post_process_bond_order_gradients`
    - :func:`core_post_process_bond_order_gradients_of_factor`
    - :func:`core_release_all_memory`
    - :func:`core_set_ewald_parameters`
    - :func:`core_update_atom_charges`
    - :func:`core_update_atom_coordinates`
    - :func:`list_atoms`
    - :func:`list_bonds`
    - :func:`list_cell`
    - :func:`list_interactions`


Full documentation of global variables in pysic_core
----------------------------------------------------
        
        
  .. data:: atoms

    type(atom)  *pointer*  *size(:)*    
    
    an array of :data:`atom` objects representing the system
    
  .. data:: atoms_created

    logical    *scalar*    

    *initial value* = .false.
    
    logical tag indicating if atom storing arrays have been created
    
  .. data:: bond_factors

    type(bond_order_parameters)  *pointer*  *size(:)*    
    
    an array of :data:`bond_order_parameters` objects representing bond order factors modifying the potentials
    
  .. data:: bond_factors_allocated

    logical    *scalar*    

    *initial value* = .false.
    
    logical tag indicating if bond order parameter storing arrays have been allocated
    
  .. data:: bond_storage_allocated

    logical    *scalar*    

    *initial value* = .false.
    
    logical tag indicating if bond order factor storing arrays have been allocated
    
  .. data:: cell

    type(supercell)    *scalar*    
    
    a :data:`supercell` object representing the simulation cell
    
  .. data:: evaluate_ewald

    logical    *scalar*    

    *initial value* = .false.
    
    switch for enabling Ewald summation of coulomb interactions
    
  .. data:: ewald_allocated

    logical    *scalar*    

    *initial value* = .false.
    
    
    
  .. data:: ewald_cutoff

    double precision    *scalar*    
    
    
    
  .. data:: ewald_epsilon

    double precision    *scalar*    
    
    
    
  .. data:: ewald_k_cutoffs

    integer    *size(3)*    
    
    
    
  .. data:: ewald_scaler

    double precision  *pointer*  *size(:)*    
    
    
    
  .. data:: ewald_sigma

    double precision    *scalar*    
    
    
    
  .. data:: group_index_save_slot

    integer  *pointer*  *size(:)*    
    
    
    
  .. data:: interactions

    type(potential)  *pointer*  *size(:)*    
    
    an array of :data:`potential` objects representing the interactions
    
  .. data:: n_bond_factors

    integer    *scalar*    

    *initial value* = 0
    
    
    
  .. data:: n_interactions

    integer    *scalar*    

    *initial value* = 0
    
    number of potentials
    
  .. data:: n_saved_bond_order_factors

    integer    *scalar*    

    *initial value* = 0
    
    number of saved bond order factors
    
  .. data:: potentials_allocated

    logical    *scalar*    

    *initial value* = .false.
    
    logical tag indicating if potential storing arrays have been allocated
    
  .. data:: saved_bond_order_factors

    double precision  *pointer*  *size(:, :)*    
    
    Array for storing calculated bond order factors. Indexing: (atom index, group_index_save_slot(group index))
    
  .. data:: saved_bond_order_gradients

    double precision  *pointer*  *size(:, :, :, :)*    
    
    Array for storing calculated bond order gradients. Indexing: (xyz, atom index, group_index_save_slot(group index), target index)
    
  .. data:: saved_bond_order_sums

    double precision  *pointer*  *size(:, :)*    
    
    Array for storing calculated bond order sums. Indexing: (atom index, group_index_save_slot(group index))
    
  .. data:: use_saved_bond_order_factors

    logical    *scalar*    

    *initial value* = .false.
    
    Logical tag which enables / disables bond order saving. If true, bond order calculation routines try to find the precalculated factors in the saved bond order arrays instead of calculating.
    
  .. data:: use_saved_bond_order_gradients

    integer  *pointer*  *size(:, :)*    
    
    Array storing the atom index of the bond gradient stored for indices (group index, target index). Since gradients are needed for all factors (N) with respect to moving all atoms (N), storing them all would require an N x N matrix. Therefore only some are stored. This array is used for searching the stroage to see if the needed gradient is there or needs to be calculated.
    

Full documentation of subroutines in pysic_core
-----------------------------------------------
        
        
            
  .. function:: core_add_bond_order_factor(n_targets, n_params, n_split, bond_name, parameters, param_split, cutoff, smooth_cut, elements, orig_elements, group_index)

    Creates one additional bond_order_factor in the core.
    The routine assumes that adequate memory has been
    allocated already using core_allocate_bond_order_factors.
    
    When the bond order parameters in the Python interface are imported
    to the Fortran core, the target specifiers (elements)
    are permutated to create all equivalent bond order parameters.
    That is, if we have parameters for Si-O, both Si-O and O-Si
    parameters are created. This is because the energy and
    force calculation loops only deal with atom pairs A-B once
    (so only A-B or B-A is considered, not both) and if, say,
    the loop only finds an O-Si pair, it is important to apply
    the Si-O parameters also on that pair.
    In some cases, such as with the tersoff factor affecting
    triplets (A-B-C), the contribution is not symmetric for all the atoms.
    Therefore it is necessary to also store the original targets of
    the potential as specified in the Python interface. These are
    to be given in the 'orig_elements' lists.
    
    called from PyInterface: :func:`add_bond_order_factor`
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets (interacting bodies)
    n_params: integer  *intent(in)*    *scalar*  
        number of parameters
    n_split: integer  *intent(in)*    *scalar*  
        number of subsets in the list of parameters, should equal n_targets
    bond_name: character(len=*)  *intent(in)*    *scalar*  
        bond order factor names
    parameters: double precision  *intent(in)*    *size(n_params)*  
        numeric parameters
    param_split: integer  *intent(in)*    *size(n_split)*  
        the numbers of parameters for 1-body, 2-body etc.
    cutoff: double precision  *intent(in)*    *scalar*  
        interaction hard cutoff
    smooth_cut: double precision  *intent(in)*    *scalar*  
        interaction soft cutoff
    elements: character(len=label_length)  *intent(in)*    *size(n_targets)*  
        atomic symbols specifying the elements the interaction acts on
    orig_elements: character(len=label_length)  *intent(in)*    *size(n_targets)*  
        original atomic symbols specifying the elements the interaction acts on
    group_index: integer  *intent(in)*    *scalar*  
        index denoting the potential to which the factor is connected
            
  .. function:: core_add_potential(n_targets, n_params, pot_name, parameters, cutoff, smooth_cut, elements, tags, indices, orig_elements, orig_tags, orig_indices, pot_index)

    Creates one additional potential in the core.
    The routine assumes that adequate memory has been
    allocated already using core_allocate_potentials.
    
    When the potentials in the Python interface are imported
    to the Fortran core, the target specifiers (elements, tags, indices)
    are permutated to create all equivalent potentials.
    That is, if we have a potential for Si-O, both Si-O and O-Si
    potentials are created. This is because the energy and
    force calculation loops only deal with atom pairs A-B once
    (so only A-B or B-A is considered, not both) and if, say,
    the loop only finds an O-Si pair, it is important to apply
    the Si-O interaction also on that pair.
    In some cases, such as with the bond-bending potential affecting
    triplets (A-B-C), the interaction is not symmetric for all the atoms.
    Therefore it is necessary to also store the original targets of
    the potential as specified in the Python interface. These are
    to be given in the 'orig_*' lists.
    
    called from PyInterface: :func:`add_potential`
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets (interacting bodies)
    n_params: integer  *intent(in)*    *scalar*  
        number of parameters
    pot_name: character(len=*)  *intent(in)*    *scalar*  
        potential names
    parameters: double precision  *intent(in)*    *size(n_params)*  
        numeric parameters
    cutoff: double precision  *intent(in)*    *scalar*  
        interaction hard cutoff
    smooth_cut: double precision  *intent(in)*    *scalar*  
        interaction soft cutoff
    elements: character(len=label_length)  *intent(in)*    *size(n_targets)*  
        atomic symbols specifying the elements the interaction acts on
    tags: integer  *intent(in)*    *size(n_targets)*  
        tags specifying the atoms the interaction acts on
    indices: integer  *intent(in)*    *size(n_targets)*  
        indices specifying the atoms the interaction acts on
    orig_elements: character(len=label_length)  *intent(in)*    *size(n_targets)*  
        original atomic symbols specifying the elements the interaction acts on
    orig_tags: integer  *intent(in)*    *size(n_targets)*  
        original tags specifying the atoms the interaction acts on
    orig_indices: integer  *intent(in)*    *size(n_targets)*  
        original indices specifying the atoms the interaction acts on
    pot_index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: core_allocate_bond_order_factors(n_bond_factors)

    Allocates pointers for storing bond order factors.
    
    called from PyInterface: :func:`allocate_bond_order_factors`
    

    Parameters:

    n_bond_factors: integer  *intent(in)*    *scalar*  
        
            
  .. function:: core_allocate_bond_order_storage(n_atoms, n_groups, n_factors)

    Allocates arrays for storing precalculated values of bond order
    factors and gradients.
    
    called from PyInterface: :func:`allocate_bond_order_factors`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    n_groups: integer  *intent(in)*    *scalar*  
        number of bond order groups
    n_factors: integer  *intent(in)*    *scalar*  
        number of bond order parameters
            
  .. function:: core_allocate_potentials(n_pots)

    Allocates pointers for storing potentials.
    
    called from PyInterface: :func:`allocate_potentials`
    

    Parameters:

    n_pots: integer  *intent(in)*    *scalar*  
        number of potentials
            
  .. function:: core_assign_bond_order_factor_indices()

    This routine finds for each atom the potentials for which the
    atom is an accepted target at the first position.
    First position here means that for instance in an A-B-C triplet.
    A is in first position.
    Being an accepted target means that the atom has the correct
    element.
    
    called from PyInterface: :func:`create_bond_order_factor_list`

            
  .. function:: core_assign_potential_indices()

    This routine finds for each atom the potentials for which the
    atom is an accepted target at the first position.
    First position here means that for instance in an A-B-C triplet.
    A is in first position.
    Being an accepted target means that the atom has the correct
    element, index or tag (one that the potential targets).
    
    called from PyInterface: :func:`create_potential_list`

            
  .. function:: core_calculate_bond_order_factors(n_atoms, group_index, total_bond_orders)

    Calculates the bond order sums of all atoms for the given group.
    
    For a factor such as
    
    .. math::
    
         b_i = f(\sum_j c_{ij})
    
    The routine calculates
    
    .. math::
    
         \sum_j c_{ij}.
    
    The full bond order factor is then obtained by applying the
    scaling function :math:`f`. This is done with
    :func:`core_post_process_bond_order_factors`.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    **total_bond_orders**: double precision  **intent(out)**    *size(n_atoms)*  
        the calculated bond order sums
            
  .. function:: core_calculate_bond_order_gradients(n_atoms, group_index, atom_index, raw_sums, total_gradient, for_factor)

    Returns the gradients of bond order factors.
    
    For a factor such as
    
    .. math::
    
         b_i = f(\sum_j c_{ij})
    
    The routine calculates
    
    .. math::
    
        \nabla_\alpha b_i = f'(\sum_j c_{ij}) \nabla_\alpha \sum_j c_{ij}.
    
    By default, the gradients of all factors :math:`i` are calculated with respect
    to moving the given atom :math:`\alpha`.
    If for_factor is .true., the gradients of the bond factor of the given
    atom are calculated with respect to moving all atoms.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom with respect to which the factors are differentiated (:math:`\alpha`), or the atoms whose factor is differentiated (:math:`i`) if for_factor is .true.
    raw_sums: double precision  *intent(in)*    *size(n_atoms)*  
        precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
    **total_gradient**: double precision  **intent(out)**    *size(3, n_atoms)*  
        the calculated bond order gradients :math:`\nabla_\alpha b_i`
    for_factor: logical  *intent(in)*    *scalar*  *optional*
        a switch for requesting the gradients for a given :math:`i` instead of a given :math:`\alpha`
            
  .. function:: core_calculate_bond_order_gradients_of_factor(n_atoms, group_index, atom_index, raw_sums, total_gradient)

    Returns the gradients of one bond order factor with respect to
    moving all atoms.
    
    This calls :func:`core_calculate_bond_order_gradients` with for_factor = .true.
    
    For a factor such as
    
    .. math::
    
         b_i = f(\sum_j c_{ij})
    
    The routine calculates
    
    .. math::
    
        \nabla_\alpha b_i = f'(\sum_j c_{ij}) \nabla_\alpha \sum_j c_{ij}.
    
    The gradients of the bond factor of the given
    atom :math:`i` are calculated with respect to moving all atoms :math:`\alpha`.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom whose factor is differentiated (:math:`i`)
    raw_sums: double precision  *intent(in)*    *size(n_atoms)*  
        precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
    **total_gradient**: double precision  **intent(out)**    *size(3, n_atoms)*  
        the calculated bond order gradients :math:`\nabla_\alpha b_i`
            
  .. function:: core_calculate_electronegativities(n_atoms, total_enegs)

    Calculates electronegativity forces acting on all atomic charges of the system.
    
    The routine calculates the electronegativities
    
    .. math::
    
       \chi_{\alpha} = -\frac{\partial V}{\partial q_\alpha}
    
    for all atoms :math:`\alpha`. This is done according to the
    the structure and potentials allocated in the core, so the
    routine does not accept arguments. Instead, the core modifying
    routines such as :func:`core_generate_atoms` must be called
    first to set up the calculation.
    
    called from PyInterface: :func:`calculate_electronegativities`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    **total_enegs**: double precision  **intent(out)**    *size(n_atoms)*  
        an array containing the calculated charge forces for all atoms
            
  .. function:: core_calculate_energy(n_atoms, total_energy)

    Calculates the total potential energy of the system.
    
    This is done according to the
    the structure and potentials allocated in the core, so the
    routine does not accept arguments. Instead, the core modifying
    routines such as :func:`core_generate_atoms` must be called
    first to set up the calculation.
    
    called from PyInterface: :func:`calculate_energy`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    **total_energy**: double precision  **intent(out)**    *scalar*  
        calculated total potential energy
            
  .. function:: core_calculate_forces(n_atoms, total_forces)

    Calculates forces acting on all atoms of the system.
    
    The routine calculates the potential gradient
    
    .. math::
    
       \mathbf{F}_\alpha = - \nabla_\alpha V
    
    for all atoms :math:`\alpha`. This is done according to the
    the structure and potentials allocated in the core, so the
    routine does not accept arguments. Instead, the core modifying
    routines such as :func:`core_generate_atoms` must be called
    first to set up the calculation.
    
    called from PyInterface: :func:`calculate_forces`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    **total_forces**: double precision  **intent(out)**    *size(3, n_atoms)*  
        an array containing the calculated forces for all atoms
            
  .. function:: core_clear_atoms()

    Deallocates the array of atoms in the core, if allocated.

            
  .. function:: core_clear_bond_order_factors()

    Deallocates pointers for bond order factors (the parameters)

            
  .. function:: core_clear_bond_order_storage()

    Deallocates pointers for bond order factors (the precalculated factor values).

            
  .. function:: core_clear_potentials()

    Deallocates pointers for potentials

            
  .. function:: core_create_cell(vectors, inverse, periodicity)

    Creates a supercell for containing the calculation geometry.
    
    called from PyInterface: :func:`create_cell`
    

    Parameters:

    vectors: double precision  *intent(in)*    *size(3, 3)*  
        A 3x3 matrix containing the vectors spanning the supercell. The first index runs over xyz and the second index runs over the three vectors.
    inverse: double precision  *intent(in)*    *size(3, 3)*  
        A 3x3 matrix containing the inverse matrix of the one given in vectors, i.e. :math:`A*B = I` for the two matrices. Since the latter represents a cell of non-zero volume, this inverse must exist. It is not tested that the given matrix actually is the inverse, the user must make sure it is.
    periodicity: logical  *intent(in)*    *size(3)*  
        A 3-element vector containing logical tags specifying if the system is periodic in the directions of the three vectors spanning the supercell.
            
  .. function:: core_create_neighbor_list(n_nbs, atom_index, neighbors, offsets)

    Assigns a precalculated neighbor list to a single atom of the given index.
    The neighbor list must be precalculated, this method only
    stores them in the core. The list must contain
    an array storing the indices of the neighboring atoms
    as well as the supercell offsets. The offsets are integer
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
    
    called from PyInterface: :func:`create_neighbor_list`
    

    Parameters:

    n_nbs: integer  *intent(in)*    *scalar*  
        number of neighbors
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom for which the neighbor list is created
    neighbors: integer  *intent(in)*    *size(n_nbs)*  
        An array containing the indices of the neighboring atoms
    offsets: integer  *intent(in)*    *size(3, n_nbs)*  
        An array containing vectors specifying the offsets of the neighbors in periodic systems.
            
  .. function:: core_empty_bond_order_gradient_storage(index)

    Clears bond order factor gradients (the precalculated gradient values)
    but does not deallocate the arrays.
    If an index is given, then only that column is emptied.
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  *optional*
        the column to be emptied
            
  .. function:: core_empty_bond_order_storage()

    Clears bond order factors (the precalculated factor values)
    but does not deallocate the arrays.

            
  .. function:: core_fill_bond_order_storage(n_atoms)

    Fills the storage for bond order factors and bond order sums.
    This is meant to be called in the beginning of force and energy
    evaluation. The routine calculates all bond order factors
    (in parallel, if run in MPI) and stores them. Then during the
    energy or force calculation, it is sufficient to just
    look up the needed values in the arrays.
    The routine does not calculate and store bond factor gradients.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
            
  .. function:: core_generate_atoms(n_atoms, masses, charges, positions, momenta, tags, elements)

    Creates the atomic particles by invoking a subroutine in the geometry module.
    
    called from PyInterface: :func:`create_atoms`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    masses: double precision  *intent(in)*    *size(n_atoms)*  
        masses of atoms
    charges: double precision  *intent(in)*    *size(n_atoms)*  
        electric charges of atoms
    positions: double precision  *intent(in)*    *size(3, n_atoms)*  
        coordinates of atoms
    momenta: double precision  *intent(in)*    *size(3, n_atoms)*  
        momenta of atoms
    tags: integer  *intent(in)*    *size(n_atoms)*  
        numeric tags for the atoms
    elements: character(len=label_length)  *intent(in)*    *size(n_atoms)*  
        atomic symbols of the atoms
            
  .. function:: core_get_bond_order_factor_of_atom(n_atoms, group_index, atom_index, bond_order_factor)

    Returns the bond order factors of the given atom for the given group.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        index for the bond order factor group
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom whose bond order factor is returned
    **bond_order_factor**: double precision  **intent(out)**    *scalar*  
        the calculated bond order factor
            
  .. function:: core_get_bond_order_factors(n_atoms, group_index, bond_order_factors)

    Returns the bond order factors of all atoms for the given group.
    The routines tries to find the values in the stored precalculated
    values first if use_saved_bond_order_factors is true, and saves
    the calculated values if it does not find them.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        index for the bond order factor group
    **bond_order_factors**: double precision  **intent(out)**    *size(n_atoms)*  
        the calculated bond order factors
            
  .. function:: core_get_bond_order_gradients(n_atoms, group_index, atom_index, slot_index, bond_order_gradients)

    Returns the gradients of the bond order factor of the given atom
    with respect to moving all atoms for the given group.
    The routine tries to find the values in the stored precalculated
    values first if use_saved_bond_order_factors is true, and saves
    the calculated values if it does not find them.
    
    The slot index is the index of the atom in the interaction being
    evaluated (so for a triplet A-B-C, A would have slot 1, B slot 2,
    and C slot 3). This is only used for storing the values.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        index for the bond order factor group
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom whose bond order factor is differentiated
    slot_index: integer  *intent(in)*    *scalar*  
        index denoting the position of the atom in an interacting group (such as A-B-C triplet)
    **bond_order_gradients**: double precision  **intent(out)**    *size(1:3, n_atoms)*  
        the calculated gradients of the bond order factor
            
  .. function:: core_get_bond_order_sums(n_atoms, group_index, bond_order_sums)

    Returns the bond order sums of all atoms for the given group.
    By 'bond order sum', we mean the summation of local terms
    without per atom scaling. E.g., for :math:`b_i = 1 + \sum c_{ij}`,
    :math:`\sum c_{ij}` is the sum.
    The routines tries to find the values in the stored precalculated
    values first if use_saved_bond_order_factors is true, and saves
    the calculated values if it does not find them.

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        index for the bond order factor group
    **bond_order_sums**: double precision  **intent(out)**    *size(n_atoms)*  
        the calculated bond order sums
            
  .. function:: core_get_cell_vectors(vectors)

    Returns the vectors defining the supercell stored in the core.
    
    called from PyInterface: :func:`get_cell_vectors`
    

    Parameters:

    **vectors**: double precision  **intent(out)**    *size(3, 3)*  
        A 3x3 matrix containing the vectors spanning the supercell. The first index runs over xyz and the second index runs over the three vectors.
            
  .. function:: core_get_ewald_energy(real_cut, reciprocal_cut, sigma, epsilon, energy)

    Debug routine for Ewald

    Parameters:

    real_cut: double precision  *intent(in)*    *scalar*  
        
    reciprocal_cut: integer  *intent(in)*    *size(3)*  
        
    sigma: double precision  *intent(in)*    *scalar*  
        
    epsilon: double precision  *intent(in)*    *scalar*  
        
    **energy**: double precision  **intent(out)**    *scalar*  
        
            
  .. function:: core_get_number_of_atoms(n_atoms)

    Returns the number of atoms in the array allocated in the core.
    
    called from PyInterface: :func:`get_number_of_atoms`
    

    Parameters:

    **n_atoms**: integer  **intent(out)**    *scalar*  
        number of atoms
            
  .. function:: core_post_process_bond_order_factors(n_atoms, group_index, raw_sums, total_bond_orders)

    Bond-order post processing, i.e., application of per-atom scaling functions.
    
    By post processing, we mean any operations done after calculating the
    sum of pair- and many-body terms. That is, if a factor is, say,
    
    .. math::
    
         b_i = f(\sum_j c_{ij}) = 1 + \sum_j c_{ij},
    
    the :math:`\sum_j c_{ij}` would have been calculated already
    (with :func:`core_calculate_bond_order_factors`)
    and the operation :math:`f(x) = 1 + x`
    remains to be carried out.
    The post processing is done per atom regardless of if the
    bond factor is of a pair or many body type.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    raw_sums: double precision  *intent(in)*    *size(n_atoms)*  
        precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
    **total_bond_orders**: double precision  **intent(out)**    *size(n_atoms)*  
        the calculated bond order factors :math:`b_i`
            
  .. function:: core_post_process_bond_order_gradients(n_atoms, group_index, raw_sums, raw_gradients, total_bond_gradients, mpi_split)

    Bond-order post processing, i.e., application of per-atom scaling functions.
    This routine does the scaling for all bond factors with the given
    bond order sums and gradients of these sums.
    
    By post processing, we mean any operations done after calculating the
    sum of pair- and many-body terms. That is, if a factor is, say,
    
    .. math::
    
         b_i = f(\sum_j c_{ij}) = 1 + \sum_j c_{ij},
    
    the :math:`\sum_j c_{ij}` would have been calculated already and the
    operation :math:`f(x) = 1 + x` remains to be carried out.
    The post processing is done per atom regardless of if the
    bond factor is of a pair or many body type.
    
    For gradients, one needs to evaluate
    
    .. math::
    
        \nabla_\alpha b_i = f'(\sum_j c_{ij}) \nabla_\alpha \sum_j c_{ij}
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    raw_sums: double precision  *intent(in)*    *size(n_atoms)*  
        precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example
    raw_gradients: double precision  *intent(in)*    *size(3, n_atoms)*  
        precalculated gradients of bond order sums, :math:`\nabla_\alpha \sum_j c_{ij}`, in the above example
    **total_bond_gradients**: double precision  **intent(out)**    *size(3, n_atoms)*  
        the calculated bond order gradients :math:`\nabla_\alpha b_i`
    mpi_split: logical  *intent(in)*    *scalar*  *optional*
        A switch for enabling MPI parallelization. By default the routine is sequential since the calculation may be called from within an already parallelized routine.
            
  .. function:: core_post_process_bond_order_gradients_of_factor(n_atoms, group_index, atom_index, raw_sum, raw_gradients, total_bond_gradients, mpi_split)

    Bond-order post processing, i.e., application of per-atom scaling functions.
    This routine does the scaling for the bond order factor of the given atom
    with respect to moving all atoms
    with the given bond order sum for the factor and
    the gradients of the sum with respect to moving all atoms.
    
    By post processing, we mean any operations done after calculating the
    sum of pair- and many-body terms. That is, if a factor is, say,
    
    .. math::
    
         b_i = f(\sum_j c_{ij}) = 1 + \sum_j c_{ij},
    
    the :math:`\sum_j c_{ij}` would have been calculated already and the operation :math:`f(x) = 1 + x`
    remains to be carried out.
    The post processing is done per atom regardless of if the
    bond factor is of a pair or many body type.
    
    For gradients, one needs to evaluate
    
    .. math::
    
        \nabla_\alpha b_i = f'(\sum_j c_{ij}) \nabla_\alpha \sum_j c_{ij}
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    atom_index: integer  *intent(in)*    *scalar*  
        the index of the atom whose factor is differentiated (:math:`i`)
    raw_sum: double precision  *intent(in)*    *scalar*  
        precalculated bond order sum for the given atom, :math:`\sum_j c_{ij}`, in the above example
    raw_gradients: double precision  *intent(in)*    *size(3, n_atoms)*  
        precalculated gradients of bond order sums, :math:`\nabla_\alpha \sum_j c_{ij}`, in the above example
    **total_bond_gradients**: double precision  **intent(out)**    *size(3, n_atoms)*  
        the calculated bond order gradients :math:`\nabla_\alpha b_i`
    mpi_split: logical  *intent(in)*    *scalar*  *optional*
        A switch for enabling MPI parallelization. By default the routine is sequential since the calculation may be called from within an already parallelized routine.
            
  .. function:: core_release_all_memory()

    Release all allocated pointer arrays in the core.

            
  .. function:: core_set_ewald_parameters(n_atoms, real_cut, reciprocal_cut, sigma, epsilon, scaler)

    Sets the parameters for Ewald summation in the core.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        
    real_cut: double precision  *intent(in)*    *scalar*  
        the real-space cutoff
    reciprocal_cut: integer  *intent(in)*    *size(3)*  
        the k-space cutoffs
    sigma: double precision  *intent(in)*    *scalar*  
        the split parameter
    epsilon: double precision  *intent(in)*    *scalar*  
        electric constant
    scaler: double precision  *intent(in)*    *size(n_atoms)*  
        scaling factors for the individual charges
            
  .. function:: core_update_atom_charges(n_atoms, charges)

    Updates the charges of atomic particles.
    
    called from PyInterface: :func:`update_atom_charges`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    charges: double precision  *intent(in)*    *size(n_atoms)*  
        new charges for the atoms
            
  .. function:: core_update_atom_coordinates(n_atoms, positions, momenta)

    Updates the positions and momenta of atomic particles.
    
    called from PyInterface: :func:`update_atom_coordinates`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    positions: double precision  *intent(in)*    *size(3, n_atoms)*  
        new coordinates for the atoms
    momenta: double precision  *intent(in)*    *size(3, n_atoms)*  
        new momenta for the atoms
            
  .. function:: list_atoms()

    Prints some information on the atoms stored in the core in stdout.

            
  .. function:: list_bonds()

    Prints some information on the bond order factors stored in the core in stdout.

            
  .. function:: list_cell()

    Prints some information on the supercell stored in the core in stdout.

            
  .. function:: list_interactions()

    Prints some information on the potentials stored in the core in stdout.
