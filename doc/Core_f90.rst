
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
    - :data:`bo_factors`
    - :data:`bo_gradients`
    - :data:`bo_scaling`
    - :data:`bo_sums`
    - :data:`bo_temp`
    - :data:`bond_factors`
    - :data:`bond_factors_allocated`
    - :data:`bond_storage_allocated`
    - :data:`cell`
    - :data:`electronegativity_evaluation_index`
    - :data:`energy_evaluation_index`
    - :data:`evaluate_ewald`
    - :data:`ewald_allocated`
    - :data:`ewald_cutoff`
    - :data:`ewald_epsilon`
    - :data:`ewald_k_cutoffs`
    - :data:`ewald_k_radius`
    - :data:`ewald_scaler`
    - :data:`ewald_sigma`
    - :data:`force_evaluation_index`
    - :data:`group_index_save_slot`
    - :data:`interactions`
    - :data:`multipliers`
    - :data:`n_bond_factors`
    - :data:`n_interactions`
    - :data:`n_multi`
    - :data:`n_nbs`
    - :data:`n_saved_bond_order_factors`
    - :data:`number_of_atoms`
    - :data:`potentials_allocated`
    - :data:`saved_bond_order_factors`
    - :data:`saved_bond_order_gradients`
    - :data:`saved_bond_order_sums`
    - :data:`saved_bond_order_virials`
    - :data:`temp_enegs`
    - :data:`temp_forces`
    - :data:`temp_gradient`
    - :data:`total_n_nbs`
    - :data:`use_saved_bond_order_factors`
    - :data:`use_saved_bond_order_gradients`

    List of subroutines in pysic_core
    ---------------------------------
        
    - :func:`core_add_bond_order_factor`
    - :func:`core_add_bond_order_forces`
    - :func:`core_add_pair_bond_order_forces`
    - :func:`core_add_potential`
    - :func:`core_allocate_bond_order_factors`
    - :func:`core_allocate_bond_order_storage`
    - :func:`core_allocate_potentials`
    - :func:`core_assign_bond_order_factor_indices`
    - :func:`core_assign_potential_indices`
    - :func:`core_build_neighbor_lists`
    - :func:`core_calculate_bond_order_factors`
    - :func:`core_calculate_bond_order_gradients`
    - :func:`core_calculate_bond_order_gradients_of_factor`
    - :func:`core_calculate_electronegativities`
    - :func:`core_calculate_energy`
    - :func:`core_calculate_forces`
    - :func:`core_calculate_pair_bond_order_factor`
    - :func:`core_calculate_pair_bond_order_gradients`
    - :func:`core_clear_atoms`
    - :func:`core_clear_bond_order_factors`
    - :func:`core_clear_bond_order_storage`
    - :func:`core_clear_ewald_arrays`
    - :func:`core_clear_potential_multipliers`
    - :func:`core_clear_potentials`
    - :func:`core_create_cell`
    - :func:`core_create_neighbor_list`
    - :func:`core_create_space_partitioning`
    - :func:`core_debug_dump`
    - :func:`core_empty_bond_order_gradient_storage`
    - :func:`core_empty_bond_order_storage`
    - :func:`core_evaluate_local_doublet`
    - :func:`core_evaluate_local_doublet_electronegativities`
    - :func:`core_evaluate_local_doublet_electronegativities_B`
    - :func:`core_evaluate_local_doublet_energy`
    - :func:`core_evaluate_local_doublet_energy_B`
    - :func:`core_evaluate_local_doublet_forces`
    - :func:`core_evaluate_local_doublet_forces_B`
    - :func:`core_evaluate_local_quadruplet`
    - :func:`core_evaluate_local_quadruplet_B`
    - :func:`core_evaluate_local_singlet`
    - :func:`core_evaluate_local_triplet`
    - :func:`core_evaluate_local_triplet_B`
    - :func:`core_fill_bond_order_storage`
    - :func:`core_generate_atoms`
    - :func:`core_get_bond_order_factor_of_atom`
    - :func:`core_get_bond_order_factors`
    - :func:`core_get_bond_order_gradients`
    - :func:`core_get_bond_order_sums`
    - :func:`core_get_cell_vectors`
    - :func:`core_get_ewald_energy`
    - :func:`core_get_neighbor_list_of_atom`
    - :func:`core_get_number_of_atoms`
    - :func:`core_get_number_of_neighbors`
    - :func:`core_loop_over_local_interactions`
    - :func:`core_post_process_bond_order_factors`
    - :func:`core_post_process_bond_order_gradients`
    - :func:`core_post_process_bond_order_gradients_of_factor`
    - :func:`core_post_process_pair_bond_order_factor`
    - :func:`core_post_process_pair_bond_order_gradients`
    - :func:`core_release_all_memory`
    - :func:`core_set_ewald_parameters`
    - :func:`core_update_atom_charges`
    - :func:`core_update_atom_coordinates`
    - :func:`expand_neighbor_storage`
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
    
  .. data:: bo_factors

    double precision  *pointer*  *size(:)*    
    
    
    
  .. data:: bo_gradients

    double precision  *pointer*  *size(:, :, :)*    
    
    
    
  .. data:: bo_scaling

    logical  *pointer*  *size(:)*    
    
    
    
  .. data:: bo_sums

    double precision  *pointer*  *size(:)*    
    
    
    
  .. data:: bo_temp

    double precision  *pointer*  *size(:)*    
    
    
    
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
    
  .. data:: electronegativity_evaluation_index

    integer    *scalar*  *parameter*  

    *initial value* = 3
    
    
    
  .. data:: energy_evaluation_index

    integer    *scalar*  *parameter*  

    *initial value* = 1
    
    
    
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
    
    
    
  .. data:: ewald_k_radius

    double precision    *scalar*    
    
    
    
  .. data:: ewald_scaler

    double precision  *pointer*  *size(:)*    
    
    
    
  .. data:: ewald_sigma

    double precision    *scalar*    
    
    
    
  .. data:: force_evaluation_index

    integer    *scalar*  *parameter*  

    *initial value* = 2
    
    
    
  .. data:: group_index_save_slot

    integer  *pointer*  *size(:)*    
    
    
    
  .. data:: interactions

    type(potential)  *pointer*  *size(:)*    
    
    an array of :data:`potential` objects representing the interactions
    
  .. data:: multipliers

    type(potential)  *allocatable*  *size(:)*    
    
    a temporary array for storing multiplying potentials before associating them with a master potential
    
  .. data:: n_bond_factors

    integer    *scalar*    

    *initial value* = 0
    
    
    
  .. data:: n_interactions

    integer    *scalar*    

    *initial value* = 0
    
    number of potentials
    
  .. data:: n_multi

    integer    *scalar*    

    *initial value* = 0
    
    number of temporary product potentials
    
  .. data:: n_nbs

    integer  *pointer*  *size(:)*    
    
    
    
  .. data:: n_saved_bond_order_factors

    integer    *scalar*    

    *initial value* = 0
    
    number of saved bond order factors
    
  .. data:: number_of_atoms

    integer    *scalar*    
    
    
    
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
    
  .. data:: saved_bond_order_virials

    double precision  *pointer*  *size(:, :, :)*    
    
    Array for storing calculated bond order virials. Indexing: (xyz, group_index_save_slot(group index), target index)
    
  .. data:: temp_enegs

    double precision  *pointer*  *size(:)*    
    
    
    
  .. data:: temp_forces

    double precision  *pointer*  *size(:, :)*    
    
    
    
  .. data:: temp_gradient

    double precision  *pointer*  *size(:, :, :)*    
    
    
    
  .. data:: total_n_nbs

    integer  *pointer*  *size(:)*    
    
    
    
  .. data:: use_saved_bond_order_factors

    logical    *scalar*    

    *initial value* = .false.
    
    Logical tag which enables / disables bond order saving. If true, bond order calculation routines try to find the precalculated factors in the saved bond order arrays instead of calculating.
    
  .. data:: use_saved_bond_order_gradients

    integer  *pointer*  *size(:, :)*    
    
    Array storing the atom index of the bond gradient stored for indices (group index, target index). Since gradients are needed for all factors (N) with respect to moving all atoms (N), storing them all would require an N x N matrix. Therefore only some are stored. This array is used for searching the stroage to see if the needed gradient is there or needs to be calculated.
    

Full documentation of subroutines in pysic_core
-----------------------------------------------
        
        
            
  .. function:: core_add_bond_order_factor(n_targets, n_params, n_split, bond_name, parameters, param_split, cutoff, smooth_cut, elements, orig_elements, group_index, success)

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
    **success**: logical  **intent(out)**    *scalar*  
        logical tag specifying if creation of the factor succeeded
            
  .. function:: core_add_bond_order_forces(group_index, atom_index, prefactor, forces, stress)


    Parameters:

    group_index: integer  *intent(in)*    *scalar*  
        
    atom_index: integer  *intent(in)*    *scalar*  
        
    prefactor: double precision  *intent(in)*    *scalar*  
        
    **forces**: double precision  **intent(inout)**    *size(:, :)*  
        
    **stress**: double precision  **intent(inout)**    *size(6)*  
        
            
  .. function:: core_add_pair_bond_order_forces(index1, index2, prefactor, separation, direction, distance, group_index, pair_bo_sums, pair_bo_factors, forces, stress)

    Evaluates the local force affecting two atoms from bond order factors.
    

    Parameters:

    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    prefactor: double precision  *intent(in)*    *scalar*  
        
    separation: double precision  *intent(in)*    *size(3)*  
        
    direction: double precision  *intent(in)*    *size(3)*  
        
    distance: double precision  *intent(in)*    *scalar*  
        
    group_index: integer  *intent(in)*    *scalar*  
        
    pair_bo_sums: double precision  *intent(in)*    *size(2)*  
        
    pair_bo_factors: double precision  *intent(in)*    *size(2)*  
        
    **forces**: double precision  **intent(inout)**    *size(:, :)*  
        calculated forces
    **stress**: double precision  **intent(inout)**    *size(6)*  
        calculated stress
            
  .. function:: core_add_potential(n_targets, n_params, pot_name, parameters, cutoff, smooth_cut, elements, tags, indices, orig_elements, orig_tags, orig_indices, pot_index, is_multiplier, success)

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
    
    If product potentials are created, all but the first one of the potentials
    are created with ``is_multiplier == .true.``. This leads to the potentials
    being stored in the global temporary array ``multipliers``. The last potential
    of a group should be created with ``is_multiplier = .false.`` and the stored
    multipliers are attached to it. The list of multipliers is not cleared automatically,
    since usually one creates copies of the same potential with permutated targets and all
    of these need the same multipiers.
    Instead the multipliers are cleared with a call of :func:`clear_potential_multipliers`.
    
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
    is_multiplier: logical  *intent(in)*    *scalar*  
        logical tag specifying if this potential should be treated as a multiplier
    **success**: logical  **intent(out)**    *scalar*  
        logical tag specifying if creation of the potential succeeded
            
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

            
  .. function:: core_build_neighbor_lists(cutoffs)

    Builds the neighbor lists in the core.
    The simulation cell must be partitioned with :func:`core_create_space_partitioning`
    before this routine can be called.
    

    Parameters:

    cutoffs: double precision  *intent(in)*    *size(:)*  
        list of cutoffs, atom by atom
            
  .. function:: core_calculate_bond_order_factors(group_index, total_bond_orders)

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

    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    **total_bond_orders**: double precision  **intent(inout)**    *size(:)*  
        the calculated bond order sums
            
  .. function:: core_calculate_bond_order_gradients(group_index, atom_index, raw_sums, total_gradient, total_virial, for_factor)

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

    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom with respect to which the factors are differentiated (:math:`\alpha`), or the atoms whose factor is differentiated (:math:`i`) if for_factor is .true.
    raw_sums: double precision  *intent(in)*    *size(:)*  
        precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
    **total_gradient**: double precision  **intent(inout)**    *size(:, :)*  
        the calculated bond order gradients :math:`\nabla_\alpha b_i`
    **total_virial**: double precision  **intent(inout)**    *size(6)*  
        the components of the virial due to the bond order gradients
    for_factor: logical  *intent(in)*    *scalar*  *optional*
        a switch for requesting the gradients for a given :math:`i` instead of a given :math:`\alpha`
            
  .. function:: core_calculate_bond_order_gradients_of_factor(group_index, atom_index, raw_sums, total_gradient, total_virial)

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

    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom whose factor is differentiated (:math:`i`)
    raw_sums: double precision  *intent(in)*    *size(:)*  
        precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
    **total_gradient**: double precision  **intent(inout)**    *size(:, :)*  
        the calculated bond order gradients :math:`\nabla_\alpha b_i`
    **total_virial**: double precision  **intent(inout)**    *size(6)*  
        the components of the virial due to the bond order gradient
            
  .. function:: core_calculate_electronegativities(total_enegs)

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

    **total_enegs**: double precision  **intent(inout)**    *size(:)*  
        an array containing the calculated charge forces for all atoms
            
  .. function:: core_calculate_energy(total_energy)

    Calculates the total potential energy of the system.
    
    This is done according to the
    the structure and potentials allocated in the core, so the
    routine does not accept arguments. Instead, the core modifying
    routines such as :func:`core_generate_atoms` must be called
    first to set up the calculation.
    
    called from PyInterface: :func:`calculate_energy`
    

    Parameters:

    **total_energy**: double precision  **intent(out)**    *scalar*  
        calculated total potential energy
            
  .. function:: core_calculate_forces(total_forces, total_stress)

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

    **total_forces**: double precision  **intent(inout)**    *size(:, :)*  
        an array containing the calculated forces for all atoms
    **total_stress**: double precision  **intent(inout)**    *size(6)*  
        as array containing the calculated stress tensor
            
  .. function:: core_calculate_pair_bond_order_factor(atom_pair, separation, distance, direction, group_index, bond_order_sum)

    Calculates the bond order sum for a given pair of atoms for the given group.
    
    For a factor such as
    
    .. math::
    
         b_ij = f(\sum_k c_{ijk})
    
    The routine calculates
    
    .. math::
    
         \sum_k c_{ijk}.
    
    The full bond order factor is then obtained by applying the
    scaling function :math:`f`. This is done with
    :func:`core_post_process_bond_order_factors`.
    

    Parameters:

    atom_pair: integer  *intent(in)*    *size(2)*  
        
    separation: double precision  *intent(in)*    *size(3)*  
        
    distance: double precision  *intent(in)*    *scalar*  
        
    direction: double precision  *intent(in)*    *size(3)*  
        
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    **bond_order_sum**: double precision  **intent(out)**    *size(2)*  
        the calculated bond order sums
            
  .. function:: core_calculate_pair_bond_order_gradients(atom_pair, separation, distance, direction, group_index, raw_sums, total_gradient, total_virial)

    Returns the gradients of a pair bond order factor.
    
    For a factor such as
    
    .. math::
    
         b_{ij} = f(\sum_k c_{ijk})
    
    The routine calculates
    
    .. math::
    
        \nabla_\alpha b_{ij} = f'(\sum_k c_{ijk}) \nabla_\alpha \sum_k c_{ijk}.
    
    By default, the gradients the factor :math:`ij` is calculated with respect
    to moving all atoms :math:`\alpha`.
    

    Parameters:

    atom_pair: integer  *intent(in)*    *size(2)*  
        
    separation: double precision  *intent(in)*    *size(3)*  
        
    distance: double precision  *intent(in)*    *scalar*  
        
    direction: double precision  *intent(in)*    *size(3)*  
        
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    raw_sums: double precision  *intent(in)*    *size(2)*  
        precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
    **total_gradient**: double precision  **intent(inout)**    *size(:, :, :)*  
        the calculated bond order gradients :math:`\nabla_\alpha b_i`
    **total_virial**: double precision  **intent(inout)**    *size(6, 2)*  
        the components of the virial due to the bond order gradient
            
  .. function:: core_clear_atoms()

    Deallocates the array of atoms in the core, if allocated.

            
  .. function:: core_clear_bond_order_factors()

    Deallocates pointers for bond order factors (the parameters)

            
  .. function:: core_clear_bond_order_storage()

    Deallocates pointers for bond order factors (the precalculated factor values).

            
  .. function:: core_clear_ewald_arrays()


            
  .. function:: core_clear_potential_multipliers()


            
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
            
  .. function:: core_create_neighbor_list(n_nbors, atom_index, neighbors, offsets)

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

    n_nbors: integer  *intent(in)*    *scalar*  
        
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom for which the neighbor list is created
    neighbors: integer  *intent(in)*    *size(n_nbors)*  
        An array containing the indices of the neighboring atoms
    offsets: integer  *intent(in)*    *size(3, n_nbors)*  
        An array containing vectors specifying the offsets of the neighbors in periodic systems.
            
  .. function:: core_create_space_partitioning(max_cutoff)

    Partitions the simulation volume in subvolumes for fast neighbor searching
    

    Parameters:

    max_cutoff: double precision  *intent(in)*    *scalar*  
        the maximum cutoff radius for neighbor search
            
  .. function:: core_debug_dump(forces)

    Write atomic coordinates and other info in a file.
    This is only for debugging.

    Parameters:

    forces: double precision  *intent(in)*    *size(:, :)*  
        
            
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

            
  .. function:: core_evaluate_local_doublet(n_atoms, atom_doublet, index1, index2, test_index1, interaction_indices, separations, directions, distances, calculation_type, energy, forces, enegs, stress, many_bodies_found)

    Evaluates the interactions affecting two atoms.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        total number of atoms in the system
    atom_doublet: type(atom)  *intent(in)*    *size(2)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2
    interaction_indices: integer  *intent()*  *pointer*  *size(:)*  
        the interactions targeting the given atoms
    separations: double precision  *intent(in)*    *size(3, 1)*  
        distance vector from 1 to 2, as an array
    directions: double precision  *intent(in)*    *size(3, 1)*  
        unit vector from 1 to 2, as an array
    distances: double precision  *intent(in)*    *size(1)*  
        distance from 1 to 2, as an array
    calculation_type: integer  *intent(in)*    *scalar*  
        the type of information requested
    **energy**: double precision  **intent(inout)**    *scalar*  
        calculated energy
    **forces**: double precision  **intent(inout)**    *size(3, n_atoms)*  
        calculated forces
    **enegs**: double precision  **intent(inout)**    *size(n_atoms)*  
        calculated electronegativities
    **stress**: double precision  **intent(inout)**    *size(6)*  
        calculated stress
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
            
  .. function:: core_evaluate_local_doublet_electronegativities(n_atoms, atom_doublet, index1, index2, test_index1, interaction_indices, separations, directions, distances, enegs, many_bodies_found)

    Evaluates the local electronegativity affecting two atoms.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        
    atom_doublet: type(atom)  *intent(in)*    *size(2)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2
    interaction_indices: integer  *intent()*  *pointer*  *size(:)*  
        the interactions targeting the given atoms
    separations: double precision  *intent(in)*    *size(3, 1)*  
        distance vector from 1 to 2, as an array
    directions: double precision  *intent(in)*    *size(3, 1)*  
        unit vector from 1 to 2, as an array
    distances: double precision  *intent(in)*    *size(1)*  
        distance from 1 to 2, as an array
    **enegs**: double precision  **intent(inout)**    *size(n_atoms)*  
        calculated electronegativities
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
            
  .. function:: core_evaluate_local_doublet_electronegativities_B(atom_doublet, index1, index2, test_index1, interaction_indices, separations, directions, distances, enegs, many_bodies_found, manybody_indices, n_manybody)

    Evaluates the local electronegativity affecting two atoms. (Rearranged internally.)
    

    Parameters:

    atom_doublet: type(atom)  *intent(in)*    *size(2)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2
    interaction_indices: integer  *intent()*  *pointer*  *size(:)*  
        the interactions targeting the given atoms
    separations: double precision  *intent(in)*    *size(3, 1)*  
        distance vector from 1 to 2, as an array
    directions: double precision  *intent(in)*    *size(3, 1)*  
        unit vector from 1 to 2, as an array
    distances: double precision  *intent(in)*    *size(1)*  
        distance from 1 to 2, as an array
    **enegs**: double precision  **intent(inout)**    *size(:)*  
        calculated electronegativities
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
    manybody_indices: integer  *intent()*  *pointer*  *size(:)*  
        
    **n_manybody**: integer  **intent(out)**    *scalar*  
        
            
  .. function:: core_evaluate_local_doublet_energy(n_atoms, atom_doublet, index1, index2, test_index1, interaction_indices, separations, directions, distances, energy, many_bodies_found)

    Evaluates the local potential affecting two atoms.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        
    atom_doublet: type(atom)  *intent(in)*    *size(2)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2
    interaction_indices: integer  *intent()*  *pointer*  *size(:)*  
        the interactions targeting the given atoms
    separations: double precision  *intent(in)*    *size(3, 1)*  
        distance vector from 1 to 2, as an array
    directions: double precision  *intent(in)*    *size(3, 1)*  
        unit vector from 1 to 2, as an array
    distances: double precision  *intent(in)*    *size(1)*  
        distance from 1 to 2, as an array
    **energy**: double precision  **intent(inout)**    *scalar*  
        calculated energy
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
            
  .. function:: core_evaluate_local_doublet_energy_B(atom_doublet, index1, index2, test_index1, interaction_indices, separations, directions, distances, energy, many_bodies_found, manybody_indices, n_manybody)

    Evaluates the local potential affecting two atoms. (Rearranged internally compared to 'A'.)
    

    Parameters:

    atom_doublet: type(atom)  *intent(in)*    *size(2)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2
    interaction_indices: integer  *intent()*  *pointer*  *size(:)*  
        the interactions targeting the given atoms
    separations: double precision  *intent(in)*    *size(3, 1)*  
        distance vector from 1 to 2, as an array
    directions: double precision  *intent(in)*    *size(3, 1)*  
        unit vector from 1 to 2, as an array
    distances: double precision  *intent(in)*    *size(1)*  
        distance from 1 to 2, as an array
    **energy**: double precision  **intent(inout)**    *scalar*  
        calculated energy
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
    manybody_indices: integer  *intent()*  *pointer*  *size(:)*  
        
    **n_manybody**: integer  **intent(out)**    *scalar*  
        
            
  .. function:: core_evaluate_local_doublet_forces(n_atoms, atom_doublet, index1, index2, test_index1, interaction_indices, separations, directions, distances, forces, stress, many_bodies_found)

    Evaluates the local force affecting two atoms.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        total number of atoms in the system
    atom_doublet: type(atom)  *intent(in)*    *size(2)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2
    interaction_indices: integer  *intent()*  *pointer*  *size(:)*  
        the interactions targeting the given atoms
    separations: double precision  *intent(in)*    *size(3, 1)*  
        distance vector from 1 to 2, as an array
    directions: double precision  *intent(in)*    *size(3, 1)*  
        unit vector from 1 to 2, as an array
    distances: double precision  *intent(in)*    *size(1)*  
        distance from 1 to 2, as an array
    **forces**: double precision  **intent(inout)**    *size(3, n_atoms)*  
        calculated forces
    **stress**: double precision  **intent(inout)**    *size(6)*  
        calculated stress
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
            
  .. function:: core_evaluate_local_doublet_forces_B(atom_doublet, index1, index2, test_index1, interaction_indices, separations, directions, distances, forces, stress, many_bodies_found, manybody_indices, n_manybody)

    Evaluates the local force affecting two atoms. (Rearranged internally.)
    

    Parameters:

    atom_doublet: type(atom)  *intent(in)*    *size(2)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the interaction targets atom1; similarly for 2
    interaction_indices: integer  *intent()*  *pointer*  *size(:)*  
        the interactions targeting the given atoms
    separations: double precision  *intent(in)*    *size(3, 1)*  
        distance vector from 1 to 2, as an array
    directions: double precision  *intent(in)*    *size(3, 1)*  
        unit vector from 1 to 2, as an array
    distances: double precision  *intent(in)*    *size(1)*  
        distance from 1 to 2, as an array
    **forces**: double precision  **intent(inout)**    *size(:, :)*  
        calculated forces
    **stress**: double precision  **intent(inout)**    *size(6)*  
        calculated stress
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
    manybody_indices: integer  *intent()*  *pointer*  *size(:)*  
        
    **n_manybody**: integer  **intent(out)**    *scalar*  
        
            
  .. function:: core_evaluate_local_quadruplet(n_atoms, atom_quadruplet, index1, index2, index3, index4, test_index1, test_index2, test_index3, interaction_indices, separations, directions, distances, calculation_type, energy, forces, enegs, stress, many_bodies_found)

    Evaluates the interactions affecting four atoms.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        total number of atoms in the system
    atom_quadruplet: type(atom)  *intent(in)*    *size(4)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    index3: integer  *intent(in)*    *scalar*  
        index of the atom 3
    index4: integer  *intent(in)*    *scalar*  
        index of the atom 4
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2, 3
    test_index2: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2, 3
    test_index3: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2, 3
    interaction_indices: integer  *intent()*  *pointer*  *size(:)*  
        the interactions targeting the given atoms
    separations: double precision  *intent(in)*    *size(3, 3)*  
        distance vector from 1 to 2, 2 to 3 and 3 to 4 as an array
    directions: double precision  *intent(in)*    *size(3, 3)*  
        unit vector from 1 to 2, 2 to 3 and 3 to 4 as an array
    distances: double precision  *intent(in)*    *size(3)*  
        distance from 1 to 2, 2 to 3 and 3 to 4 as an array
    calculation_type: integer  *intent(in)*    *scalar*  
        the type of information requested
    **energy**: double precision  **intent(out)**    *scalar*  
        calculated energy
    **forces**: double precision  **intent(out)**    *size(3, n_atoms)*  
        calculated forces
    **enegs**: double precision  **intent(out)**    *size(n_atoms)*  
        calculated electronegativities
    **stress**: double precision  **intent(out)**    *size(6)*  
        calculated stress
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
            
  .. function:: core_evaluate_local_quadruplet_B(atom_quadruplet, index1, index2, index3, index4, test_index1, test_index2, test_index3, separations, directions, distances, calculation_type, energy, forces, enegs, stress, many_bodies_found, manybody_indices, n_manybody)

    Evaluates the interactions affecting four atoms. (Rearranged internally.)
    

    Parameters:

    atom_quadruplet: type(atom)  *intent(in)*    *size(4)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    index3: integer  *intent(in)*    *scalar*  
        index of the atom 3
    index4: integer  *intent(in)*    *scalar*  
        index of the atom 4
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2, 3
    test_index2: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2, 3
    test_index3: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2, 3
    separations: double precision  *intent(in)*    *size(3, 3)*  
        distance vector from 1 to 2, 2 to 3 and 3 to 4 as an array
    directions: double precision  *intent(in)*    *size(3, 3)*  
        unit vector from 1 to 2, 2 to 3 and 3 to 4 as an array
    distances: double precision  *intent(in)*    *size(3)*  
        distance from 1 to 2, 2 to 3 and 3 to 4 as an array
    calculation_type: integer  *intent(in)*    *scalar*  
        the type of information requested
    **energy**: double precision  **intent(inout)**    *scalar*  
        calculated energy
    **forces**: double precision  **intent(inout)**    *size(:, :)*  
        calculated forces
    **enegs**: double precision  **intent(inout)**    *size(:)*  
        calculated electronegativities
    **stress**: double precision  **intent(inout)**    *size(6)*  
        calculated stress
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
    manybody_indices: integer  *intent()*  *pointer*  *size(:)*  
        
    n_manybody: integer  *intent(in)*    *scalar*  
        
            
  .. function:: core_evaluate_local_singlet(index1, atom_singlet, interaction_indices, calculation_type, energy, forces, stress, enegs)

    Evaluates the local potential affecting a single atom
    

    Parameters:

    index1: integer  *intent(in)*    *scalar*  
        index of the atom
    atom_singlet: type(atom)  *intent(in)*    *scalar*  
        the atom that is targeted
    interaction_indices: integer  *intent()*  *pointer*  *size(:)*  
        the interactions targeting the given atom
    calculation_type: integer  *intent(in)*    *scalar*  
        specifies if we are evaluating the energy, forces, or electronegativities
    **energy**: double precision  **intent(inout)**    *scalar*  
        calculated energy
    **forces**: double precision  **intent(inout)**    *size(:, :)*  
        calculated forces
    **stress**: double precision  **intent(inout)**    *size(6)*  
        calculated stress
    **enegs**: double precision  **intent(inout)**    *size(:)*  
        calculated electronegativities
            
  .. function:: core_evaluate_local_triplet(n_atoms, atom_triplet, index1, index2, index3, test_index1, test_index2, interaction_indices, separations, directions, distances, calculation_type, energy, forces, enegs, stress, many_bodies_found)

    Evaluates the interactions affecting three atoms.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        total number of atoms in the system
    atom_triplet: type(atom)  *intent(in)*    *size(3)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    index3: integer  *intent(in)*    *scalar*  
        index of the atom 3
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2, 3
    test_index2: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2, 3
    interaction_indices: integer  *intent()*  *pointer*  *size(:)*  
        the interactions targeting the given atoms
    separations: double precision  *intent(in)*    *size(3, 2)*  
        distance vector from 1 to 2 and 2 to 3 as an array
    directions: double precision  *intent(in)*    *size(3, 2)*  
        unit vector from 1 to 2 and 2 to 3 as an array
    distances: double precision  *intent(in)*    *size(2)*  
        distance from 1 to 2 and 2 to 3 as an array
    calculation_type: integer  *intent(in)*    *scalar*  
        the type of information requested
    **energy**: double precision  **intent(out)**    *scalar*  
        calculated energy
    **forces**: double precision  **intent(out)**    *size(3, n_atoms)*  
        calculated forces
    **enegs**: double precision  **intent(out)**    *size(n_atoms)*  
        calculated electronegativities
    **stress**: double precision  **intent(out)**    *size(6)*  
        calculated stress
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
            
  .. function:: core_evaluate_local_triplet_B(atom_triplet, index1, index2, index3, test_index1, test_index2, separations, directions, distances, calculation_type, energy, forces, enegs, stress, many_bodies_found, manybody_indices, n_manybody)

    Evaluates the interactions affecting three atoms. (Rearranged internally.)
    

    Parameters:

    atom_triplet: type(atom)  *intent(in)*    *size(3)*  
        the atoms that are targeted
    index1: integer  *intent(in)*    *scalar*  
        index of the atom 1
    index2: integer  *intent(in)*    *scalar*  
        index of the atom 2
    index3: integer  *intent(in)*    *scalar*  
        index of the atom 3
    test_index1: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2, 3
    test_index2: integer  *intent(in)*    *scalar*  
        if 1, test if the ineraction targets atom1; similarly for 2, 3
    separations: double precision  *intent(in)*    *size(3, 2)*  
        distance vector from 1 to 2 and 2 to 3 as an array
    directions: double precision  *intent(in)*    *size(3, 2)*  
        unit vector from 1 to 2 and 2 to 3 as an array
    distances: double precision  *intent(in)*    *size(2)*  
        distance from 1 to 2 and 2 to 3 as an array
    calculation_type: integer  *intent(in)*    *scalar*  
        the type of information requested
    **energy**: double precision  **intent(inout)**    *scalar*  
        calculated energy
    **forces**: double precision  **intent(inout)**    *size(:, :)*  
        calculated forces
    **enegs**: double precision  **intent(inout)**    *size(:)*  
        calculated electronegativities
    **stress**: double precision  **intent(inout)**    *size(6)*  
        calculated stress
    **many_bodies_found**: logical  **intent(out)**    *scalar*  
        returns true if the loop finds an interaction with 3 or more targets
    manybody_indices: integer  *intent()*  *pointer*  *size(:)*  
        
    n_manybody: integer  *intent(in)*    *scalar*  
        
            
  .. function:: core_fill_bond_order_storage()

    Fills the storage for bond order factors and bond order sums.
    This is meant to be called in the beginning of force and energy
    evaluation. The routine calculates all bond order factors
    (in parallel, if run in MPI) and stores them. Then during the
    energy or force calculation, it is sufficient to just
    look up the needed values in the arrays.
    The routine does not calculate and store bond factor gradients.
    

            
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
            
  .. function:: core_get_bond_order_factor_of_atom(group_index, atom_index, bond_order_factor)

    Returns the bond order factors of the given atom for the given group.
    

    Parameters:

    group_index: integer  *intent(in)*    *scalar*  
        index for the bond order factor group
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom whose bond order factor is returned
    **bond_order_factor**: double precision  **intent(inout)**    *scalar*  
        the calculated bond order factor
            
  .. function:: core_get_bond_order_factors(group_index, bond_order_factors)

    Returns the bond order factors of all atoms for the given group.
    The routines tries to find the values in the stored precalculated
    values first if use_saved_bond_order_factors is true, and saves
    the calculated values if it does not find them.
    

    Parameters:

    group_index: integer  *intent(in)*    *scalar*  
        index for the bond order factor group
    **bond_order_factors**: double precision  **intent(inout)**    *size(:)*  
        the calculated bond order factors
            
  .. function:: core_get_bond_order_gradients(group_index, atom_index, slot_index, bond_order_gradients, bond_order_virial)

    Returns the gradients of the bond order factor of the given atom
    with respect to moving all atoms, for the given group.
    The routine tries to find the values in the stored precalculated
    values first if use_saved_bond_order_factors is true, and saves
    the calculated values if it does not find them.
    
    The slot index is the index of the atom in the interaction being
    evaluated (so for a triplet A-B-C, A would have slot 1, B slot 2,
    and C slot 3). This is only used for storing the values.
    

    Parameters:

    group_index: integer  *intent(in)*    *scalar*  
        index for the bond order factor group
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom whose bond order factor is differentiated
    slot_index: integer  *intent(in)*    *scalar*  
        index denoting the position of the atom in an interacting group (such as A-B-C triplet)
    **bond_order_gradients**: double precision  **intent(inout)**    *size(:, :)*  
        the calculated gradients of the bond order factor
    **bond_order_virial**: double precision  **intent(inout)**    *size(6)*  
        the components of the virial due to the bond order factors
            
  .. function:: core_get_bond_order_sums(group_index, bond_order_sums)

    Returns the bond order sums of all atoms for the given group.
    By 'bond order sum', we mean the summation of local terms
    without per atom scaling. E.g., for :math:`b_i = 1 + \sum c_{ij}`,
    :math:`\sum c_{ij}` is the sum.
    The routines tries to find the values in the stored precalculated
    values first if use_saved_bond_order_factors is true, and saves
    the calculated values if it does not find them.
    

    Parameters:

    group_index: integer  *intent(in)*    *scalar*  
        index for the bond order factor group
    **bond_order_sums**: double precision  **intent(inout)**    *size(:)*  
        the calculated bond order sums
            
  .. function:: core_get_cell_vectors(vectors)

    Returns the vectors defining the supercell stored in the core.
    
    called from PyInterface: :func:`get_cell_vectors`
    

    Parameters:

    **vectors**: double precision  **intent(out)**    *size(3, 3)*  
        A 3x3 matrix containing the vectors spanning the supercell. The first index runs over xyz and the second index runs over the three vectors.
            
  .. function:: core_get_ewald_energy(real_cut, k_cut, reciprocal_cut, sigma, epsilon, energy)

    Debug routine for Ewald

    Parameters:

    real_cut: double precision  *intent(in)*    *scalar*  
        
    k_cut: double precision  *intent(in)*    *scalar*  
        
    reciprocal_cut: integer  *intent(in)*    *size(3)*  
        
    sigma: double precision  *intent(in)*    *scalar*  
        
    epsilon: double precision  *intent(in)*    *scalar*  
        
    **energy**: double precision  **intent(out)**    *scalar*  
        
            
  .. function:: core_get_neighbor_list_of_atom(atom_index, n_neighbors, neighbors, offsets)

    Returns the list of neighbros for an atom
    

    Parameters:

    atom_index: integer  *intent(in)*    *scalar*  
        the index of the atom whose neighbors are returned
    n_neighbors: integer  *intent(in)*    *scalar*  
        the number of neighbors
    **neighbors**: integer  **intent(out)**    *size(n_neighbors)*  
        the indices of the neighboring atoms
    **offsets**: integer  **intent(out)**    *size(3, n_neighbors)*  
        the offsets for periodic boundaries
            
  .. function:: core_get_number_of_atoms(n_atoms)

    Returns the number of atoms in the array allocated in the core.
    
    called from PyInterface: :func:`get_number_of_atoms`
    

    Parameters:

    **n_atoms**: integer  **intent(out)**    *scalar*  
        number of atoms
            
  .. function:: core_get_number_of_neighbors(atom_index, n_neighbors)

    Returns the number of neighbors for an atom
    

    Parameters:

    atom_index: integer  *intent(in)*    *scalar*  
        the index of the atoms
    **n_neighbors**: integer  **intent(out)**    *scalar*  
        the number of neighbors
            
  .. function:: core_loop_over_local_interactions(calculation_type, total_energy, total_forces, total_enegs, total_stress)

    Loops over atoms, atomic pairs, atomic triplets, and atomic quadruplets
    and calculates the contributions from local potentials to energy, forces,
    or electronegativities. This routine is called from the routines
    
     - :func:`core_calculate_energy`
     - :func:`core_calculate_forces`
     - :func:`core_calculate_electronegaivities`
    

    Parameters:

    calculation_type: integer  *intent(in)*    *scalar*  
        index to specify if the loop calculates energies, forces, or e-negativities
    **total_energy**: double precision  **intent(inout)**    *scalar*  
        calculated energy
    **total_forces**: double precision  **intent(inout)**    *size(:, :)*  
        calculated forces
    **total_enegs**: double precision  **intent(inout)**    *size(:)*  
        calculated electronegativities
    **total_stress**: double precision  **intent(inout)**    *size(6)*  
        calculated stress
            
  .. function:: core_post_process_bond_order_factors(group_index, raw_sums, total_bond_orders)

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

    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    raw_sums: double precision  *intent(in)*    *size(:)*  
        precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
    **total_bond_orders**: double precision  **intent(inout)**    *size(:)*  
        the calculated bond order factors :math:`b_i`
            
  .. function:: core_post_process_bond_order_gradients(group_index, raw_sums, raw_gradients, total_bond_gradients, mpi_split)

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

    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    raw_sums: double precision  *intent(in)*    *size(:)*  
        precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example
    raw_gradients: double precision  *intent(in)*    *size(:, :)*  
        precalculated gradients of bond order sums, :math:`\nabla_\alpha \sum_j c_{ij}`, in the above example
    **total_bond_gradients**: double precision  **intent(inout)**    *size(:, :)*  
        the calculated bond order gradients :math:`\nabla_\alpha b_i`
    mpi_split: logical  *intent(in)*    *scalar*  *optional*
        A switch for enabling MPI parallelization. By default the routine is sequential since the calculation may be called from within an already parallelized routine.
            
  .. function:: core_post_process_bond_order_gradients_of_factor(group_index, atom_index, raw_sum, raw_gradients, total_bond_gradients, raw_virial, total_virial, mpi_split)

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

    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    atom_index: integer  *intent(in)*    *scalar*  
        the index of the atom whose factor is differentiated (:math:`i`)
    raw_sum: double precision  *intent(in)*    *scalar*  
        precalculated bond order sum for the given atom, :math:`\sum_j c_{ij}`, in the above example
    raw_gradients: double precision  *intent(in)*    *size(:, :)*  
        precalculated gradients of bond order sums, :math:`\nabla_\alpha \sum_j c_{ij}`, in the above example
    **total_bond_gradients**: double precision  **intent(inout)**    *size(:, :)*  
        the calculated bond order gradients :math:`\nabla_\alpha b_i`
    raw_virial: double precision  *intent(in)*    *size(6)*  
        the precalculated virial due to the bond order gradient
    **total_virial**: double precision  **intent(inout)**    *size(6)*  
        the scaled  virial due to the bond order gradient
    mpi_split: logical  *intent(in)*    *scalar*  *optional*
        A switch for enabling MPI parallelization. By default the routine is sequential since the calculation may be called from within an already parallelized routine.
            
  .. function:: core_post_process_pair_bond_order_factor(atom1, group_index, raw_sum, total_bond_order)

    Bond-order post processing, i.e., application of per-pair scaling functions.
    
    By post processing, we mean any operations done after calculating the
    sum of pair- and many-body terms. That is, if a factor is, say,
    
    .. math::
    
         b_{ij} = f(\sum_k c_{ijk}) = 1 + \sum_k c_{ijk},
    
    the :math:`\sum_k c_{ijk}` would have been calculated already
    (with :func:`core_calculate_pair_bond_order_factor`)
    and the operation :math:`f(x) = 1 + x`
    remains to be carried out.
    

    Parameters:

    atom1: type(atom)  *intent(in)*    *scalar*  
        the central atom of the pair bond order factor
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    raw_sum: double precision  *intent(in)*    *scalar*  
        precalculated bond order sum, :math:`\sum_k c_{ijk}`, in the above example.
    **total_bond_order**: double precision  **intent(out)**    *scalar*  
        the calculated bond order factor :math:`b_{ij}`
            
  .. function:: core_post_process_pair_bond_order_gradients(group_index, atom1, raw_sum, raw_gradients, total_bond_gradients, raw_virial, total_virial, mpi_split)

    Bond-order post processing, i.e., application of per-pair scaling functions.
    This routine does the scaling for the bond order factor of the given pair
    with respect to moving all atoms
    with the given bond order sum for the factor and
    the gradients of the sum with respect to moving all atoms.
    
    By post processing, we mean any operations done after calculating the
    sum of pair- and many-body terms. That is, if a factor is, say,
    
    .. math::
    
         b_{ij} = f(\sum_k c_{ijk}) = 1 + \sum_k c_{ijk},
    
    the :math:`\sum_k c_{ijk}` would have been calculated already and the operation :math:`f(x) = 1 + x`
    remains to be carried out.
    The post processing is done per pair.
    
    For gradients, one needs to evaluate
    
    .. math::
    
        \nabla_\alpha b_{ij} = f'(\sum_k c_{ijk}) \nabla_\alpha \sum_k c_{ijk}
    

    Parameters:

    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    atom1: type(atom)  *intent(in)*    *scalar*  
        the central atom of the pair bond order factor
    raw_sum: double precision  *intent(in)*    *scalar*  
        precalculated bond order sum for the given atom, :math:`\sum_j c_{ij}`, in the above example
    raw_gradients: double precision  *intent(in)*    *size(:, :)*  
        precalculated gradients of bond order sums, :math:`\nabla_\alpha \sum_j c_{ij}`, in the above example
    **total_bond_gradients**: double precision  **intent(out)**    *size(:, :)*  
        the calculated bond order gradients :math:`\nabla_\alpha b_i`
    raw_virial: double precision  *intent(in)*    *size(6)*  
        the precalculated virial due to the bond order gradient
    **total_virial**: double precision  **intent(out)**    *size(6)*  
        the scaled  virial due to the bond order gradient
    mpi_split: logical  *intent(in)*    *scalar*  *optional*
        A switch for enabling MPI parallelization. By default the routine is sequential since the calculation may be called from within an already parallelized routine.
            
  .. function:: core_release_all_memory()

    Release all allocated pointer arrays in the core.

            
  .. function:: core_set_ewald_parameters(real_cut, k_radius, reciprocal_cut, sigma, epsilon, scaler)

    Sets the parameters for Ewald summation in the core.
    

    Parameters:

    real_cut: double precision  *intent(in)*    *scalar*  
        the real-space cutoff
    k_radius: double precision  *intent(in)*    *scalar*  
        the k-space cutoff (in inverse length)
    reciprocal_cut: integer  *intent(in)*    *size(3)*  
        the k-space cutoffs (in numbers of k-space cells)
    sigma: double precision  *intent(in)*    *scalar*  
        the split parameter
    epsilon: double precision  *intent(in)*    *scalar*  
        electric constant
    scaler: double precision  *intent(in)*    *size(:)*  
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
            
  .. function:: expand_neighbor_storage(nbors_and_offsets, length, new_length, n_atoms)

    Expands the allocated memory for storing neighbor lists

    Parameters:

    nbors_and_offsets: integer  *intent()*  *pointer*  *size(:, :, :)*  
        
    length: integer  *intent(in)*    *scalar*  
        
    new_length: integer  *intent(in)*    *scalar*  
        
    n_atoms: integer  *intent(in)*    *scalar*  
        
            
  .. function:: list_atoms()

    Prints some information on the atoms stored in the core in stdout.

            
  .. function:: list_bonds()

    Prints some information on the bond order factors stored in the core in stdout.

            
  .. function:: list_cell()

    Prints some information on the supercell stored in the core in stdout.

            
  .. function:: list_interactions()

    Prints some information on the potentials stored in the core in stdout.
