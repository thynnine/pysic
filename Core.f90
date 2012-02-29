module pysic_core
  use geometry
  use potentials
  use mpi
  implicit none
  
  type(atom), pointer :: atoms(:)
  type(supercell) :: cell
  type(potential), pointer :: interactions(:)
  type(bond_order_parameters), pointer :: bond_factors(:)
  integer :: n_interactions = 0, n_bond_factors = 0
  logical :: &
       atoms_created = .false., &
       potentials_allocated = .false., &
       bond_factors_allocated = .false., &
       bond_storage_allocated = .false.

  ! storage for bond order factors during force and energy evaluation
  double precision, pointer :: saved_bond_order_sums(:,:)
  double precision, pointer :: saved_bond_order_factors(:,:)
  double precision, pointer :: saved_bond_order_gradients(:,:,:,:)
  integer :: n_saved_bond_order_factors = 0
  integer, pointer :: group_index_save_slot(:)
  logical :: use_saved_bond_order_factors = .false.
  integer, pointer :: use_saved_bond_order_gradients(:,:)

contains

  subroutine core_release_all_memory()
    implicit none

    call core_clear_atoms()
    call core_clear_potentials()
    call core_clear_bond_order_factors()
    call core_clear_bond_order_storage()

  end subroutine core_release_all_memory



  ! Creates the atomic particles by invoking a subroutine in the geometry module
  subroutine core_generate_atoms(n_atoms,masses,charges,positions,momenta,tags,elements)
    implicit none
    integer, intent(in) :: n_atoms, tags(n_atoms)
    double precision, intent(in) :: masses(n_atoms), charges(n_atoms), positions(3,n_atoms), &
         momenta(3,n_atoms)
    character(len=label_length), intent(in) :: elements(n_atoms)

    call generate_atoms(n_atoms,masses,charges,positions,momenta,tags,elements,atoms)
    atoms_created = .true.

  end subroutine core_generate_atoms



  subroutine core_update_atom_coordinates(n_atoms,positions,momenta)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(in) :: positions(3,n_atoms), momenta(3,n_atoms)

    call update_atomic_positions(n_atoms,positions,momenta,atoms)

  end subroutine core_update_atom_coordinates
 


  subroutine core_get_number_of_atoms(n_atoms)
    implicit none
    integer, intent(out) :: n_atoms
    
    if(atoms_created)then
       n_atoms = size(atoms)
    else
       n_atoms = 0
    end if

  end subroutine core_get_number_of_atoms



  subroutine core_clear_atoms()
    implicit none
    integer :: i

    if(atoms_created)then
       do i = 1, size(atoms)

          if(atoms(i)%neighbor_list%max_length > 0)then
             deallocate(atoms(i)%neighbor_list%neighbors)
             deallocate(atoms(i)%neighbor_list%pbc_offsets)
          end if
          if(atoms(i)%potentials_listed)then
             deallocate(atoms(i)%potential_indices)
          end if
          if(atoms(i)%bond_order_factors_listed)then
             deallocate(atoms(i)%bond_indices)
          end if

       end do
       deallocate(atoms)
       atoms_created = .false.
    end if

  end subroutine core_clear_atoms


  ! Creates a supercell for containing the calculation geometry
  subroutine core_create_cell(vectors,inverse,periodicity)
    implicit none
    double precision, intent(in) :: vectors(3,3), inverse(3,3)
    logical, intent(in) :: periodicity(3)

    call generate_supercell(vectors,inverse,periodicity,cell)

  end subroutine core_create_cell


  subroutine core_get_cell_vectors(vectors)
    implicit none
    double precision, intent(out) :: vectors(3,3)
    
    vectors = cell%vectors

  end subroutine core_get_cell_vectors


  ! Assigns a precalculated neighbor list to a single atom of the given index
  subroutine core_create_neighbor_list(n_nbs,atom_index,neighbors,offsets)
    implicit none
    integer, intent(in) :: n_nbs, atom_index
    integer, intent(in) :: neighbors(n_nbs), offsets(3,n_nbs)

    call assign_neighbor_list(n_nbs,atoms(atom_index)%neighbor_list,neighbors,offsets)

  end subroutine core_create_neighbor_list


  subroutine core_clear_potentials()
    implicit none
    integer :: i

    if(n_interactions > 0)then
       do i = 1, n_interactions
          deallocate(interactions(i)%parameters)
          deallocate(interactions(i)%apply_elements)
          deallocate(interactions(i)%apply_tags)
          deallocate(interactions(i)%apply_indices)
          deallocate(interactions(i)%original_elements)
          deallocate(interactions(i)%original_tags)
          deallocate(interactions(i)%original_indices)
          deallocate(interactions(i)%derived_parameters)  
       end do
    end if
    if(potentials_allocated)then
       deallocate(interactions)
    end if
    n_interactions = 0
    potentials_allocated = .false.

  end subroutine core_clear_potentials
  

  
  subroutine core_clear_bond_order_factors()
    implicit none
    integer :: i

    if(n_bond_factors > 0)then
       do i = 1, n_bond_factors
          deallocate(bond_factors(i)%parameters)
          deallocate(bond_factors(i)%n_params)
          deallocate(bond_factors(i)%apply_elements)
          deallocate(bond_factors(i)%original_elements)
          deallocate(bond_factors(i)%derived_parameters)  
       end do
    end if
    if(bond_factors_allocated)then
       deallocate(bond_factors)
    end if
    n_bond_factors = 0
    bond_factors_allocated = .false.

  end subroutine core_clear_bond_order_factors


  subroutine core_clear_bond_order_storage()
    implicit none

    if(bond_storage_allocated)then
       deallocate(saved_bond_order_sums)
       deallocate(saved_bond_order_factors)
       deallocate(saved_bond_order_gradients)
       deallocate(use_saved_bond_order_gradients)
       deallocate(group_index_save_slot)
    else
       nullify(saved_bond_order_sums)
       nullify(saved_bond_order_factors)
       nullify(saved_bond_order_gradients)
       nullify(use_saved_bond_order_gradients)
       nullify(group_index_save_slot)
    end if
    n_saved_bond_order_factors = 0
    use_saved_bond_order_factors = .false.
    bond_storage_allocated = .false.

  end subroutine core_clear_bond_order_storage


  subroutine core_empty_bond_order_storage()
    implicit none
    
    if(bond_storage_allocated)then
       saved_bond_order_sums = 0.d0
       saved_bond_order_factors = 0.d0
       call core_empty_bond_order_gradient_storage()
       group_index_save_slot = -1
       n_saved_bond_order_factors = 0
    end if
    
  end subroutine core_empty_bond_order_storage

  subroutine core_empty_bond_order_gradient_storage(index)
    implicit none
    integer, optional, intent(in) :: index

    if(present(index))then
       if(bond_storage_allocated)then
          use_saved_bond_order_gradients(:,index) = -1
       end if
    else
       if(bond_storage_allocated)then
          saved_bond_order_gradients = 0.d0
          use_saved_bond_order_gradients = -1
       end if
    end if

  end subroutine core_empty_bond_order_gradient_storage

  subroutine core_allocate_bond_order_storage(n_atoms,n_groups,n_factors)
    implicit none
    integer, intent(in) :: n_atoms, n_groups, n_factors

    call core_clear_bond_order_storage()
    allocate(saved_bond_order_sums(n_atoms,n_factors))
    allocate(saved_bond_order_factors(n_atoms,n_factors))
    allocate(saved_bond_order_gradients(3,n_atoms,n_factors,3))
    allocate(use_saved_bond_order_gradients(n_factors,3))
    allocate(group_index_save_slot(0:n_groups))
    bond_storage_allocated = .true.
    call core_empty_bond_order_storage()


  end subroutine core_allocate_bond_order_storage
  

  subroutine core_fill_bond_order_storage(n_atoms)
    implicit none
    integer, intent(in) :: n_atoms
    integer :: i, group_index
    double precision :: dummy_factors(n_atoms)

    do i = 1, n_interactions
       group_index = interactions(i)%pot_index
       if(group_index > -1)then
          call core_get_bond_order_sums(n_atoms,group_index,dummy_factors)
       end if
    end do

  end subroutine core_fill_bond_order_storage


  subroutine core_get_bond_order_gradients(n_atoms,group_index,atom_index,slot_index,bond_order_gradients)
    implicit none
    integer, intent(in) :: n_atoms, group_index, atom_index, slot_index
    double precision, intent(out) :: bond_order_gradients(1:3,n_atoms)
    double precision :: bond_order_sums(n_atoms)
    logical :: found_grads
    integer :: save_slot

    if(use_saved_bond_order_factors)then
       found_grads = .false.
       save_slot = group_index_save_slot(group_index)       
       if(save_slot > 0)then
          if(use_saved_bond_order_gradients(save_slot,slot_index) == atom_index)then
             !write(*,*) "access saved grads: ", save_slot, slot_index, atom_index
             found_grads = .true.
             bond_order_gradients(1:3,1:n_atoms) = saved_bond_order_gradients(1:3,1:n_atoms,save_slot,slot_index)
          end if
       end if
       if(.not.found_grads)then
          !write(*,*) "calculate grads: ", save_slot, slot_index, atom_index
          call core_get_bond_order_sums(n_atoms,group_index,bond_order_sums)
          call core_calculate_bond_order_gradients_of_factor(n_atoms,&
               group_index,&
               atom_index,&
               bond_order_sums,&
               bond_order_gradients)
          saved_bond_order_gradients(1:3,1:n_atoms,save_slot,slot_index) = bond_order_gradients(1:3,1:n_atoms)
          use_saved_bond_order_gradients(save_slot,slot_index) = atom_index
       end if
    else
       call core_get_bond_order_sums(n_atoms,group_index,bond_order_sums)
       call core_calculate_bond_order_gradients_of_factor(n_atoms,&
            group_index,&
            atom_index,&
            bond_order_sums,&
            bond_order_gradients)
    end if

  end subroutine core_get_bond_order_gradients


  subroutine core_get_bond_order_sums(n_atoms,group_index,bond_order_sums)
    implicit none
    integer, intent(in) :: n_atoms, group_index
    double precision, intent(out) :: bond_order_sums(n_atoms)
    double precision :: bond_order_factors(n_atoms)
    integer :: save_slot

    if(use_saved_bond_order_factors)then
       save_slot = group_index_save_slot(group_index)
       if(save_slot > 0)then
          bond_order_sums(1:n_atoms) = saved_bond_order_sums(1:n_atoms,save_slot)
       else
          call core_calculate_bond_order_factors(n_atoms,group_index,bond_order_sums)
          call core_post_process_bond_order_factors(n_atoms,group_index,bond_order_sums,bond_order_factors)
          n_saved_bond_order_factors = n_saved_bond_order_factors + 1
          group_index_save_slot(group_index) = n_saved_bond_order_factors
          saved_bond_order_sums(1:n_atoms,n_saved_bond_order_factors) = bond_order_sums(1:n_atoms)
          saved_bond_order_factors(1:n_atoms,n_saved_bond_order_factors) = bond_order_factors(1:n_atoms)
       end if
    else
       call core_calculate_bond_order_factors(n_atoms,group_index,bond_order_sums)       
    end if

  end subroutine core_get_bond_order_sums



  subroutine core_get_bond_order_factors(n_atoms,group_index,bond_order_factors)
    implicit none
    integer, intent(in) :: n_atoms, group_index
    double precision, intent(out) :: bond_order_factors(n_atoms)
    double precision :: bond_order_sums(n_atoms)
    integer :: save_slot

    if(use_saved_bond_order_factors)then
       save_slot = group_index_save_slot(group_index)
       if(save_slot > 0)then
          bond_order_factors(1:n_atoms) = saved_bond_order_factors(1:n_atoms,save_slot)
       else
          call core_calculate_bond_order_factors(n_atoms,group_index,bond_order_sums)
          call core_post_process_bond_order_factors(n_atoms,group_index,bond_order_sums,bond_order_factors)
          n_saved_bond_order_factors = n_saved_bond_order_factors + 1
          group_index_save_slot(group_index) = n_saved_bond_order_factors
          saved_bond_order_sums(1:n_atoms,n_saved_bond_order_factors) = bond_order_sums(1:n_atoms)
          saved_bond_order_factors(1:n_atoms,n_saved_bond_order_factors) = bond_order_factors(1:n_atoms)
       end if
    else
       call core_calculate_bond_order_factors(n_atoms,group_index,bond_order_sums)       
       call core_post_process_bond_order_factors(n_atoms,group_index,bond_order_sums,bond_order_factors)
    end if

  end subroutine core_get_bond_order_factors


  subroutine core_get_bond_order_factor_of_atom(n_atoms,group_index,atom_index,bond_order_factor)
    implicit none
    integer, intent(in) :: n_atoms, group_index, atom_index
    double precision, intent(out) :: bond_order_factor
    double precision :: bond_order_factors(n_atoms)

    call core_get_bond_order_factors(n_atoms,group_index,bond_order_factors)
    bond_order_factor = bond_order_factors(atom_index)

  end subroutine core_get_bond_order_factor_of_atom


  subroutine core_allocate_potentials(n_pots)
    implicit none
    integer, intent(in) :: n_pots
    integer :: i

    call core_clear_potentials()
    allocate(interactions(n_pots))
    do i = 1, n_pots
       nullify(interactions(i)%parameters)
       nullify(interactions(i)%apply_elements)
       nullify(interactions(i)%apply_tags)
       nullify(interactions(i)%apply_indices)
       nullify(interactions(i)%original_elements)
       nullify(interactions(i)%original_tags)
       nullify(interactions(i)%original_indices)
       nullify(interactions(i)%derived_parameters) 
    end do
    potentials_allocated = .true.

  end subroutine core_allocate_potentials


  subroutine core_allocate_bond_order_factors(n_bond_factors)
    implicit none
    integer, intent(in) :: n_bond_factors
    integer :: i

    call core_clear_bond_order_factors()
    allocate(bond_factors(n_bond_factors))
    do i = 1, n_bond_factors
       nullify(bond_factors(i)%parameters)
       nullify(bond_factors(i)%n_params)
       nullify(bond_factors(i)%apply_elements)
       nullify(bond_factors(i)%original_elements)
       nullify(bond_factors(i)%derived_parameters)
    end do
    bond_factors_allocated = .true.

  end subroutine core_allocate_bond_order_factors



  subroutine core_add_potential(n_targets,n_params,pot_name,parameters,cutoff,smooth_cut,&
       elements,tags,indices,orig_elements,orig_tags,orig_indices,pot_index)
    implicit none
    integer, intent(in) :: n_targets, n_params, pot_index
    character(len=*), intent(in) :: pot_name
    double precision, intent(in) :: parameters(n_params)
    double precision, intent(in) :: cutoff, smooth_cut
    character(len=label_length), intent(in) :: elements(n_targets)
    integer, intent(in) :: tags(n_targets), indices(n_targets)
    character(len=label_length), intent(in) :: orig_elements(n_targets)
    integer, intent(in) :: orig_tags(n_targets), orig_indices(n_targets)
    type(potential) :: new_interaction

    n_interactions = n_interactions + 1
    call create_potential(n_targets,n_params,&
         pot_name,parameters,cutoff,smooth_cut,&
         elements,tags,indices,&
         orig_elements,orig_tags,orig_indices,pot_index,&
         new_interaction)
    interactions(n_interactions) = new_interaction

  end subroutine core_add_potential



  subroutine core_add_bond_order_factor(n_targets,n_params,n_split,bond_name,parameters,param_split,&
       cutoff,smooth_cut,elements,orig_elements,group_index)
    implicit none
    integer, intent(in) :: n_targets, n_params, n_split, group_index
    integer, intent(in) :: param_split(n_split)
    character(len=*), intent(in) :: bond_name
    double precision, intent(in) :: parameters(n_params)
    double precision, intent(in) :: cutoff, smooth_cut
    character(len=label_length), intent(in) :: elements(n_targets)
    character(len=label_length), intent(in) :: orig_elements(n_targets)
    type(bond_order_parameters) :: new_bond_factor

    n_bond_factors = n_bond_factors + 1
    call create_bond_order_factor(n_targets,n_params,n_split,&
         bond_name,parameters,param_split,cutoff,smooth_cut,&
         elements,orig_elements,group_index,&
         new_bond_factor)
    bond_factors(n_bond_factors) = new_bond_factor

  end subroutine core_add_bond_order_factor


  subroutine core_assign_bond_order_factor_indices()
    implicit none
    logical :: affects(n_bond_factors)
    integer :: i, j, k, total
    integer, allocatable :: bond_indices(:)

    do i = 1, size(atoms)       
       total = 0
       do j = 1, n_bond_factors
          call bond_order_factor_affects_atom(bond_factors(j),atoms(i),affects(j),1)
          if(affects(j))then
             total = total+1
          end if
       end do

       allocate(bond_indices(total))
       
       k = 0
       do j = 1, n_bond_factors
          if(affects(j))then
             k = k+1
             bond_indices(k) = j
          end if
       end do

       ! pointer allocations are handled in the subroutine
       call assign_bond_order_factor_indices(total,atoms(i),bond_indices)

       deallocate(bond_indices)
    end do

  end subroutine core_assign_bond_order_factor_indices


  subroutine core_assign_potential_indices()
    implicit none
    logical :: affects(n_interactions)
    integer :: i, j, k, total
    integer, allocatable :: pot_indices(:)
    
    do i = 1, size(atoms)       
       total = 0
       do j = 1, n_interactions
          call potential_affects_atom(interactions(j),atoms(i),affects(j),1)
          if(affects(j))then
             total = total+1
          end if
       end do

       allocate(pot_indices(total))
       
       k = 0
       do j = 1, n_interactions
          if(affects(j))then
             k = k+1
             pot_indices(k) = j
          end if
       end do

       ! pointer allocations are handled in the subroutine
       call assign_potential_indices(total,atoms(i),pot_indices)

       deallocate(pot_indices)
    end do

  end subroutine core_assign_potential_indices


  ! calculates the number of neighbors for each atom
  subroutine core_calculate_coordinations(n_atoms,cutoffs,total_coordinations)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(in) :: cutoffs(2)
    double precision, intent(out) :: total_coordinations(n_atoms)
    double precision :: coordinations(n_atoms), distance, separation(3), proximity_factor
    integer :: index1, index2, j
    type(atom) :: atom1, atom2
    type(neighbor_list) :: nbors1

    coordinations = 0.d0
    total_coordinations = 0.d0

    ! loop over atoms
    do index1 = 1, size(atoms)
       
       ! in MPI, only consider the atoms of this cpu
       if(is_my_atom(index1))then

          atom1 = atoms(index1)
          nbors1 = atom1%neighbor_list

          ! loop over neighbors
          do j = 1, nbors1%n_neighbors
             
             index2 = nbors1%neighbors(j)
             if(index2 > index1)then ! do not double count

                atom2 = atoms(index2)
                call separation_vector(atom1%position, &
                     atom2%position, &
                     nbors1%pbc_offsets(1:3,j), &
                     cell, &
                     separation)
                distance = .norm.(separation)

                call smoothening_factor(distance,&
                     cutoffs(2),cutoffs(1),&
                     proximity_factor)

                coordinations(index1) = coordinations(index1) + proximity_factor
                coordinations(index2) = coordinations(index2) + proximity_factor

             end if ! index2 > index1

          end do ! j = 1, nbors1%n_neighbors

       end if ! is_my_atom(index1)

    end do ! index1 = 1, size(atoms)


#ifdef MPI
    call mpi_allreduce(coordinations,total_coordinations,size(coordinations),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
#else
    total_coordinations = coordinations
#endif

  end subroutine core_calculate_coordinations


  ! Returns the gradients of all bond order factors with respect to
  ! moving a given atom.
  subroutine core_calculate_bond_order_gradients(n_atoms,group_index,&
       atom_index,raw_sums,total_gradient)
    implicit none
    integer, intent(in) :: n_atoms, group_index, atom_index
    double precision, intent(out) :: total_gradient(3,n_atoms)
    double precision, intent(in) :: raw_sums(n_atoms)
    double precision :: gradient(3,n_atoms)
    type(atom) :: atom1, atom2, atom3
    type(atom) :: atom_list(3)
    type(neighbor_list) :: nbors1, nbors2
    type(bond_order_parameters) :: bond_params(2)
    integer, pointer :: bond_indices(:), bond_indices2(:)
    integer :: index1, index2, index3, k1, k2, j, l, n_targets
    double precision :: separations(3,2), distances(2), directions(3,2), tmp_grad(3,3,3)
    logical :: is_active, is_in_group, many_bodies_found, separation3_unknown

    gradient = 0.d0
    total_gradient = 0.d0
    
    ! target atom
    index1 = atom_index
    atom1 = atoms(atom_index)
    nbors1 = atom1%neighbor_list
    bond_indices => atom1%bond_indices
          
    ! loop over neighbors
    do j = 1, nbors1%n_neighbors
       
       ! neighboring atom
       index2 = nbors1%neighbors(j)
          
       atom2 = atoms(index2)
       call separation_vector(atom1%position, &
            atom2%position, &
            nbors1%pbc_offsets(1:3,j), &
            cell, &
            separations(1:3,1))
       distances(1) = .norm.(separations(1:3,1))
       if(distances(1) == 0.d0)then
          directions(1:3,1) = (/ 0.d0, 0.d0, 0.d0 /)
       else
          directions(1:3,1) = separations(1:3,1) / distances(1)
       end if
          
       many_bodies_found = .false.
       ! 2-body bond order factors
       do k1 = 1, size(bond_indices)
                   
          bond_params(1) = bond_factors(bond_indices(k1))
          call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,n_targets)
          call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2)
          call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group)
                   
          if( is_active .and. is_in_group )then !.and. bond_params(1)%cutoff > distances(1) )then
             if( n_targets == 2 )then

                call evaluate_bond_order_gradient(2,&
                     separations(1:3,1),&
                     distances(1),&
                     bond_params(1),&
                     tmp_grad(1:3,1:2,1:2))
                ! store the gradients of the atom1 and atom2 terms
                ! with respect to movind atom1 (the target atom)
                gradient(1:3,index1) = gradient(1:3,index1) + tmp_grad(1:3,1,1)
                gradient(1:3,index2) = gradient(1:3,index2) + tmp_grad(1:3,2,1)

             else if( n_targets > 2 )then
                   
                many_bodies_found = .true.
                   
             end if
          end if
             
       end do ! k1 = 1, size(bond_indices)


       if(many_bodies_found)then

          nbors2 = atom2%neighbor_list

          ! loop over neighbors of atom 1
          do l = 1, nbors1%n_neighbors
             index3 = nbors1%neighbors(l)

             ! Since we first loop over the neighbors of atom1 to get atom2 candidates
             ! and then again to find atom3 candidates, we will find the same triplet
             ! atom2-atom1-atom3 = atom3-atom1-atom2 twice.
             ! In order to filter out the duplicate, only consider the case where
             ! index of atom3 is bigger than the index of atom2. This condition is
             ! bound to be true for one and false for the other triplet.
             if(index3 > index2)then

                atom3 = atoms(index3)
                separation3_unknown = .true.
                atom_list = (/ atom2, atom1, atom3 /)
                
                ! search for the first bond params containing the parameters for atom1-atom2
                do k1 = 1, size(bond_indices)
                   
                   bond_params(1) = bond_factors(bond_indices(k1))
                   call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                        n_targets)
                   call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2)
                   call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group)
                   
                   if( is_active .and. is_in_group .and. n_targets == 3 )then
                      call bond_order_factor_affects_atom(bond_params(1),atom3,is_active,3)
                      
                      if( is_active )then
                         
                         
                         ! search for the second bond params containing the parameters for atom1-atom3
                         do k2 = 1, size(bond_indices)
                            
                            bond_params(2) = bond_factors(bond_indices(k2))
                            call get_number_of_targets_of_bond_order_factor_index(bond_params(2)%type_index,&
                                 n_targets)
                            call bond_order_factor_affects_atom(bond_params(2),atom2,is_active,2)
                            call bond_order_factor_is_in_group(bond_params(2),group_index,is_in_group)
                            
                            if( is_active .and. is_in_group .and. n_targets == 3 )then
                               call bond_order_factor_affects_atom(bond_params(2),atom3,is_active,3)
                               
                               if( is_active )then
                                  
                                  if( separation3_unknown )then
                                     call separation_vector(atom1%position, &
                                          atom3%position, &
                                          nbors1%pbc_offsets(1:3,j), &
                                          cell, &
                                          separations(1:3,2))
                                     separation3_unknown = .false.
                                     distances(2) = .norm.(separations(1:3,2))
                                     if(distances(2) == 0.d0)then
                                        directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                     else
                                        directions(1:3,2) = separations(1:3,2) / distances(2)
                                     end if
                                  end if
                                  
                                  call evaluate_bond_order_gradient(3,&
                                       separations(1:3,1:2),&
                                       distances(1:2),&
                                       bond_params(1:2),&
                                       tmp_grad(1:3,1:3,1:3),&
                                       atom_list)
                                  ! store the gradients of the atom1, atom2, and atom3 terms
                                  ! with respect to movind atom1 (the target atom)
                                  ! Note that here atom1 is the middle atom, so we take index 2
                                  ! in the third column of tmp_grad.
                                  gradient(1:3,index2) = gradient(1:3,index2) + tmp_grad(1:3,1,2)
                                  gradient(1:3,index1) = gradient(1:3,index1) + tmp_grad(1:3,2,2)
                                  gradient(1:3,index3) = gradient(1:3,index3) + tmp_grad(1:3,3,2)
                                  
                               end if ! is_active
                            end if ! is_active and is_in_group and n_targets == 3
                            
                         end do ! k2
                         
                      end if ! is_active
                   end if ! is_active and is_in_group and n_targets == 3
                   
                end do ! k1
             end if ! index3 /= index2
          end do ! l = 1, nbors1%n_neighbors


          ! Next we try to find ordered triplets atom1 -- atom2 -- atom3
          ! Therefore we need separations a2--a1 and a2--a3.
          separations(1:3,1) = -separations(1:3,1)

          ! loop over neighbors of atom 2
          do l = 1, nbors2%n_neighbors
             index3 = nbors2%neighbors(l)

             if(index3 /= index1)then
                atom3 = atoms(index3)
                separation3_unknown = .true.
                atom_list = (/ atom1, atom2, atom3 /)
                
                ! search for the first bond params containing the parameters for atom2-atom1
                do k1 = 1, size(bond_indices)
                   
                   bond_params(1) = bond_factors(bond_indices(k1))
                   call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,n_targets)
                   call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2)
                   call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group)
                   
                   if( is_active .and. is_in_group .and. n_targets == 3 )then
                      call bond_order_factor_affects_atom(bond_params(1),atom3,is_active,3)
                      
                      if( is_active )then
                         
                         ! search for the second bond params containing the parameters for atom2-atom3
                         do k2 = 1, size(bond_indices)
                            
                            bond_params(2) = bond_factors(bond_indices(k2))
                            call get_number_of_targets_of_bond_order_factor_index(bond_params(2)%type_index,&
                                 n_targets)
                            call bond_order_factor_affects_atom(bond_params(2),atom2,is_active,2)
                            call bond_order_factor_is_in_group(bond_params(2),group_index,is_in_group)
                            
                            if( is_active .and. is_in_group .and. n_targets == 3 )then
                               call bond_order_factor_affects_atom(bond_params(2),atom3,is_active,3)
                               
                               if( is_active )then
                                  
                                  if( separation3_unknown )then
                                     call separation_vector(atom2%position, &
                                          atom3%position, &
                                          nbors2%pbc_offsets(1:3,l), &
                                          cell, &
                                          separations(1:3,2))
                                     separation3_unknown = .false.
                                     distances(2) = .norm.(separations(1:3,2))
                                     if(distances(2) == 0.d0)then
                                        directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                     else
                                        directions(1:3,2) = separations(1:3,2) / distances(2)
                                     end if
                                  end if
                                  
                                  call evaluate_bond_order_gradient(3,&
                                       separations(1:3,1:2),&
                                       distances(1:2),&
                                       bond_params(1:2),&
                                       tmp_grad(1:3,1:3,1:3),&
                                       atom_list)
                                  ! store the gradients of the atom1, atom2, and atom3 terms
                                  ! with respect to movind atom1 (the target atom)
                                  ! Note that here atom1 is the first atom, so we take index 1
                                  ! in the third column of tmp_grad.
                                  gradient(1:3,index1) = gradient(1:3,index1) + tmp_grad(1:3,1,1)
                                  gradient(1:3,index2) = gradient(1:3,index2) + tmp_grad(1:3,2,1)
                                  gradient(1:3,index3) = gradient(1:3,index3) + tmp_grad(1:3,3,1)
                                  
                               end if ! is_active
                            end if ! is_active and is_in_group and n_targets == 3
                            
                         end do ! k2
                         
                      end if ! is_active
                   end if ! is_active and is_in_group and n_targets == 3
                   
                end do ! k1
                
             end if ! index3 /= index1
          end do ! l = 1, nbors2%n_neighbors
             
       end if ! many_bodies_found
    end do ! j = nbors1%n_neighbors

    call core_post_process_bond_order_gradients(n_atoms,group_index,raw_sums,&
         gradient,total_gradient)

  end subroutine core_calculate_bond_order_gradients


  ! Returns the gradients of one bond order factor with respect to
  ! moving all atoms.
  subroutine core_calculate_bond_order_gradients_of_factor(n_atoms,group_index,&
       atom_index,raw_sums,total_gradient)
    implicit none
    integer, intent(in) :: n_atoms, group_index, atom_index
    double precision, intent(out) :: total_gradient(3,n_atoms)
    double precision, intent(in) :: raw_sums(n_atoms)
    double precision :: gradient(3,n_atoms)
    type(atom) :: atom1, atom2, atom3
    type(atom) :: atom_list(3)
    type(neighbor_list) :: nbors1, nbors2
    type(bond_order_parameters) :: bond_params(2)
    integer, pointer :: bond_indices(:), bond_indices2(:)
    integer :: index1, index2, index3, k1, k2, j, l, n_targets
    double precision :: separations(3,2), distances(2), directions(3,2), tmp_grad(3,3,3)
    logical :: is_active, is_in_group, many_bodies_found, separation3_unknown

    gradient = 0.d0
    total_gradient = 0.d0
    
    ! target atom
    index1 = atom_index
    atom1 = atoms(atom_index)
    nbors1 = atom1%neighbor_list
    bond_indices => atom1%bond_indices
          
    ! loop over neighbors
    do j = 1, nbors1%n_neighbors
       
       ! neighboring atom
       index2 = nbors1%neighbors(j)
          
       atom2 = atoms(index2)
       call separation_vector(atom1%position, &
            atom2%position, &
            nbors1%pbc_offsets(1:3,j), &
            cell, &
            separations(1:3,1))
       distances(1) = .norm.(separations(1:3,1))
       if(distances(1) == 0.d0)then
          directions(1:3,1) = (/ 0.d0, 0.d0, 0.d0 /)
       else
          directions(1:3,1) = separations(1:3,1) / distances(1)
       end if
          
       many_bodies_found = .false.
       ! 2-body bond order factors
       do k1 = 1, size(bond_indices)
                   
          bond_params(1) = bond_factors(bond_indices(k1))
          call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,n_targets)
          call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2)
          call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group)
                   
          if( is_active .and. is_in_group )then !.and. bond_params(1)%cutoff > distances(1) )then
             if( n_targets == 2 )then

                call evaluate_bond_order_gradient(2,&
                     separations(1:3,1),&
                     distances(1),&
                     bond_params(1),&
                     tmp_grad(1:3,1:2,1:2))
                ! store the gradients of the atom1 term
                ! with respect to moving atom1 and atom2
                gradient(1:3,index1) = gradient(1:3,index1) + tmp_grad(1:3,1,1)
                gradient(1:3,index2) = gradient(1:3,index2) + tmp_grad(1:3,1,2)

             else if( n_targets > 2 )then
                   
                many_bodies_found = .true.
                   
             end if
          end if
             
       end do ! k1 = 1, size(bond_indices)


       if(many_bodies_found)then

          nbors2 = atom2%neighbor_list

          ! loop over neighbors of atom 1
          do l = 1, nbors1%n_neighbors
             index3 = nbors1%neighbors(l)

             ! Since we first loop over the neighbors of atom1 to get atom2 candidates
             ! and then again to find atom3 candidates, we will find the same triplet
             ! atom2-atom1-atom3 = atom3-atom1-atom2 twice.
             ! In order to filter out the duplicate, only consider the case where
             ! index of atom3 is bigger than the index of atom2. This condition is
             ! bound to be true for one and false for the other triplet.
             if(index3 > index2)then

                atom3 = atoms(index3)
                separation3_unknown = .true.
                atom_list = (/ atom2, atom1, atom3 /)
                
                ! search for the first bond params containing the parameters for atom1-atom2
                do k1 = 1, size(bond_indices)
                   
                   bond_params(1) = bond_factors(bond_indices(k1))
                   call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                        n_targets)
                   call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2)
                   call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group)
                   
                   if( is_active .and. is_in_group .and. n_targets == 3 )then
                      call bond_order_factor_affects_atom(bond_params(1),atom3,is_active,3)
                      
                      if( is_active )then
                         
                         
                         ! search for the second bond params containing the parameters for atom1-atom3
                         do k2 = 1, size(bond_indices)
                            
                            bond_params(2) = bond_factors(bond_indices(k2))
                            call get_number_of_targets_of_bond_order_factor_index(bond_params(2)%type_index,&
                                 n_targets)
                            call bond_order_factor_affects_atom(bond_params(2),atom2,is_active,2)
                            call bond_order_factor_is_in_group(bond_params(2),group_index,is_in_group)
                            
                            if( is_active .and. is_in_group .and. n_targets == 3 )then
                               call bond_order_factor_affects_atom(bond_params(2),atom3,is_active,3)
                               
                               if( is_active )then
                                  
                                  if( separation3_unknown )then
                                     call separation_vector(atom1%position, &
                                          atom3%position, &
                                          nbors1%pbc_offsets(1:3,j), &
                                          cell, &
                                          separations(1:3,2))
                                     separation3_unknown = .false.
                                     distances(2) = .norm.(separations(1:3,2))
                                     if(distances(2) == 0.d0)then
                                        directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                     else
                                        directions(1:3,2) = separations(1:3,2) / distances(2)
                                     end if
                                  end if
                                  
                                  call evaluate_bond_order_gradient(3,&
                                       separations(1:3,1:2),&
                                       distances(1:2),&
                                       bond_params(1:2),&
                                       tmp_grad(1:3,1:3,1:3),&
                                       atom_list)
                                  ! store the gradients of the atom1 terms
                                  ! with respect to moving atom1, atom2, atom3 
                                  ! Note that here atom1 is the middle atom
                                  gradient(1:3,index2) = gradient(1:3,index2) + tmp_grad(1:3,2,1)
                                  gradient(1:3,index1) = gradient(1:3,index1) + tmp_grad(1:3,2,2)
                                  gradient(1:3,index3) = gradient(1:3,index3) + tmp_grad(1:3,2,3)
                                  
                               end if ! is_active
                            end if ! is_active and is_in_group and n_targets == 3
                            
                         end do ! k2
                         
                      end if ! is_active
                   end if ! is_active and is_in_group and n_targets == 3
                   
                end do ! k1
             end if ! index3 /= index2
          end do ! l = 1, nbors1%n_neighbors


          ! Next we try to find ordered triplets atom1 -- atom2 -- atom3
          ! Therefore we need separations a2--a1 and a2--a3.
          separations(1:3,1) = -separations(1:3,1)

          ! loop over neighbors of atom 2
          do l = 1, nbors2%n_neighbors
             index3 = nbors2%neighbors(l)

             if(index3 /= index1)then
                atom3 = atoms(index3)
                separation3_unknown = .true.
                atom_list = (/ atom1, atom2, atom3 /)
                
                ! search for the first bond params containing the parameters for atom2-atom1
                do k1 = 1, size(bond_indices)
                   
                   bond_params(1) = bond_factors(bond_indices(k1))
                   call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,n_targets)
                   call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2)
                   call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group)
                   
                   if( is_active .and. is_in_group .and. n_targets == 3 )then
                      call bond_order_factor_affects_atom(bond_params(1),atom3,is_active,3)
                      
                      if( is_active )then
                         
                         ! search for the second bond params containing the parameters for atom2-atom3
                         do k2 = 1, size(bond_indices)
                            
                            bond_params(2) = bond_factors(bond_indices(k2))
                            call get_number_of_targets_of_bond_order_factor_index(bond_params(2)%type_index,&
                                 n_targets)
                            call bond_order_factor_affects_atom(bond_params(2),atom2,is_active,2)
                            call bond_order_factor_is_in_group(bond_params(2),group_index,is_in_group)
                            
                            if( is_active .and. is_in_group .and. n_targets == 3 )then
                               call bond_order_factor_affects_atom(bond_params(2),atom3,is_active,3)
                               
                               if( is_active )then
                                  
                                  if( separation3_unknown )then
                                     call separation_vector(atom2%position, &
                                          atom3%position, &
                                          nbors2%pbc_offsets(1:3,l), &
                                          cell, &
                                          separations(1:3,2))
                                     separation3_unknown = .false.
                                     distances(2) = .norm.(separations(1:3,2))
                                     if(distances(2) == 0.d0)then
                                        directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                     else
                                        directions(1:3,2) = separations(1:3,2) / distances(2)
                                     end if
                                  end if
                                  
                                  call evaluate_bond_order_gradient(3,&
                                       separations(1:3,1:2),&
                                       distances(1:2),&
                                       bond_params(1:2),&
                                       tmp_grad(1:3,1:3,1:3),&
                                       atom_list)
                                  ! store the gradients of the atom1 term
                                  ! with respect to moving atom1, atom2, atom3
                                  ! Note that here atom1 is the first atom
                                  gradient(1:3,index1) = gradient(1:3,index1) + tmp_grad(1:3,1,1)
                                  gradient(1:3,index2) = gradient(1:3,index2) + tmp_grad(1:3,1,2)
                                  gradient(1:3,index3) = gradient(1:3,index3) + tmp_grad(1:3,1,3)
                                  
                               end if ! is_active
                            end if ! is_active and is_in_group and n_targets == 3
                            
                         end do ! k2
                         
                      end if ! is_active
                   end if ! is_active and is_in_group and n_targets == 3
                   
                end do ! k1
                
             end if ! index3 /= index1
          end do ! l = 1, nbors2%n_neighbors
             
       end if ! many_bodies_found
    end do ! j = nbors1%n_neighbors

    call core_post_process_bond_order_gradients_of_factor(n_atoms,group_index,raw_sums(index1),&
         gradient,total_gradient)

  end subroutine core_calculate_bond_order_gradients_of_factor




  subroutine core_calculate_bond_order_factors(n_atoms,group_index,total_bond_orders)
    implicit none
    integer, intent(in) :: n_atoms, group_index
    double precision, intent(out) :: total_bond_orders(n_atoms)
    integer :: index1, index2, index3, k1, k2, j, l, n_targets
    double precision :: separations(3,2), distances(2), directions(3,2)
    double precision :: tmp_factor(3), bond_orders(n_atoms)
    type(atom) :: atom1, atom2, atom3
    type(atom) :: atom_list(3)
    type(neighbor_list) :: nbors1, nbors2
    type(bond_order_parameters) :: bond_params(2)
    integer, pointer :: bond_indices(:), bond_indices2(:)
    logical :: is_active, is_in_group, many_bodies_found, separation3_unknown

    bond_orders = 0.d0
    total_bond_orders = 0.d0

    ! loop over atoms
    do index1 = 1, size(atoms)
       
       ! in MPI, only consider the atoms allocated to this particular cpu
       if(is_my_atom(index1))then

          ! target atom
          atom1 = atoms(index1)
          nbors1 = atom1%neighbor_list
          bond_indices => atom1%bond_indices

          ! loop over neighbors
          do j = 1, nbors1%n_neighbors

             ! neighboring atom
             index2 = nbors1%neighbors(j)
             if(index2 > index1)then ! do not double count

                atom2 = atoms(index2)
                call separation_vector(atom1%position, &
                     atom2%position, &
                     nbors1%pbc_offsets(1:3,j), &
                     cell, &
                     separations(1:3,1))
                distances(1) = .norm.(separations(1:3,1))
                if(distances(1) == 0.d0)then
                   directions(1:3,1) = (/ 0.d0, 0.d0, 0.d0 /)
                else
                   directions(1:3,1) = separations(1:3,1) / distances(1)
                end if
                
                many_bodies_found = .false.
                ! 2-body bond order factors
                do k1 = 1, size(bond_indices)
                   
                   bond_params(1) = bond_factors(bond_indices(k1))
                   call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,n_targets)
                   call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2)
                   call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group)
                   
                   if( is_active .and. is_in_group )then !.and. bond_params(1)%cutoff > distances(1) )then
                      if( n_targets == 2 )then
                         
                         call evaluate_bond_order_factor(2,separations(1:3,1),distances(1),bond_params(1),tmp_factor(1:2))

                         bond_orders(index1) = bond_orders(index1) + tmp_factor(1)
                         bond_orders(index2) = bond_orders(index2) + tmp_factor(2)
                         
                      else if( n_targets > 2 )then
                         
                         many_bodies_found = .true.
                         
                      end if
                   end if
                   
                end do ! k1 = 1, size(bond_indices)

                if(many_bodies_found)then
                   
                   nbors2 = atom2%neighbor_list
                   !bond_indices2 => atom2%bond_indices
                                      
                   ! loop over neighbors of atom 1
                   do l = 1, nbors1%n_neighbors
                      index3 = nbors1%neighbors(l)
                      atom3 = atoms(index3)
                      separation3_unknown = .true.
                      atom_list = (/ atom2, atom1, atom3 /)
                      
                      ! the condition for finding each triplet once is such that
                      ! index 3 must be higher than the index of the atom whose
                      ! neighbors are NOT currently searched
                      if(index3 > index2)then

                         ! search for the first bond params containing the parameters for atom1-atom2
                         do k1 = 1, size(bond_indices)
                         
                            bond_params(1) = bond_factors(bond_indices(k1))
                            call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                                 n_targets)
                            call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2)
                            call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group)

                            if( is_active .and. is_in_group .and. n_targets == 3 )then
                               call bond_order_factor_affects_atom(bond_params(1),atom3,is_active,3)
                            
                               if( is_active )then


                         ! search for the second bond params containing the parameters for atom1-atom3
                         do k2 = 1, size(bond_indices)
                         
                            bond_params(2) = bond_factors(bond_indices(k2))
                            call get_number_of_targets_of_bond_order_factor_index(bond_params(2)%type_index,&
                                 n_targets)
                            call bond_order_factor_affects_atom(bond_params(2),atom2,is_active,2)
                            call bond_order_factor_is_in_group(bond_params(2),group_index,is_in_group)
                         
                            if( is_active .and. is_in_group .and. n_targets == 3 )then
                               call bond_order_factor_affects_atom(bond_params(2),atom3,is_active,3)
                            
                               if( is_active )then

                                  if( separation3_unknown )then
                                     call separation_vector(atom1%position, &
                                          atom3%position, &
                                          nbors1%pbc_offsets(1:3,j), &
                                          cell, &
                                          separations(1:3,2))
                                     separation3_unknown = .false.
                                     distances(2) = .norm.(separations(1:3,2))
                                     if(distances(2) == 0.d0)then
                                        directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                     else
                                        directions(1:3,2) = separations(1:3,2) / distances(2)
                                     end if
                                  end if
                                        
                                  call evaluate_bond_order_factor(3,separations(1:3,1:2),&
                                       distances(1:2),bond_params(1:2),tmp_factor(1:3),atom_list)

                                  ! bond factor:
                                  bond_orders(index2) = bond_orders(index2) + tmp_factor(1)
                                  bond_orders(index1) = bond_orders(index1) + tmp_factor(2)
                                  bond_orders(index3) = bond_orders(index3) + tmp_factor(3)
                                        
                               end if ! is_active
                            end if ! is_active and is_in_group and n_targets == 3
                            
                         end do ! k2
   
                               end if ! is_active
                            end if ! is_active and is_in_group and n_targets == 3
                            
                         end do ! k1
                      end if ! index3 > index2

                   end do ! l = 1, nbors1%n_neighbors


                   ! Next we try to find ordered triplets atom1 -- atom2 -- atom3
                   ! Therefore we need separations a2--a1 and a2--a3.
                   separations(1:3,1) = -separations(1:3,1)

                   ! loop over neighbors of atom 2
                   do l = 1, nbors2%n_neighbors
                      index3 = nbors2%neighbors(l)
                      atom3 = atoms(index3)
                      separation3_unknown = .true.
                      atom_list = (/ atom1, atom2, atom3 /)
                      
                      if(index3 > index1)then

                         ! search for the first bond params containing the parameters for atom2-atom1
                         do k1 = 1, size(bond_indices)
                         
                            bond_params(1) = bond_factors(bond_indices(k1))
                            call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,n_targets)
                            call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2)
                            call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group)
                         
                            if( is_active .and. is_in_group .and. n_targets == 3 )then
                               call bond_order_factor_affects_atom(bond_params(1),atom3,is_active,3)
                            
                               if( is_active )then
                               
                         ! search for the second bond params containing the parameters for atom2-atom3
                         do k2 = 1, size(bond_indices)
                         
                            bond_params(2) = bond_factors(bond_indices(k2))
                            call get_number_of_targets_of_bond_order_factor_index(bond_params(2)%type_index,&
                                 n_targets)
                            call bond_order_factor_affects_atom(bond_params(2),atom2,is_active,2)
                            call bond_order_factor_is_in_group(bond_params(2),group_index,is_in_group)
                         
                            if( is_active .and. is_in_group .and. n_targets == 3 )then
                               call bond_order_factor_affects_atom(bond_params(2),atom3,is_active,3)
                            
                               if( is_active )then

                                  if( separation3_unknown )then
                                     call separation_vector(atom2%position, &
                                          atom3%position, &
                                          nbors2%pbc_offsets(1:3,l), &
                                          cell, &
                                          separations(1:3,2))
                                     separation3_unknown = .false.
                                     distances(2) = .norm.(separations(1:3,2))
                                     if(distances(2) == 0.d0)then
                                        directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                     else
                                        directions(1:3,2) = separations(1:3,2) / distances(2)
                                     end if
                                  end if
                                        
                                  call evaluate_bond_order_factor(3,separations(1:3,1:2),&
                                       distances(1:2),bond_params(1:2),tmp_factor(1:3),atom_list)

                                  ! bond factor:
                                  bond_orders(index1) = bond_orders(index1) + tmp_factor(1)
                                  bond_orders(index2) = bond_orders(index2) + tmp_factor(2)
                                  bond_orders(index3) = bond_orders(index3) + tmp_factor(3)
                                     
                               end if ! is_active
                            end if ! is_active and is_in_group and n_targets == 3
                            
                         end do ! k2
     
                               end if ! is_active
                            end if ! is_active and is_in_group and n_targets == 3
                            
                         end do ! k1
                      end if ! index3 > index2

                   end do ! l = 1, nbors2%n_neighbors

                end if ! many_bodies_found
             end if
          end do

       end if

    end do

#ifdef MPI
    call mpi_allreduce(bond_orders,total_bond_orders,size(bond_orders),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
    bond_orders = total_bond_orders
#else
    total_bond_orders = bond_orders
#endif

  end subroutine core_calculate_bond_order_factors


    ! bond-order post processing 
    !
    ! By post processing, we mean any operations done after calculating the
    ! sum of pair- and many-body terms. That is, if a factor is, say,
    !      b_i = 1 + sum_j c_ij,
    ! the 'sum_j c_ij' would have been calculated already and the operation '1 +'
    ! remains to be carried out.
    ! The post processing is done per atom regardless of if the
    ! bond factor is of a pair or many body type.
  subroutine core_post_process_bond_order_factors(n_atoms,group_index,raw_sums,total_bond_orders)
    implicit none
    integer, intent(in) :: n_atoms, group_index
    double precision, intent(out) :: total_bond_orders(n_atoms)
    integer :: index1, index2
    double precision, intent(in) :: raw_sums(n_atoms)
    double precision :: bond_orders(n_atoms)
    type(atom) :: atom1
    type(bond_order_parameters) :: bond_params
    integer, pointer :: bond_indices(:)
    integer :: post_process

    do index1 = 1, size(atoms)
       
       ! in MPI, only consider the atoms allocated to this particular cpu
       if(is_my_atom(index1))then

          ! target atom
          atom1 = atoms(index1)
          bond_indices => atom1%bond_indices

          ! Check all the bond factors that affect the atom and see
          ! if any of them require post processing.
          ! If such a factor is found, the corresponding parameters
          ! are saved to be used for the post processing.
          ! Note that only one set of post processing parameters will
          ! be used: the first found. This is because there may well
          ! be several factors acting on the same type of atom and
          ! we do not want to apply the post processing several times.
          post_process = -1
          do index2 = 1, size(bond_indices)
             bond_params = bond_factors(bond_indices(index2))
             if( bond_params%includes_post_processing )then
                if( bond_params%original_elements(1) == atom1%element )then
                   post_process = bond_indices(index2)
                   exit
                end if
             end if
          end do

          if( post_process > 0 )then
             call post_process_bond_order_factor(raw_sums(index1),&
                  bond_factors( post_process ), &
                  bond_orders(index1) )
          else
             bond_orders(index1) = raw_sums(index1)
          end if

       else
          bond_orders(index1) = 0.d0
       end if

    end do

#ifdef MPI
    call mpi_allreduce(bond_orders,total_bond_orders,size(bond_orders),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
#else
    total_bond_orders = bond_orders
#endif

  end subroutine core_post_process_bond_order_factors

  subroutine core_post_process_bond_order_gradients(n_atoms,group_index,raw_sums,&
       raw_gradients,total_bond_gradients,mpi_split)
    implicit none
    integer, intent(in) :: n_atoms, group_index
    double precision, intent(out) :: total_bond_gradients(3,n_atoms)
    double precision, intent(in) :: raw_sums(n_atoms), raw_gradients(3,n_atoms)
    logical, optional, intent(in) :: mpi_split
    integer :: index1, index2
    double precision :: bond_gradients(3,n_atoms)
    type(atom) :: atom1
    type(neighbor_list) :: nbors1
    type(bond_order_parameters) :: bond_params
    integer, pointer :: bond_indices(:)
    integer :: post_process
    logical :: evaluate

    do index1 = 1, size(atoms)
       
       evaluate = .true.
       ! in MPI, only consider the atoms allocated to this particular cpu
       if(present(mpi_split))then
          if(is_my_atom(index1) .or. .not.mpi_split)then
             evaluate = .true.
          else
             evaluate = .false.
          end if
       end if
       if(evaluate)then
          ! target atom
          atom1 = atoms(index1)
          bond_indices => atom1%bond_indices

          ! Check all the bond factors that affect the atom and see
          ! if any of them require post processing.
          ! If such a factor is found, the corresponding parameters
          ! are saved to be used for the post processing.
          ! Note that only one set of post processing parameters will
          ! be used: the first found. This is because there may well
          ! be several factors acting on the same type of atom and
          ! we do not want to apply the post processing several times.
          post_process = -1
          do index2 = 1, size(bond_indices)
             bond_params = bond_factors(bond_indices(index2))
             if( bond_params%includes_post_processing )then
                if( bond_params%original_elements(1) == atom1%element )then
                   post_process = bond_indices(index2)
                   exit
                end if
             end if
          end do

          if( post_process > 0 )then
             call post_process_bond_order_gradient(raw_sums(index1),&
                  raw_gradients(1:3,index1),&
                  bond_factors( post_process ), &
                  bond_gradients(1:3,index1) )
          else
             bond_gradients(1:3,index1) = raw_gradients(1:3,index1)
          end if

       else
          bond_gradients(1:3,index1) = 0.d0
       end if

    end do

#ifdef MPI
    if(present(mpi_split))then
       if(mpi_split)then
          call mpi_allreduce(bond_gradients,total_bond_gradients,size(bond_gradients),mpi_double_precision,&
               mpi_sum,mpi_comm_world,mpistat)
       else
          total_bond_gradients = bond_gradients
       end if
    else
       total_bond_gradients = bond_gradients
    end if
#else
    total_bond_gradients = bond_gradients
#endif

  end subroutine core_post_process_bond_order_gradients



  subroutine core_post_process_bond_order_gradients_of_factor(n_atoms,group_index,raw_sum,&
       raw_gradients,total_bond_gradients,mpi_split)
    implicit none
    integer, intent(in) :: n_atoms, group_index
    double precision, intent(out) :: total_bond_gradients(3,n_atoms)
    double precision, intent(in) :: raw_sum, raw_gradients(3,n_atoms)
    logical, optional, intent(in) :: mpi_split
    integer :: index1, index2
    double precision :: bond_gradients(3,n_atoms)
    type(atom) :: atom1
    type(neighbor_list) :: nbors1
    type(bond_order_parameters) :: bond_params
    integer, pointer :: bond_indices(:)
    integer :: post_process
    logical :: evaluate

    do index1 = 1, size(atoms)
       
       evaluate = .true.
       ! in MPI, only consider the atoms allocated to this particular cpu
       if(present(mpi_split))then
          if(is_my_atom(index1) .or. .not.mpi_split)then
             evaluate = .true.
          else
             evaluate = .false.
          end if
       end if
       if(evaluate)then
          ! target atom
          atom1 = atoms(index1)
          bond_indices => atom1%bond_indices

          ! Check all the bond factors that affect the atom and see
          ! if any of them require post processing.
          ! If such a factor is found, the corresponding parameters
          ! are saved to be used for the post processing.
          ! Note that only one set of post processing parameters will
          ! be used: the first found. This is because there may well
          ! be several factors acting on the same type of atom and
          ! we do not want to apply the post processing several times.
          post_process = -1
          do index2 = 1, size(bond_indices)
             bond_params = bond_factors(bond_indices(index2))
             if( bond_params%includes_post_processing )then
                if( bond_params%original_elements(1) == atom1%element )then
                   post_process = bond_indices(index2)
                   exit
                end if
             end if
          end do

          if( post_process > 0 )then
             call post_process_bond_order_gradient(raw_sum,&
                  raw_gradients(1:3,index1),&
                  bond_factors( post_process ), &
                  bond_gradients(1:3,index1) )
          else
             bond_gradients(1:3,index1) = raw_gradients(1:3,index1)
          end if

       else
          bond_gradients(1:3,index1) = 0.d0
       end if

    end do

#ifdef MPI
    if(present(mpi_split))then
       if(mpi_split)then
          call mpi_allreduce(bond_gradients,total_bond_gradients,size(bond_gradients),mpi_double_precision,&
               mpi_sum,mpi_comm_world,mpistat)
       else
          total_bond_gradients = bond_gradients
       end if
    else
       total_bond_gradients = bond_gradients
    end if
#else
    total_bond_gradients = bond_gradients
#endif

  end subroutine core_post_process_bond_order_gradients_of_factor



  subroutine core_calculate_forces(n_atoms,total_forces)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(out) :: total_forces(3,n_atoms)
    integer :: j, k, l, n_targets, index1, index2, index3, &
         indexB1, indexB2, indexB3, iB, jB, kB
    double precision :: forces(3,n_atoms), tmp_forces(3,3), &
         separations(3,2), distances(2), directions(3,2), &
         dummy_sep(3,0), dummy_dist(0), &
         cut_factors(2), cut_gradients(3,2), tmp_energy, stopwatch_0, stopwatch_1
    double precision :: bo_factors(n_atoms), bo_sums(n_atoms), bo_gradients(3,n_atoms,3)
    type(atom) :: atom1, atom2, atom3, atomB1, atomB2, atomB3
    type(atom) :: atom_list(3), atom_listB(3)
    type(neighbor_list) :: nbors1, nbors2
    type(potential) :: interaction
    integer, pointer :: interaction_indices(:)
    logical :: is_active, many_bodies_found, separation3_unknown

    forces = 0.d0
    total_forces = 0.d0

    call start_timer()

    bo_factors = 1.d0
    bo_sums = 0.d0
    bo_gradients = 0.d0
    use_saved_bond_order_factors = .true.
    call core_fill_bond_order_storage(n_atoms)

    ! loop over atoms
    do index1 = 1, size(atoms)

       ! in MPI, only consider the atoms allocated to this particular cpu
       if(is_my_atom(index1))then
          
          call core_empty_bond_order_gradient_storage()

          atom1 = atoms(index1)
          nbors1 = atom1%neighbor_list
          interaction_indices => atom1%potential_indices
          
          ! 1-body interactions
          do k = 1, size(interaction_indices)
             
             interaction = interactions(interaction_indices(k))
             call get_number_of_targets_of_potential_index(interaction%type_index,n_targets)
             
             if( n_targets == 1 )then

                call evaluate_forces(1,dummy_sep,dummy_dist,interaction,tmp_forces(1:3,1),atoms(index1:index1))                
                if(interaction%pot_index > -1)then
                   call core_get_bond_order_factors(n_atoms,&
                        interaction%pot_index,&
                        bo_factors)
                   call core_get_bond_order_gradients(n_atoms,&
                        interaction%pot_index,&
                        index1,& ! atom index
                        1, & ! slot_index
                        bo_gradients(1:3,1:n_atoms,1))

                   ! add the bond order gradient terms
                   call evaluate_energy(1,dummy_sep,dummy_dist,interaction,tmp_energy,atoms(index1:index1))
                   forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) - tmp_energy*bo_gradients(1:3,1:n_atoms,1)

                else
                   bo_factors = 1.d0
                   bo_sums = 0.d0
                   bo_gradients = 0.d0
                end if

                forces(1:3,index1) = forces(1:3,index1) + tmp_forces(1:3,1)*bo_factors(index1)
                
             end if
          end do
          
          ! loop over neighbors
          do j = 1, nbors1%n_neighbors
             
             ! Note that we loop over the neighbors in the outer loop and
             ! over the interactions in the inner loop. This is to avoid calculating
             ! the interatomic distances repeatedly for multiple potentials affecting
             ! the same pair of atoms.
             
             index2 = nbors1%neighbors(j)
             if(index2 > index1)then ! do not double count
                
                call core_empty_bond_order_gradient_storage(2)
                atom2 = atoms(index2)
                atom_list(1) = atom1
                atom_list(2) = atom2
                call separation_vector(atom1%position, &
                     atom2%position, &
                     nbors1%pbc_offsets(1:3,j), &
                     cell, &
                     separations(1:3,1))
                distances(1) = .norm.(separations(1:3,1))
                if(distances(1) == 0.d0)then
                   directions(1:3,1) = (/ 0.d0, 0.d0, 0.d0 /)
                else
                   directions(1:3,1) = separations(1:3,1) / distances(1)
                end if
                
                many_bodies_found = .false.
                ! 2-body interactions
                do k = 1, size(interaction_indices)
                   
                   interaction = interactions(interaction_indices(k))
                   call get_number_of_targets_of_potential_index(interaction%type_index,n_targets)
                   call potential_affects_atom(interaction,atom2,is_active,2)
                   
                   if( is_active .and. interaction%cutoff > distances(1) )then
                      if( n_targets == 2 )then

                         if((interaction%pot_index > -1) .or. &
                            interaction%smoothened)then
                            call evaluate_energy(2,separations(1:3,1),distances(1),&
                                 interaction,tmp_energy,atom_list(1:2))
                         else
                            tmp_energy = 0.d0
                         end if

                         if(interaction%smoothened)then
                            call smoothening_factor(distances(1),&
                                 interaction%cutoff,interaction%soft_cutoff,&
                                 cut_factors(1))
                            call smoothening_gradient(directions(1:3,1),distances(1),&
                                 interaction%cutoff,interaction%soft_cutoff,&
                                 cut_gradients(1:3,1))
                         else
                            cut_factors(1) = 1.d0
                            cut_gradients(1:3,1) = 0.d0
                         end if

                         if(interaction%pot_index > -1)then
                            call core_get_bond_order_factors(n_atoms,&
                                 interaction%pot_index,&
                                 bo_factors)
                            call core_get_bond_order_gradients(n_atoms,&
                                 interaction%pot_index,&
                                 index1,& ! atom index
                                 1, & ! slot_index
                                 bo_gradients(1:3,1:n_atoms,1))
                            call core_get_bond_order_gradients(n_atoms,&
                                 interaction%pot_index,&
                                 index2,& ! atom index
                                 2, & ! slot_index
                                 bo_gradients(1:3,1:n_atoms,2))

                            forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) &
                                 - tmp_energy*cut_factors(1)*&
                                 (bo_gradients(1:3,1:n_atoms,1)+bo_gradients(1:3,1:n_atoms,2))*0.5d0

                         else
                            bo_factors = 1.d0
                            bo_sums = 0.d0
                            bo_gradients = 0.d0
                         end if

                         ! the effect of smooth cutoff:
                         ! -D (f V) = - (D f) V - f (D V)
                         call evaluate_forces(2,separations(1:3,1),distances(1),interaction,tmp_forces(1:3,1:2))
                         forces(1:3,index1) = forces(1:3,index1) + &
                              ( tmp_forces(1:3,1) * cut_factors(1) + &
                              tmp_energy * cut_gradients(1:3,1) ) * &
                              ( bo_factors(index1) +  bo_factors(index2) ) * 0.5d0
                         forces(1:3,index2) = forces(1:3,index2) + &
                              ( tmp_forces(1:3,2) * cut_factors(1) - &
                              tmp_energy * cut_gradients(1:3,1) ) * &
                              ( bo_factors(index1) +  bo_factors(index2) ) * 0.5d0
                         
                      else if( n_targets > 2)then
                         many_bodies_found = .true.
                      end if
                   end if

                end do
                
                if(many_bodies_found)then

                   nbors2 = atom2%neighbor_list
                   
                   ! loop over neighbors of atom 1
                   do l = 1, nbors1%n_neighbors
                      index3 = nbors1%neighbors(l)
                      atom3 = atoms(index3)
                      separation3_unknown = .true.
                      atom_list = (/ atom2, atom1, atom3 /)
                      
                      ! the condition for finding each triplet once is such that
                      ! index 3 must be higher than the index of the atom whose
                      ! neighbors are NOT currently searched
                      if(index3 > index2)then
                         !call core_empty_bond_order_gradient_storage(3)
                         do k = 1, size(interaction_indices)
                            
                            interaction = interactions(interaction_indices(k))
                            call get_number_of_targets_of_potential_index(interaction%type_index,n_targets)
                            call potential_affects_atom(interaction,atom2,is_active,2)
                            
                            if( is_active .and. n_targets == 3 )then
                               call potential_affects_atom(interaction,atom3,is_active,3)                            
                               if( is_active )then
                                  
                                  ! The ordered triplet found is atom2 -- atom1 -- atom3
                                  ! Calculate the separations and distances between the particles
                                  ! starting from atom1: a1--a2, a1--a3
                                  ! (atom2 -- atom1 is already known though from 2-body calculation)
                                  
                                  if( separation3_unknown )then
                                     call separation_vector(atom1%position, &
                                          atom3%position, &
                                          nbors1%pbc_offsets(1:3,j), &
                                          cell, &
                                          separations(1:3,2))
                                     separation3_unknown = .false.
                                     distances(2) = .norm.(separations(1:3,2))
                                     if(distances(2) == 0.d0)then
                                        directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                     else
                                        directions(1:3,2) = separations(1:3,2) / distances(2)
                                     end if
                                  end if
                                  
                                  if( interaction%cutoff > distances(2) )then
                                     
                                     if((interaction%pot_index > -1) .or. &
                                          interaction%smoothened)then
                                        call evaluate_energy(3,separations(1:3,1:2),distances(1:2),&
                                             interaction,tmp_energy,atom_list)
                                     else
                                        tmp_energy = 0.d0
                                     end if

                                     if(interaction%smoothened)then
                                        call smoothening_factor(distances(1),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_factors(1))
                                        call smoothening_gradient(directions(1:3,1),distances(1),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_gradients(1:3,1))
                                        call smoothening_factor(distances(2),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_factors(2))
                                        call smoothening_gradient(directions(1:3,2),distances(2),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_gradients(1:3,2))
                                     else
                                        cut_factors(1:2) = 1.d0
                                        cut_gradients(1:3,1:2) = 0.d0
                                     end if
                                     
                                     if(interaction%pot_index > -1)then
                                        call core_get_bond_order_factors(n_atoms,&
                                             interaction%pot_index,&
                                             bo_factors)
                                        call core_get_bond_order_gradients(n_atoms,&
                                             interaction%pot_index,&
                                             index1,& ! atom index
                                             1, & ! slot_index
                                             bo_gradients(1:3,1:n_atoms,1))
                                        call core_get_bond_order_gradients(n_atoms,&
                                             interaction%pot_index,&
                                             index2,& ! atom index
                                             2, & ! slot_index
                                             bo_gradients(1:3,1:n_atoms,2))
                                        call core_get_bond_order_gradients(n_atoms,&
                                             interaction%pot_index,&
                                             index3,& ! atom index
                                             3, & ! slot_index
                                             bo_gradients(1:3,1:n_atoms,3))
                                        
                                        forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) &
                                             - tmp_energy*cut_factors(1)*cut_factors(2)*&
                                             ( bo_gradients(1:3,1:n_atoms,1) &
                                             + bo_gradients(1:3,1:n_atoms,2) &
                                             + bo_gradients(1:3,1:n_atoms,3) )/3.d0
                                        
                                     else
                                        bo_factors = 1.d0
                                        bo_sums = 0.d0
                                        bo_gradients = 0.d0
                                     end if

                                     call evaluate_forces(3,separations(1:3,1:2),distances(1:2),interaction,&
                                          tmp_forces(1:3,1:3),atom_list)
                                     
                                     ! the influence of smooth cutoffs:
                                     ! -D (f1 f2 V) = - (D f1) f2 V - f1 (D f2) V - f1 f2 (D V)
                                     ! f1: a1 -- a2, f2: a1 -- a3
                                     
                                     ! force on atom 2:
                                     forces(1:3,index2) = forces(1:3,index2) + &
                                          ( tmp_forces(1:3,1)*cut_factors(1)*cut_factors(2) - &
                                          cut_gradients(1:3,1)*cut_factors(2)*tmp_energy ) * &
                                          ( bo_factors(index1) &
                                          + bo_factors(index2) &
                                          + bo_factors(index3) )/3.d0
                                     
                                     ! force on atom 1:
                                     forces(1:3,index1) = forces(1:3,index1) + &
                                          ( tmp_forces(1:3,2)*cut_factors(1)*cut_factors(2) + &
                                          (cut_gradients(1:3,1)*cut_factors(2) + &
                                          cut_gradients(1:3,2)*cut_factors(1)) * tmp_energy ) * &
                                          ( bo_factors(index1) &
                                          + bo_factors(index2) &
                                          + bo_factors(index3) )/3.d0
                                     
                                     ! force on atom 3:
                                     forces(1:3,index3) = forces(1:3,index3) + &
                                          ( tmp_forces(1:3,3)*cut_factors(1)*cut_factors(2) - &
                                          cut_gradients(1:3,2)*cut_factors(1)*tmp_energy ) * &
                                          ( bo_factors(index1) &
                                          + bo_factors(index2) &
                                          + bo_factors(index3) )/3.d0
                                     
                                  end if
                                  
                               end if ! is_active
                            end if ! is_active .and. n_targets == 3
                            
                         end do ! k
                      end if ! index3 > index2
                      
                   end do
                   
                   ! Next we try to find ordered triplets atom1 -- atom2 -- atom3
                   ! Therefore we need separations a2--a1 and a2--a3.
                   separations(1:3,1) = -separations(1:3,1)

                   ! loop over neighbors of atom 2
                   do l = 1, nbors2%n_neighbors
                      index3 = nbors2%neighbors(l)
                      atom3 = atoms(index3)
                      separation3_unknown = .true.
                      atom_list = (/ atom1, atom2, atom3 /)
                      
                      if(index3 > index1)then
                         !call core_empty_bond_order_gradient_storage(3)
                         do k = 1, size(interaction_indices)
                            
                            interaction = interactions(interaction_indices(k))
                            call get_number_of_targets_of_potential_index(interaction%type_index,n_targets)
                            call potential_affects_atom(interaction,atom2,is_active,2)
                            
                            if( is_active .and. n_targets == 3 )then
                               call potential_affects_atom(interaction,atom3,is_active,3)
                               if( is_active )then
                                  
                                  if( separation3_unknown )then
                                     call separation_vector(atom2%position, &
                                          atom3%position, &
                                          nbors2%pbc_offsets(1:3,l), &
                                          cell, &
                                          separations(1:3,2))
                                     separation3_unknown = .false.
                                     distances(2) = .norm.(separations(1:3,2))
                                     if(distances(2) == 0.d0)then
                                        directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                     else
                                        directions(1:3,2) = separations(1:3,2)/distances(2)
                                     end if
                                  end if
                                  
                                  if( interaction%cutoff > distances(2) )then
                                     
                                     if((interaction%pot_index > -1) .or. &
                                          interaction%smoothened)then
                                        call evaluate_energy(3,separations(1:3,1:2),distances(1:2),&
                                             interaction,tmp_energy,atom_list)
                                     else
                                        tmp_energy = 0.d0
                                     end if

                                     if(interaction%smoothened)then
                                        call smoothening_factor(distances(1),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_factors(1))
                                        call smoothening_gradient(directions(1:3,1),distances(1),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_gradients(1:3,1))
                                        call smoothening_factor(distances(2),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_factors(2))
                                        call smoothening_gradient(directions(1:3,2),distances(2),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_gradients(1:3,2))
                                     else
                                        cut_factors(1:2) = 1.d0
                                        cut_gradients(1:3,1:2) = 0.d0
                                     end if
                                     
                                     if(interaction%pot_index > -1)then
                                        call core_get_bond_order_factors(n_atoms,&
                                             interaction%pot_index,&
                                             bo_factors)
                                        call core_get_bond_order_gradients(n_atoms,&
                                             interaction%pot_index,&
                                             index1,& ! atom index
                                             1, & ! slot_index
                                             bo_gradients(1:3,1:n_atoms,1))
                                        call core_get_bond_order_gradients(n_atoms,&
                                             interaction%pot_index,&
                                             index2,& ! atom index
                                             2, & ! slot_index
                                             bo_gradients(1:3,1:n_atoms,2))
                                        call core_get_bond_order_gradients(n_atoms,&
                                             interaction%pot_index,&
                                             index3,& ! atom index
                                             3, & ! slot_index
                                             bo_gradients(1:3,1:n_atoms,3))
                                        
                                        forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) &
                                             - tmp_energy*cut_factors(1)*cut_factors(2)*&
                                             ( bo_gradients(1:3,1:n_atoms,1) &
                                             + bo_gradients(1:3,1:n_atoms,2) &
                                             + bo_gradients(1:3,1:n_atoms,3) )/3.d0
                                        
                                     else
                                        bo_factors = 1.d0
                                        bo_sums = 0.d0
                                        bo_gradients = 0.d0
                                     end if

                                     call evaluate_forces(3,separations(1:3,1:2),distances(1:2),interaction,&
                                          tmp_forces(1:3,1:3),atom_list)
                                     
                                     ! the influence of smooth cutoffs:
                                     ! -D (f1 f2 V) = - (D f1) f2 V - f1 (D f2) V - f1 f2 (D V)
                                     ! f1: a1 -- a2, f2: a2 -- a3
                                     
                                     ! force on atom 1:
                                     forces(1:3,index1) = forces(1:3,index1) + &
                                          ( tmp_forces(1:3,1)*cut_factors(1)*cut_factors(2) + &
                                          cut_gradients(1:3,1)*cut_factors(2)*tmp_energy ) * &
                                          ( bo_factors(index1) &
                                          + bo_factors(index2) &
                                          + bo_factors(index3) )/3.d0
                                     
                                     ! force on atom 2:
                                     forces(1:3,index2) = forces(1:3,index2) + &
                                          ( tmp_forces(1:3,2)*cut_factors(1)*cut_factors(2) - &
                                          cut_gradients(1:3,1)*cut_factors(2)*tmp_energy + &
                                          cut_gradients(1:3,2)*cut_factors(1)*tmp_energy ) * &
                                          ( bo_factors(index1) &
                                          + bo_factors(index2) &
                                          + bo_factors(index3) )/3.d0
                                     
                                     ! force on atom 3:
                                     forces(1:3,index3) = forces(1:3,index3) + &
                                          ( tmp_forces(1:3,3)*cut_factors(1)*cut_factors(2) - &
                                          cut_gradients(1:3,2)*cut_factors(1)*tmp_energy ) * &
                                          ( bo_factors(index1) &
                                          + bo_factors(index2) &
                                          + bo_factors(index3) )/3.d0
                                     
                                  end if ! cutoff
                                  
                               end if ! is_active
                            end if ! is_active .and. n_targets == 3
                            
                         end do ! k
                      end if ! index3 > index1
                      
                   end do ! l
                   
                end if ! many-bodies_found
                
             end if ! index2 < i
             
          end do ! j = 1, size(nbors%neighbors)

       end if ! is_my_atom
    end do ! i = 1, size(atoms)

    call timer(stopwatch_0)

#ifdef MPI
    call record_load(stopwatch_0)
    call balance_loads()
#endif

#ifdef MPI
    call mpi_allreduce(forces,total_forces,size(forces),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
#else
    total_forces = forces
#endif

  end subroutine core_calculate_forces


  subroutine core_calculate_energy(n_atoms,total_energy)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(out) :: total_energy
    integer :: j, k, l, n_targets, index1, index2, index3
    double precision :: energy, tmp_energy, &
         separations(3,2), distances(2), &
         dummy_sep(3,0), dummy_dist(0), &
         cut_factors(2)
    double precision :: bo_factors(n_atoms)
    type(atom) :: atom1, atom2, atom3
    type(atom) :: atom_list(3)
    type(neighbor_list) :: nbors1, nbors2
    type(potential) :: interaction
    integer, pointer :: interaction_indices(:)
    logical :: is_active, many_bodies_found, separation3_unknown

    energy = 0.d0
    total_energy = 0.d0

    bo_factors = 1.d0
    use_saved_bond_order_factors = .true.
    call core_fill_bond_order_storage(n_atoms)

    do index1 = 1, n_atoms

       ! in MPI, only consider the atoms allocated to this particular cpu
       if(is_my_atom(index1))then
          
          atom1 = atoms(index1)
          nbors1 = atom1%neighbor_list
          interaction_indices => atom1%potential_indices

          ! 1-body interactions
          do k = 1, size(interaction_indices)

             interaction = interactions(interaction_indices(k))
             call get_number_of_targets_of_potential_index(interaction%type_index,n_targets)

             if( n_targets == 1 )then

                if(interaction%pot_index > -1)then
                   call core_get_bond_order_factors(n_atoms,&
                        interaction%pot_index,&
                        bo_factors)
                else
                   bo_factors = 1.d0
                end if

                call evaluate_energy(1,dummy_sep,dummy_dist,interaction,tmp_energy,atoms(index1:index1))
                energy = energy + tmp_energy*bo_factors(index1)

             end if
          end do

          do j = 1, nbors1%n_neighbors

             index2 = nbors1%neighbors(j)
             if(index2 > index1)then

                atom2 = atoms(index2)
                call separation_vector(atom1%position, &
                     atom2%position, &
                     nbors1%pbc_offsets(1:3,j), &
                     cell, &
                     separations(1:3,1))
                distances(1) = .norm.(separations(1:3,1))

                many_bodies_found = .false.
                ! 2-body interactions
                do k = 1, size(interaction_indices)

                   interaction = interactions(interaction_indices(k))
                   call get_number_of_targets_of_potential_index(interaction%type_index,n_targets)
                   call potential_affects_atom(interaction,atom2,is_active,2)

                   if( is_active .and. interaction%cutoff > distances(1) )then  
                      if( n_targets == 2 )then
                         if(interaction%pot_index > -1)then
                            call core_get_bond_order_factors(size(atoms),&
                                 interaction%pot_index,&
                                 bo_factors)
                         else
                            bo_factors = 1.d0
                         end if
                         if(interaction%smoothened)then
                            call smoothening_factor(distances(1),&
                                 interaction%cutoff,interaction%soft_cutoff,&
                                 cut_factors(1))
                         else
                            cut_factors(1) = 1.d0
                         end if
                         call evaluate_energy(2,separations(1:3,1),distances(1),interaction,tmp_energy)
                         energy = energy + tmp_energy*cut_factors(1)*&
                              (bo_factors(index1)+bo_factors(index2))*0.5d0
                      else if(n_targets > 2)then
                         many_bodies_found = .true.
                      end if

                   end if ! is_active
                end do ! k

                if(many_bodies_found)then

                   nbors2 = atom2%neighbor_list

                   ! loop over neighbors atom 1
                   do l = 1, nbors1%n_neighbors
                      index3 = nbors1%neighbors(l)
                      atom3 = atoms(index3)
                      separation3_unknown = .true.
                      atom_list = (/ atom2, atom1, atom3 /)

                      ! the condition for finding each triplet once is such that
                      ! index 3 must be higher than the index of the atom whose
                      ! neighbors are NOT currently searched
                      if(index3 > index2)then
                         do k = 1, size(interaction_indices)

                            interaction = interactions(interaction_indices(k))
                            call get_number_of_targets_of_potential_index(interaction%type_index,n_targets)
                            call potential_affects_atom(interaction,atom2,is_active,2)

                            if( is_active .and. n_targets == 3 .and. interaction%cutoff > distances(1) )then
                               call potential_affects_atom(interaction,atom3,is_active,3)                            
                               if( is_active )then

                                  ! The ordered triplet found is atom2 -- atom1 -- atom3
                                  ! Calculate the separations and distances between the particles
                                  ! starting from atom1: a1--a2, a1--a3
                                  ! (atom2 -- atom1 is already known though from 2-body calculation)

                                  if( separation3_unknown )then
                                     call separation_vector(atom1%position, &
                                          atom3%position, &
                                          nbors1%pbc_offsets(1:3,j), &
                                          cell, &
                                          separations(1:3,2))
                                     separation3_unknown = .false.
                                     distances(2) = .norm.(separations(1:3,2))
                                  end if

                                  if( interaction%cutoff > distances(2) )then

                                     if(interaction%pot_index > -1)then
                                        call core_get_bond_order_factors(n_atoms,&
                                             interaction%pot_index,&
                                             bo_factors)
                                     else
                                        bo_factors = 1.d0
                                     end if
                                     if(interaction%smoothened)then
                                        call smoothening_factor(distances(1),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_factors(1))
                                        call smoothening_factor(distances(2),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_factors(2))
                                     else
                                        cut_factors(1:2) = 1.d0
                                     end if

                                     call evaluate_energy(3,separations(1:3,1:2),distances(1:2),&
                                          interaction,tmp_energy,atom_list)
                                     energy = energy + tmp_energy*cut_factors(1)*cut_factors(2)*&
                                          (bo_factors(index1)+bo_factors(index2)+bo_factors(index3))/3.d0

                                  end if

                               end if
                            end if

                         end do ! k
                      end if ! index3 > index2

                   end do

                   ! Next we try to find ordered triplets atom1 -- atom2 -- atom3
                   ! Therefore we need separations a2--a1 and a2--a3.
                   separations(1:3,1) = -separations(1:3,1)

                   ! loop over neighbors of atom 2
                   do l = 1, nbors2%n_neighbors
                      index3 = nbors2%neighbors(l)
                      atom3 = atoms(index3)
                      separation3_unknown = .true.
                      atom_list = (/ atom1, atom2, atom3 /)

                      if(index3 > index1)then
                         do k = 1, size(interaction_indices)

                            interaction = interactions(interaction_indices(k))
                            call get_number_of_targets_of_potential_index(interaction%type_index,n_targets)
                            call potential_affects_atom(interaction,atom2,is_active,2)

                            if( is_active .and. n_targets == 3 )then
                               call potential_affects_atom(interaction,atom3,is_active,3)
                               if( is_active )then

                                  if( separation3_unknown )then
                                     call separation_vector(atom2%position, &
                                          atom3%position, &
                                          nbors1%pbc_offsets(1:3,j), &
                                          cell, &
                                          separations(1:3,2))
                                     separation3_unknown = .false.
                                     distances(2) = .norm.(separations(1:3,2))
                                  end if

                                  if( interaction%cutoff > distances(2) )then

                                     if(interaction%pot_index > -1)then
                                        call core_get_bond_order_factors(n_atoms,&
                                             interaction%pot_index,&
                                             bo_factors)
                                     else
                                        bo_factors = 1.d0
                                     end if
                                     if(interaction%smoothened)then
                                        call smoothening_factor(distances(1),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_factors(1))
                                        call smoothening_factor(distances(2),&
                                             interaction%cutoff,interaction%soft_cutoff,&
                                             cut_factors(2))
                                     else
                                        cut_factors(1:2) = 1.d0
                                     end if

                                     call evaluate_energy(3,separations(1:3,1:2),distances(1:2),&
                                          interaction,tmp_energy,atom_list)
                                     energy = energy + tmp_energy*cut_factors(1)*cut_factors(2)*&
                                          (bo_factors(index1)+bo_factors(index2)+bo_factors(index3))/3.d0

                                  end if

                               end if
                            end if

                         end do ! k
                      end if ! index3 > index1

                   end do

                end if ! many-bodies_found

             end if ! index2 > index1

          end do ! j

       end if ! is_my_atom
    end do ! index1

#ifdef MPI
    call mpi_allreduce(energy,total_energy,1,mpi_double_precision,mpi_sum,&
        mpi_comm_world,mpistat)
#else
    total_energy = energy
#endif

    use_saved_bond_order_factors = .false.
    call core_empty_bond_order_storage()

  end subroutine core_calculate_energy


  subroutine list_atoms()
    implicit none
    integer :: i

    if(associated(atoms))then
       write(*,*) "element, position, charge"    
       do i = 1, size(atoms)
          write(*,'(I5,A,F10.4,F10.4,F10.4,F10.4)') i, atoms(i)%element, &
               atoms(i)%position(1), atoms(i)%position(2), atoms(i)%position(3), atoms(i)%charge
       end do
    else
       write(*,*) "no atoms defined"
    end if

  end subroutine list_atoms

  subroutine list_cell()
    implicit none
    integer :: i

    write(*,*) "cell"
    do i = 1, 3
        write(*,'(F10.4,F10.4,F10.4)') cell%vectors(1:3,i)
    end do
    write(*,*) "pbc: ", cell%periodic(1:3)

  end subroutine list_cell

  subroutine list_interactions()
    implicit none
    integer :: i
    
    write(*,*) "interactions"
    do i = 1, n_interactions
       write(*,'(A,I5,F10.4)') "type, cutoff ", interactions(i)%type_index, interactions(i)%cutoff
       write(*,*) "params ", interactions(i)%parameters
       if(interactions(i)%filter_elements)then
          write(*,*) "symbols ",interactions(i)%apply_elements
       end if
       if(interactions(i)%filter_tags)then
          write(*,*) "tags ", interactions(i)%apply_tags
       end if
       if(interactions(i)%filter_indices)then
          write(*,*) "indices ", interactions(i)%apply_indices
       end if
       if(interactions(i)%pot_index >= 0)then
          write(*,*) "bond order index ", interactions(i)%pot_index
       end if
       write(*,*) ""
    end do

  end subroutine list_interactions



  subroutine list_bonds()
    implicit none
    integer :: i,j

    write(*,*) "bond order factors"
    do i = 1, n_bond_factors
       write(*,'(A,I5,F10.4)') "type, cutoff ", bond_factors(i)%type_index, bond_factors(i)%cutoff
       write(*,*) "params "
       do j = 1, size(bond_factors(i)%n_params)
          write(*,*) j, " : ", bond_factors(i)%parameters(:,j)
       end do
       write(*,*) "symbols ",bond_factors(i)%apply_elements
       write(*,*) "bond order index ", bond_factors(i)%group_index
    end do

  end subroutine list_bonds


  


end module pysic_core
