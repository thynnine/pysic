!
! Core, true to its name, is the heart of the Fortran core
! of Pysic. It contains the data structures defining the simulation
! geometry and interactions and defines the central routines for
! calculating the total energy of and the forces acting on the system.
!
! Many of the routines in :ref:`pysic_interface` which `f2py`_ interfaces 
! to Python are simply calling routines here.
!
! 
! .. _f2py: http://www.scipy.org/F2py
!
module pysic_core
  use geometry
  use potentials
  use mpi
  implicit none
  
  ! storage for structure

  ! *atoms an array of :data:`atom` objects representing the system
  type(atom), pointer :: atoms(:)
  ! *cell a :data:`supercell` object representing the simulation cell
  type(supercell) :: cell
  ! *interactions an array of :data:`potential` objects representing the interactions
  type(potential), pointer :: interactions(:)
  ! *multipliers a temporary array for storing multiplying potentials before associating them with a master potential
  type(potential), allocatable :: multipliers(:)
  ! *bond_factors an array of :data:`bond_order_parameters` objects representing bond order factors modifying the potentials
  type(bond_order_parameters), pointer :: bond_factors(:)
  ! *n_interactions number of potentials
  ! *n_bond_order_factors number of bond order factors
  ! *n_multi number of temporary product potentials
  integer :: n_interactions = 0, n_bond_factors = 0, n_multi = 0
  ! logical tags monitoring the allocation and deallocation
  ! of the corresponding pointers
  ! *atoms_created logical tag indicating if atom storing arrays have been created
  ! *potentials_allocated logical tag indicating if potential storing arrays have been allocated
  ! *bond_factors_allocated logical tag indicating if bond order parameter storing arrays have been allocated
  ! *bond_storage_allocated logical tag indicating if bond order factor storing arrays have been allocated
  logical :: &
       atoms_created = .false., &
       potentials_allocated = .false., &
       bond_factors_allocated = .false., &
       bond_storage_allocated = .false.

  ! storage for bond order factors during force and energy evaluation

  ! arrays for storing calculated bond order sums, factors and gradients
  ! 

  ! *saved_bond_order_sums Array for storing calculated bond order sums. Indexing: (atom index, group_index_save_slot(group index))
  double precision, pointer :: saved_bond_order_sums(:,:) 
  ! *saved_bond_order_factors Array for storing calculated bond order factors. Indexing: (atom index, group_index_save_slot(group index))
  double precision, pointer :: saved_bond_order_factors(:,:)
  ! *saved_bond_order_gradients Array for storing calculated bond order gradients. Indexing: (xyz, atom index, group_index_save_slot(group index), target index)
  double precision, pointer :: saved_bond_order_gradients(:,:,:,:)
  ! *saved_bond_order_virials Array for storing calculated bond order virials. Indexing: (xyz, group_index_save_slot(group index), target index)
  double precision, pointer :: saved_bond_order_virials(:,:,:)
  ! *n_saved_bond_order_factors number of saved bond order factors
  integer :: n_saved_bond_order_factors = 0
  ! * group_index_save_slot A list joining group indices and bond factor save slots: Group indices are indices for potentials, but not every potential needs to have a bond order factor. Therefore, the saved bond order arrays should have less columns than the number of groups. This array changes the group index into the column index of the saved bond order arrays.
  integer, pointer :: group_index_save_slot(:)
  ! *use_saved_bond_order_factors Logical tag which enables / disables bond order saving. If true, bond order calculation routines try to find the precalculated factors in the saved bond order arrays instead of calculating.
  logical :: use_saved_bond_order_factors = .false.
  ! *use_saved_bond_order_gradients Array storing the atom index of the bond gradient stored for indices (group index, target index). Since gradients are needed for all factors (N) with respect to moving all atoms (N), storing them all would require an N x N matrix. Therefore only some are stored. This array is used for searching the stroage to see if the needed gradient is there or needs to be calculated.
  integer, pointer :: use_saved_bond_order_gradients(:,:)

  ! temporary arrays needed in force and energy evaluation
  integer :: number_of_atoms
  double precision, pointer :: bo_factors(:), bo_gradients(:,:,:), bo_sums(:), bo_temp(:), &
       temp_gradient(:,:,:), temp_forces(:,:), temp_enegs(:)
  integer, pointer :: n_nbs(:), total_n_nbs(:)
  logical, pointer :: bo_scaling(:)
  


  ! ToDo: allow other methods than Ewald for 1/r-summations and make a general framework

  ! *evaluate_ewald switch for enabling Ewald summation of coulomb interactions
  ! *ewald_k_cutoffs the number of reciprocal cells used in the k-space sum of Coulomb interactions in the Ewald method
  ! *ewald_cutoff the real-space cutoff for Ewald summation
  ! *ewald_sigma the splitting parameter :math:`\sigma` for Ewald summation
  ! *ewald_epsilon the electric constant :math:`\varepsilon_0` used in Ewald summation
  ! *ewald_scaler scaling factors for individual charges in Ewald summation
  ! *ewald_allocated logical tag for tracking allocation of the arrays
  logical :: evaluate_ewald = .false.
  integer :: ewald_k_cutoffs(3)
  double precision :: ewald_cutoff, ewald_k_radius, ewald_sigma, ewald_epsilon
  double precision, pointer :: ewald_scaler(:)
  logical :: ewald_allocated = .false.


  ! indices for specifying the type of quantity evaluated during local structure loops
  integer, parameter :: energy_evaluation_index = 1
  integer, parameter :: force_evaluation_index = 2
  integer, parameter :: electronegativity_evaluation_index = 3
  

contains

! !!!: core_release_all_memory

  ! Release all allocated pointer arrays in the core.
  subroutine core_release_all_memory()
    implicit none

    call core_clear_atoms()
    call core_clear_potentials()
    call core_clear_bond_order_factors()
    call core_clear_bond_order_storage()
    call core_clear_ewald_arrays()

  end subroutine core_release_all_memory

! !!!: core_generate_atoms

  ! Creates the atomic particles by invoking a subroutine in the geometry module.
  !
  ! called from PyInterface: :func:`create_atoms`
  !
  ! *n_atoms number of atoms
  ! *masses masses of atoms
  ! *charges electric charges of atoms
  ! *positions coordinates of atoms
  ! *momenta momenta of atoms
  ! *tags numeric tags for the atoms
  ! *elements atomic symbols of the atoms
  subroutine core_generate_atoms(n_atoms,masses,charges,positions,momenta,tags,elements)
    implicit none
    integer, intent(in) :: n_atoms, tags(n_atoms)
    double precision, intent(in) :: masses(n_atoms), charges(n_atoms), positions(3,n_atoms), &
         momenta(3,n_atoms)
    character(len=label_length), intent(in) :: elements(n_atoms)

    call generate_atoms(n_atoms,masses,charges,positions,momenta,tags,elements,atoms) ! in Geometry.f90
    number_of_atoms = n_atoms
    atoms_created = .true.

  end subroutine core_generate_atoms


! !!!: core_update_atom_coordinates

  ! Updates the positions and momenta of atomic particles.
  !
  ! called from PyInterface: :func:`update_atom_coordinates`
  !
  ! *n_atoms number of atoms
  ! *positions new coordinates for the atoms
  ! *momenta new momenta for the atoms
  subroutine core_update_atom_coordinates(n_atoms,positions,momenta)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(in) :: positions(3,n_atoms), momenta(3,n_atoms)

    call update_atomic_positions(n_atoms,positions,momenta,atoms) ! in Geometry.f90

  end subroutine core_update_atom_coordinates
 

! !!!: core_update_atom_charges

  ! Updates the charges of atomic particles.
  !
  ! called from PyInterface: :func:`update_atom_charges`
  !
  ! *n_atoms number of atoms
  ! *charges new charges for the atoms
  subroutine core_update_atom_charges(n_atoms,charges)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(in) :: charges(n_atoms)

    call update_atomic_charges(n_atoms,charges,atoms) ! in Geometry.f90

  end subroutine core_update_atom_charges


! !!!: core_get_number_of_atoms

  ! Returns the number of atoms in the array allocated in the core.
  !
  ! called from PyInterface: :func:`get_number_of_atoms`
  !
  ! *n_atoms number of atoms
  subroutine core_get_number_of_atoms(n_atoms)
    implicit none
    integer, intent(out) :: n_atoms
    
    if(atoms_created)then
       n_atoms = size(atoms)
    else
       n_atoms = 0
    end if

  end subroutine core_get_number_of_atoms


! !!!: core_clear_atoms

  ! Deallocates the array of atoms in the core, if allocated.
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


! !!!: core_create_cell

  ! Creates a supercell for containing the calculation geometry.
  !
  ! called from PyInterface: :func:`create_cell`
  !
  ! *vectors A 3x3 matrix containing the vectors spanning the supercell. The first index runs over xyz and the second index runs over the three vectors.
  ! *inverse A 3x3 matrix containing the inverse matrix of the one given in vectors, i.e. :math:`A*B = I` for the two matrices. Since the latter represents a cell of non-zero volume, this inverse must exist. It is not tested that the given matrix actually is the inverse, the user must make sure it is.
  ! *periodicity A 3-element vector containing logical tags specifying if the system is periodic in the directions of the three vectors spanning the supercell.
  subroutine core_create_cell(vectors,inverse,periodicity)
    implicit none
    double precision, intent(in) :: vectors(3,3), inverse(3,3)
    logical, intent(in) :: periodicity(3)

    call generate_supercell(vectors,inverse,periodicity,cell) ! in Geometry.f90

  end subroutine core_create_cell


! !!!: core_get_cell_vectors

  ! Returns the vectors defining the supercell stored in the core.
  !
  ! called from PyInterface: :func:`get_cell_vectors`
  ! 
  ! *vectors A 3x3 matrix containing the vectors spanning the supercell. The first index runs over xyz and the second index runs over the three vectors.
  subroutine core_get_cell_vectors(vectors)
    implicit none
    double precision, intent(out) :: vectors(3,3)
    
    vectors = cell%vectors

  end subroutine core_get_cell_vectors


! !!!: core_create_neighbor_list

  ! Assigns a precalculated neighbor list to a single atom of the given index.
  ! The neighbor list must be precalculated, this method only
  ! stores them in the core. The list must contain 
  ! an array storing the indices of the neighboring atoms
  ! as well as the supercell offsets. The offsets are integer
  ! triplets showing how many times must the supercell vectors
  ! be added to the position of the neighbor to find the
  ! neighboring image in a periodic system.
  ! For example, let the supercell be::
  !
  !  [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]],
  !
  ! i.e., a unit cube, with periodic boundaries.
  ! Now, if we have particles with coordinates::
  !
  !  a = [1.5, 0.5, 0.5]
  !  b = [0.4, 1.6, 3.3]
  ! 
  ! the closest separation vector :math:`\mathbf{r}_b-\mathbf{r}_a` between the particles is::
  !
  !   [-.1, .1, -.2]
  !
  ! obtained if we add the vector of periodicity::
  !
  !   [1.0, -1.0, -3.0]
  !
  ! to the coordinates of particle b. The offset vector
  ! (for particle b, when listing neighbors of a) is then::
  !
  !   [1, -1, -3]
  !
  ! Note that if the system is small, one atom can in
  ! principle appear several times in the neighbor list with
  ! different offsets.
  !
  ! called from PyInterface: :func:`create_neighbor_list`
  !
  ! *n_nbs number of neighbors
  ! *atom_index index of the atom for which the neighbor list is created
  ! *neighbors An array containing the indices of the neighboring atoms
  ! *offsets An array containing vectors specifying the offsets of the neighbors in periodic systems. 
  subroutine core_create_neighbor_list(n_nbors,atom_index,neighbors,offsets)
    implicit none
    integer, intent(in) :: n_nbors, atom_index
    integer, intent(in) :: neighbors(n_nbors), offsets(3,n_nbors)

    call assign_neighbor_list(n_nbors,atoms(atom_index)%neighbor_list,neighbors,offsets)

  end subroutine core_create_neighbor_list


! !!!: core_clear_potentials

  ! Deallocates pointers for potentials
  subroutine core_clear_potentials()
    implicit none
    integer :: i, j

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
          deallocate(interactions(i)%table)
          deallocate(interactions(i)%multipliers)
       end do
    end if
    if(potentials_allocated)then
       deallocate(interactions)
    end if
    n_interactions = 0
    potentials_allocated = .false.

  end subroutine core_clear_potentials
  

  subroutine core_clear_potential_multipliers()
    implicit none

    if(n_multi > 0)then
       n_multi = 0
       deallocate(multipliers)             
    end if
    
  end subroutine core_clear_potential_multipliers


! !!!: core_clear_bond_order_factors

  ! Deallocates pointers for bond order factors (the parameters)
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



! !!!: core_clear_bond_order_storage

  ! Deallocates pointers for bond order factors (the precalculated factor values).
  subroutine core_clear_bond_order_storage()
    implicit none

    if(bond_storage_allocated)then
       deallocate(saved_bond_order_sums)
       deallocate(saved_bond_order_factors)
       ! saved gradients
       !deallocate(saved_bond_order_gradients)
       !deallocate(saved_bond_order_virials)
       deallocate(use_saved_bond_order_gradients)
       deallocate(group_index_save_slot)

       deallocate(bo_factors)
       deallocate(bo_sums)
       deallocate(bo_gradients)
       deallocate(temp_gradient)
       deallocate(bo_temp)
       deallocate(bo_scaling)
       deallocate(temp_forces)
       deallocate(temp_enegs)
       deallocate(n_nbs)
       deallocate(total_n_nbs)

    else
       nullify(saved_bond_order_sums)
       nullify(saved_bond_order_factors)
       ! saved gradients
       !nullify(saved_bond_order_gradients)
       !nullify(saved_bond_order_virials)
       nullify(use_saved_bond_order_gradients)
       nullify(group_index_save_slot)

       nullify(bo_factors)
       nullify(bo_sums)
       nullify(bo_gradients)
       nullify(temp_gradient)
       nullify(bo_temp)
       nullify(bo_scaling)
       nullify(temp_forces)
       nullify(temp_enegs)
       nullify(n_nbs)
       nullify(total_n_nbs)

    end if
    n_saved_bond_order_factors = 0
    use_saved_bond_order_factors = .false.
    bond_storage_allocated = .false.

  end subroutine core_clear_bond_order_storage


  subroutine core_clear_ewald_arrays()
    implicit none

    if(ewald_allocated)then
       deallocate(ewald_scaler)
    else
       nullify(ewald_scaler)
    end if
    ewald_allocated = .false.

    call deallocate_ewald_arrays()

  end subroutine core_clear_ewald_arrays

! !!!: core_empty_bond_order_storage

  ! Clears bond order factors (the precalculated factor values) 
  ! but does not deallocate the arrays.
  subroutine core_empty_bond_order_storage()
    implicit none
    
    if(bond_storage_allocated)then
       saved_bond_order_sums = 0.d0
       saved_bond_order_factors = 0.d0
       ! saved gradients
       !call core_empty_bond_order_gradient_storage()
       group_index_save_slot = -1
       n_saved_bond_order_factors = 0
    end if
    
  end subroutine core_empty_bond_order_storage



! !!!: core_empty_bond_order_gradient_storage

  ! Clears bond order factor gradients (the precalculated gradient values) 
  ! but does not deallocate the arrays.
  ! If an index is given, then only that column is emptied.
  !
  ! *index the column to be emptied
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
          saved_bond_order_virials = 0.d0
          use_saved_bond_order_gradients = -1
       end if
    end if

  end subroutine core_empty_bond_order_gradient_storage



! !!!: core_allocate_bond_order_storage

  ! Allocates arrays for storing precalculated values of bond order
  ! factors and gradients.
  !
  ! called from PyInterface: :func:`allocate_bond_order_factors`
  !
  ! *n_atoms number of atoms
  ! *n_groups number of bond order groups
  ! *n_factors number of bond order parameters 
  subroutine core_allocate_bond_order_storage(n_atoms,n_groups,n_factors)
    implicit none
    integer, intent(in) :: n_atoms, n_groups, n_factors

    call core_clear_bond_order_storage()
    allocate(saved_bond_order_sums(n_atoms,n_factors))
    allocate(saved_bond_order_factors(n_atoms,n_factors))
    ! saved gradients
    !allocate(saved_bond_order_gradients(3,n_atoms,n_factors,4))
    !allocate(saved_bond_order_virials(6,n_factors,4))
    allocate(use_saved_bond_order_gradients(n_factors,4))
    allocate(group_index_save_slot(0:n_groups))

       allocate(bo_factors(n_atoms))
       allocate(bo_sums(n_atoms))
       allocate(bo_gradients(3,n_atoms,4))
       allocate(bo_temp(n_atoms))
       allocate(bo_scaling(n_atoms))
       allocate(temp_gradient(3,n_atoms,2))
       allocate(temp_forces(3,n_atoms))
       allocate(temp_enegs(n_atoms))
       allocate(n_nbs(n_atoms))
       allocate(total_n_nbs(n_atoms))

       call allocate_ewald_arrays(n_atoms)

    bond_storage_allocated = .true.
    call core_empty_bond_order_storage()

  end subroutine core_allocate_bond_order_storage
  


! !!!: core_fill_bond_order_storage

  ! Fills the storage for bond order factors and bond order sums.
  ! This is meant to be called in the beginning of force and energy
  ! evaluation. The routine calculates all bond order factors
  ! (in parallel, if run in MPI) and stores them. Then during the 
  ! energy or force calculation, it is sufficient to just
  ! look up the needed values in the arrays.
  ! The routine does not calculate and store bond factor gradients.
  !
  ! *n_atoms number of atoms
  subroutine core_fill_bond_order_storage()
    implicit none
    integer :: n_atoms
    integer :: i, group_index

    n_atoms = size(atoms)

    do i = 1, n_interactions
       group_index = interactions(i)%pot_index
       ! If group index is non-negative, a bond order factor is
       ! connected to the potential and needs to be calculated.
       ! The routine called will check if the factors are already calculated 
       ! for this group index and calculates and stores them if not.
       if(group_index > -1)then
          call core_get_bond_order_sums(group_index,bo_factors)
       end if
    end do

  end subroutine core_fill_bond_order_storage




! !!!: core_get_bond_order_gradients

  ! Returns the gradients of the bond order factor of the given atom
  ! with respect to moving all atoms, for the given group.
  ! The routine tries to find the values in the stored precalculated
  ! values first if use_saved_bond_order_factors is true, and saves
  ! the calculated values if it does not find them.
  !
  ! The slot index is the index of the atom in the interaction being
  ! evaluated (so for a triplet A-B-C, A would have slot 1, B slot 2,
  ! and C slot 3). This is only used for storing the values.
  !
  ! *group_index index for the bond order factor group
  ! *atom_index index of the atom whose bond order factor is differentiated
  ! *slot_index index denoting the position of the atom in an interacting group (such as A-B-C triplet)
  ! *bond_order_gradients the calculated gradients of the bond order factor
  ! *bond_order_virial the components of the virial due to the bond order factors
  subroutine core_get_bond_order_gradients(group_index,atom_index,slot_index,&
       bond_order_gradients,bond_order_virial)
    implicit none
    integer, intent(in) ::group_index, atom_index, slot_index
    double precision, intent(inout) :: bond_order_gradients(:,:), bond_order_virial(6)

    call core_get_bond_order_sums(group_index,bo_sums)
    call core_calculate_bond_order_gradients_of_factor(group_index,&
         atom_index,&
         bo_sums,&
         bond_order_gradients,&
         bond_order_virial)

  end subroutine core_get_bond_order_gradients
!!$
!!$
!!$
!!$! !!!: core_get_bond_order_gradients
!!$
!!$  ! Returns the gradients of the bond order factor of the given atom
!!$  ! with respect to moving all atoms, for the given group.
!!$  ! The routine tries to find the values in the stored precalculated
!!$  ! values first if use_saved_bond_order_factors is true, and saves
!!$  ! the calculated values if it does not find them.
!!$  !
!!$  ! The slot index is the index of the atom in the interaction being
!!$  ! evaluated (so for a triplet A-B-C, A would have slot 1, B slot 2,
!!$  ! and C slot 3). This is only used for storing the values.
!!$  !
!!$  ! *group_index index for the bond order factor group
!!$  ! *atom_index index of the atom whose bond order factor is differentiated
!!$  ! *slot_index index denoting the position of the atom in an interacting group (such as A-B-C triplet)
!!$  ! *bond_order_gradients the calculated gradients of the bond order factor
!!$  ! *bond_order_virial the components of the virial due to the bond order factors
!!$  subroutine core_get_bond_order_gradients(group_index,atom_index,slot_index,&
!!$       bond_order_gradients,bond_order_virial)
!!$    implicit none
!!$    integer, intent(in) :: group_index, atom_index, slot_index
!!$    double precision, intent(inout) :: bond_order_gradients(:,:), bond_order_virial(6)
!!$    logical :: found_grads
!!$    integer :: save_slot, n_atoms
!!$
!!$    n_atoms = size(atoms)
!!$    if(use_saved_bond_order_factors)then
!!$       found_grads = .false.
!!$       save_slot = group_index_save_slot(group_index)       
!!$       if(save_slot > 0)then
!!$          if(use_saved_bond_order_gradients(save_slot,slot_index) == atom_index)then
!!$             found_grads = .true.
!!$             bond_order_gradients(1:3,1:n_atoms) = saved_bond_order_gradients(1:3,1:n_atoms,save_slot,slot_index)
!!$             bond_order_virial(1:6) = saved_bond_order_virials(1:6,save_slot,slot_index)
!!$          end if
!!$       end if
!!$       if(.not.found_grads)then
!!$          call core_get_bond_order_sums(group_index,bo_sums)
!!$          call core_calculate_bond_order_gradients_of_factor(&
!!$               group_index,&
!!$               atom_index,&
!!$               bo_sums,&
!!$               bond_order_gradients,&
!!$               bond_order_virial)
!!$          saved_bond_order_gradients(1:3,1:n_atoms,save_slot,slot_index) = bond_order_gradients(1:3,1:n_atoms)
!!$          saved_bond_order_virials(1:6,save_slot,slot_index) = bond_order_virial(1:6)
!!$          use_saved_bond_order_gradients(save_slot,slot_index) = atom_index
!!$       end if
!!$    else
!!$       call core_get_bond_order_sums(group_index,bo_sums)
!!$       call core_calculate_bond_order_gradients_of_factor(&
!!$            group_index,&
!!$            atom_index,&
!!$            bo_sums,&
!!$            bond_order_gradients,&
!!$            bond_order_virial)
!!$    end if
!!$
!!$  end subroutine core_get_bond_order_gradients


! !!!: core_get_bond_order_sums

  ! Returns the bond order sums of all atoms for the given group.
  ! By 'bond order sum', we mean the summation of local terms
  ! without per atom scaling. E.g., for :math:`b_i = 1 + \sum c_{ij}`,
  ! :math:`\sum c_{ij}` is the sum.
  ! The routines tries to find the values in the stored precalculated
  ! values first if use_saved_bond_order_factors is true, and saves
  ! the calculated values if it does not find them.
  !
  ! *group_index index for the bond order factor group
  ! *bond_order_sums the calculated bond order sums
  subroutine core_get_bond_order_sums(group_index,bond_order_sums)
    implicit none
    integer, intent(in) :: group_index
    double precision, intent(inout) :: bond_order_sums(:)
    integer :: save_slot, n_atoms

    n_atoms = size(atoms)
    if(use_saved_bond_order_factors)then
       save_slot = group_index_save_slot(group_index)
       if(save_slot > 0)then
          bond_order_sums(1:n_atoms) = saved_bond_order_sums(1:n_atoms,save_slot)
       else
          call core_calculate_bond_order_factors(group_index,bond_order_sums)
          call core_post_process_bond_order_factors(group_index,bond_order_sums,bo_factors)
          n_saved_bond_order_factors = n_saved_bond_order_factors + 1
          group_index_save_slot(group_index) = n_saved_bond_order_factors
          saved_bond_order_sums(1:n_atoms,n_saved_bond_order_factors) = bond_order_sums(1:n_atoms)
          saved_bond_order_factors(1:n_atoms,n_saved_bond_order_factors) = bo_factors(1:n_atoms)
       end if
    else
       call core_calculate_bond_order_factors(group_index,bond_order_sums)       
    end if

  end subroutine core_get_bond_order_sums




!!$! !!!: core_get_bond_order_sums
!!$
!!$  ! Returns the bond order sums of all atoms for the given group.
!!$  ! By 'bond order sum', we mean the summation of local terms
!!$  ! without per atom scaling. E.g., for :math:`b_i = 1 + \sum c_{ij}`,
!!$  ! :math:`\sum c_{ij}` is the sum.
!!$  ! The routines tries to find the values in the stored precalculated
!!$  ! values first if use_saved_bond_order_factors is true, and saves
!!$  ! the calculated values if it does not find them.
!!$  ! *n_atoms number of atoms
!!$  ! *group_index index for the bond order factor group
!!$  ! *bond_order_sums the calculated bond order sums
!!$  subroutine core_get_bond_order_sums(n_atoms,group_index,bond_order_sums)
!!$    implicit none
!!$    integer, intent(in) :: n_atoms, group_index
!!$    double precision, intent(out) :: bond_order_sums(n_atoms)
!!$    double precision :: bond_order_factors(n_atoms)
!!$    integer :: save_slot
!!$
!!$    if(use_saved_bond_order_factors)then
!!$       save_slot = group_index_save_slot(group_index)
!!$       if(save_slot > 0)then
!!$          bond_order_sums(1:n_atoms) = saved_bond_order_sums(1:n_atoms,save_slot)
!!$       else
!!$          call core_calculate_bond_order_factors(n_atoms,group_index,bond_order_sums)
!!$          call core_post_process_bond_order_factors(n_atoms,group_index,bond_order_sums,bond_order_factors)
!!$          n_saved_bond_order_factors = n_saved_bond_order_factors + 1
!!$          group_index_save_slot(group_index) = n_saved_bond_order_factors
!!$          saved_bond_order_sums(1:n_atoms,n_saved_bond_order_factors) = bond_order_sums(1:n_atoms)
!!$          saved_bond_order_factors(1:n_atoms,n_saved_bond_order_factors) = bond_order_factors(1:n_atoms)
!!$       end if
!!$    else
!!$       call core_calculate_bond_order_factors(n_atoms,group_index,bond_order_sums)       
!!$    end if
!!$
!!$  end subroutine core_get_bond_order_sums



! !!!: core_get_bond_order_factors

  ! Returns the bond order factors of all atoms for the given group.
  ! The routines tries to find the values in the stored precalculated
  ! values first if use_saved_bond_order_factors is true, and saves
  ! the calculated values if it does not find them.
  !
  ! *group_index index for the bond order factor group
  ! *bond_order_factors the calculated bond order factors
  subroutine core_get_bond_order_factors(group_index,bond_order_factors)
    implicit none
    integer, intent(in) :: group_index
    double precision, intent(inout) :: bond_order_factors(:)
    integer :: save_slot, n_atoms

    n_atoms = size(atoms)
    if(use_saved_bond_order_factors)then
       save_slot = group_index_save_slot(group_index)
       if(save_slot > 0)then
          bond_order_factors(1:n_atoms) = saved_bond_order_factors(1:n_atoms,save_slot)
       else
          call core_calculate_bond_order_factors(group_index,bo_sums)
          call core_post_process_bond_order_factors(group_index,bo_sums,bond_order_factors)
          n_saved_bond_order_factors = n_saved_bond_order_factors + 1
          group_index_save_slot(group_index) = n_saved_bond_order_factors
          saved_bond_order_sums(1:n_atoms,n_saved_bond_order_factors) = bo_sums(1:n_atoms)
          saved_bond_order_factors(1:n_atoms,n_saved_bond_order_factors) = bond_order_factors(1:n_atoms)
       end if
    else
       call core_calculate_bond_order_factors(group_index,bo_sums)       
       call core_post_process_bond_order_factors(group_index,bo_sums,bond_order_factors)
    end if

  end subroutine core_get_bond_order_factors


! !!!: core_get_bond_order_factor_of_atom

  ! Returns the bond order factors of the given atom for the given group.
  !
  ! *n_atoms number of atoms
  ! *group_index index for the bond order factor group
  ! *atom_index index of the atom whose bond order factor is returned
  ! *bond_order_factor the calculated bond order factor
  subroutine core_get_bond_order_factor_of_atom(group_index,atom_index,bond_order_factor)
    implicit none
    integer, intent(in) :: group_index, atom_index
    double precision, intent(inout) :: bond_order_factor

    call core_get_bond_order_factors(group_index,bo_factors)
    bond_order_factor = bo_factors(atom_index)

  end subroutine core_get_bond_order_factor_of_atom


! !!!: core_allocate_potentials

  ! Allocates pointers for storing potentials.
  !
  ! called from PyInterface: :func:`allocate_potentials`
  !
  ! *n_pots number of potentials
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


! !!!: core_allocate_bond_order_factors

  ! Allocates pointers for storing bond order factors.
  !
  ! called from PyInterface: :func:`allocate_bond_order_factors`
  !
  ! *n_bond_order_factors number of bond order parameters
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



! !!!: core_add_potential

  ! Creates one additional potential in the core.
  ! The routine assumes that adequate memory has been
  ! allocated already using core_allocate_potentials.
  !
  ! When the potentials in the Python interface are imported
  ! to the Fortran core, the target specifiers (elements, tags, indices)
  ! are permutated to create all equivalent potentials.
  ! That is, if we have a potential for Si-O, both Si-O and O-Si
  ! potentials are created. This is because the energy and
  ! force calculation loops only deal with atom pairs A-B once
  ! (so only A-B or B-A is considered, not both) and if, say,
  ! the loop only finds an O-Si pair, it is important to apply
  ! the Si-O interaction also on that pair.
  ! In some cases, such as with the bond-bending potential affecting
  ! triplets (A-B-C), the interaction is not symmetric for all the atoms.
  ! Therefore it is necessary to also store the original targets of
  ! the potential as specified in the Python interface. These are
  ! to be given in the 'orig_*' lists.
  !
  ! If product potentials are created, all but the first one of the potentials
  ! are created with ``is_multiplier == .true.``. This leads to the potentials
  ! being stored in the global temporary array ``multipliers``. The last potential
  ! of a group should be created with ``is_multiplier = .false.`` and the stored
  ! multipliers are attached to it. The list of multipliers is not cleared automatically,
  ! since usually one creates copies of the same potential with permutated targets and all
  ! of these need the same multipiers.
  ! Instead the multipliers are cleared with a call of :func:`clear_potential_multipliers`.
  !
  ! called from PyInterface: :func:`add_potential`
  !
  ! *n_targets number of targets (interacting bodies)
  ! *n_params number of parameters
  ! *pot_name potential names
  ! *parameters numeric parameters
  ! *cutoff interaction hard cutoff
  ! *smooth_cut interaction soft cutoff
  ! *elements atomic symbols specifying the elements the interaction acts on
  ! *tags tags specifying the atoms the interaction acts on
  ! *indices indices specifying the atoms the interaction acts on
  ! *orig_elements original atomic symbols specifying the elements the interaction acts on
  ! *orig_tags original tags specifying the atoms the interaction acts on
  ! *orig_indices original indices specifying the atoms the interaction acts on
  ! *pot_index index of the potential
  ! *success logical tag specifying if creation of the potential succeeded
  ! *is_multiplier logical tag specifying if this potential should be treated as a multiplier
  subroutine core_add_potential(n_targets,n_params,pot_name,parameters,cutoff,smooth_cut,&
       elements,tags,indices,orig_elements,orig_tags,orig_indices,pot_index,is_multiplier,&
       success)
    implicit none
    integer, intent(in) :: n_targets, n_params, pot_index
    character(len=*), intent(in) :: pot_name
    double precision, intent(in) :: parameters(n_params)
    double precision, intent(in) :: cutoff, smooth_cut
    character(len=label_length), intent(in) :: elements(n_targets)
    integer, intent(in) :: tags(n_targets), indices(n_targets)
    character(len=label_length), intent(in) :: orig_elements(n_targets)
    integer, intent(in) :: orig_tags(n_targets), orig_indices(n_targets)
    logical, intent(in) :: is_multiplier
    logical, intent(out) :: success
    type(potential) :: new_interaction, dummy_multiplier(0)
    type(potential), allocatable :: old_multipliers(:)
    integer :: i

    if(is_multiplier)then
       call create_potential(n_targets,n_params,&
            pot_name,parameters,cutoff,smooth_cut,&
            elements,tags,indices,&
            orig_elements,orig_tags,orig_indices,pot_index,&
            0,dummy_multiplier,&
            new_interaction,success) ! in Potentials.f90

       if(success)then
          if(n_multi > 0)then
             allocate(old_multipliers(n_multi))
             old_multipliers(1:n_multi) = multipliers(1:n_multi)
             deallocate(multipliers)
             allocate(multipliers(n_multi+1))
             multipliers(1:n_multi) = old_multipliers(1:n_multi)
             deallocate(old_multipliers)
          else
             allocate(multipliers(1))
          end if
          n_multi = n_multi + 1
          multipliers(n_multi) = new_interaction
       end if
    else
       call create_potential(n_targets,n_params,&
            pot_name,parameters,cutoff,smooth_cut,&
            elements,tags,indices,&
            orig_elements,orig_tags,orig_indices,pot_index,&
            n_multi,multipliers,&
            new_interaction,success) ! in Potentials.f90

       if(success)then
          n_interactions = n_interactions + 1
          interactions(n_interactions) = new_interaction
       end if
    end if

  end subroutine core_add_potential



! !!!: core_add_bond_order_factor

  ! Creates one additional bond_order_factor in the core.
  ! The routine assumes that adequate memory has been
  ! allocated already using core_allocate_bond_order_factors.
  !
  ! When the bond order parameters in the Python interface are imported
  ! to the Fortran core, the target specifiers (elements)
  ! are permutated to create all equivalent bond order parameters.
  ! That is, if we have parameters for Si-O, both Si-O and O-Si
  ! parameters are created. This is because the energy and
  ! force calculation loops only deal with atom pairs A-B once
  ! (so only A-B or B-A is considered, not both) and if, say,
  ! the loop only finds an O-Si pair, it is important to apply
  ! the Si-O parameters also on that pair.
  ! In some cases, such as with the tersoff factor affecting
  ! triplets (A-B-C), the contribution is not symmetric for all the atoms.
  ! Therefore it is necessary to also store the original targets of
  ! the potential as specified in the Python interface. These are
  ! to be given in the 'orig_elements' lists.
  !
  ! called from PyInterface: :func:`add_bond_order_factor`
  !
  ! *n_targets number of targets (interacting bodies)
  ! *n_params number of parameters
  ! *n_split number of subsets in the list of parameters, should equal n_targets
  ! *bond_name bond order factor names
  ! *parameters numeric parameters
  ! *param_split the numbers of parameters for 1-body, 2-body etc.
  ! *cutoff interaction hard cutoff
  ! *smooth_cut interaction soft cutoff
  ! *elements atomic symbols specifying the elements the interaction acts on
  ! *orig_elements original atomic symbols specifying the elements the interaction acts on
  ! *group_index index denoting the potential to which the factor is connected
  ! *success logical tag specifying if creation of the factor succeeded
  subroutine core_add_bond_order_factor(n_targets,n_params,n_split,bond_name,parameters,param_split,&
       cutoff,smooth_cut,elements,orig_elements,group_index,success)
    implicit none
    integer, intent(in) :: n_targets, n_params, n_split, group_index
    integer, intent(in) :: param_split(n_split)
    character(len=*), intent(in) :: bond_name
    double precision, intent(in) :: parameters(n_params)
    double precision, intent(in) :: cutoff, smooth_cut
    character(len=label_length), intent(in) :: elements(n_targets)
    character(len=label_length), intent(in) :: orig_elements(n_targets)
    logical, intent(out) :: success
    type(bond_order_parameters) :: new_bond_factor

    n_bond_factors = n_bond_factors + 1
    call create_bond_order_factor(n_targets,n_params,n_split,&
         bond_name,parameters,param_split,cutoff,smooth_cut,&
         elements,orig_elements,group_index,&
         new_bond_factor,success)
    bond_factors(n_bond_factors) = new_bond_factor

  end subroutine core_add_bond_order_factor


! !!!: core_assign_bond_order_factor_indices

  ! This routine finds for each atom the potentials for which the
  ! atom is an accepted target at the first position.
  ! First position here means that for instance in an A-B-C triplet.
  ! A is in first position.
  ! Being an accepted target means that the atom has the correct
  ! element.
  !
  ! called from PyInterface: :func:`create_bond_order_factor_list`
  subroutine core_assign_bond_order_factor_indices()
    implicit none
    logical :: affects(n_bond_factors)
    integer :: i, j, k, total
    integer, allocatable :: bond_indices(:)
    double precision :: test_cut, max_cut

    do i = 1, size(atoms)       
       total = 0
       max_cut = 0.d0
       test_cut = 0.d0
       do j = 1, n_bond_factors
          call bond_order_factor_affects_atom(bond_factors(j),atoms(i),affects(j),1) ! in Potentials.f90
          if(affects(j))then
             total = total+1
             test_cut = bond_factors(j)%cutoff
             if(test_cut > max_cut)then
                max_cut = test_cut
             end if
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
       call assign_bond_order_factor_indices(total,atoms(i),bond_indices) ! in Potentials.f90
       call assign_max_bond_order_factor_cutoff(atoms(i),max_cut)

       deallocate(bond_indices)
    end do

  end subroutine core_assign_bond_order_factor_indices


! !!!: core_assign_potential_indices

  ! This routine finds for each atom the potentials for which the
  ! atom is an accepted target at the first position.
  ! First position here means that for instance in an A-B-C triplet.
  ! A is in first position.
  ! Being an accepted target means that the atom has the correct
  ! element, index or tag (one that the potential targets).
  !
  ! called from PyInterface: :func:`create_potential_list`
  subroutine core_assign_potential_indices()
    implicit none
    logical :: affects(n_interactions)
    integer :: i, j, k, total
    integer, allocatable :: pot_indices(:)
    double precision :: test_cut, max_cut
    
    do i = 1, size(atoms)       
       total = 0
       max_cut = 0.d0
       test_cut = 0.d0
       do j = 1, n_interactions
          call potential_affects_atom(interactions(j),atoms(i),affects(j),1) ! in Potentials.f90
          if(affects(j))then
             total = total+1
             test_cut = interactions(j)%cutoff
             if(test_cut > max_cut)then
                max_cut = test_cut
             end if
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
       call assign_potential_indices(total,atoms(i),pot_indices) ! in Potentials.f90
       call assign_max_potential_cutoff(atoms(i),max_cut)

       deallocate(pot_indices)
    end do

  end subroutine core_assign_potential_indices




! !!!: core_calculate_bond_order_gradients

  ! Returns the gradients of bond order factors.
  !
  ! For a factor such as
  !
  ! .. math::
  !
  !      b_i = f(\sum_j c_{ij})
  !
  ! The routine calculates
  !
  ! .. math::
  !
  !     \nabla_\alpha b_i = f'(\sum_j c_{ij}) \nabla_\alpha \sum_j c_{ij}.
  !
  ! By default, the gradients of all factors :math:`i` are calculated with respect
  ! to moving the given atom :math:`\alpha`.
  ! If for_factor is .true., the gradients of the bond factor of the given
  ! atom are calculated with respect to moving all atoms.
  !
  ! *group_index an index denoting the potential to which the factor is connected
  ! *atom_index index of the atom with respect to which the factors are differentiated (:math:`\alpha`), or the atoms whose factor is differentiated (:math:`i`) if for_factor is .true.
  ! *raw_sums precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
  ! *total_gradient the calculated bond order gradients :math:`\nabla_\alpha b_i`
  ! *total_virial the components of the virial due to the bond order gradients
  ! *for_factor a switch for requesting the gradients for a given :math:`i` instead of a given :math:`\alpha`
  subroutine core_calculate_bond_order_gradients(group_index,&
       atom_index,raw_sums,total_gradient,total_virial,for_factor)
    implicit none
    integer, intent(in) :: group_index, atom_index
    double precision, intent(inout) :: total_gradient(:,:), total_virial(6)
    double precision, intent(in) :: raw_sums(:)
    logical, optional, intent(in) :: for_factor
    double precision :: virial(6)
    type(atom) :: atom1, atom2, atom3
    type(atom) :: atom_list(3)
    type(neighbor_list) :: nbors1, nbors2
    type(bond_order_parameters) :: bond_params(2)
    integer, pointer :: bond_indices(:), bond_indices2(:)
    integer :: index1, index2, index3, k1, k2, j, l, n_targets, n_atoms
    double precision, save :: separations(3,2), distances(2), directions(3,2), tmp_grad(3,3,3)
    logical :: is_active, is_in_group, many_bodies_found, separation3_unknown, gradients_for_factor

    n_atoms = size(atoms)
    temp_gradient(1:3,1:n_atoms,1) = 0.d0
    total_gradient = 0.d0
    virial = 0.d0
    total_virial = 0.d0
    bo_scaling = .false.
    
    ! If for_factor is given and true, we return
    ! gradients for the given factor.
    ! Otherwise, we return gradients for all factors
    ! with respect to moving the given atom.
    gradients_for_factor = .false.
    if(present(for_factor))then
       if(for_factor)then
          gradients_for_factor = .true.
       end if
    end if

    ! target atom, given as argument 
    ! (the atom moved in differentiation or 
    ! the atom whose factor is differentiated)
    index1 = atom_index
    atom1 = atoms(atom_index)
    nbors1 = atom1%neighbor_list
    bond_indices => atom1%bond_indices
          
    ! loop over neighbors
    do j = 1, nbors1%n_neighbors
       
       ! neighboring atom
       index2 = nbors1%neighbors(j)          
       atom2 = atoms(index2)
       atom_list(1) = atom1
       atom_list(2) = atom2

       ! calculate atom1-atom2 separation vector
       ! and distance
       call separation_vector(atom1%position, &
            atom2%position, &
            nbors1%pbc_offsets(1:3,j), &
            cell, &
            separations(1:3,1)) ! in Geometry.f90
       distances(1) = .norm.(separations(1:3,1))
       if(distances(1) == 0.d0)then
          directions(1:3,1) = (/ 0.d0, 0.d0, 0.d0 /)
       else
          directions(1:3,1) = separations(1:3,1) / distances(1)
       end if
          
       many_bodies_found = .false.

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! apply 2-body bond order factors !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       if(distances(1) < atom1%max_bond_radius)then

          ! loop over bond order factors affecting atom1
          do k1 = 1, size(bond_indices)
          
             bond_params(1) = bond_factors(bond_indices(k1))

             ! filter the bond indices before applying by checking:
             ! number of targets,
             ! is atom2 affected by the factor,
             ! is the factor in the correct group
             call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group) ! in Potentials.f90
             if( is_in_group .and. bond_params(1)%cutoff > distances(1) )then
                call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2) ! in Potentials.f90
                if( is_active )then
                   call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                        n_targets) ! in Potentials.f90
                   if( n_targets == 2 )then
                      
                      ! evaluate the atom1-atom2 term in the bond gradient sum
                      call evaluate_bond_order_gradient(2,&
                           separations(1:3,1),&
                           distances(1),&
                           bond_params(1),&
                           tmp_grad(1:3,1:2,1:2),&
                           atom_list(1:2)) ! in Potentials.f90
                      
                      if(gradients_for_factor)then
                         ! store the gradients of the atom1 term (the target atom)
                         ! with respect to moving atom1 and atom2
                         temp_gradient(1:3,index1,1) = temp_gradient(1:3,index1,1) + tmp_grad(1:3,1,1)
                         temp_gradient(1:3,index2,1) = temp_gradient(1:3,index2,1) + tmp_grad(1:3,1,2)
                         bo_scaling(index1) = .true.
                         bo_scaling(index2) = .true.

                         ! virial: r_a*f_b (xx, yy, zz, yz, xz, xy)
                         virial(1) = virial(1) + separations(1,1) * tmp_grad(1,1,2)
                         virial(2) = virial(2) + separations(2,1) * tmp_grad(2,1,2)
                         virial(3) = virial(3) + separations(3,1) * tmp_grad(3,1,2)
                         virial(4) = virial(4) + separations(2,1) * tmp_grad(3,1,2)
                         virial(5) = virial(5) + separations(1,1) * tmp_grad(3,1,2)
                         virial(6) = virial(6) + separations(1,1) * tmp_grad(2,1,2)
                      else
                         ! store the gradients of the atom1 and atom2 terms
                         ! with respect to movind atom1 (the target atom)
                         temp_gradient(1:3,index1,1) = temp_gradient(1:3,index1,1) + tmp_grad(1:3,1,1)
                         temp_gradient(1:3,index2,1) = temp_gradient(1:3,index2,1) + tmp_grad(1:3,2,1)
                         bo_scaling(index1) = .true.
                         bo_scaling(index2) = .true.
                         
                      end if
                      
                   else if( n_targets > 2 )then
                      
                      ! If the number of targets is greater than 2,
                      ! we have found a many-body bond order factor.
                      ! Make a note that we must also evaluate the many-body factors.
                      many_bodies_found = .true.
                      
                   end if ! n_targets == 2
                end if ! is_active
             end if ! is_in_group
             
          end do ! k1 = 1, size(bond_indices)
       end if ! distance < max_cut

       ! Only do the 3-body loop if we found many-body factors 
       ! during 2-body evaluation.
       ! 
       ! In the 3-body loop, we search the neighbors of both atom1
       ! and atom2 to find the triplets 
       ! atom1-atom2-atom3 and
       ! atom2-atom1-atom3
       ! These are considered to be different, since the middle atom of
       ! the triplet is different (atom2 vs. atom1).       
       if(many_bodies_found)then
             
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! apply 3-body bond order factors !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! neighbors of atom2
          nbors2 = atom2%neighbor_list
          
          ! First we try to find ordered triplets atom2 -- atom1 -- atom3
          ! Therefore we need separations a2--a1 and a1--a3.
          separations(1:3,1) = -separations(1:3,1)

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

                ! third atom of the triplet
                atom3 = atoms(index3)
                ! atom3 is new so we don't know the separation from atom1
                separation3_unknown = .true.
                ! The list of atoms is passed to bond factor evaluation routine
                ! for further filtering.
                ! This is triplet atom2 - atom1 - atom3, since we loop over
                ! neighbors of atom1.
                atom_list = (/ atom2, atom1, atom3 /)

                ! Since a triplet contains both pairs 
                ! atom1-atom2 and
                ! atom1-atom3,
                ! We will need bond order parameters for both.
                ! Therefore we must loop over the bond factors of atom1
                ! twice to find the suitable parameters.
                
                ! search for the first bond params containing the parameters for atom1-atom2
                do k1 = 1, size(bond_indices)
                   
                   bond_params(1) = bond_factors(bond_indices(k1))

                   ! filter the parameters by:
                   ! number of targets, atom2 and atom3 being targets, group index                   
                   call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group) ! in Potentials.f90
                   
                   if( is_in_group )then
                      call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2) ! in Potentials.f90
                      if( is_active )then
                         call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                              n_targets) ! in Potentials.f90
                         if( n_targets == 3 )then
                            call bond_order_factor_affects_atom(bond_params(1),&
                                 atom3,is_active,3) ! in Potentials.f90
                            if( is_active )then
                                                  
                ! search for the second bond params containing the parameters for atom1-atom3
                do k2 = 1, size(bond_indices)

                   bond_params(2) = bond_factors(bond_indices(k2))
                   call bond_order_factor_is_in_group(bond_params(2),group_index,is_in_group) ! in Potentials.f90
                   if( is_in_group )then 
                      call bond_order_factor_affects_atom(bond_params(2),atom2,is_active,2) ! in Potentials.f90
                      if( is_active )then
                         call get_number_of_targets_of_bond_order_factor_index(bond_params(2)%type_index,&
                              n_targets) ! in Potentials.f90
                         if( n_targets == 3 )then
                            call bond_order_factor_affects_atom(bond_params(2),&
                                 atom3,is_active,3) ! in Potentials.f90
                            if( is_active )then

                               ! When we loop over the bond factors
                               ! we may need the atom1-atom3 distances
                               ! repeatedly. We only calculate it the first
                               ! time.
                               ! This could be done once before the loop as well,
                               ! as is done for 2-body terms,
                               ! but then we could end up calculating many
                               ! separations and distances for atom pairs
                               ! that do not have an interaction.
                               if( separation3_unknown )then
                                  call separation_vector(atom1%position, &
                                       atom3%position, &
                                       nbors1%pbc_offsets(1:3,j), &
                                       cell, &
                                       separations(1:3,2)) ! in Geometry.f90
                                  separation3_unknown = .false.
                                  distances(2) = .norm.(separations(1:3,2))
                                  if(distances(2) == 0.d0)then
                                     directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                  else
                                     directions(1:3,2) = separations(1:3,2) / distances(2)
                                  end if
                               end if
                               
                               ! Evaluate the atom2-atom1-atom3 triplet contribution
                               ! in the bond gradient sum.
                               ! The gradient is stored in tmp_grad whose indices are
                               ! (xyz, factor index, moved atom index)
                               call evaluate_bond_order_gradient(3,&
                                    separations(1:3,1:2),&
                                    distances(1:2),&
                                    bond_params(1:2),&
                                    tmp_grad(1:3,1:3,1:3),&
                                    atom_list) ! in Potentials.f90

                               if(gradients_for_factor)then
                                  ! store the gradients of the atom1 term (the target atom)
                                  ! with respect to moving atom1, atom2, atom3
                                  ! Note that here atom1 is the middle atom, so we take index 2
                                  ! in the second column of tmp_grad.
                                  ! The triplet is atom2-atom1-atom3, so index2 gets the first
                                  ! entry in the third column of tmp_grad.
                                  temp_gradient(1:3,index2,1) = temp_gradient(1:3,index2,1) + tmp_grad(1:3,2,1)
                                  temp_gradient(1:3,index1,1) = temp_gradient(1:3,index1,1) + tmp_grad(1:3,2,2)
                                  temp_gradient(1:3,index3,1) = temp_gradient(1:3,index3,1) + tmp_grad(1:3,2,3)
                                  bo_scaling(index1) = .true.
                                  bo_scaling(index2) = .true.
                                  bo_scaling(index3) = .true.

                                  ! virial: r_a*f_b (xx, yy, zz, yz, xz, xy) 
                                  ! r from atom1 to other atoms - now atom1 is center atom
                                  virial(1) = virial(1) - separations(1,1) * tmp_grad(1,2,1) & ! r_12*f2
                                                        + separations(1,2) * tmp_grad(1,2,3)   ! r_13*f3
                                  virial(2) = virial(2) - separations(2,1) * tmp_grad(2,2,1) & ! r_12*f2
                                                        + separations(2,2) * tmp_grad(2,2,3)   ! r_13*f3
                                  virial(3) = virial(3) - separations(3,1) * tmp_grad(3,2,1) & ! r_12*f2
                                                        + separations(3,2) * tmp_grad(3,2,3)   ! r_13*f3
                                  virial(4) = virial(4) - separations(2,1) * tmp_grad(3,2,1) & ! r_12*f2
                                                        + separations(2,2) * tmp_grad(3,2,3)   ! r_13*f3
                                  virial(5) = virial(5) - separations(1,1) * tmp_grad(3,2,1) & ! r_12*f2
                                                        + separations(1,2) * tmp_grad(3,2,3)   ! r_13*f3
                                  virial(6) = virial(6) - separations(1,1) * tmp_grad(2,2,1) & ! r_12*f2
                                                        + separations(1,2) * tmp_grad(2,2,3)   ! r_13*f3

                               else
                                  ! store the gradients of the atom1, atom2, and atom3 terms
                                  ! with respect to moving atom1 (the target atom)
                                  ! Note that here atom1 is the middle atom, so we take index 2
                                  ! in the third column of tmp_grad.
                                  ! The triplet is atom2-atom1-atom3, so index2 gets the first
                                  ! entry in the second column of tmp_grad.
                                  temp_gradient(1:3,index2,1) = temp_gradient(1:3,index2,1) + tmp_grad(1:3,1,2)
                                  temp_gradient(1:3,index1,1) = temp_gradient(1:3,index1,1) + tmp_grad(1:3,2,2)
                                  temp_gradient(1:3,index3,1) = temp_gradient(1:3,index3,1) + tmp_grad(1:3,3,2)
                                  bo_scaling(index1) = .true.
                                  bo_scaling(index2) = .true.
                                  bo_scaling(index3) = .true.

                               end if

                            end if ! is_active
                         end if ! n_targets == 3
                      end if ! is_active
                   end if ! is_in_group
                end do ! k2
                         
                            end if ! is_active
                         end if ! n_targets == 3
                      end if ! is_active
                   end if ! is_in_group
                end do ! k1

             end if ! index3 > index2
          end do ! l = 1, nbors1%n_neighbors


          ! Next we try to find ordered triplets atom1 -- atom2 -- atom3
          ! Therefore we need separations a1--a2 and a2--a3.
          separations(1:3,1) = -separations(1:3,1)

          ! loop over neighbors of atom 2
          do l = 1, nbors2%n_neighbors
             index3 = nbors2%neighbors(l)

             ! In the similar loop above, we filter by index3 > index2 to
             ! prevent double counting.
             ! Here we are looping over neighbors of atom2, so the double counting
             ! problem of looping over atom1 neighbors twice is not present.
             ! Also, unlike in the calculation of the bond factors or forces,
             ! we cannot find the atom1-atom2-atom3 triplet by starting from
             ! atom3, so the double counting is not happening that way either.
             ! (This is because we are not looping over all atoms to get all
             ! gradients, but just calculate the gradients of a certain factor
             ! or gradients with respect to moving a certain particle!)
             ! What must be filtered is the case atom1-atom2-atom1, which will
             ! happen when atom2 finds atom1 as its neighbor.
             if(index3 /= index1)then

                ! Third atom of the triplet.
                atom3 = atoms(index3)
                separation3_unknown = .true.
                atom_list = (/ atom1, atom2, atom3 /)
                
                ! Since a triplet contains both pairs 
                ! atom1-atom2 and
                ! atom2-atom3,
                ! We will need bond order parameters for both.
                ! Therefore we must loop over the bond factors of atom1
                ! twice to find the suitable parameters.

                ! search for the first bond params containing the parameters for atom2-atom1
                do k1 = 1, size(bond_indices)
                   
                   bond_params(1) = bond_factors(bond_indices(k1))

                   ! filter the parameters by:
                   ! number of targets, atom2 and atom3 being targets, group index   
                   call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group) ! in Potentials.f90
                   if( is_in_group )then 
                      call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2) ! in Potentials.f90
                      if( is_active )then
                         call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,n_targets) ! in Potentials.f90
                         if( n_targets == 3 )then
                            call bond_order_factor_affects_atom(bond_params(1),atom3,is_active,3) ! in Potentials.f90
                            if( is_active )then
                         
                 ! search for the second bond params containing the parameters for atom2-atom3
                 do k2 = 1, size(bond_indices)
                            
                    bond_params(2) = bond_factors(bond_indices(k2))
                    call bond_order_factor_is_in_group(bond_params(2),group_index,is_in_group) ! in Potentials.f90
                    if( is_in_group )then 
                       call bond_order_factor_affects_atom(bond_params(2),atom2,is_active,2) ! in Potentials.f90
                       if( is_active )then 
                          call get_number_of_targets_of_bond_order_factor_index(bond_params(2)%type_index,&
                               n_targets) ! in Potentials.f90
                          if( n_targets == 3 )then
                             call bond_order_factor_affects_atom(bond_params(2),atom3,is_active,3) ! in Potentials.f90
                             if( is_active )then
                                  
                                if( separation3_unknown )then
                                   call separation_vector(atom2%position, &
                                        atom3%position, &
                                        nbors2%pbc_offsets(1:3,l), &
                                        cell, &
                                        separations(1:3,2)) ! in Geometry.f90
                                   separation3_unknown = .false.
                                   distances(2) = .norm.(separations(1:3,2))
                                   if(distances(2) == 0.d0)then
                                      directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                   else
                                      directions(1:3,2) = separations(1:3,2) / distances(2)
                                   end if
                                end if
                                
                               ! Evaluate the atom1-atom2-atom3 triplet contribution
                               ! in the bond gradient sum.
                               ! The gradient is stored in tmp_grad whose indices are
                               ! (xyz, factor index, moved atom index)
                                call evaluate_bond_order_gradient(3,&
                                     separations(1:3,1:2),&
                                     distances(1:2),&
                                     bond_params(1:2),&
                                     tmp_grad(1:3,1:3,1:3),&
                                     atom_list) ! in Potentials.f90

                               if(gradients_for_factor)then
                                  ! store the gradients of the atom1 term (the target atom)
                                  ! with respect to moving atom1, atom2, atom3
                                  ! Note that here atom1 is the first atom, so we take index 1
                                  ! in the second column of tmp_grad.
                                  ! The triplet is atom1-atom2-atom3, so index1 gets the first
                                  ! entry in the third column of tmp_grad.
                                  temp_gradient(1:3,index1,1) = temp_gradient(1:3,index1,1) + tmp_grad(1:3,1,1)
                                  temp_gradient(1:3,index2,1) = temp_gradient(1:3,index2,1) + tmp_grad(1:3,1,2)
                                  temp_gradient(1:3,index3,1) = temp_gradient(1:3,index3,1) + tmp_grad(1:3,1,3)
                                  bo_scaling(index1) = .true.
                                  bo_scaling(index2) = .true.
                                  bo_scaling(index3) = .true.

                                  ! virial: r_a*f_b (xx, yy, zz, yz, xz, xy) 
                                  ! r from atom1 to other atoms - now atom1 is first atom
                                  virial(1) = virial(1) + separations(1,1) * tmp_grad(1,2,1) & ! r_12*f2
                                       + (separations(1,1)+separations(1,2)) * tmp_grad(1,2,3)   ! r_13*f3
                                  virial(2) = virial(2) + separations(2,1) * tmp_grad(2,2,1) & ! r_12*f2
                                       + (separations(2,1)+separations(2,2)) * tmp_grad(2,2,3)   ! r_13*f3
                                  virial(3) = virial(3) + separations(3,1) * tmp_grad(3,2,1) & ! r_12*f2
                                       + (separations(3,1)+separations(3,2)) * tmp_grad(3,2,3)   ! r_13*f3
                                  virial(4) = virial(4) + separations(2,1) * tmp_grad(3,2,1) & ! r_12*f2
                                       + (separations(2,1)+separations(2,2)) * tmp_grad(3,2,3)   ! r_13*f3
                                  virial(5) = virial(5) + separations(1,1) * tmp_grad(3,2,1) & ! r_12*f2
                                       + (separations(1,1)+separations(1,2)) * tmp_grad(3,2,3)   ! r_13*f3
                                  virial(6) = virial(6) + separations(1,1) * tmp_grad(2,2,1) & ! r_12*f2
                                       + (separations(1,1)+separations(1,2)) * tmp_grad(2,2,3)   ! r_13*f3
                               else
                                  ! store the gradients of the atom1, atom2, and atom3 terms
                                  ! with respect to movind atom1 (the target atom)
                                  ! Note that here atom1 is the first atom, so we take index 1
                                  ! in the third column of tmp_grad.
                                  ! The triplet is atom1-atom2-atom3, so index1 gets the first
                                  ! entry in the second column of tmp_grad.
                                  temp_gradient(1:3,index1,1) = temp_gradient(1:3,index1,1) + tmp_grad(1:3,1,1)
                                  temp_gradient(1:3,index2,1) = temp_gradient(1:3,index2,1) + tmp_grad(1:3,2,1)
                                  temp_gradient(1:3,index3,1) = temp_gradient(1:3,index3,1) + tmp_grad(1:3,3,1)
                                  bo_scaling(index1) = .true.
                                  bo_scaling(index2) = .true.
                                  bo_scaling(index3) = .true.

                               end if


                            end if ! is_active
                         end if ! n_targets == 3
                      end if ! is_active
                   end if ! is_in_group
                end do ! k2
                         
                            end if ! is_active
                         end if ! n_targets == 3
                      end if ! is_active
                   end if ! is_in_group
                end do ! k1
                
             end if ! index3 /= index1
          end do ! l = 1, nbors2%n_neighbors
             
       end if ! many_bodies_found
    end do ! j = nbors1%n_neighbors

    ! Above we calculated the gradient for the bond order sums, e.g.,
    ! for b_i = f( \sum_ij c_ij ) we now have
    ! \nabla_a \sum_c_ij.
    ! To get \nabla_a b_i, we need to evaluate
    ! \nabla_a b_i = f'( \sum c_ij ) * \nabla_a \sum_c_ij.

    if(gradients_for_factor)then
       call core_post_process_bond_order_gradients_of_factor(group_index,index1,raw_sums(index1),&
            temp_gradient(1:3,1:n_atoms,1),total_gradient,virial,total_virial)
    else
       call core_post_process_bond_order_gradients(group_index,raw_sums,&
            temp_gradient(1:3,1:n_atoms,1),total_gradient)
    end if

  end subroutine core_calculate_bond_order_gradients




! !!!: core_calculate_pair_bond_order_gradients

  ! Returns the gradients of a pair bond order factor.
  !
  ! For a factor such as
  !
  ! .. math::
  !
  !      b_{ij} = f(\sum_k c_{ijk})
  !
  ! The routine calculates
  !
  ! .. math::
  !
  !     \nabla_\alpha b_{ij} = f'(\sum_k c_{ijk}) \nabla_\alpha \sum_k c_{ijk}.
  !
  ! By default, the gradients the factor :math:`ij` is calculated with respect
  ! to moving all atoms :math:`\alpha`.
  !
  ! *group_index an index denoting the potential to which the factor is connected
  ! *atom_index index of the atom with respect to which the factors are differentiated (:math:`\alpha`)
  ! *raw_sums precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
  ! *total_gradient the calculated bond order gradients :math:`\nabla_\alpha b_i`
  ! *total_virial the components of the virial due to the bond order gradient
  subroutine core_calculate_pair_bond_order_gradients(atom_pair,separation,distance,direction,group_index,&
       raw_sums,total_gradient,total_virial)
    implicit none
    integer, intent(in) :: atom_pair(2), group_index
    double precision, intent(in) :: separation(3), distance, direction(3)
    double precision, intent(inout) :: total_gradient(:,:,:), total_virial(6,2)
    double precision, intent(in) :: raw_sums(2)
    double precision :: virial(6,2)
    type(atom) :: atom1, atom2, atom3
    type(atom) :: atom_list(3)
    type(neighbor_list) :: nbors1, nbors2
    type(bond_order_parameters) :: bond_params(1)
    integer, pointer :: bond_indices(:), bond_indices2(:)
    integer :: index1, index2, index3, k1, j, l, n_targets, offset(3), n_atoms
    double precision :: separations(3,2), distances(2), directions(3,2), tmp_grad(3,3)
    logical :: is_active, is_in_group

    n_atoms = size(atoms)
    temp_gradient = 0.d0
    total_gradient = 0.d0
    virial = 0.d0
    total_virial = 0.d0
    
  
    ! target atom 1
    index1 = atom_pair(1)
    atom1 = atoms(index1)
    nbors1 = atom1%neighbor_list
    bond_indices => atom1%bond_indices

    ! target atom 2
    index2 = atom_pair(2)
    atom2 = atoms(index2)
    nbors2 = atom2%neighbor_list
    bond_indices2 => atom2%bond_indices


    ! In the 3-body loop, we search the neighbors of both atom1
    ! and atom2 to find the triplets 
    ! atom1-atom2-atom3 and
    ! atom2-atom1-atom3
    ! These are considered to be different, since the middle atom of
    ! the triplet is different (atom2 vs. atom1).
    
    ! First we try to find ordered triplets atom2 -- atom1 -- atom3
    ! Therefore we need separations a2--a1 and a1--a3.

    ! atom2 - atom1 vector
    separations(1:3,1) = -separation(1:3)
    distances(1) = distance
    directions(1:3,1) = -direction(1:3)
        

    ! loop over neighbors of atom 1
    do j = 1, nbors1%n_neighbors
       
       ! neighboring atom
       index3 = nbors1%neighbors(j)
       offset(1:3) = nbors1%pbc_offsets(1:3,j)

       ! To prevent the case index2 - index1 - index2
       if(index2 /= index3)then

          atom3 = atoms(index3)
          atom_list = (/ atom2, atom1, atom3 /)

          ! calculate atom1-atom3 separation vector
          ! and distance
          call separation_vector(atom1%position, &
               atom3%position, &
               offset(1:3), &
               cell, &
               separations(1:3,2)) ! in Geometry.f90
          distances(2) = .norm.(separations(1:3,2))
          if(distances(2) == 0.d0)then
             directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
          else
             directions(1:3,2) = separations(1:3,2) / distances(2)
          end if
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! apply 3-body bond order factors !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if(distances(2) < atom1%max_bond_radius)then

             ! search for the first bond params containing the parameters for (atom2-atom1) -- atom3
             do k1 = 1, size(bond_indices)
                
                bond_params(1) = bond_factors(bond_indices(k1))
                
                ! filter the parameters by:
                ! number of targets, atom2 and atom3 being targets, group index                   
                call bond_order_factor_is_in_group(bond_params(1),&
                     group_index,is_in_group) ! in Potentials.f90
                            
                if( is_in_group )then
                   call bond_order_factor_affects_atom(bond_params(1),&
                        atom2,is_active,2) ! in Potentials.f90
                   if( is_active )then
                      call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                           n_targets) ! in Potentials.f90
                      if( n_targets == 3 )then
                         call bond_order_factor_affects_atom(bond_params(1),&
                              atom3,is_active,3) ! in Potentials.f90
                         if( is_active )then
                            
                            
                            ! Evaluate the atom2-atom1-atom3 triplet contribution
                            ! in the bond gradient sum.
                            ! The gradient is stored in tmp_grad whose indices are
                            ! (xyz, factor index, moved atom index)
                            call evaluate_pair_bond_order_gradient(3,&
                                 separations(1:3,1:2),&
                                 distances(1:2),&
                                 bond_params(1),&
                                 tmp_grad(1:3,1:3),&
                                 atom_list) ! in Potentials.f90
                            
                            ! store the gradients of the pair bond factor
                            ! with respect to moving atom1, atom2, atom3
                            ! The triplet is atom2-atom1-atom3, so index2 gets the first
                            ! entry in the second column of tmp_grad.
                            temp_gradient(1:3,index2,1) = temp_gradient(1:3,index2,1) + tmp_grad(1:3,1)
                            temp_gradient(1:3,index1,1) = temp_gradient(1:3,index1,1) + tmp_grad(1:3,2)
                            temp_gradient(1:3,index3,1) = temp_gradient(1:3,index3,1) + tmp_grad(1:3,3)
                            
                            ! virial: r_a*f_b (xx, yy, zz, yz, xz, xy) 
                            ! r from atom1 to other atoms - now atom1 is center atom
                            virial(1,1) = virial(1,1) - separations(1,1) * tmp_grad(1,1) & ! r_12*f2
                                 + separations(1,2) * tmp_grad(1,3)   ! r_13*f3
                            virial(2,1) = virial(2,1) - separations(2,1) * tmp_grad(2,1) & ! r_12*f2
                                 + separations(2,2) * tmp_grad(2,3)   ! r_13*f3
                            virial(3,1) = virial(3,1) - separations(3,1) * tmp_grad(3,1) & ! r_12*f2
                                 + separations(3,2) * tmp_grad(3,3)   ! r_13*f3
                            virial(4,1) = virial(4,1) - separations(2,1) * tmp_grad(3,1) & ! r_12*f2
                                 + separations(2,2) * tmp_grad(3,3)   ! r_13*f3
                            virial(5,1) = virial(5,1) - separations(1,1) * tmp_grad(3,1) & ! r_12*f2
                                 + separations(1,2) * tmp_grad(3,3)   ! r_13*f3
                            virial(6,1) = virial(6,1) - separations(1,1) * tmp_grad(2,1) & ! r_12*f2
                                 + separations(1,2) * tmp_grad(2,3)   ! r_13*f3
                            
                         end if ! is_active
                      end if ! n_targets == 3
                   end if ! is_active
                end if ! is_in_group
             end do ! k1
          end if ! distance < max_cut
       end if ! index3 /= index2

    end do ! j = 1, nbors1%n_neighbors


    ! Next we try to find ordered triplets atom1 -- atom2 -- atom3
    ! Therefore we need separations a1--a2 and a2--a3.
    separations(1:3,1) = -separations(1:3,1)
    directions(1:3,1) = -directions(1:3,1)

    ! loop over neighbors of atom 2
    do j = 1, nbors2%n_neighbors
       
       ! neighboring atom
       index3 = nbors2%neighbors(j)
       offset(1:3) = nbors2%pbc_offsets(1:3,j)

       ! To prevent the case index1 - index2 - index1
       if(index1 /= index3)then

          atom3 = atoms(index3)
          atom_list = (/ atom1, atom2, atom3 /)

          ! calculate atom2-atom3 separation vector
          ! and distance
          call separation_vector(atom2%position, &
               atom3%position, &
               offset(1:3), &
               cell, &
               separations(1:3,2)) ! in Geometry.f90
          distances(2) = .norm.(separations(1:3,2))
          if(distances(2) == 0.d0)then
             directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
          else
             directions(1:3,2) = separations(1:3,2) / distances(2)
          end if
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! apply 3-body bond order factors !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if(distances(2) < atom2%max_bond_radius)then

             ! search for the bond params containing the parameters for (atom1-atom2) -- atom3
             do k1 = 1, size(bond_indices2)
                
                bond_params(1) = bond_factors(bond_indices2(k1))
                
                ! filter the parameters by:
                ! number of targets, atom2 and atom3 being targets, group index                   
                call bond_order_factor_is_in_group(bond_params(1),&
                     group_index,is_in_group) ! in Potentials.f90
                
                if( is_in_group )then
                   call bond_order_factor_affects_atom(bond_params(1),&
                        atom1,is_active,2) ! in Potentials.f90
                   if( is_active )then
                      call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                           n_targets) ! in Potentials.f90
                      if( n_targets == 3 )then
                         call bond_order_factor_affects_atom(bond_params(1),&
                              atom3,is_active,3) ! in Potentials.f90
                         if( is_active )then            
                            
                            ! Evaluate the atom2-atom1-atom3 triplet contribution
                            ! in the bond gradient sum.
                            ! The gradient is stored in tmp_grad whose indices are
                            ! (xyz, factor index, moved atom index)
                            call evaluate_pair_bond_order_gradient(3,&
                                 separations(1:3,1:2),&
                                 distances(1:2),&
                                 bond_params(1),&
                                 tmp_grad(1:3,1:3),&
                                 atom_list) ! in Potentials.f90
                            
                            ! store the gradients of the pair bond factor
                            ! with respect to moving atom1, atom2, atom3
                            ! The triplet is atom1-atom2-atom3, so index1 gets the first
                            ! entry in the second column of tmp_grad.
                            temp_gradient(1:3,index1,2) = temp_gradient(1:3,index1,2) + tmp_grad(1:3,1)
                            temp_gradient(1:3,index2,2) = temp_gradient(1:3,index2,2) + tmp_grad(1:3,2)
                            temp_gradient(1:3,index3,2) = temp_gradient(1:3,index3,2) + tmp_grad(1:3,3)
                            
                            ! virial: r_a*f_b (xx, yy, zz, yz, xz, xy) 
                            ! r from atom1 to other atoms - now atom1 is center atom
                            virial(1,2) = virial(1,2) - separations(1,1) * tmp_grad(1,1) & ! r_12*f2
                                 + separations(1,2) * tmp_grad(1,3)   ! r_13*f3
                            virial(2,2) = virial(2,2) - separations(2,1) * tmp_grad(2,1) & ! r_12*f2
                                 + separations(2,2) * tmp_grad(2,3)   ! r_13*f3
                            virial(3,2) = virial(3,2) - separations(3,1) * tmp_grad(3,1) & ! r_12*f2
                                 + separations(3,2) * tmp_grad(3,3)   ! r_13*f3
                            virial(4,2) = virial(4,2) - separations(2,1) * tmp_grad(3,1) & ! r_12*f2
                                 + separations(2,2) * tmp_grad(3,3)   ! r_13*f3
                            virial(5,2) = virial(5,2) - separations(1,1) * tmp_grad(3,1) & ! r_12*f2
                                 + separations(1,2) * tmp_grad(3,3)   ! r_13*f3
                            virial(6,2) = virial(6,2) - separations(1,1) * tmp_grad(2,1) & ! r_12*f2
                                 + separations(1,2) * tmp_grad(2,3)   ! r_13*f3
                            
                         end if ! is_active
                      end if ! n_targets == 3
                   end if ! is_active
                end if ! is_in_group
             end do ! k1
          end if ! distance < max_cut
       end if ! index3 /= index1

    end do ! j = 1, nbors2%n_neighbors






    ! Above we calculated the gradient for the bond order sums, e.g.,
    ! for b_ij = f( \sum_k c_ijk ) we now have
    ! \nabla_a \sum_k c_ijk.
    ! To get \nabla_a b_ij, we need to evaluate
    ! \nabla_a b_ij = f'( \sum c_ijk ) * \nabla_a \sum_c_ijk.

    call core_post_process_pair_bond_order_gradients(group_index,atom1,raw_sums(1),&
         temp_gradient(1:3,1:n_atoms,1),total_gradient(1:3,1:n_atoms,1),virial(1:6,1),total_virial(1:6,1))
    call core_post_process_pair_bond_order_gradients(group_index,atom2,raw_sums(2),&
         temp_gradient(1:3,1:n_atoms,2),total_gradient(1:3,1:n_atoms,2),virial(1:6,2),total_virial(1:6,2))


  end subroutine core_calculate_pair_bond_order_gradients






! !!!: core_calculate_bond_order_gradients_of_factor

  ! Returns the gradients of one bond order factor with respect to
  ! moving all atoms.
  ! 
  ! This calls :func:`core_calculate_bond_order_gradients` with for_factor = .true.
  !
  ! For a factor such as
  !
  ! .. math::
  !
  !      b_i = f(\sum_j c_{ij})
  !
  ! The routine calculates
  !
  ! .. math::
  !
  !     \nabla_\alpha b_i = f'(\sum_j c_{ij}) \nabla_\alpha \sum_j c_{ij}.
  !
  ! The gradients of the bond factor of the given
  ! atom :math:`i` are calculated with respect to moving all atoms :math:`\alpha`.
  !
  ! *group_index an index denoting the potential to which the factor is connected
  ! *atom_index index of the atom whose factor is differentiated (:math:`i`)
  ! *raw_sums precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
  ! *total_gradient the calculated bond order gradients :math:`\nabla_\alpha b_i`
  ! *total_virial the components of the virial due to the bond order gradient
  subroutine core_calculate_bond_order_gradients_of_factor(group_index,&
       atom_index,raw_sums,total_gradient,total_virial)
    implicit none
    integer, intent(in) :: group_index, atom_index
    double precision, intent(inout) :: total_gradient(:,:),total_virial(6)
    double precision, intent(in) :: raw_sums(:)

    call core_calculate_bond_order_gradients(group_index,&
       atom_index,raw_sums,total_gradient,total_virial,.true.)

  end subroutine core_calculate_bond_order_gradients_of_factor



! !!!: core_calculate_bond_order_factors

  ! Calculates the bond order sums of all atoms for the given group.
  !
  ! For a factor such as
  !
  ! .. math::
  !
  !      b_i = f(\sum_j c_{ij})
  !
  ! The routine calculates
  !
  ! .. math::
  !
  !      \sum_j c_{ij}.
  !
  ! The full bond order factor is then obtained by applying the
  ! scaling function :math:`f`. This is done with
  ! :func:`core_post_process_bond_order_factors`.
  !
  ! *group_index an index denoting the potential to which the factor is connected
  ! *total_bond_orders the calculated bond order sums
  subroutine core_calculate_bond_order_factors(group_index,total_bond_orders)
    implicit none
    integer, intent(in) :: group_index
    double precision, intent(inout) :: total_bond_orders(:)
    integer :: index1, index2, index3, k1, k2, j, l, n_targets
    double precision :: separations(3,2), distances(2), directions(3,2)
    double precision :: tmp_factor(3)
    type(atom) :: atom1, atom2, atom3
    type(atom) :: atom_list(3)
    type(neighbor_list) :: nbors1, nbors2
    type(bond_order_parameters) :: bond_params(2)
    integer, pointer :: bond_indices(:), bond_indices2(:)
    logical :: is_active, is_in_group, many_bodies_found, separation3_unknown
    integer :: offset(3), n_atoms

    n_atoms = size(atoms)
    bo_temp = 0.d0
    total_bond_orders = 0.d0

    ! loop over atoms
    do index1 = 1, n_atoms
       
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
             offset(1:3) = nbors1%pbc_offsets(1:3,j)

             ! Since we loop over the neighbors of all atoms, we will find the pair
             ! atom1-atom2 = atom2-atom1 twice.
             ! To prevent the double counting, we filter by index2 > index1.
             if(pick(index1,index2,offset))then

                atom2 = atoms(index2)
                atom_list(1) = atom1
                atom_list(2) = atom2

                ! calculate atom1-atom2 separation vector
                ! and distance
                call separation_vector(atom1%position, &
                     atom2%position, &
                     nbors1%pbc_offsets(1:3,j), &
                     cell, &
                     separations(1:3,1)) ! in Geometry.f90
                distances(1) = .norm.(separations(1:3,1))
                if(distances(1) == 0.d0)then
                   directions(1:3,1) = (/ 0.d0, 0.d0, 0.d0 /)
                else
                   directions(1:3,1) = separations(1:3,1) / distances(1)
                end if

                many_bodies_found = .false.

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! apply 2-body bond order factors !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if(distances(1) < atom1%max_bond_radius)then

                   ! loop over bond order factors affecting atom1
                   do k1 = 1, size(bond_indices)
                      
                      bond_params(1) = bond_factors(bond_indices(k1))
                      
                      ! filter the bond indices before applying by checking:
                      ! number of targets,
                      ! is atom2 affected by the factor,
                      ! is the factor in the correct group
                      if( bond_params(1)%n_level == 1)then
                         call bond_order_factor_is_in_group(bond_params(1),group_index,is_in_group) ! in Potentials.f90
                         if( is_in_group .and. bond_params(1)%cutoff > distances(1) )then
                            call bond_order_factor_affects_atom(bond_params(1),atom2,is_active,2) ! in Potentials.f90
                            if( is_active )then 
                               call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                                    n_targets) ! in Potentials.f90
                               if( n_targets == 2 )then
                                  
                                  ! evaluate the atom1-atom2 term in the bond order factor sum
                                  call evaluate_bond_order_factor(2,separations(1:3,1),&
                                       distances(1),&
                                       bond_params(1),&
                                       tmp_factor(1:2),&
                                       atom_list(1:2)) ! in Potentials.f90
                                  
                                  bo_temp(index1) = bo_temp(index1) + tmp_factor(1)
                                  bo_temp(index2) = bo_temp(index2) + tmp_factor(2)
                                  
                               else if( n_targets > 2 )then
                                  
                                  ! If the number of targets is greater than 2,
                                  ! we have found a many-body bond order factor.
                                  ! Make a note that we must also evaluate the many-body factors.
                                  many_bodies_found = .true.
                                  
                               end if ! n_targets == 2
                            end if ! is_active
                         end if ! is_in_group
                      end if ! n_level == 1

                   end do ! k1 = 1, size(bond_indices)
                end if ! distance < max_cut

                ! Only do the 3-body loop if we found many-body factors 
                ! during 2-body evaluation.
                ! 
                ! In the 3-body loop, we search the neighbors of both atom1
                ! and atom2 to find the triplets 
                ! atom1-atom2-atom3 and
                ! atom2-atom1-atom3
                ! These are considered to be different, since the middle atom of
                ! the triplet is different (atom2 vs. atom1).
                !
                ! We loop over all atoms to get atom1, 
                ! then over the neighbors of the atom1 to get atom2,
                ! then again over the neighbors of both atom1 and atom2 to get atom3.
                ! We want to find every triplet A-B-C, B-A-C, A-C-B exactly once filtering by the
                ! ordering of the indices of the atoms.
                ! Triplets A-B-C are considered equal to C-B-A and should be only found once.
                !
                ! Consider the indices A: 1, B: 2, C: 3. For other orderings, we can just
                ! permutate the names A, B and C so this is not affecting the generality
                ! of the argument.
                !
                ! We already filter by index1 < index2 when searching for atom2.
                ! Therefore the possible ways to get atom1 and atom2 for these orderings are:
                !
                !  A B C   atom1 atom2  or  atom1 atom2  or  atom1 atom2 
                !  1 2 3   A : 1 B : 2      A : 1 C : 3      B : 2 C : 3
                !
                ! If we filter atom3 by index3 > index2 when searching atom1 neighbors
                ! and index3 > index1 when searching atom2 neighbors (i.e., index3
                ! greater than the index of the atom whose neighbors are not searched),
                ! we get:
                !
                !  A B C   atom1 atom2  atom3 as atom1 nbor / atom2 nbor
                !  1 2 3   A : 1 B : 2  C : 3 -> found B-A-C
                !                       C : 3 -> found A-B-C
                !          A : 1 C : 3  B : 2 -> B < C (2 < 3) so ignored
                !                       B : 2 -> found A-C-B
                !          B : 2 C : 3  A : 1 -> A < B (1 < 2) so ignored
                !                       A : 1 -> A < C (1 < 3) so ignored

                if(many_bodies_found)then

                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ! apply 3-body bond order factors !
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                   ! neighbors of atom2
                   nbors2 = atom2%neighbor_list                   
                   bond_indices2 => atom2%bond_indices

                   ! First we try to find ordered triplets atom2 -- atom1 -- atom3
                   ! Therefore we need separations a2--a1 and a1--a3.
                   separations(1:3,1) = -separations(1:3,1)

                   ! loop over neighbors of atom 1
                   do l = 1, nbors1%n_neighbors
                      index3 = nbors1%neighbors(l)

                      ! the condition for finding each triplet once is such that
                      ! index 3 must be higher than the index of the atom whose
                      ! neighbors are NOT currently searched
                      ! For the offset check, we want atom2->atom3 offset, which is
                      ! (atom1->atom3) - (atom1->atom2), the latter being stored in offset
                      if(pick(index2,index3,nbors1%pbc_offsets(1:3,l)-offset))then

                         ! third atom of the triplet
                         atom3 = atoms(index3)
                         ! atom3 is new so we don't know the separation from atom1
                         separation3_unknown = .true.
                         ! The list of atoms is passed to bond factor evaluation routine
                         ! for further filtering.
                         ! This is triplet atom2 - atom1 - atom3, since we loop over
                         ! neighbors of atom1.
                         atom_list = (/ atom2, atom1, atom3 /)

                         ! search for the first bond params containing the parameters for atom1-atom2
                         do k1 = 1, size(bond_indices)
                            
                            bond_params(1) = bond_factors(bond_indices(k1))
                            
                            ! filter the parameters by:
                            ! number of targets, atom2 and atom3 being targets, group index                   
                            call bond_order_factor_is_in_group(bond_params(1),&
                                 group_index,is_in_group) ! in Potentials.f90
                            
                            if( is_in_group .and. bond_params(1)%n_level == 1 )then
                               call bond_order_factor_affects_atom(bond_params(1),&
                                    atom2,is_active,2) ! in Potentials.f90
                               if( is_active )then
                                  call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                                       n_targets) ! in Potentials.f90
                                  if( n_targets == 3 )then
                                     call bond_order_factor_affects_atom(bond_params(1),&
                                          atom3,is_active,3) ! in Potentials.f90
                                     if( is_active )then
                                                  
                         ! search for the second bond params containing the 
                         ! parameters for atom1-atom3
                         do k2 = 1, size(bond_indices)

                            bond_params(2) = bond_factors(bond_indices(k2))
                            call bond_order_factor_is_in_group(bond_params(2),&
                                 group_index,is_in_group) ! in Potentials.f90
                            if( is_in_group .and. bond_params(2)%n_level == 1 )then 
                               call bond_order_factor_affects_atom(bond_params(2),&
                                    atom3,is_active,2) ! in Potentials.f90
                               if( is_active )then
                                  call get_number_of_targets_of_bond_order_factor_index(&
                                       bond_params(2)%type_index,&
                                       n_targets) ! in Potentials.f90
                                  if( n_targets == 3 )then
                                     call bond_order_factor_affects_atom(bond_params(2),&
                                          atom2,is_active,3) ! in Potentials.f90
                                     if( is_active )then
                                        
                                        ! When we loop over the bond factors
                                        ! we may need the atom1-atom3 distance
                                        ! repeatedly. We only calculate it the first
                                        ! time.
                                        if( separation3_unknown )then
                                           call separation_vector(atom1%position, &
                                                atom3%position, &
                                                nbors1%pbc_offsets(1:3,j), &
                                                cell, &
                                                separations(1:3,2)) ! in Geometry.f90
                                           separation3_unknown = .false.
                                           distances(2) = .norm.(separations(1:3,2))
                                           if(distances(2) == 0.d0)then
                                              directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                           else
                                              directions(1:3,2) = separations(1:3,2) / distances(2)
                                           end if
                                        end if
                                        
                                        ! Evaluate the atom2-atom1-atom3 triplet contribution
                                        ! in the bond order factor sum.
                                        call evaluate_bond_order_factor(3,separations(1:3,1:2),&
                                             distances(1:2),bond_params(1:2),tmp_factor(1:3),atom_list)
                                        
                                        bo_temp(index2) = bo_temp(index2) + tmp_factor(1)
                                        bo_temp(index1) = bo_temp(index1) + tmp_factor(2)
                                        bo_temp(index3) = bo_temp(index3) + tmp_factor(3)
                                        
                                     end if ! is_active
                                  end if ! n_targets == 3
                               end if ! is_active
                            end if ! is_in_group
                         end do ! k2
                         
                                     end if ! is_active
                                  end if ! n_targets == 3
                               end if ! is_active
                            end if ! is_in_group
                         end do ! k1

                      end if ! index3 > index2

                   end do ! l = 1, nbors1%n_neighbors


                   ! Next we try to find ordered triplets atom1 -- atom2 -- atom3
                   ! Therefore we need separations a2--a1 and a2--a3.
                   separations(1:3,1) = -separations(1:3,1)
                   directions(1:3,1) = -directions(1:3,1)

                   ! loop over neighbors of atom 2
                   do l = 1, nbors2%n_neighbors
                      index3 = nbors2%neighbors(l)
                      
                      ! the condition for finding each triplet once is such that
                      ! index 3 must be higher than the index of the atom whose
                      ! neighbors are NOT currently searched
                      ! For the offset check we need atom1->atom3 which equals
                      ! (atom1->atom2) + (atom2->atom3), the former being stored in offset
                      if(pick(index1,index3,nbors2%pbc_offsets(1:3,l)+offset))then

                         ! Third atom of the triplet.
                         atom3 = atoms(index3)
                         separation3_unknown = .true.
                         atom_list = (/ atom1, atom2, atom3 /)

                         ! search for the first bond params containing the parameters for atom2-atom1
                         do k1 = 1, size(bond_indices2)
                         
                            bond_params(1) = bond_factors(bond_indices2(k1))
                            
                            ! filter the parameters by:
                            ! number of targets, atom2 and atom3 being targets, group index
                            call bond_order_factor_is_in_group(bond_params(1),&
                                 group_index,is_in_group) ! in Potentials.f90                         
                            if( is_in_group .and. bond_params(1)%n_level == 1 )then 
                               call bond_order_factor_affects_atom(bond_params(1),&
                                    atom1,is_active,2) ! in Potentials.f90
                               if( is_active )then 
                                  call get_number_of_targets_of_bond_order_factor_index(&
                                       bond_params(1)%type_index,&
                                       n_targets) ! in Potentials.f90
                                  if( n_targets == 3 )then
                                     call bond_order_factor_affects_atom(bond_params(1),&
                                          atom3,is_active,3) ! in Potentials.f90                            
                                     if( is_active )then
                               
                         ! search for the second bond params containing the parameters for atom2-atom3
                         do k2 = 1, size(bond_indices2)
                         
                            bond_params(2) = bond_factors(bond_indices2(k2))
                            call bond_order_factor_is_in_group(bond_params(2),&
                                 group_index,is_in_group) ! in Potentials.f90                         
                            if( is_in_group .and. bond_params(2)%n_level == 1 )then 
                               call bond_order_factor_affects_atom(bond_params(2),&
                                    atom3,is_active,2) ! in Potentials.f90
                               if( is_active )then 
                                  call get_number_of_targets_of_bond_order_factor_index(&
                                       bond_params(2)%type_index,&
                                       n_targets) ! in Potentials.f90
                                  if( n_targets == 3 )then
                                     call bond_order_factor_affects_atom(bond_params(2),&
                                          atom1,is_active,3) ! in Potentials.f90
                                     if( is_active )then

                                        if( separation3_unknown )then
                                           call separation_vector(atom2%position, &
                                                atom3%position, &
                                                nbors2%pbc_offsets(1:3,l), &
                                                cell, &
                                                separations(1:3,2))  ! in Geometry.f90
                                           separation3_unknown = .false.
                                           distances(2) = .norm.(separations(1:3,2))
                                           if(distances(2) == 0.d0)then
                                              directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                                           else
                                              directions(1:3,2) = separations(1:3,2) / distances(2)
                                           end if
                                        end if
                                        
                                        ! Evaluate the atom1-atom2-atom3 triplet contribution
                                        ! in the bond order factor sum.
                                        call evaluate_bond_order_factor(3,separations(1:3,1:2),&
                                             distances(1:2),bond_params(1:2),tmp_factor(1:3),atom_list)

                                        bo_temp(index1) = bo_temp(index1) + tmp_factor(1)
                                        bo_temp(index2) = bo_temp(index2) + tmp_factor(2)
                                        bo_temp(index3) = bo_temp(index3) + tmp_factor(3)
                                     
                                     end if ! is_active
                                  end if ! n_targets == 3
                               end if ! is_active
                            end if ! is_in_group
                         end do ! k2
                         
                                     end if ! is_active
                                  end if ! n_targets == 3
                               end if ! is_active
                            end if ! is_in_group
                         end do ! k1

                      end if ! index3 > index2

                   end do ! l = 1, nbors2%n_neighbors

                end if ! many_bodies_found

             end if ! index2 > index1
          end do ! j = 1, nbors1%n_neighbors

       end if ! is_my_atom

    end do ! index1 = 1, size(atoms) 

#ifdef MPI
    ! sum the contributions from all the processors
    call mpi_allreduce(bo_temp,total_bond_orders,size(total_bond_orders),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
#else
    total_bond_orders = bo_temp
#endif

  end subroutine core_calculate_bond_order_factors











! !!!: core_calculate_pair_bond_order_factor

  ! Calculates the bond order sum for a given pair of atoms for the given group.
  !
  ! For a factor such as
  !
  ! .. math::
  !
  !      b_ij = f(\sum_k c_{ijk})
  !
  ! The routine calculates
  !
  ! .. math::
  !
  !      \sum_k c_{ijk}.
  !
  ! The full bond order factor is then obtained by applying the
  ! scaling function :math:`f`. This is done with
  ! :func:`core_post_process_bond_order_factors`.
  !
  ! *n_atoms number of atoms
  ! *group_index an index denoting the potential to which the factor is connected
  ! *bond_order_sum the calculated bond order sums
  subroutine core_calculate_pair_bond_order_factor(atom_pair,separation,distance,direction,group_index,bond_order_sum)
    implicit none
    integer, intent(in) :: atom_pair(2)
    integer, intent(in) :: group_index
    double precision, intent(in) :: separation(3), distance, direction(3)
    double precision, intent(out) :: bond_order_sum(2)
    integer :: index1, index2, index3, k1, j, n_targets
    double precision :: separations(3,2), distances(2), directions(3,2)
    double precision :: tmp_factor
    type(atom) :: atom1, atom2, atom3
    type(atom) :: atom_list(3)
    type(neighbor_list) :: nbors1, nbors2
    type(bond_order_parameters) :: bond_params(1)
    integer, pointer :: bond_indices(:), bond_indices2(:)
    logical :: is_active, is_in_group
    integer :: offset(3)

    bond_order_sum = 0.d0

    ! target atom 1
    index1 = atom_pair(1)
    atom1 = atoms(index1)
    nbors1 = atom1%neighbor_list
    bond_indices => atom1%bond_indices

    ! target atom 2
    index2 = atom_pair(2)
    atom2 = atoms(index2)
    nbors2 = atom2%neighbor_list
    bond_indices2 => atom2%bond_indices


    ! In the 3-body loop, we search the neighbors of both atom1
    ! and atom2 to find the triplets 
    ! atom1-atom2-atom3 and
    ! atom2-atom1-atom3
    ! These are considered to be different, since the middle atom of
    ! the triplet is different (atom2 vs. atom1).
    
    ! First we try to find ordered triplets atom2 -- atom1 -- atom3
    ! Therefore we need separations a2--a1 and a1--a3.

    ! atom2 - atom1 vector
    separations(1:3,1) = -separation(1:3)
    distances(1) = distance
    directions(1:3,1) = -direction(1:3)


    ! loop over neighbors of atom 1
    do j = 1, nbors1%n_neighbors
       
       ! neighboring atom
       index3 = nbors1%neighbors(j)
       offset(1:3) = nbors1%pbc_offsets(1:3,j)

       ! To prevent the case index2 - index1 - index2
       if(index2 /= index3)then

          atom3 = atoms(index3)
          atom_list = (/ atom2, atom1, atom3 /)

          ! calculate atom1-atom3 separation vector
          ! and distance
          call separation_vector(atom1%position, &
               atom3%position, &
               offset(1:3), &
               cell, &
               separations(1:3,2)) ! in Geometry.f90
          distances(2) = .norm.(separations(1:3,2))
          if(distances(2) == 0.d0)then
             directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
          else
             directions(1:3,2) = separations(1:3,2) / distances(2)
          end if
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! apply 3-body bond order factors !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if(distances(2) < atom1%max_bond_radius)then

             ! search for the first bond params containing the parameters for (atom2-atom1) -- atom3
             do k1 = 1, size(bond_indices)
                            
                bond_params(1) = bond_factors(bond_indices(k1))
                            
                ! filter the parameters by:
                ! number of targets, atom2 and atom3 being targets, group index                   
                call bond_order_factor_is_in_group(bond_params(1),&
                     group_index,is_in_group) ! in Potentials.f90
                
                !write(*,*) "group, level", is_in_group, bond_params(1)%n_level
                
                if( is_in_group .and. bond_params(1)%n_level == 2 )then
                   call bond_order_factor_affects_atom(bond_params(1),&
                        atom2,is_active,2) ! in Potentials.f90
                   if( is_active )then
                      call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                           n_targets) ! in Potentials.f90
                      if( n_targets == 3 )then
                         call bond_order_factor_affects_atom(bond_params(1),&
                              atom3,is_active,3) ! in Potentials.f90
                         if( is_active )then
                            
                            ! Evaluate the atom2-atom1-atom3 triplet contribution
                            ! in the bond order factor sum.
                            call evaluate_pair_bond_order_factor(3,separations(1:3,1:2),&
                                 distances(1:2),bond_params(1),tmp_factor,atom_list)
                            
                            bond_order_sum(1) = bond_order_sum(1) + tmp_factor
                            
                         end if ! is_active
                      end if ! n_targets == 3
                   end if ! is_active
                end if ! is_in_group
             end do ! k1
          end if ! distance < max_cut
       end if ! index3 /= index2

    end do ! j = 1, nbors1%n_neighbors

    
    ! Next we try to find ordered triplets atom1 -- atom2 -- atom3
    ! Therefore we need separations a1--a2 and a2--a3.
    separations(1:3,1) = -separations(1:3,1)
    directions(1:3,1) = -directions(1:3,1)

    ! loop over neighbors of atom 2
    do j = 1, nbors2%n_neighbors
       
       ! neighboring atom
       index3 = nbors2%neighbors(j)
       offset(1:3) = nbors2%pbc_offsets(1:3,j)

       ! To prevent the case index1 - index2 - index1
       if(index1 /= index3)then

          atom3 = atoms(index3)
          atom_list = (/ atom1, atom2, atom3 /)

          ! calculate atom1-atom3 separation vector
          ! and distance
          call separation_vector(atom2%position, &
               atom3%position, &
               offset(1:3), &
               cell, &
               separations(1:3,2)) ! in Geometry.f90
          distances(2) = .norm.(separations(1:3,2))
          if(distances(2) == 0.d0)then
             directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
          else
             directions(1:3,2) = separations(1:3,2) / distances(2)
          end if
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! apply 3-body bond order factors !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if(distances(2) < atom2%max_bond_radius)then

             ! search for the bond params containing the parameters for (atom1-atom2) -- atom3
             do k1 = 1, size(bond_indices2)
                
                bond_params(1) = bond_factors(bond_indices2(k1))
                
                ! filter the parameters by:
                ! number of targets, atom2 and atom3 being targets, group index                   
                call bond_order_factor_is_in_group(bond_params(1),&
                     group_index,is_in_group) ! in Potentials.f90
                
                if( is_in_group .and. bond_params(1)%n_level == 2 )then
                   call bond_order_factor_affects_atom(bond_params(1),&
                        atom1,is_active,2) ! in Potentials.f90
                   if( is_active )then
                      call get_number_of_targets_of_bond_order_factor_index(bond_params(1)%type_index,&
                           n_targets) ! in Potentials.f90
                      if( n_targets == 3 )then
                         call bond_order_factor_affects_atom(bond_params(1),&
                              atom3,is_active,3) ! in Potentials.f90
                         if( is_active )then
                            
                            ! Evaluate the atom2-atom1-atom3 triplet contribution
                            ! in the bond order factor sum.
                            call evaluate_pair_bond_order_factor(3,separations(1:3,1:2),&
                                 distances(1:2),bond_params(1),tmp_factor,atom_list)
                            
                            bond_order_sum(2) = bond_order_sum(2) + tmp_factor
                            
                         end if ! is_active
                      end if ! n_targets == 3
                   end if ! is_active
                end if ! is_in_group
             end do ! k1
          end if ! distance < max_cut

       end if ! index3 /= index1

    end do ! j = 1, nbors2%n_neighbors


  end subroutine core_calculate_pair_bond_order_factor












! !!!: core_post_process_bond_order_factors

  ! Bond-order post processing, i.e., application of per-atom scaling functions.
  !
  ! By post processing, we mean any operations done after calculating the
  ! sum of pair- and many-body terms. That is, if a factor is, say,
  !
  ! .. math::
  !
  !      b_i = f(\sum_j c_{ij}) = 1 + \sum_j c_{ij},
  !
  ! the :math:`\sum_j c_{ij}` would have been calculated already 
  ! (with :func:`core_calculate_bond_order_factors`)
  ! and the operation :math:`f(x) = 1 + x`
  ! remains to be carried out.
  ! The post processing is done per atom regardless of if the
  ! bond factor is of a pair or many body type.
  !
  ! *n_atoms number of atoms
  ! *group_index an index denoting the potential to which the factor is connected
  ! *raw_sums precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
  ! *total_bond_orders the calculated bond order factors :math:`b_i`
  subroutine core_post_process_bond_order_factors(group_index,raw_sums,total_bond_orders)
    implicit none
    integer, intent(in) :: group_index
    double precision, intent(inout) :: total_bond_orders(:)
    integer :: index1, index2
    double precision, intent(in) :: raw_sums(:)
    type(atom) :: atom1
    type(bond_order_parameters) :: bond_params
    integer, pointer :: bond_indices(:)
    integer :: post_process

    ! loop over all atoms
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
             if( bond_params%includes_post_processing .and. bond_params%n_level == 1 )then
                if( bond_params%original_elements(1) == atom1%element )then
                   post_process = bond_indices(index2)
                   exit
                end if
             end if
          end do

          if( post_process > 0 )then
             call post_process_bond_order_factor(raw_sums(index1),&
                  bond_factors( post_process ), &
                  bo_temp(index1) ) ! in Potentials.f90
          else
             bo_temp(index1) = raw_sums(index1)
          end if

       else
          bo_temp(index1) = 0.d0
       end if

    end do

#ifdef MPI
    ! gather the bond order factors from all cpus
    call mpi_allreduce(bo_temp,total_bond_orders,size(total_bond_orders),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
#else
    total_bond_orders = bo_temp
#endif

  end subroutine core_post_process_bond_order_factors







! !!!: core_post_process_pair_bond_order_factors

  ! Bond-order post processing, i.e., application of per-atom scaling functions.
  !
  ! By post processing, we mean any operations done after calculating the
  ! sum of pair- and many-body terms. That is, if a factor is, say,
  !
  ! .. math::
  !
  !      b_{ij} = f(\sum_k c_{ijk}) = 1 + \sum_k c_{ijk},
  !
  ! the :math:`\sum_k c_{ijk}` would have been calculated already 
  ! (with :func:`core_calculate_pair_bond_order_factor`)
  ! and the operation :math:`f(x) = 1 + x`
  ! remains to be carried out.
  !
  ! *n_atoms number of atoms
  ! *group_index an index denoting the potential to which the factor is connected
  ! *raw_sums precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example.
  ! *total_bond_orders the calculated bond order factors :math:`b_i`
  subroutine core_post_process_pair_bond_order_factor(atom1,group_index,raw_sum,total_bond_order)
    implicit none
    type(atom), intent(in) :: atom1
    integer, intent(in) :: group_index
    double precision, intent(out) :: total_bond_order
    integer :: index1, index2
    double precision, intent(in) :: raw_sum
    type(bond_order_parameters) :: bond_params
    integer, pointer :: bond_indices(:)
    integer :: post_process

    ! target atom bond parameters
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
       if( bond_params%includes_post_processing .and. bond_params%n_level == 2 )then
          if( bond_params%original_elements(1) == atom1%element )then
             post_process = bond_indices(index2)
             exit
          end if
       end if
    end do
          
    if( post_process > 0 )then
       call post_process_bond_order_factor(raw_sum,&
            bond_factors( post_process ), &
            total_bond_order ) ! in Potentials.f90
    else
       total_bond_order = raw_sum
    end if
    

  end subroutine core_post_process_pair_bond_order_factor





! !!!: core_post_process_bond_order_gradients

  ! Bond-order post processing, i.e., application of per-atom scaling functions.
  ! This routine does the scaling for all bond factors with the given
  ! bond order sums and gradients of these sums.
  !
  ! By post processing, we mean any operations done after calculating the
  ! sum of pair- and many-body terms. That is, if a factor is, say,
  !
  ! .. math::
  !
  !      b_i = f(\sum_j c_{ij}) = 1 + \sum_j c_{ij},
  !
  ! the :math:`\sum_j c_{ij}` would have been calculated already and the 
  ! operation :math:`f(x) = 1 + x` remains to be carried out.
  ! The post processing is done per atom regardless of if the
  ! bond factor is of a pair or many body type.
  !
  ! For gradients, one needs to evaluate
  !
  ! .. math::
  !
  !     \nabla_\alpha b_i = f'(\sum_j c_{ij}) \nabla_\alpha \sum_j c_{ij}
  !
  ! *group_index an index denoting the potential to which the factor is connected
  ! *raw_sums precalculated bond order sums, :math:`\sum_j c_{ij}`, in the above example
  ! *raw_gradients precalculated gradients of bond order sums, :math:`\nabla_\alpha \sum_j c_{ij}`, in the above example
  ! *total_bond_gradients the calculated bond order gradients :math:`\nabla_\alpha b_i`
  ! *mpi_split A switch for enabling MPI parallelization. By default the routine is sequential since the calculation may be called from within an already parallelized routine.
  subroutine core_post_process_bond_order_gradients(group_index,raw_sums,&
       raw_gradients,total_bond_gradients,mpi_split)
    implicit none
    integer, intent(in) :: group_index
    double precision, intent(inout) :: total_bond_gradients(:,:)
    double precision, intent(in) :: raw_sums(:), raw_gradients(:,:)
    logical, optional, intent(in) :: mpi_split
    integer :: index1, index2
    type(atom) :: atom1
    type(neighbor_list) :: nbors1
    type(bond_order_parameters) :: bond_params
    integer, pointer :: bond_indices(:)
    integer :: post_process, n_atoms
    logical :: evaluate, is_in_group

    n_atoms = size(atoms)
    ! loop over oll bond order factors
    do index1 = 1, size(atoms)
       
       evaluate = .true.
       ! In MPI, only consider the atoms allocated to this particular cpu.
       ! However, the routine may be called from within an already
       ! parallel part of the code (force calculation).
       ! Therefore, parallellization is made optional.
       if(present(mpi_split))then
          if(is_my_atom(index1) .or. .not.mpi_split)then
             evaluate = .true.
          else
             evaluate = .false.
          end if
       end if
       if(evaluate)then
          ! target atom whose factor is differentiated
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
             call bond_order_factor_is_in_group(bond_params,group_index,is_in_group) ! in Potentials.f90
             if(is_in_group)then
                if( bond_params%includes_post_processing .and. bond_params%n_level == 1 )then
                   if( bond_params%original_elements(1) == atom1%element )then
                      ! store the index of the parameters 
                      post_process = bond_indices(index2)
                      exit
                   end if
                end if
             end if
          end do

          ! if scaling parameters were found, do the post processing
          if( post_process > 0 )then
             call post_process_bond_order_gradient(raw_sums(index1),&
                  raw_gradients(1:3,index1),&
                  bond_factors( post_process ), &
                  temp_gradient(1:3,index1,1) ) ! in Potentials.f90
          else
             temp_gradient(1:3,index1,1) = raw_gradients(1:3,index1)
          end if

       else
          temp_gradient(1:3,index1,1) = 0.d0
       end if

    end do

#ifdef MPI
    ! Gather data from all cpus if needed.
    if(present(mpi_split))then
       if(mpi_split)then
          call mpi_allreduce(temp_gradient(1:3,1:n_atoms,1),total_bond_gradients(1:3,1:n_atoms),&
               size(total_bond_gradients),mpi_double_precision,&
               mpi_sum,mpi_comm_world,mpistat)
       else
          total_bond_gradients(1:3,1:n_atoms) = temp_gradient(1:3,1:n_atoms,1)
       end if
    else
       total_bond_gradients(1:3,1:n_atoms) = temp_gradient(1:3,1:n_atoms,1)
    end if
#else
    total_bond_gradients = bond_gradients
#endif

  end subroutine core_post_process_bond_order_gradients



! !!!: core_post_process_bond_order_gradients_of_factor

  ! Bond-order post processing, i.e., application of per-atom scaling functions.
  ! This routine does the scaling for the bond order factor of the given atom
  ! with respect to moving all atoms
  ! with the given bond order sum for the factor and 
  ! the gradients of the sum with respect to moving all atoms.
  !
  ! By post processing, we mean any operations done after calculating the
  ! sum of pair- and many-body terms. That is, if a factor is, say,
  !
  ! .. math::
  !
  !      b_i = f(\sum_j c_{ij}) = 1 + \sum_j c_{ij},
  !
  ! the :math:`\sum_j c_{ij}` would have been calculated already and the operation :math:`f(x) = 1 + x`
  ! remains to be carried out.
  ! The post processing is done per atom regardless of if the
  ! bond factor is of a pair or many body type.
  !
  ! For gradients, one needs to evaluate
  !
  ! .. math::
  !
  !     \nabla_\alpha b_i = f'(\sum_j c_{ij}) \nabla_\alpha \sum_j c_{ij}
  !
  ! *group_index an index denoting the potential to which the factor is connected
  ! *atom_index the index of the atom whose factor is differentiated (:math:`i`)
  ! *raw_sum precalculated bond order sum for the given atom, :math:`\sum_j c_{ij}`, in the above example
  ! *raw_gradients precalculated gradients of bond order sums, :math:`\nabla_\alpha \sum_j c_{ij}`, in the above example
  ! *total_bond_gradients the calculated bond order gradients :math:`\nabla_\alpha b_i`
  ! *raw_virial the precalculated virial due to the bond order gradient
  ! *total_virial the scaled  virial due to the bond order gradient
  ! *mpi_split A switch for enabling MPI parallelization. By default the routine is sequential since the calculation may be called from within an already parallelized routine.
  subroutine core_post_process_bond_order_gradients_of_factor(group_index,atom_index,raw_sum,&
       raw_gradients,total_bond_gradients,raw_virial,total_virial,mpi_split)
    implicit none
    integer, intent(in) :: group_index, atom_index
    double precision, intent(inout) :: total_bond_gradients(:,:), total_virial(6)
    double precision, intent(in) :: raw_sum, raw_gradients(:,:), raw_virial(6)
    logical, optional, intent(in) :: mpi_split
    integer :: index1, index2, n_atoms
    type(atom) :: atom1
    type(neighbor_list) :: nbors1
    type(bond_order_parameters) :: bond_params
    integer, pointer :: bond_indices(:)
    integer :: post_process
    logical :: evaluate, is_in_group

    n_atoms = size(atoms)
    ! target atom whose factor is differentiated
    atom1 = atoms(atom_index)
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
       call bond_order_factor_is_in_group(bond_params,group_index,is_in_group) ! in Potentials.f90
       if(is_in_group)then
          if( bond_params%includes_post_processing .and. bond_params%n_level == 1 )then
             if( bond_params%original_elements(1) == atom1%element )then
                ! store the index of the parameters 
                post_process = bond_indices(index2)
                exit
             end if
          end if
       end if
    end do

    ! if scaling parameters were found, do the post processing
    if( post_process > 0 )then

       ! loop over moved atoms
       do index1 = 1, size(atoms)
          
          if(bo_scaling(index1))then
             evaluate = .true.
          else
             evaluate = .false.
          end if
          ! in MPI, only consider the atoms allocated to this particular cpu
          if(evaluate .and. present(mpi_split))then
             evaluate = .false.
             if(is_my_atom(index1) .or. .not.mpi_split)then
                evaluate = .true.
             end if
          end if
          if(evaluate)then
             
             call post_process_bond_order_gradient(raw_sum,&
                  raw_gradients(1:3,index1),&
                  bond_factors( post_process ), &
                  temp_gradient(1:3,index1,1) ) ! in Potentials.f90
          else
             temp_gradient(1:3,index1,1) = 0.d0
          end if ! evaluate
          
       end do ! index1
       
       ! Post process the virial in two batches
       ! The scaling function is the same as for the
       ! gradient itself
       call post_process_bond_order_gradient(raw_sum, &
            raw_virial(1:3), &
            bond_factors( post_process), &
            total_virial(1:3) ) ! in Potentials.f90
       call post_process_bond_order_gradient(raw_sum, &
            raw_virial(4:6), &
            bond_factors( post_process), &
            total_virial(4:6) ) ! in Potentials.f90

    else
       ! there is no post processing for this factor
       temp_gradient(1:3,1:n_atoms,1) = raw_gradients(1:3,1:n_atoms)
       total_virial(1:6) = raw_virial(1:6)
    end if ! post_process

#ifdef MPI
    ! collect data from all cpus if needed
    if(present(mpi_split))then
       if(mpi_split .and. post_process > 0)then ! if we did not scale, the values are already good
          call mpi_allreduce(temp_gradient(1:3,1:n_atoms,1),total_bond_gradients(1:3,1:n_atoms),&
               size(total_bond_gradients),mpi_double_precision,&
               mpi_sum,mpi_comm_world,mpistat)
       else
          total_bond_gradients(1:3,1:n_atoms) = temp_gradient(1:3,1:n_atoms,1)
       end if
    else
       total_bond_gradients(1:3,1:n_atoms) = temp_gradient(1:3,1:n_atoms,1)
    end if
#else
    total_bond_gradients(1:3,1:n_atoms) = temp_gradient(1:3,1:n_atoms,1)
#endif

  end subroutine core_post_process_bond_order_gradients_of_factor





! !!!: core_post_process_pair_bond_order_gradients

  ! Bond-order post processing, i.e., application of per-pair scaling functions.
  ! This routine does the scaling for the bond order factor of the given pair
  ! with respect to moving all atoms
  ! with the given bond order sum for the factor and 
  ! the gradients of the sum with respect to moving all atoms.
  !
  ! By post processing, we mean any operations done after calculating the
  ! sum of pair- and many-body terms. That is, if a factor is, say,
  !
  ! .. math::
  !
  !      b_{ij} = f(\sum_k c_{ijk}) = 1 + \sum_k c_{ijk},
  !
  ! the :math:`\sum_k c_{ijk}` would have been calculated already and the operation :math:`f(x) = 1 + x`
  ! remains to be carried out.
  ! The post processing is done per pair.
  !
  ! For gradients, one needs to evaluate
  !
  ! .. math::
  !
  !     \nabla_\alpha b_{ij} = f'(\sum_k c_{ijk}) \nabla_\alpha \sum_k c_{ijk}
  !
  ! *group_index an index denoting the potential to which the factor is connected
  ! *atom_index the index of the atom whose factor is differentiated (:math:`i`)
  ! *raw_sum precalculated bond order sum for the given atom, :math:`\sum_j c_{ij}`, in the above example
  ! *raw_gradients precalculated gradients of bond order sums, :math:`\nabla_\alpha \sum_j c_{ij}`, in the above example
  ! *total_bond_gradients the calculated bond order gradients :math:`\nabla_\alpha b_i`
  ! *raw_virial the precalculated virial due to the bond order gradient
  ! *total_virial the scaled  virial due to the bond order gradient
  ! *mpi_split A switch for enabling MPI parallelization. By default the routine is sequential since the calculation may be called from within an already parallelized routine.
  subroutine core_post_process_pair_bond_order_gradients(group_index,atom1,raw_sum,&
       raw_gradients,total_bond_gradients,raw_virial,total_virial,mpi_split)
    implicit none
    integer, intent(in) :: group_index
    double precision, intent(out) :: total_bond_gradients(:,:), total_virial(6)
    double precision, intent(in) :: raw_sum, raw_gradients(:,:), raw_virial(6)
    logical, optional, intent(in) :: mpi_split
    integer :: index1, index2
    type(atom), intent(in) :: atom1
    type(neighbor_list) :: nbors1
    type(bond_order_parameters) :: bond_params
    integer, pointer :: bond_indices(:)
    integer :: post_process, n_atoms
    logical :: evaluate, is_in_group

    n_atoms = size(atoms)

    ! target atom whose factor is differentiated
    !atom1 = atoms(atom_index)
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
       call bond_order_factor_is_in_group(bond_params,group_index,is_in_group) ! in Potentials.f90
       if(is_in_group)then
          if( bond_params%includes_post_processing .and. bond_params%n_level == 2 )then
             if( bond_params%original_elements(1) == atom1%element )then
                ! store the index of the parameters 
                post_process = bond_indices(index2)
                exit
             end if
          end if
       end if
    end do

    ! if scaling parameters were found, do the post processing
    if( post_process > 0 )then

       ! loop over moved atoms
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
             
             call post_process_bond_order_gradient(raw_sum,&
                  raw_gradients(1:3,index1),&
                  bond_factors( post_process ), &
                  temp_gradient(1:3,index1,1) ) ! in Potentials.f90

          else
             temp_gradient(1:3,index1,1) = 0.d0
          end if ! evaluate
          
       end do ! index1
       
       ! Post process the virial in two batches
       ! The scaling function is the same as for the
       ! gradient itself
       call post_process_bond_order_gradient(raw_sum, &
            raw_virial(1:3), &
            bond_factors( post_process), &
            total_virial(1:3) ) ! in Potentials.f90
       call post_process_bond_order_gradient(raw_sum, &
            raw_virial(4:6), &
            bond_factors( post_process), &
            total_virial(4:6) ) ! in Potentials.f90

    else
       ! there is no post processing for this factor
       temp_gradient(1:3,1:n_atoms,1) = raw_gradients(1:3,1:n_atoms)
       total_virial(1:6) = raw_virial(1:6)
    end if ! post_process

#ifdef MPI
    ! collect data from all cpus if needed
    if(present(mpi_split))then
       if(mpi_split .and. post_process > 0)then ! if we did not scale, the values are already good
          call mpi_allreduce(temp_gradient(1:3,1:n_atoms,1),total_bond_gradients(1:3,1:n_atoms),&
               size(total_bond_gradients),mpi_double_precision,&
               mpi_sum,mpi_comm_world,mpistat)
       else
          total_bond_gradients(1:3,1:n_atoms) = temp_gradient(1:3,1:n_atoms,1)
       end if
    else
       total_bond_gradients(1:3,1:n_atoms) = temp_gradient(1:3,1:n_atoms,1)
    end if
#else
    total_bond_gradients(1:3,1:n_atoms) = temp_gradient(1:3,1:n_atoms,1)
#endif

  end subroutine core_post_process_pair_bond_order_gradients






! !!!: core_loop_over_local_interactions

  ! Loops over atoms, atomic pairs, atomic triplets, and atomic quadruplets
  ! and calculates the contributions from local potentials to energy, forces, 
  ! or electronegativities. This routine is called from the routines
  !
  !  - :func:`core_calculate_energy`
  !  - :func:`core_calculate_forces`
  !  - :func:`core_calculate_electronegaivities`
  !
  ! *calculation_type index to specify if the loop calculates energies, forces, or e-negativities
  ! *total_energy calculated energy
  ! *total_forces calculated forces
  ! *total_enegs calculated electronegativities
  ! *total_stress calculated stress
  subroutine core_loop_over_local_interactions(calculation_type,&
       total_energy,total_forces,total_enegs,total_stress)
    implicit none
    integer, intent(in) :: calculation_type
    double precision, intent(inout) :: total_energy, total_forces(:,:), &
         total_enegs(:), total_stress(6)
    integer :: j, k, l, m, n_targets, index1, index2, index3, index4
    double precision,save :: energy,  &
         stress(6), &
         separations(3,3), distances(3), directions(3,3), &
         stopwatch_0, stopwatch_1
    type(atom) :: atom1, atom2, atom3, atom4
    type(atom), save :: atom_list(4)
    type(neighbor_list) :: nbors1, nbors2, nbors3
    integer, pointer :: interaction_indices(:)
    logical :: is_active, many_bodies_found
    double precision :: max_cutoff, sigma, t00, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, ta, tb
    integer, save :: reci_length(3), offset(3), tripleoffset(3), quadoffset(3), n_atoms

    t00 = 0.0
    t0 = 0.0
    t1 = 0.0
    t2 = 0.0
    t3 = 0.0
    t4 = 0.0
    t5 = 0.0
    t6 = 0.0
    t7 = 0.0
    t8 = 0.0
    t9 = 0.0
    ta = 0.0
    tb = 0.0

    n_atoms = size(atoms)

    energy = 0.d0
    total_energy = 0.d0

    temp_forces = 0.d0
    total_forces = 0.d0

    temp_enegs = 0.d0
    total_enegs = 0.d0

    stress = 0.d0
    total_stress = 0.d0

    ! For MPI load balancing, the execution time of each cpu
    ! is recorded. After the forces have been calculated, the
    ! workload of all cpus are examined and load is transferred
    ! between the cpus in order to make the workloads as equal
    ! as possible.
    call start_timer()
    call mpi_wall_clock(ta)
    t00 = ta

    ! Before starting the calculation proper,
    ! all per-atom bond order factors are calculated and
    ! stored in arrays.
    ! Thus, they need not be recalculated during the
    ! evaluation loops.
    use_saved_bond_order_factors = .true.
    call core_fill_bond_order_storage()

    call mpi_wall_clock(t0)
    t1 = t1+t0-t00

    ! loop over atoms
    do index1 = 1, n_atoms

       !write(*,*) "debug: atom ", index1, "/", n_atoms


       ! in MPI, only consider the atoms allocated to this particular cpu
       if(is_my_atom(index1))then

          ! Bond order gradients are not stored since there are potentially
          ! so many. Some most recent ones are saved, though.
          ! At the start of the first atom loop, we clear the storage.
          ! call core_empty_bond_order_gradient_storage()

          atom1 = atoms(index1)
          nbors1 = atom1%neighbor_list
          interaction_indices => atom1%potential_indices

          !*********************!
          ! 1-body interactions !
          !*********************!

          call mpi_wall_clock(t00)
          call core_evaluate_local_singlet(index1, &
               atom1,&
               interaction_indices,&
               calculation_type,energy,temp_forces,stress,temp_enegs)
          call mpi_wall_clock(t0)
          t2 = t2+t0-t00

          ! loop over neighbors
          do j = 1, nbors1%n_neighbors


             ! Note that we loop over the neighbors in the outer loop and
             ! over the interactions in the inner loop. This is to avoid calculating
             ! the interatomic distances repeatedly for multiple potentials affecting
             ! the same pair of atoms.

             ! neighboring atom
             index2 = nbors1%neighbors(j)
             offset(1:3) = nbors1%pbc_offsets(1:3,j) ! offset atom1 -> atom2

             ! Since we loop over the neighbors of all atoms, we will find the pair
             ! atom1-atom2 = atom2-atom1 twice.
             ! To prevent the double counting, we filter by index2 > index1.
             if(pick(index1,index2,offset))then

                call mpi_wall_clock(t00)
                ! Empty bond gradient storage for atom2 slot (since we have a new atom2)
                ! call core_empty_bond_order_gradient_storage(2)

                call mpi_wall_clock(t0)
                t7 = t7+t0-t00
                
                atom2 = atoms(index2)
                atom_list(1) = atom1
                atom_list(2) = atom2
                
                call mpi_wall_clock(t00)
                ! calculate atom1-atom2 separation vector
                ! and distance
                call separation_vector(atom1%position, &
                     atom2%position, &
                     nbors1%pbc_offsets(1:3,j), &
                     cell, &
                     separations(1:3,1)) ! in Geometry.f90
                distances(1) = .norm.(separations(1:3,1))
                if(distances(1) == 0.d0)then
                   directions(1:3,1) = (/ 0.d0, 0.d0, 0.d0 /)
                else
                   directions(1:3,1) = separations(1:3,1) / distances(1)
                end if
                call mpi_wall_clock(t0)
                t3 = t3+t0-t00

                !*********************!
                ! 2-body interactions !
                !*********************!

                t00 = t0
                if(distances(1) < atom1%max_potential_radius)then
                   ! differentiate between energy, force, and electronegativity evaluation
                   select case(calculation_type)
                   case(energy_evaluation_index)
                      call core_evaluate_local_doublet_energy_B(atom_list(1:2), &
                           index1, index2, &
                           2, & ! test for atom2, since interaction indices filters for atom1 already
                           interaction_indices, &
                           separations(1:3,1), directions(1:3,1), distances(1), &
                           energy, &
                           many_bodies_found)
                   case(force_evaluation_index)
                      if(.true.)then
                         call core_evaluate_local_doublet_forces_B(atom_list(1:2), &
                              index1, index2, &
                              2, & ! test for atom2, since interaction indices filters for atom1 already
                              interaction_indices, &
                              separations(1:3,1), directions(1:3,1), distances(1), &
                              temp_forces,stress, &
                              many_bodies_found)
                      end if
                   case(electronegativity_evaluation_index)
                      call core_evaluate_local_doublet_electronegativities_B(atom_list(1:2), &
                           index1, index2, &
                           2, & ! test for atom2, since interaction indices filters for atom1 already
                           interaction_indices, &
                           separations(1:3,1), directions(1:3,1), distances(1), &
                           temp_enegs, &
                           many_bodies_found)
                   end select
                end if
                call mpi_wall_clock(t0)
                t4 = t4+t0-t00

                ! Only do the 3-body loop if we found many-body potentials 
                ! during 2-body evaluation.
                ! 
                ! In the 3-body loop, we search the neighbors of both atom1
                ! and atom2 to find the triplets 
                ! atom1-atom2-atom3 and
                ! atom2-atom1-atom3
                ! These are considered to be different, since the middle atom of
                ! the triplet is different (atom2 vs. atom1).
                !
                ! We loop over all atoms to get atom1, 
                ! then over the neighbors of the atom1 to get atom2,
                ! then again over the neighbors of both atom1 and atom2 to get atom3.
                ! We want to find every triplet A-B-C, B-A-C, A-C-B exactly once filtering by the
                ! ordering of the indices of the atoms.
                ! Triplets A-B-C are considered equal to C-B-A and should be only found once.
                !
                ! Consider the indices A: 1, B: 2, C: 3. For other orderings, we can just
                ! permutate the names A, B and C so this is not affecting the generality
                ! of the argument.
                !
                ! We already filter by index1 < index2 when searching for atom2.
                ! Therefore the possible ways to get atom1 and atom2 for these orderings are:
                !
                !  A B C   atom1 atom2  or  atom1 atom2  or  atom1 atom2 
                !  1 2 3   A : 1 B : 2      A : 1 C : 3      B : 2 C : 3
                !
                ! If we filter atom3 by index3 > index2 when searching atom1 neighbors
                ! and index3 > index1 when searching atom2 neighbors (i.e., index3
                ! greater than the index of the atom whose neighbors are not searched),
                ! we get:
                !
                !  A B C   atom1 atom2  atom3 as atom1 nbor / atom2 nbor
                !  1 2 3   A : 1 B : 2  C : 3 -> found B-A-C
                !                       C : 3 -> found A-B-C
                !          A : 1 C : 3  B : 2 -> B < C (2 < 3) so ignored
                !                       B : 2 -> found A-C-B
                !          B : 2 C : 3  A : 1 -> A < B (1 < 2) so ignored
                !                       A : 1 -> A < C (1 < 3) so ignored

                if(many_bodies_found)then

                   many_bodies_found = .false.

                   !*********************!
                   ! 3-body interactions !
                   !*********************!

                   ! neighbors of atom2
                   nbors2 = atom2%neighbor_list

                   ! First we try to find ordered triplets atom2 -- atom1 -- atom3
                   ! Therefore we need separations a2--a1 and a1--a3.
                   separations(1:3,1) = -separations(1:3,1)
                   directions(1:3,1) = -directions(1:3,1)

                   ! loop over neighbors atom 1
                   do l = 1, nbors1%n_neighbors
                      index3 = nbors1%neighbors(l)

                      ! the condition for finding each triplet once is such that
                      ! index 3 must be higher than the index of the atom whose
                      ! neighbors are NOT currently searched
                      ! For the offset check we need atom2->atom3 which equals
                      ! (atom1->atom3) - (atom1->atom2), the latter being stored in offset
                      tripleoffset = nbors1%pbc_offsets(1:3,l)-offset ! offset atom2 -> atom3
                      if(pick(index2,index3,tripleoffset))then

                         ! third atom of the triplet
                         atom3 = atoms(index3)
                         ! The list of atoms is passed to force evaluation routine
                         ! for further filtering.
                         ! This is triplet atom2 - atom1 - atom3, since we loop over
                         ! neighbors of atom1.
                         atom_list(1:3) = (/ atom2, atom1, atom3 /)

                         ! Calculate the separations and distances between the particles
                         ! starting from atom2: a2--a1, a1--a3
                         ! (atom2 -- atom1 is already known though from 2-body calculation)
                         call separation_vector(atom1%position, &
                              atom3%position, &
                              nbors1%pbc_offsets(1:3,j), &
                              cell, &
                              separations(1:3,2)) ! in Geometry.f90
                         distances(2) = .norm.(separations(1:3,2))
                         if(distances(2) == 0.d0)then
                            directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                         else
                            directions(1:3,2) = separations(1:3,2) / distances(2)
                         end if

                         if(distances(2) < atom1%max_potential_radius)then
                            call core_evaluate_local_triplet_B(atom_list(1:3), &
                                 index2, index1, index3, & ! atom2 - atom1 - atom3
                                 1, 3, & ! atom1 is 2nd in the triplet, so test for 1st and 3rd
                                 interaction_indices, &
                                 separations(1:3,1:2), directions(1:3,1:2), distances(1:2), &
                                 calculation_type, energy, temp_forces, temp_enegs, stress, &
                                 many_bodies_found)
                         end if

                         if(many_bodies_found)then

                            many_bodies_found = .false.

                            !*********************!
                            ! 4-body interactions !
                            !*********************!

                            ! We search for atomic chain quadruplets A-B-C-D where the
                            ! triplet A-B-C or B-C-D is the triplet considered above.
                            ! Similarly to the triplets, the different quadruplets are found
                            ! by searching the neighbors of the end atoms of the triplets.

                            ! We have the triplet  atom2 -- atom1 -- atom3, so searching for neighbors of
                            ! atom2 and atom3 gives quadruplets
                            ! atom4 -- atom2 -- atom1 -- atom3 and
                            ! atom2 -- atom1 -- atom3 -- atom4

                            ! first, atom4 -- atom2 -- atom1 -- atom3:

                            ! neighbors of atom3 (not yet needed, but later)
                            nbors3 = atom3%neighbor_list

                            ! we need separations a4-a2, a2-a1, a1-a3, the latter two are already
                            ! known but in the wrong place
                            separations(1:3,2:3) = separations(1:3,1:2)
                            directions(1:3,2:3) = directions(1:3,1:2)
                            distances(2:3) = distances(1:2)

                            ! loop over neighbors atom 2
                            do m = 1, nbors2%n_neighbors
                               index4 = nbors2%neighbors(m)

                               ! the condition for finding each triplet once is such that
                               ! index 4 must be higher than the index of the atom whose
                               ! neighbors are NOT currently searched
                               ! For the offset check we need atom3->atom4 which equals
                               ! (atom2->atom4) - (atom2->atom3), the latter being stored in tripleoffset
                               quadoffset = nbors2%pbc_offsets(1:3,m)-tripleoffset ! offset atom3 -> atom4
                               if(pick(index3,index4,quadoffset) .and. &
                                  index4 /= index1)then

                                  ! fourth atom of the quadruplet
                                  atom4 = atoms(index4)
                                  ! The list of atoms is passed to force evaluation routine
                                  ! for further filtering.
                                  atom_list(1:4) = (/ atom4, atom2, atom1, atom3 /)

                                  call separation_vector(atom4%position, &
                                       atom2%position, &
                                       -nbors2%pbc_offsets(1:3,m), &
                                       cell, &
                                       separations(1:3,1)) ! in Geometry.f90
                                  distances(1) = .norm.(separations(1:3,1))
                                  if(distances(1) == 0.d0)then
                                     directions(1:3,1) = (/ 0.d0, 0.d0, 0.d0 /)
                                  else
                                     directions(1:3,1) = separations(1:3,1) / distances(1)
                                  end if

                                  if(distances(1) < atom2%max_potential_radius)then
                                     call core_evaluate_local_quadruplet_B(atom_list(1:4), &
                                       index4, index2, index1, index3, & ! atom4 - atom2 - atom1 - atom3
                                       1, 2, 4, & ! atom1 is 3rd in the quadruplet
                                       interaction_indices, &
                                       separations(1:3,1:3), directions(1:3,1:3), distances(1:3), &
                                       calculation_type, energy, temp_forces, temp_enegs, stress, &
                                       many_bodies_found)
                                  end if

                               end if ! index4 > index3
                            end do ! m = 1, nbors2%n_neighbors

                            ! move the a2-a1, a1-a3 separations back to their original place
                            separations(1:3,1:2) = separations(1:3,2:3)
                            directions(1:3,1:2) = directions(1:3,2:3)
                            distances(1:2) = distances(2:3)



                            ! second, atom2 -- atom1 -- atom3 -- atom4:

                            ! we need separations a2-a1, a1-a3, a3-a4, the first two are already
                            ! known and now in the right place as well

                            ! loop over neighbors atom 3
                            do m = 1, nbors3%n_neighbors
                               index4 = nbors3%neighbors(m)

                               ! the condition for finding each triplet once is such that
                               ! index 4 must be higher than the index of the atom whose
                               ! neighbors are NOT currently searched
                               ! For the offset check we need atom2->atom4 which equals
                               ! (atom3->atom4) + (atom2->atom3), the latter being stored in tripleoffset
                               quadoffset = nbors3%pbc_offsets(1:3,m)+tripleoffset ! offset atom2 -> atom4
                               if(pick(index2,index4,quadoffset) .and. &
                                    index4 /= index1)then

                                  ! fourth atom of the quadruplet
                                  atom4 = atoms(index4)
                                  ! The list of atoms is passed to force evaluation routine
                                  ! for further filtering.
                                  atom_list(1:4) = (/ atom2, atom1, atom3, atom4 /)

                                  call separation_vector(atom3%position, &
                                       atom4%position, &
                                       nbors3%pbc_offsets(1:3,m), &
                                       cell, &
                                       separations(1:3,3)) ! in Geometry.f90
                                  distances(3) = .norm.(separations(1:3,3))
                                  if(distances(3) == 0.d0)then
                                     directions(1:3,3) = (/ 0.d0, 0.d0, 0.d0 /)
                                  else
                                     directions(1:3,3) = separations(1:3,3) / distances(3)
                                  end if

                                  if(distances(3) < atom3%max_potential_radius)then
                                     call core_evaluate_local_quadruplet_B(atom_list(1:4), &
                                          index2, index1, index3, index4, & ! atom2 - atom1 - atom3 - atom4
                                          1, 3, 4, & ! atom1 is 2nd in the quadruplet
                                          interaction_indices, &
                                          separations(1:3,1:3), directions(1:3,1:3), distances(1:3), &
                                          calculation_type, energy, temp_forces, temp_enegs, stress, &
                                          many_bodies_found)
                                  end if

                               end if ! index4 > index3
                            end do ! m = 1, nbors3%n_neighbors

                         end if ! many_bodies_found

                      end if ! index3 > index2

                   end do ! l = 1, nbors1%n_neighbors

                   many_bodies_found = .false.

                   ! Next we try to find ordered triplets atom1 -- atom2 -- atom3
                   ! Therefore we need separations a1--a2 and a2--a3.
                   separations(1:3,1) = -separations(1:3,1)
                   directions(1:3,1) = -directions(1:3,1)

                   ! loop over neighbors of atom 2
                   do l = 1, nbors2%n_neighbors
                      index3 = nbors2%neighbors(l)

                      ! the condition for finding each triplet once is such that
                      ! index 3 must be higher than the index of the atom whose
                      ! neighbors are NOT currently searched
                      ! For the offset check we need atom1->atom3 which equals
                      ! (atom1->atom2) + (atom2->atom3), the former being stored in offset
                      tripleoffset = nbors2%pbc_offsets(1:3,l)+offset ! offset atom1 -> atom3
                      if(pick(index1,index3,tripleoffset))then

                         ! third atom of the triplet
                         atom3 = atoms(index3)
                         ! The list of atoms is passed to force evaluation routine
                         ! for further filtering.
                         ! This is triplet atom1 - atom2 - atom3, since we loop over
                         ! neighbors of atom2.
                         atom_list(1:3) = (/ atom1, atom2, atom3 /)

                         ! Calculate the separations and distances between the particles
                         ! starting from atom1: a1--a2, a2--a3
                         ! (atom1 -- atom2 is already known though from 2-body calculation)
                         call separation_vector(atom2%position, &
                              atom3%position, &
                              nbors2%pbc_offsets(1:3,l), &
                              cell, &
                              separations(1:3,2)) ! in Geometry.f90
                         distances(2) = .norm.(separations(1:3,2))
                         if(distances(2) == 0.d0)then
                            directions(1:3,2) = (/ 0.d0, 0.d0, 0.d0 /)
                         else
                            directions(1:3,2) = separations(1:3,2)/distances(2)
                         end if

                         if(distances(2) < atom2%max_potential_radius)then
                            call core_evaluate_local_triplet_B(atom_list(1:3), &
                                 index1, index2, index3, & ! atom1 - atom2 - atom3
                                 2, 3, & ! atom1 is 1st in the triplet, so test for 2nd and 3rd
                                 interaction_indices, &
                                 separations(1:3,1:2), directions(1:3,1:2), distances(1:2), &
                                 calculation_type, energy, temp_forces, temp_enegs, stress, &
                                 many_bodies_found)
                         end if

                         if(many_bodies_found)then

                            many_bodies_found = .false.

                            !*********************!
                            ! 4-body interactions !
                            !*********************!

                            ! We search for atomic chain quadruplets A-B-C-D where the
                            ! triplet A-B-C or B-C-D is the triplet considered above.
                            ! Similarly to the triplets, the different quadruplets are found
                            ! by searching the neighbors of the end atoms of the triplets.

                            ! We have the triplet  atom1 -- atom2 -- atom3, so searching for neighbors of
                            ! atom1 and atom3 gives quadruplets
                            ! atom4 -- atom1 -- atom2 -- atom3 and
                            ! atom1 -- atom2 -- atom3 -- atom4

                            ! first, atom4 -- atom1 -- atom2 -- atom3:

                            ! neighbors of atom3
                            nbors3 = atom3%neighbor_list

                            ! we need separations a4-a1, a1-a2, a2-a3, the latter two are already
                            ! known but in the wrong place
                            separations(1:3,2:3) = separations(1:3,1:2)
                            directions(1:3,2:3) = directions(1:3,1:2)
                            distances(2:3) = distances(1:2)

                            ! loop over neighbors atom 1
                            do m = 1, nbors1%n_neighbors
                               index4 = nbors1%neighbors(m)

                               ! the condition for finding each triplet once is such that
                               ! index 4 must be higher than the index of the atom whose
                               ! neighbors are NOT currently searched
                               ! For the offset check we need atom3->atom4 which equals
                               ! (atom1->atom4) - (atom1->atom3), the latter being stored in tripleoffset
                               quadoffset = nbors1%pbc_offsets(1:3,m)-tripleoffset ! offset atom3 -> atom4
                               if(pick(index3,index4,quadoffset) .and. &
                                    index4 /= index2)then

                                  ! fourth atom of the quadruplet
                                  atom4 = atoms(index4)
                                  ! The list of atoms is passed to force evaluation routine
                                  ! for further filtering.
                                  atom_list(1:4) = (/ atom4, atom1, atom2, atom3 /)

                                  call separation_vector(atom4%position, &
                                       atom1%position, &
                                       -nbors1%pbc_offsets(1:3,m), &
                                       cell, &
                                       separations(1:3,1)) ! in Geometry.f90
                                  distances(1) = .norm.(separations(1:3,1))
                                  if(distances(1) == 0.d0)then
                                     directions(1:3,1) = (/ 0.d0, 0.d0, 0.d0 /)
                                  else
                                     directions(1:3,1) = separations(1:3,1) / distances(1)
                                  end if


                                  if(distances(1) < atom1%max_potential_radius)then
                                     call core_evaluate_local_quadruplet_B(atom_list(1:4), &
                                          index4, index1, index2, index3, & ! atom4 - atom1 - atom2 - atom3
                                          1, 3, 4, & ! atom1 is 2nd in the quadruplet
                                          interaction_indices, &
                                          separations(1:3,1:3), directions(1:3,1:3), distances(1:3), &
                                          calculation_type, energy, temp_forces, temp_enegs, stress, &
                                          many_bodies_found)
                                  end if

                               end if ! index4 > index3
                            end do ! m = 1, nbors2%n_neighbors

                            ! move the a1-a2, a2-a3 separations back to their original place
                            separations(1:3,1:2) = separations(1:3,2:3)
                            directions(1:3,1:2) = directions(1:3,2:3)
                            distances(1:2) = distances(2:3)



                            ! second, atom1 -- atom2 -- atom3 -- atom4:

                            ! we need separations a1-a2, a2-a3, a3-a4, the first two are already
                            ! known and now in the right place as well

                            ! loop over neighbors atom 3
                            do m = 1, nbors3%n_neighbors
                               index4 = nbors3%neighbors(m)

                               ! the condition for finding each triplet once is such that
                               ! index 4 must be higher than the index of the atom whose
                               ! neighbors are NOT currently searched
                               ! For the offset check we need atom1->atom4 which equals
                               ! (atom3->atom4) + (atom1->atom3), the latter being stored in tripleoffset
                               quadoffset = nbors3%pbc_offsets(1:3,m)+tripleoffset ! offset atom1 -> atom4
                               if(pick(index1,index4,quadoffset) .and. &
                                    index4 /= index2)then

                                  ! fourth atom of the quadruplet
                                  atom4 = atoms(index4)
                                  ! The list of atoms is passed to force evaluation routine
                                  ! for further filtering.
                                  atom_list(1:4) = (/ atom1, atom2, atom3, atom4 /)

                                  call separation_vector(atom3%position, &
                                       atom4%position, &
                                       nbors3%pbc_offsets(1:3,m), &
                                       cell, &
                                       separations(1:3,3)) ! in Geometry.f90
                                  distances(3) = .norm.(separations(1:3,3))
                                  if(distances(3) == 0.d0)then
                                     directions(1:3,3) = (/ 0.d0, 0.d0, 0.d0 /)
                                  else
                                     directions(1:3,3) = separations(1:3,3) / distances(3)
                                  end if


                                  if(distances(3) < atom3%max_potential_radius)then
                                     call core_evaluate_local_quadruplet_B(atom_list(1:4), &
                                          index1, index2, index3, index4, & ! atom1 - atom2 - atom3 - atom4
                                          2, 3, 4, & ! atom1 is 1st in the quadruplet
                                          interaction_indices, &
                                          separations(1:3,1:3), directions(1:3,1:3), distances(1:3), &
                                          calculation_type, energy, temp_forces, temp_enegs, stress, &
                                          many_bodies_found)
                                  end if

                               end if ! index4 > index3
                            end do ! m = 1, nbors3%n_neighbors

                         end if ! many_bodies_found

                      end if ! index3 > index1

                   end do ! l

                end if ! many-bodies_found

             end if ! index2 > index1

          end do ! j

       end if ! is_my_atom
    end do ! index1

    ! Stop the load timer
    call timer(stopwatch_0)

    ! differentiate between energy, force, and electronegativity evaluation
    select case(calculation_type)
    case(energy_evaluation_index)

       !*****************************!
       ! energy collecting and ewald !
       !*****************************!

       call mpi_wall_clock(t00)
#ifdef MPI
       ! In MPI, calculate the loads for all cpus and try to balance the loads
       call record_load(stopwatch_0)
       call balance_loads()

       ! collect data from all cpus in MPI
       call mpi_allreduce(energy,total_energy,1,mpi_double_precision,mpi_sum,&
            mpi_comm_world,mpistat)
#else
       total_energy = energy
#endif
       call mpi_wall_clock(t0)
       t5 = t5+t0-t00
       t00 = t0

       if(evaluate_ewald)then
          call calculate_ewald_energy(atoms,cell,ewald_cutoff,ewald_k_radius,ewald_k_cutoffs,ewald_sigma,&
               ewald_epsilon,ewald_scaler,.false.,energy)

          total_energy = total_energy + energy
       call mpi_wall_clock(t0)
       t6 = t6+t0-t00

       end if

    case(force_evaluation_index)

       !****************************!
       ! force collecting and ewald !
       !****************************!

       call mpi_wall_clock(t00)
#ifdef MPI
       ! In MPI, calculate the loads for all cpus and try to balance the loads
       call record_load(stopwatch_0)
       call balance_loads()

       ! collect data from all cpus in MPI
       call mpi_allreduce(temp_forces,total_forces,size(temp_forces),mpi_double_precision,&
            mpi_sum,mpi_comm_world,mpistat)
       call mpi_allreduce(stress,total_stress,size(stress),mpi_double_precision,&
            mpi_sum,mpi_comm_world,mpistat)
#else
       total_forces = temp_forces
       total_stress = stress
#endif
       call mpi_wall_clock(t0)
       t5 = t5+t0-t00
       t00 = t0

       if(evaluate_ewald)then
          call calculate_ewald_forces(atoms,cell,ewald_cutoff,ewald_k_radius,ewald_k_cutoffs,ewald_sigma,&
               ewald_epsilon,ewald_scaler,.false.,temp_forces,stress)
          total_forces = total_forces + temp_forces
          total_stress = total_stress + stress
       end if
       call mpi_wall_clock(t0)
       t6 = t6+t0-t00

    case(electronegativity_evaluation_index)

       !****************************************!
       ! electronegativity collecting and ewald !
       !****************************************!

#ifdef MPI
       ! In MPI, calculate the loads for all cpus and try to balance the loads
       call record_load(stopwatch_0)
       call balance_loads()

       ! collect data from all cpus in MPI
       call mpi_allreduce(temp_enegs,total_enegs,size(temp_enegs),mpi_double_precision,&
            mpi_sum,mpi_comm_world,mpistat)
#else
       total_enegs = temp_enegs
#endif

       ! ewald summation
       if(evaluate_ewald)then
          temp_enegs = 0.d0
          call calculate_ewald_electronegativities(atoms,cell,ewald_cutoff,ewald_k_cutoffs,ewald_sigma,&
               ewald_epsilon,ewald_scaler,.false.,temp_enegs)
          total_enegs = total_enegs + temp_enegs
       end if

    end select

    ! Empty the bond order factor storage and stop searching them from memory.
    ! This is done so that if the geometry changes due to atoms moving, for instance,
    ! then the obsolete factors are not used in error.
    use_saved_bond_order_factors = .false.
    call core_empty_bond_order_storage()

    call mpi_wall_clock(tb)
    t9 = t1+t2+t3+t4+t5+t6+t7
    t8 = tb-ta
    write(*,*) ""
    write(*,'(A)') "evaluation timing"
    write(*,'(A, F7.3, F7.1, A)') "bond_factors: ", t1,  100*t1/t8, " %"
    write(*,'(A, F7.3, F7.1, A)') "singlets:     ", t2,  100*t2/t8, " %"
    write(*,'(A, F7.3, F7.1, A)') "distances:    ", t3,  100*t3/t8, " %"
    write(*,'(A, F7.3, F7.1, A)') "doublets:     ", t4,  100*t4/t8, " %"
    write(*,'(A, F7.3, F7.1, A)') "mpi:          ", t5,  100*t5/t8, " %"
    write(*,'(A, F7.3, F7.1, A)') "ewald:        ", t6,  100*t6/t8, " %"
    write(*,'(A, F7.3, F7.1, A)') "clear:        ", t7,  100*t7/t8, " %"
    write(*,'(A, F7.3)') "sum:          ", t9
    write(*,'(A, F7.3)') "total:        ", t8
    write(*,*) ""

  end subroutine core_loop_over_local_interactions

  ! Evaluates the local potential affecting a single atom
  !
  ! *index1 index of the atom
  ! *atom_singlet the atom that is targeted
  ! *interaction_indices the interactions targeting the given atom
  ! *calculation_type specifies if we are evaluating the energy, forces, or electronegativities
  ! *energy calculated energy
  ! *forces calculated forces
  ! *enegs calculated electronegativities
  subroutine core_evaluate_local_singlet(index1, &
       atom_singlet, &
       interaction_indices, &
       calculation_type,energy,forces,stress,enegs)
    implicit none
    integer, intent(in) :: calculation_type, index1
    integer, pointer :: interaction_indices(:)
    double precision, intent(inout) :: energy, forces(:,:), enegs(:), stress(6)
    type(atom), intent(in) :: atom_singlet
    integer :: k, n_targets, n_atoms, index2
    type(potential) :: interaction
    double precision, save :: &
         tmp_energy, tmp_forces(3,1), tmp_enegs(1), &
         dummy_sep(3,0), dummy_dist(0), bo_virial(6)

    n_atoms = size(atoms)

    ! loop over potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))

       ! filter the potentials according to number of targets
       call get_number_of_targets_of_potential_index(interaction%type_index,n_targets)
       if( n_targets == 1 )then

          ! differentiate between energy, force, and electronegativity evaluation
          select case(calculation_type)
          case(energy_evaluation_index)

             !***********************!
             ! 1-body energy (atom1) !
             !***********************!

             ! evaluate the 1-body energy involving atom1
             call evaluate_energy(1,interaction%n_product,dummy_sep,dummy_dist,&
                  interaction,tmp_energy,atoms(index1:index1)) ! in Potentials.f90


             ! If there is a bond order factor associated with the potential,
             ! we add the contribution is brings:
             !
             ! V = \sum_i b_i v_i
             if(interaction%pot_index > -1)then
                call core_get_bond_order_factors(interaction%pot_index,&
                     bo_factors)
                     
                 energy = energy + tmp_energy*bo_factors(index1)
             else
                !bo_factors = 1.d0                
                 energy = energy + tmp_energy
             end if

          case(force_evaluation_index)

             !**********************!
             ! 1-body force (atom1) !
             !**********************!

             ! evaluate the 1-body force involving atom1
             call evaluate_forces(1,interaction%n_product,dummy_sep,dummy_dist,&
                  interaction,tmp_forces(1:3,1),atoms(index1:index1)) ! in Potentials.f90

             ! If there is a bond order factor associated with the potential,
             ! we add the contribution is brings:
             !
             ! V = \sum_i b_i v_i
             ! F_a = - \nabla_a V 
             !     = - \sum_i (\nabla_a b_i) v_i + b_i (\nabla_a v_i)
             !     = - \sum_i (\nabla_a b_i) v_i + b_i f_a,i
             !
             if(interaction%pot_index > -1)then
                call core_get_bond_order_factors(interaction%pot_index,&
                     bo_factors)
                call core_get_bond_order_gradients(interaction%pot_index,&
                     index1,& ! atom index
                     1, & ! slot_index
                     bo_gradients(1:3,1:n_atoms,1),&
                     bo_virial)

                ! Add the bond order gradient terms involving the atom1 self energy for all atoms.
                ! That is, add the (\nabla_a b_i) v_i term with the given i (atom1) for all a.
                call evaluate_energy(1,interaction%n_product,dummy_sep,dummy_dist,interaction,&
                     tmp_energy,atoms(index1:index1))  ! in Potentials.f90
!!$                forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) - &
!!$                     tmp_energy*bo_gradients(1:3,1:n_atoms,1)
                do index2 = 1, n_atoms
                   if(bo_scaling(index2))then
                      forces(1:3,index2) = forces(1:3,index2) - tmp_energy*bo_gradients(1:3,index2,1)
                   end if
                end do

                ! Add the contribution of bond order gradients to the stress tensor
                stress(1:6) = stress(1:6) - tmp_energy*bo_virial(1:6)

                 ! Add the force due to potential gradient             
                 forces(1:3,index1) = forces(1:3,index1) + &
                      tmp_forces(1:3,1)*bo_factors(index1)

             else

                 ! Add the force due to potential gradient             
                 forces(1:3,index1) = forces(1:3,index1) + &
                      tmp_forces(1:3,1)
             end if


          case(electronegativity_evaluation_index)

             !**********************************!
             ! 1-body electronegativity (atom1) !
             !**********************************!

             ! evaluate the 1-body energy involving atom1
             call evaluate_electronegativity(1,interaction%n_product,dummy_sep,dummy_dist,&
                  interaction,tmp_enegs(1),atoms(index1:index1)) ! in Potentials.f90

             ! If there is a bond order factor associated with the potential,
             ! we add the contribution is brings:
             if(interaction%pot_index > -1)then
                call core_get_bond_order_factors(interaction%pot_index,&
                     bo_factors)

                 ! Add the electronegativity due to potential gradient
                 enegs(index1) = enegs(index1) + tmp_enegs(1)*bo_factors(index1)
             else

                 ! Add the electronegativity due to potential gradient
                 enegs(index1) = enegs(index1) + tmp_enegs(1)
             end if


          end select

       end if
    end do


  end subroutine core_evaluate_local_singlet


 


  subroutine core_evaluate_local_doublet_energy_B(atom_doublet, &
       index1, index2, &
       test_index1, &
       interaction_indices, &
       separations, directions, distances, &
       energy, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: index1, index2, test_index1
    integer, pointer :: interaction_indices(:)
    double precision, intent(inout) :: energy
    type(atom), intent(in) :: atom_doublet(2)
    double precision, intent(in) :: separations(3,1), directions(3,1), distances(1)
    logical, intent(out) :: many_bodies_found
    
    type(atom) :: atom1, atom2
    integer :: k, n_targets, index_pair(2)
    type(potential) :: interaction
    double precision :: &
         tmp_energy, &
         cut_factors(1), pair_bo_factors(2), pair_bo_sums(2)
    logical :: is_active

    many_bodies_found = .false.
    atom1 = atom_doublet(1)
    atom2 = atom_doublet(2)
    index_pair = (/ index1, index2 /)

    ! loop over potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is it a 2-body potential
       if( interaction%cutoff > distances(1) )then 
          call potential_affects_atom(interaction,atom_doublet(test_index1),is_active,2) ! in Potentials.f90
          if( is_active )then 
             call get_number_of_targets_of_potential_index(interaction%type_index,&
                  n_targets) ! in Potentials.f90
             if( n_targets == 2 )then
                
                
                !*****************************!
                ! 2-body energy (atom1-atom2) !
                !*****************************!

                ! If there is a bond order factor associated with the potential,
                ! we add the contribution is brings:
                !
                ! V = \sum_ij b_ij v_ij
                ! b_ij = (b_i + b_j) / 2 + B_ij
                if(interaction%pot_index > -1)then
                   ! get b_i (for all i, they have been precalculated)
                   ! call core_get_bond_order_factors(n_atoms,&
                   !     interaction%pot_index,&
                   !     bo_factors)


                   ! get B_ij for this ij (not precalculated)
                   ! this is 
                   ! B_ij = (B*_ij + B*_ji) / 2
                   call core_calculate_pair_bond_order_factor(index_pair, &
                        separations(1:3,1), distances(1), directions(1:3,1), &
                        interaction%pot_index, &
                        pair_bo_sums(1:2))
                   call core_post_process_pair_bond_order_factor(atom1, &
                        interaction%pot_index, &
                        pair_bo_sums(1), &
                        pair_bo_factors(1))
                   call core_post_process_pair_bond_order_factor(atom2, &
                        interaction%pot_index, &
                        pair_bo_sums(2), &
                        pair_bo_factors(2))

                else

                end if

                ! If a smooth cutoff is present, we add the
                ! contribution it brings:
                ! 
                ! V = \sum_ij v_ij f(r_ij)
                if(interaction%smoothened)then
                   ! get f(r_ij)
                   call smoothening_factor(distances(1),&
                        interaction%cutoff,interaction%soft_cutoff,&
                        cut_factors(1)) ! in Potentials.f90
                else
                   cut_factors(1) = 1.d0
                end if

                ! evaluate the 2-body energy involving atom1-atom2 interaction
                call evaluate_energy(2,interaction%n_product,separations(1:3,1),distances(1),&
                     interaction,tmp_energy,atom_doublet)  ! in Potentials.f90

                ! add the term: b_ij v_ij f(rij)
                if(interaction%pot_index > -1)then
                    energy = energy + tmp_energy*cut_factors(1)*&
                         (saved_bond_order_factors(index1,interaction%pot_index)+&
                         saved_bond_order_factors(index2,interaction%pot_index)+&
                         pair_bo_factors(1)+pair_bo_factors(2))*0.5d0
                else
                    energy = energy + tmp_energy*cut_factors(1)
                end if

          else if(n_targets > 2)then

             ! If the number of targets is greater than 2,
             ! we have found a many-body potential.
             ! Make a note that we must also evaluate the many-body terms.
             many_bodies_found = .true.

          end if ! n_targets == 2

       end if ! is_active
       end if ! cutoff
    end do ! k

  end subroutine core_evaluate_local_doublet_energy_B



  subroutine core_evaluate_local_doublet_energy(n_atoms, &
       atom_doublet, &
       index1, index2, &
       test_index1, &
       interaction_indices, &
       separations, directions, distances, &
       energy, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: n_atoms, index1, index2, test_index1
    integer, pointer :: interaction_indices(:)
    double precision, intent(inout) :: energy
    type(atom), intent(in) :: atom_doublet(2)
    double precision, intent(in) :: separations(3,1), directions(3,1), distances(1)
    logical, intent(out) :: many_bodies_found
    
    type(atom) :: atom1, atom2
    integer :: k, n_targets, index_pair(2)
    type(potential) :: interaction
    double precision :: bo_factors(n_atoms), bo_sums(n_atoms), &
         tmp_energy, &
         cut_factors(1), pair_bo_factors(2), pair_bo_sums(2)
    logical :: is_active

    many_bodies_found = .false.
    atom1 = atom_doublet(1)
    atom2 = atom_doublet(2)
    index_pair = (/ index1, index2 /)

    ! loop over potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is it a 2-body potential
       if( interaction%cutoff > distances(1) )then 
       call potential_affects_atom(interaction,atom_doublet(test_index1),is_active,2) ! in Potentials.f90
       if( is_active )then 
          call get_number_of_targets_of_potential_index(interaction%type_index,&
               n_targets) ! in Potentials.f90
          if( n_targets == 2 )then


                !*****************************!
                ! 2-body energy (atom1-atom2) !
                !*****************************!

                ! If there is a bond order factor associated with the potential,
                ! we add the contribution is brings:
                !
                ! V = \sum_ij b_ij v_ij
                ! b_ij = (b_i + b_j) / 2 + B_ij
                if(interaction%pot_index > -1)then
                   ! get b_i (for all i, they have been precalculated)
                   call core_get_bond_order_factors(interaction%pot_index,&
                        bo_factors)

                   ! get B_ij for this ij (not precalculated)
                   ! this is 
                   ! B_ij = (B*_ij + B*_ji) / 2
                   call core_calculate_pair_bond_order_factor(index_pair, &
                        separations(1:3,1), distances(1), directions(1:3,1), &
                        interaction%pot_index, &
                        pair_bo_sums(1:2))
                   call core_post_process_pair_bond_order_factor(atom1, &
                        interaction%pot_index, &
                        pair_bo_sums(1), &
                        pair_bo_factors(1))
                   call core_post_process_pair_bond_order_factor(atom2, &
                        interaction%pot_index, &
                        pair_bo_sums(2), &
                        pair_bo_factors(2))

                else
                   !bo_factors = 1.d0
                end if

                ! If a smooth cutoff is present, we add the
                ! contribution it brings:
                ! 
                ! V = \sum_ij v_ij f(r_ij)
                if(interaction%smoothened)then
                   ! get f(r_ij)
                   call smoothening_factor(distances(1),&
                        interaction%cutoff,interaction%soft_cutoff,&
                        cut_factors(1)) ! in Potentials.f90
                else
                   cut_factors(1) = 1.d0
                end if

                ! evaluate the 2-body energy involving atom1-atom2 interaction
                call evaluate_energy(2,interaction%n_product,separations(1:3,1),distances(1),&
                     interaction,tmp_energy,atom_doublet)  ! in Potentials.f90

                ! add the term: b_ij v_ij f(rij)
                if(interaction%pot_index > -1)then
                    energy = energy + tmp_energy*cut_factors(1)*&
                         (bo_factors(index1)+bo_factors(index2)+pair_bo_factors(1)+pair_bo_factors(2))*0.5d0
                else
                    energy = energy + tmp_energy*cut_factors(1)
                end if

          else if(n_targets > 2)then

             ! If the number of targets is greater than 2,
             ! we have found a many-body potential.
             ! Make a note that we must also evaluate the many-body terms.
             many_bodies_found = .true.

          end if ! n_targets == 2

       end if ! is_active
       end if ! cutoff
    end do ! k

  end subroutine core_evaluate_local_doublet_energy



  subroutine core_evaluate_local_doublet_forces(n_atoms, &
       atom_doublet, &
       index1, index2, &
       test_index1, &
       interaction_indices, &
       separations, directions, distances, &
       forces,stress, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: n_atoms, index1, index2, test_index1
    integer, pointer :: interaction_indices(:)
    double precision, intent(inout) :: forces(3,n_atoms), stress(6)
    type(atom), intent(in) :: atom_doublet(2)
    double precision, intent(in) :: separations(3,1), directions(3,1), distances(1)
    logical, intent(out) :: many_bodies_found
    
    type(atom) :: atom1, atom2
    integer :: k, n_targets, index_pair(2)
    type(potential) :: interaction
    double precision :: bo_factors(n_atoms), bo_sums(n_atoms), bo_gradients(3,n_atoms,4), &
         tmp_energy, tmp_forces(3,2), tmp_enegs(2), &
         cut_factors(1), cut_gradients(3,1), &
         pair_forces(3,2), bo_virial(6,2), &
         pair_bo_sums(2), pair_bo_factors(2), &
         pair_bo_virial(6,2)
    logical :: is_active

    many_bodies_found = .false.
    atom1 = atom_doublet(1)
    atom2 = atom_doublet(2)
    index_pair = (/ index1, index2 /)
    
    ! loop over potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is it a 2-body potential
       if( interaction%cutoff > distances(1) )then 
       call potential_affects_atom(interaction,atom_doublet(test_index1),is_active,2) ! in Potentials.f90
       if( is_active )then 
          call get_number_of_targets_of_potential_index(interaction%type_index,&
               n_targets) ! in Potentials.f90
          if( n_targets == 2 )then

             !****************************!
             ! 2-body force (atom1-atom2) !
             !****************************!
             
             ! We will need the energy contribution from atom1-atom2
             ! interaction if smooth cutoffs or bond factors are used,
             ! since we are mulplying the potential.
             if((interaction%pot_index > -1) .or. &
                  interaction%smoothened)then
                call evaluate_energy(2,interaction%n_product,separations(1:3,1),distances(1),&
                     interaction,tmp_energy,atom_doublet) ! in Potentials.f90
             else
                tmp_energy = 0.d0
             end if
             
             ! If a smooth cutoff is present, we add the
             ! contribution it brings:
             ! 
             ! V = \sum_ij v_ij f(r_ij)
             ! F_a = - \nabla_a V 
             !     = - \sum_ij v_ij f'(r_ij) (\nabla_a r_ij) + (\nabla_a v_ij) f(r_ij) 
             !     = - \sum_ij v_ij f'(r_ij) (\nabla_a r_ij) + f_a,ij f(r_ij)
             !
             if(interaction%smoothened)then
                ! get f(r_ij)
                call smoothening_factor(distances(1),&
                     interaction%cutoff,interaction%soft_cutoff,&
                     cut_factors(1)) ! in Potentials.f90
                ! get f'(r_ij) (\nabla_a r_ij)
                call smoothening_gradient(directions(1:3,1),distances(1),&
                     interaction%cutoff,interaction%soft_cutoff,&
                     cut_gradients(1:3,1)) ! in Potentials.f90
             else
                cut_factors(1) = 1.d0
                cut_gradients(1:3,1) = 0.d0
             end if
             
             ! If there is a bond order factor associated with the potential,
             ! we add the contribution is brings:
             !
             ! V = \sum_ij b_ij v_ij
             ! b_ij = (b_i + b_j) / 2 + B_ij
             ! F_a = - \nabla_a V 
             !     = - \sum_ij (\nabla_a b_ij) v_ij + b_ij (\nabla_a v_ij)
             !     = - \sum_ij (\nabla_a b_ij) v_ij + b_ij f_a,ij
             !
             if(interaction%pot_index > -1)then
                ! get b_i (for all i, they have been precalculated)
                call core_get_bond_order_factors(interaction%pot_index,&
                     bo_factors)
                ! get (\nabla_a b_i) (for all a)
                call core_get_bond_order_gradients(interaction%pot_index,&
                     index1,& ! atom index
                     1, & ! slot_index
                     bo_gradients(1:3,1:n_atoms,1), &
                     bo_virial(1:6,1))
                ! get (\nabla_a b_j) (for all a)
                call core_get_bond_order_gradients(interaction%pot_index,&
                     index2,& ! atom index
                     2, & ! slot_index
                     bo_gradients(1:3,1:n_atoms,2), &
                     bo_virial(1:6,2))
                
                
                ! get B_ij and \nabla_a B_ij for this ij (not precalculated)
                ! this is 
                ! B_ij = (B*_ij + B*_ji) / 2
                call core_calculate_pair_bond_order_factor(index_pair, &
                     separations(1:3,1), distances(1), directions(1:3,1), &
                     interaction%pot_index, &
                     pair_bo_sums(1:2))
                call core_post_process_pair_bond_order_factor(atom1, &
                     interaction%pot_index, &
                     pair_bo_sums(1), &
                     pair_bo_factors(1))
                call core_post_process_pair_bond_order_factor(atom2, &
                     interaction%pot_index, &
                     pair_bo_sums(2), &
                     pair_bo_factors(2))
                
                call core_calculate_pair_bond_order_gradients(index_pair, &
                     separations(1:3,1), distances(1), directions(1:3,1), &
                     interaction%pot_index, &                        
                     pair_bo_sums(1:2), &
                     bo_gradients(1:3,1:n_atoms,3:4), &
                     pair_bo_virial(1:6,1:2))
                
                
                ! Add the bond order gradient terms involving the 
                ! atom1-atom2 energy for all atoms.
                ! That is, add the (\nabla_a b_ij) v_ij term with 
                ! the given ij (atom1,atom2) for all a.
                forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) &
                     - tmp_energy*cut_factors(1)*&
                     (bo_gradients(1:3,1:n_atoms,1) + bo_gradients(1:3,1:n_atoms,2) + &
                     bo_gradients(1:3,1:n_atoms,3) + bo_gradients(1:3,1:n_atoms,4))*0.5d0
                
                ! Add the contribution of bond order gradients to the stress tensor
                stress(1:6) = stress(1:6) - tmp_energy*cut_factors(1)* &
                     (bo_virial(1:6,1) + bo_virial(1:6,2) + pair_bo_virial(1:6,1) + pair_bo_virial(1:6,2))*0.5d0
                
             else
                !bo_factors = 1.d0
                !bo_sums = 0.d0
                !bo_gradients = 0.d0
             end if
             
             ! evaluate the 2-body force involving atom1-atom2 interaction
             call evaluate_forces(2,interaction%n_product,separations(1:3,1),distances(1),&
                  interaction,tmp_forces(1:3,1:2),atom_doublet) ! in Potentials.f90
             
             if(interaction%pot_index > -1)then
                ! force on atom 1:                 
                pair_forces(1:3,1) = ( tmp_forces(1:3,1) * cut_factors(1) + &
                     tmp_energy * cut_gradients(1:3,1) ) * &
                     ( bo_factors(index1) + bo_factors(index2) + &
                     pair_bo_factors(1) + pair_bo_factors(2) ) * 0.5d0
                
                ! force on atom 2:
                pair_forces(1:3,2) = ( tmp_forces(1:3,2) * cut_factors(1) - &
                     tmp_energy * cut_gradients(1:3,1) ) * &
                     ( bo_factors(index1) + bo_factors(index2) + &
                     pair_bo_factors(1) + pair_bo_factors(2) ) * 0.5d0
             else
                pair_forces(1:3,1) = ( tmp_forces(1:3,1) * cut_factors(1) + &
                     tmp_energy * cut_gradients(1:3,1) )
                pair_forces(1:3,2) = ( tmp_forces(1:3,2) * cut_factors(1) - &
                     tmp_energy * cut_gradients(1:3,1) )
             end if
             
             forces(1:3,index1) = forces(1:3,index1) + pair_forces(1:3,1)
             forces(1:3,index2) = forces(1:3,index2) + pair_forces(1:3,2)
             
             !***************!
             ! stress tensor !
             !***************!
             
             ! s_xx, s_yy, s_zz, s_yz, s_xz, s_xy:
             stress(1) = stress(1) + separations(1,1) * pair_forces(1,2)
             stress(2) = stress(2) + separations(2,1) * pair_forces(2,2)
             stress(3) = stress(3) + separations(3,1) * pair_forces(3,2)
             stress(4) = stress(4) + separations(2,1) * pair_forces(3,2)
             stress(5) = stress(5) + separations(1,1) * pair_forces(3,2)
             stress(6) = stress(6) + separations(1,1) * pair_forces(2,2)
             
          else if(n_targets > 2)then

             ! If the number of targets is greater than 2,
             ! we have found a many-body potential.
             ! Make a note that we must also evaluate the many-body terms.
             many_bodies_found = .true.

          end if ! n_targets == 2

       end if ! is_active
       end if ! cutoff
    end do ! k

  end subroutine core_evaluate_local_doublet_forces




  subroutine core_evaluate_local_doublet_forces_B( &
       atom_doublet, &
       index1, index2, &
       test_index1, &
       interaction_indices, &
       separations, directions, distances, &
       forces,stress, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: index1, index2, test_index1
    integer, pointer :: interaction_indices(:)
    double precision, intent(inout) :: forces(:,:), stress(6)
    type(atom), intent(in) :: atom_doublet(2)
    double precision, intent(in) :: separations(3,1), directions(3,1), distances(1)
    logical, intent(out) :: many_bodies_found
    
    type(atom) :: atom1, atom2
    integer :: k, n_targets, index_pair(2), n_atoms
    type(potential) :: interaction
    double precision :: &
         tmp_energy, tmp_forces(3,2), tmp_enegs(2), &
         cut_factors(1), cut_gradients(3,1), &
         pair_forces(3,2), bo_virial(6,2), &
         pair_bo_sums(2), pair_bo_factors(2), &
         pair_bo_virial(6,2)
    logical :: is_active

    n_atoms = size(atoms)
    many_bodies_found = .false.
    atom1 = atom_doublet(1)
    atom2 = atom_doublet(2)
    index_pair = (/ index1, index2 /)
    
    ! loop over potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is it a 2-body potential
       if( interaction%cutoff > distances(1) )then 
       call potential_affects_atom(interaction,atom_doublet(test_index1),is_active,2) ! in Potentials.f90
       if( is_active )then 
          call get_number_of_targets_of_potential_index(interaction%type_index,&
               n_targets) ! in Potentials.f90
          if( n_targets == 2 )then

             !****************************!
             ! 2-body force (atom1-atom2) !
             !****************************!
             
             ! We will need the energy contribution from atom1-atom2
             ! interaction if smooth cutoffs or bond factors are used,
             ! since we are mulplying the potential.
             if((interaction%pot_index > -1) .or. &
                  interaction%smoothened)then
                call evaluate_energy(2,interaction%n_product,separations(1:3,1),distances(1),&
                     interaction,tmp_energy,atom_doublet) ! in Potentials.f90
             else
                tmp_energy = 0.d0
             end if
             
             ! If a smooth cutoff is present, we add the
             ! contribution it brings:
             ! 
             ! V = \sum_ij v_ij f(r_ij)
             ! F_a = - \nabla_a V 
             !     = - \sum_ij v_ij f'(r_ij) (\nabla_a r_ij) + (\nabla_a v_ij) f(r_ij) 
             !     = - \sum_ij v_ij f'(r_ij) (\nabla_a r_ij) + f_a,ij f(r_ij)
             !
             if(interaction%smoothened)then
                ! get f(r_ij)
                call smoothening_factor(distances(1),&
                     interaction%cutoff,interaction%soft_cutoff,&
                     cut_factors(1)) ! in Potentials.f90
                ! get f'(r_ij) (\nabla_a r_ij)
                call smoothening_gradient(directions(1:3,1),distances(1),&
                     interaction%cutoff,interaction%soft_cutoff,&
                     cut_gradients(1:3,1)) ! in Potentials.f90
             else
                cut_factors(1) = 1.d0
                cut_gradients(1:3,1) = 0.d0
             end if
             
             ! If there is a bond order factor associated with the potential,
             ! we add the contribution is brings:
             !
             ! V = \sum_ij b_ij v_ij
             ! b_ij = (b_i + b_j) / 2 + B_ij
             ! F_a = - \nabla_a V 
             !     = - \sum_ij (\nabla_a b_ij) v_ij + b_ij (\nabla_a v_ij)
             !     = - \sum_ij (\nabla_a b_ij) v_ij + b_ij f_a,ij
             !
             if(interaction%pot_index > -1)then
                ! get b_i (for all i, they have been precalculated)
                call core_get_bond_order_factors(interaction%pot_index,&
                     bo_factors)
                ! get (\nabla_a b_i) (for all a)
                call core_get_bond_order_gradients(interaction%pot_index,&
                     index1,& ! atom index
                     1, & ! slot_index
                     bo_gradients(1:3,1:n_atoms,1), &
                     bo_virial(1:6,1))
                ! get (\nabla_a b_j) (for all a)
                call core_get_bond_order_gradients(interaction%pot_index,&
                     index2,& ! atom index
                     2, & ! slot_index
                     bo_gradients(1:3,1:n_atoms,2), &
                     bo_virial(1:6,2))
                
                
                ! get B_ij and \nabla_a B_ij for this ij (not precalculated)
                ! this is 
                ! B_ij = (B*_ij + B*_ji) / 2
                call core_calculate_pair_bond_order_factor(index_pair, &
                     separations(1:3,1), distances(1), directions(1:3,1), &
                     interaction%pot_index, &
                     pair_bo_sums(1:2))
                call core_post_process_pair_bond_order_factor(atom1, &
                     interaction%pot_index, &
                     pair_bo_sums(1), &
                     pair_bo_factors(1))
                call core_post_process_pair_bond_order_factor(atom2, &
                     interaction%pot_index, &
                     pair_bo_sums(2), &
                     pair_bo_factors(2))
                
                call core_calculate_pair_bond_order_gradients(index_pair, &
                     separations(1:3,1), distances(1), directions(1:3,1), &
                     interaction%pot_index, &                        
                     pair_bo_sums(1:2), &
                     bo_gradients(1:3,1:n_atoms,3:4), &
                     pair_bo_virial(1:6,1:2))
                
                
                ! Add the bond order gradient terms involving the 
                ! atom1-atom2 energy for all atoms.
                ! That is, add the (\nabla_a b_ij) v_ij term with 
                ! the given ij (atom1,atom2) for all a.
                forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) &
                     - tmp_energy*cut_factors(1)*&
                     (bo_gradients(1:3,1:n_atoms,1) + bo_gradients(1:3,1:n_atoms,2) + &
                     bo_gradients(1:3,1:n_atoms,3) + bo_gradients(1:3,1:n_atoms,4))*0.5d0
                
                ! Add the contribution of bond order gradients to the stress tensor
                stress(1:6) = stress(1:6) - tmp_energy*cut_factors(1)* &
                     (bo_virial(1:6,1) + bo_virial(1:6,2) + pair_bo_virial(1:6,1) + pair_bo_virial(1:6,2))*0.5d0
                
             else
                !bo_factors = 1.d0
                !bo_sums = 0.d0
                !bo_gradients = 0.d0
             end if
             
             ! evaluate the 2-body force involving atom1-atom2 interaction
             call evaluate_forces(2,interaction%n_product,separations(1:3,1),distances(1),&
                  interaction,tmp_forces(1:3,1:2),atom_doublet) ! in Potentials.f90
             
             if(interaction%pot_index > -1)then
                ! force on atom 1:                 
                pair_forces(1:3,1) = ( tmp_forces(1:3,1) * cut_factors(1) + &
                     tmp_energy * cut_gradients(1:3,1) ) * &
                     ( bo_factors(index1) + bo_factors(index2) + &
                     pair_bo_factors(1) + pair_bo_factors(2) ) * 0.5d0
                
                ! force on atom 2:
                pair_forces(1:3,2) = ( tmp_forces(1:3,2) * cut_factors(1) - &
                     tmp_energy * cut_gradients(1:3,1) ) * &
                     ( bo_factors(index1) + bo_factors(index2) + &
                     pair_bo_factors(1) + pair_bo_factors(2) ) * 0.5d0
             else
                pair_forces(1:3,1) = ( tmp_forces(1:3,1) * cut_factors(1) + &
                     tmp_energy * cut_gradients(1:3,1) )
                pair_forces(1:3,2) = ( tmp_forces(1:3,2) * cut_factors(1) - &
                     tmp_energy * cut_gradients(1:3,1) )
             end if
             
             forces(1:3,index1) = forces(1:3,index1) + pair_forces(1:3,1)
             forces(1:3,index2) = forces(1:3,index2) + pair_forces(1:3,2)
             
             !***************!
             ! stress tensor !
             !***************!
             
             ! s_xx, s_yy, s_zz, s_yz, s_xz, s_xy:
             stress(1) = stress(1) + separations(1,1) * pair_forces(1,2)
             stress(2) = stress(2) + separations(2,1) * pair_forces(2,2)
             stress(3) = stress(3) + separations(3,1) * pair_forces(3,2)
             stress(4) = stress(4) + separations(2,1) * pair_forces(3,2)
             stress(5) = stress(5) + separations(1,1) * pair_forces(3,2)
             stress(6) = stress(6) + separations(1,1) * pair_forces(2,2)
             
          else if(n_targets > 2)then

             ! If the number of targets is greater than 2,
             ! we have found a many-body potential.
             ! Make a note that we must also evaluate the many-body terms.
             many_bodies_found = .true.

          end if ! n_targets == 2

       end if ! is_active
       end if ! cutoff
    end do ! k

  end subroutine core_evaluate_local_doublet_forces_B



  subroutine core_evaluate_local_doublet_electronegativities_B(atom_doublet, &
       index1, index2, &
       test_index1, &
       interaction_indices, &
       separations, directions, distances, &
       enegs, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: index1, index2, test_index1
    integer, pointer :: interaction_indices(:)
    double precision, intent(inout) :: enegs(:)
    type(atom), intent(in) :: atom_doublet(2)
    double precision, intent(in) :: separations(3,1), directions(3,1), distances(1)
    logical, intent(out) :: many_bodies_found
    
    type(atom) :: atom1, atom2
    integer :: k, n_targets, index_pair(2), n_atoms
    type(potential) :: interaction
    double precision :: tmp_energy, tmp_forces(3,2), tmp_enegs(2), &
         cut_factors(1), &
         pair_forces(3,2), bo_virial(6,2), &
         pair_bo_sums(2), pair_bo_factors(2)
    logical :: is_active

    n_atoms = size(atoms)
    many_bodies_found = .false.
    atom1 = atom_doublet(1)
    atom2 = atom_doublet(2)
    index_pair = (/ index1, index2 /)

    ! loop over potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is it a 2-body potential
       if( interaction%cutoff > distances(1) )then 
          call potential_affects_atom(interaction,atom_doublet(test_index1),is_active,2) ! in Potentials.f90
          if( is_active )then 
             call get_number_of_targets_of_potential_index(interaction%type_index,&
                  n_targets) ! in Potentials.f90
             if( n_targets == 2 )then


                !****************************************!
                ! 2-body electronegativity (atom1-atom2) !
                !****************************************!

                ! If a smooth cutoff is present, we add the
                ! contribution it brings:
                if(interaction%smoothened)then
                   ! get f(r_ij)
                   call smoothening_factor(distances(1),&
                        interaction%cutoff,interaction%soft_cutoff,&
                        cut_factors(1)) ! in Potentials.f90
                else
                   cut_factors(1) = 1.d0
                end if

                ! If there is a bond order factor associated with the potential,
                ! we add the contribution is brings:
                if(interaction%pot_index > -1)then
                   ! get b_i (for all i, they have been precalculated)
                   call core_get_bond_order_factors(interaction%pot_index,&
                        bo_factors)

                   ! get B_ij and \nabla_a B_ij for this ij (not precalculated)
                   ! this is 
                   ! B_ij = (B*_ij + B*_ji) / 2
                   call core_calculate_pair_bond_order_factor(index_pair, &
                        separations(1:3,1), distances(1), directions(1:3,1), &
                        interaction%pot_index, &
                        pair_bo_sums(1:2))
                   call core_post_process_pair_bond_order_factor(atom1, &
                        interaction%pot_index, &
                        pair_bo_sums(1), &
                        pair_bo_factors(1))
                   call core_post_process_pair_bond_order_factor(atom2, &
                        interaction%pot_index, &
                        pair_bo_sums(2), &
                        pair_bo_factors(2))

                else

                end if

                ! evaluate the 2-body e-neg involving atom1-atom2 interaction
                call evaluate_electronegativity(2,interaction%n_product,separations(1:3,1),distances(1),&
                     interaction,tmp_enegs(1:2),atom_doublet) ! in Potentials.f90


                if(interaction%pot_index > -1)then
                    ! e-neg on atom 1:
                    enegs(index1) = enegs(index1) + &
                     ( tmp_enegs(1) * cut_factors(1) ) * &
                     ( bo_factors(index1) +  bo_factors(index2) + &
                     pair_bo_factors(1) + pair_bo_factors(2) ) * 0.5d0

                    ! e-neg on atom 2:
                    enegs(index2) = enegs(index2) + &
                     ( tmp_enegs(2) * cut_factors(1) ) * &
                     ( bo_factors(index1) +  bo_factors(index2) + &
                     pair_bo_factors(1) + pair_bo_factors(2)  ) * 0.5d0
                else
                    ! e-neg on atom 1:
                    enegs(index1) = enegs(index1) + &
                     ( tmp_enegs(1) * cut_factors(1) ) 

                    ! e-neg on atom 2:
                    enegs(index2) = enegs(index2) + &
                     ( tmp_enegs(2) * cut_factors(1) ) 
                
                end if

          else if(n_targets > 2)then

             ! If the number of targets is greater than 2,
             ! we have found a many-body potential.
             ! Make a note that we must also evaluate the many-body terms.
             many_bodies_found = .true.

          end if ! n_targets == 2

       end if ! is_active
       end if ! cutoff
    end do ! k

  end subroutine core_evaluate_local_doublet_electronegativities_B



  subroutine core_evaluate_local_doublet_electronegativities(n_atoms, &
       atom_doublet, &
       index1, index2, &
       test_index1, &
       interaction_indices, &
       separations, directions, distances, &
       enegs, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: n_atoms, index1, index2, test_index1
    integer, pointer :: interaction_indices(:)
    double precision, intent(inout) :: enegs(n_atoms)
    type(atom), intent(in) :: atom_doublet(2)
    double precision, intent(in) :: separations(3,1), directions(3,1), distances(1)
    logical, intent(out) :: many_bodies_found
    
    type(atom) :: atom1, atom2
    integer :: k, n_targets, index_pair(2)
    type(potential) :: interaction
    double precision :: bo_factors(n_atoms), bo_sums(n_atoms), bo_gradients(3,n_atoms,3), &
         tmp_energy, tmp_forces(3,2), tmp_enegs(2), &
         cut_factors(1), &
         pair_forces(3,2), bo_virial(6,2), &
         pair_bo_sums(2), pair_bo_factors(2)
    logical :: is_active

    many_bodies_found = .false.
    atom1 = atom_doublet(1)
    atom2 = atom_doublet(2)
    index_pair = (/ index1, index2 /)

    ! loop over potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is it a 2-body potential
       if( interaction%cutoff > distances(1) )then 
       call potential_affects_atom(interaction,atom_doublet(test_index1),is_active,2) ! in Potentials.f90
       if( is_active )then 
          call get_number_of_targets_of_potential_index(interaction%type_index,&
               n_targets) ! in Potentials.f90
          if( n_targets == 2 )then


                !****************************************!
                ! 2-body electronegativity (atom1-atom2) !
                !****************************************!

                ! If a smooth cutoff is present, we add the
                ! contribution it brings:
                if(interaction%smoothened)then
                   ! get f(r_ij)
                   call smoothening_factor(distances(1),&
                        interaction%cutoff,interaction%soft_cutoff,&
                        cut_factors(1)) ! in Potentials.f90
                else
                   cut_factors(1) = 1.d0
                end if

                ! If there is a bond order factor associated with the potential,
                ! we add the contribution is brings:
                if(interaction%pot_index > -1)then
                   ! get b_i (for all i, they have been precalculated)
                   call core_get_bond_order_factors(interaction%pot_index,&
                        bo_factors)

                   ! get B_ij and \nabla_a B_ij for this ij (not precalculated)
                   ! this is 
                   ! B_ij = (B*_ij + B*_ji) / 2
                   call core_calculate_pair_bond_order_factor(index_pair, &
                        separations(1:3,1), distances(1), directions(1:3,1), &
                        interaction%pot_index, &
                        pair_bo_sums(1:2))
                   call core_post_process_pair_bond_order_factor(atom1, &
                        interaction%pot_index, &
                        pair_bo_sums(1), &
                        pair_bo_factors(1))
                   call core_post_process_pair_bond_order_factor(atom2, &
                        interaction%pot_index, &
                        pair_bo_sums(2), &
                        pair_bo_factors(2))

                else
                   !bo_factors = 1.d0
                end if

                ! evaluate the 2-body e-neg involving atom1-atom2 interaction
                call evaluate_electronegativity(2,interaction%n_product,separations(1:3,1),distances(1),&
                     interaction,tmp_enegs(1:2),atom_doublet) ! in Potentials.f90


                if(interaction%pot_index > -1)then
                    ! e-neg on atom 1:
                    enegs(index1) = enegs(index1) + &
                     ( tmp_enegs(1) * cut_factors(1) ) * &
                     ( bo_factors(index1) +  bo_factors(index2) + &
                     pair_bo_factors(1) + pair_bo_factors(2) ) * 0.5d0

                    ! e-neg on atom 2:
                    enegs(index2) = enegs(index2) + &
                     ( tmp_enegs(2) * cut_factors(1) ) * &
                     ( bo_factors(index1) +  bo_factors(index2) + &
                     pair_bo_factors(1) + pair_bo_factors(2)  ) * 0.5d0
                else
                    ! e-neg on atom 1:
                    enegs(index1) = enegs(index1) + &
                     ( tmp_enegs(1) * cut_factors(1) ) 

                    ! e-neg on atom 2:
                    enegs(index2) = enegs(index2) + &
                     ( tmp_enegs(2) * cut_factors(1) ) 
                
                end if

          else if(n_targets > 2)then

             ! If the number of targets is greater than 2,
             ! we have found a many-body potential.
             ! Make a note that we must also evaluate the many-body terms.
             many_bodies_found = .true.

          end if ! n_targets == 2

       end if ! is_active
       end if ! cutoff
    end do ! k

  end subroutine core_evaluate_local_doublet_electronegativities



  subroutine core_evaluate_local_doublet(n_atoms, &
       atom_doublet, &
       index1, index2, &
       test_index1, &
       interaction_indices, &
       separations, directions, distances, &
       calculation_type,energy,forces,enegs,stress, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: calculation_type, n_atoms, index1, index2, test_index1
    integer, pointer :: interaction_indices(:)
    double precision, intent(inout) :: energy, forces(3,n_atoms), enegs(n_atoms), stress(6)
    type(atom), intent(in) :: atom_doublet(2)
    double precision, intent(in) :: separations(3,1), directions(3,1), distances(1)
    logical, intent(out) :: many_bodies_found
    
    type(atom) :: atom1, atom2
    integer :: k, n_targets, index_pair(2)
    type(potential) :: interaction
    double precision :: bo_factors(n_atoms), bo_sums(n_atoms), bo_gradients(3,n_atoms,3), &
         tmp_energy, tmp_forces(3,2), tmp_enegs(2), &
         cut_factors(1), cut_gradients(3,1), &
         pair_forces(3,2), bo_virial(6,2), &
         pair_bo_sums(2), pair_bo_factors(2), &
         pair_bo_gradients(3,n_atoms,2), pair_bo_virial(6,2)
    logical :: is_active

    many_bodies_found = .false.
    atom1 = atom_doublet(1)
    atom2 = atom_doublet(2)
    index_pair = (/ index1, index2 /)

    ! loop over potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is it a 2-body potential
       if( interaction%cutoff > distances(1) )then 
       call potential_affects_atom(interaction,atom_doublet(test_index1),is_active,2) ! in Potentials.f90
       if( is_active )then 
          call get_number_of_targets_of_potential_index(interaction%type_index,&
               n_targets) ! in Potentials.f90
          if( n_targets == 2 )then

             ! differentiate between energy, force, and electronegativity evaluation
             select case(calculation_type)
             case(energy_evaluation_index)

                !*****************************!
                ! 2-body energy (atom1-atom2) !
                !*****************************!

                ! If there is a bond order factor associated with the potential,
                ! we add the contribution is brings:
                !
                ! V = \sum_ij b_ij v_ij
                ! b_ij = (b_i + b_j) / 2 + B_ij
                if(interaction%pot_index > -1)then
                   ! get b_i (for all i, they have been precalculated)
                   call core_get_bond_order_factors(interaction%pot_index,&
                        bo_factors)

                   ! get B_ij for this ij (not precalculated)
                   ! this is 
                   ! B_ij = (B*_ij + B*_ji) / 2
                   call core_calculate_pair_bond_order_factor(index_pair, &
                        separations(1:3,1), distances(1), directions(1:3,1), &
                        interaction%pot_index, &
                        pair_bo_sums(1:2))
                   call core_post_process_pair_bond_order_factor(atom1, &
                        interaction%pot_index, &
                        pair_bo_sums(1), &
                        pair_bo_factors(1))
                   call core_post_process_pair_bond_order_factor(atom2, &
                        interaction%pot_index, &
                        pair_bo_sums(2), &
                        pair_bo_factors(2))

                else
                   !bo_factors = 1.d0
                end if

                ! If a smooth cutoff is present, we add the
                ! contribution it brings:
                ! 
                ! V = \sum_ij v_ij f(r_ij)
                if(interaction%smoothened)then
                   ! get f(r_ij)
                   call smoothening_factor(distances(1),&
                        interaction%cutoff,interaction%soft_cutoff,&
                        cut_factors(1)) ! in Potentials.f90
                else
                   cut_factors(1) = 1.d0
                end if

                ! evaluate the 2-body energy involving atom1-atom2 interaction
                call evaluate_energy(2,interaction%n_product,separations(1:3,1),distances(1),&
                     interaction,tmp_energy,atom_doublet)  ! in Potentials.f90

                ! add the term: b_ij v_ij f(rij)
                if(interaction%pot_index > -1)then
                    energy = energy + tmp_energy*cut_factors(1)*&
                         (bo_factors(index1)+bo_factors(index2)+pair_bo_factors(1)+pair_bo_factors(2))*0.5d0
                else
                    energy = energy + tmp_energy*cut_factors(1)
                end if

             case(force_evaluation_index)

                !****************************!
                ! 2-body force (atom1-atom2) !
                !****************************!

                ! We will need the energy contribution from atom1-atom2
                ! interaction if smooth cutoffs or bond factors are used,
                ! since we are mulplying the potential.
                if((interaction%pot_index > -1) .or. &
                     interaction%smoothened)then
                   call evaluate_energy(2,interaction%n_product,separations(1:3,1),distances(1),&
                        interaction,tmp_energy,atom_doublet) ! in Potentials.f90
                else
                   tmp_energy = 0.d0
                end if

                ! If a smooth cutoff is present, we add the
                ! contribution it brings:
                ! 
                ! V = \sum_ij v_ij f(r_ij)
                ! F_a = - \nabla_a V 
                !     = - \sum_ij v_ij f'(r_ij) (\nabla_a r_ij) + (\nabla_a v_ij) f(r_ij) 
                !     = - \sum_ij v_ij f'(r_ij) (\nabla_a r_ij) + f_a,ij f(r_ij)
                !
                if(interaction%smoothened)then
                   ! get f(r_ij)
                   call smoothening_factor(distances(1),&
                        interaction%cutoff,interaction%soft_cutoff,&
                        cut_factors(1)) ! in Potentials.f90
                   ! get f'(r_ij) (\nabla_a r_ij)
                   call smoothening_gradient(directions(1:3,1),distances(1),&
                        interaction%cutoff,interaction%soft_cutoff,&
                        cut_gradients(1:3,1)) ! in Potentials.f90
                else
                   cut_factors(1) = 1.d0
                   cut_gradients(1:3,1) = 0.d0
                end if

                ! If there is a bond order factor associated with the potential,
                ! we add the contribution is brings:
                !
                ! V = \sum_ij b_ij v_ij
                ! b_ij = (b_i + b_j) / 2 + B_ij
                ! F_a = - \nabla_a V 
                !     = - \sum_ij (\nabla_a b_ij) v_ij + b_ij (\nabla_a v_ij)
                !     = - \sum_ij (\nabla_a b_ij) v_ij + b_ij f_a,ij
                !
                if(interaction%pot_index > -1)then
                   ! get b_i (for all i, they have been precalculated)
                   call core_get_bond_order_factors(interaction%pot_index,&
                        bo_factors)
                   ! get (\nabla_a b_i) (for all a)
                   call core_get_bond_order_gradients(interaction%pot_index,&
                        index1,& ! atom index
                        1, & ! slot_index
                        bo_gradients(1:3,1:n_atoms,1), &
                        bo_virial(1:6,1))
                   ! get (\nabla_a b_j) (for all a)
                   call core_get_bond_order_gradients(interaction%pot_index,&
                        index2,& ! atom index
                        2, & ! slot_index
                        bo_gradients(1:3,1:n_atoms,2), &
                        bo_virial(1:6,2))


                   ! get B_ij and \nabla_a B_ij for this ij (not precalculated)
                   ! this is 
                   ! B_ij = (B*_ij + B*_ji) / 2
                   call core_calculate_pair_bond_order_factor(index_pair, &
                        separations(1:3,1), distances(1), directions(1:3,1), &
                        interaction%pot_index, &
                        pair_bo_sums(1:2))
                   call core_post_process_pair_bond_order_factor(atom1, &
                        interaction%pot_index, &
                        pair_bo_sums(1), &
                        pair_bo_factors(1))
                   call core_post_process_pair_bond_order_factor(atom2, &
                        interaction%pot_index, &
                        pair_bo_sums(2), &
                        pair_bo_factors(2))

                   call core_calculate_pair_bond_order_gradients(index_pair, &
                        separations(1:3,1), distances(1), directions(1:3,1), &
                        interaction%pot_index, &                        
                        pair_bo_sums(1:2), &
                        pair_bo_gradients(1:3,1:n_atoms,1:2), &
                        pair_bo_virial(1:6,1:2))


                   ! Add the bond order gradient terms involving the 
                   ! atom1-atom2 energy for all atoms.
                   ! That is, add the (\nabla_a b_ij) v_ij term with 
                   ! the given ij (atom1,atom2) for all a.
                   forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) &
                        - tmp_energy*cut_factors(1)*&
                        (bo_gradients(1:3,1:n_atoms,1) + bo_gradients(1:3,1:n_atoms,2) + &
                        pair_bo_gradients(1:3,1:n_atoms,1) + pair_bo_gradients(1:3,1:n_atoms,2))*0.5d0

                   ! Add the contribution of bond order gradients to the stress tensor
                   stress(1:6) = stress(1:6) - tmp_energy*cut_factors(1)* &
                        (bo_virial(1:6,1) + bo_virial(1:6,2) + pair_bo_virial(1:6,1) + pair_bo_virial(1:6,2))*0.5d0

                else
                   !bo_factors = 1.d0
                   !bo_sums = 0.d0
                   !bo_gradients = 0.d0
                end if

                ! evaluate the 2-body force involving atom1-atom2 interaction
                call evaluate_forces(2,interaction%n_product,separations(1:3,1),distances(1),&
                     interaction,tmp_forces(1:3,1:2),atom_doublet) ! in Potentials.f90

                if(interaction%pot_index > -1)then
                    ! force on atom 1:                 
                    pair_forces(1:3,1) = ( tmp_forces(1:3,1) * cut_factors(1) + &
                         tmp_energy * cut_gradients(1:3,1) ) * &
                         ( bo_factors(index1) + bo_factors(index2) + &
                         pair_bo_factors(1) + pair_bo_factors(2) ) * 0.5d0

                    ! force on atom 2:
                    pair_forces(1:3,2) = ( tmp_forces(1:3,2) * cut_factors(1) - &
                         tmp_energy * cut_gradients(1:3,1) ) * &
                         ( bo_factors(index1) + bo_factors(index2) + &
                         pair_bo_factors(1) + pair_bo_factors(2) ) * 0.5d0
                else
                    pair_forces(1:3,1) = ( tmp_forces(1:3,1) * cut_factors(1) + &
                         tmp_energy * cut_gradients(1:3,1) )
                    pair_forces(1:3,2) = ( tmp_forces(1:3,2) * cut_factors(1) - &
                         tmp_energy * cut_gradients(1:3,1) )
                end if

                forces(1:3,index1) = forces(1:3,index1) + pair_forces(1:3,1)
                forces(1:3,index2) = forces(1:3,index2) + pair_forces(1:3,2)

                !***************!
                ! stress tensor !
                !***************!
                
                ! s_xx, s_yy, s_zz, s_yz, s_xz, s_xy:
                stress(1) = stress(1) + separations(1,1) * pair_forces(1,2)
                stress(2) = stress(2) + separations(2,1) * pair_forces(2,2)
                stress(3) = stress(3) + separations(3,1) * pair_forces(3,2)
                stress(4) = stress(4) + separations(2,1) * pair_forces(3,2)
                stress(5) = stress(5) + separations(1,1) * pair_forces(3,2)
                stress(6) = stress(6) + separations(1,1) * pair_forces(2,2)

             case(electronegativity_evaluation_index)

                !****************************************!
                ! 2-body electronegativity (atom1-atom2) !
                !****************************************!

                ! If a smooth cutoff is present, we add the
                ! contribution it brings:
                if(interaction%smoothened)then
                   ! get f(r_ij)
                   call smoothening_factor(distances(1),&
                        interaction%cutoff,interaction%soft_cutoff,&
                        cut_factors(1)) ! in Potentials.f90
                else
                   cut_factors(1) = 1.d0
                end if

                ! If there is a bond order factor associated with the potential,
                ! we add the contribution is brings:
                if(interaction%pot_index > -1)then
                   ! get b_i (for all i, they have been precalculated)
                   call core_get_bond_order_factors(interaction%pot_index,&
                        bo_factors)

                   ! get B_ij and \nabla_a B_ij for this ij (not precalculated)
                   ! this is 
                   ! B_ij = (B*_ij + B*_ji) / 2
                   call core_calculate_pair_bond_order_factor(index_pair, &
                        separations(1:3,1), distances(1), directions(1:3,1), &
                        interaction%pot_index, &
                        pair_bo_sums(1:2))
                   call core_post_process_pair_bond_order_factor(atom1, &
                        interaction%pot_index, &
                        pair_bo_sums(1), &
                        pair_bo_factors(1))
                   call core_post_process_pair_bond_order_factor(atom2, &
                        interaction%pot_index, &
                        pair_bo_sums(2), &
                        pair_bo_factors(2))

                else
                   !bo_factors = 1.d0
                end if

                ! evaluate the 2-body e-neg involving atom1-atom2 interaction
                call evaluate_electronegativity(2,interaction%n_product,separations(1:3,1),distances(1),&
                     interaction,tmp_enegs(1:2),atom_doublet) ! in Potentials.f90


                if(interaction%pot_index > -1)then
                    ! e-neg on atom 1:
                    enegs(index1) = enegs(index1) + &
                     ( tmp_enegs(1) * cut_factors(1) ) * &
                     ( bo_factors(index1) +  bo_factors(index2) + &
                     pair_bo_factors(1) + pair_bo_factors(2) ) * 0.5d0

                    ! e-neg on atom 2:
                    enegs(index2) = enegs(index2) + &
                     ( tmp_enegs(2) * cut_factors(1) ) * &
                     ( bo_factors(index1) +  bo_factors(index2) + &
                     pair_bo_factors(1) + pair_bo_factors(2)  ) * 0.5d0
                else
                    ! e-neg on atom 1:
                    enegs(index1) = enegs(index1) + &
                     ( tmp_enegs(1) * cut_factors(1) ) 

                    ! e-neg on atom 2:
                    enegs(index2) = enegs(index2) + &
                     ( tmp_enegs(2) * cut_factors(1) ) 
                
                end if
             end select

          else if(n_targets > 2)then

             ! If the number of targets is greater than 2,
             ! we have found a many-body potential.
             ! Make a note that we must also evaluate the many-body terms.
             many_bodies_found = .true.

          end if ! n_targets == 2

       end if ! is_active
       end if ! cutoff
    end do ! k

  end subroutine core_evaluate_local_doublet




  subroutine core_evaluate_local_triplet(n_atoms, &
       atom_triplet, &
       index1, index2, index3, &
       test_index1, test_index2, &
       interaction_indices, &
       separations, directions, distances, &
       calculation_type,energy,forces,enegs,stress, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: calculation_type, n_atoms, index1, index2, index3, test_index1, test_index2
    integer, pointer :: interaction_indices(:)
    double precision, intent(out) :: energy, forces(3,n_atoms), enegs(n_atoms), stress(6)
    double precision, intent(in) :: separations(3,2), directions(3,2), distances(2)
    logical, intent(out) :: many_bodies_found
    type(atom), intent(in) :: atom_triplet(3)

    type(atom) :: atom1, atom2, atom3
    integer :: k, n_targets
    type(potential) :: interaction
    double precision :: bo_factors(n_atoms), bo_sums(n_atoms), bo_gradients(3,n_atoms,3), &
         tmp_energy, tmp_forces(3,3), tmp_enegs(3), &
         cut_factors(2), cut_gradients(3,2), &
         triplet_forces(3,3), bo_virial(6,3)
    logical :: is_active

    many_bodies_found = .false.
    atom1 = atom_triplet(1)
    atom2 = atom_triplet(2)
    atom3 = atom_triplet(3)

    ! loop over the potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))
       call get_number_of_targets_of_potential_index(interaction%type_index,&
            n_targets) ! in Potentials.f90
       call potential_affects_atom(interaction,atom_triplet(test_index1),is_active,2) ! in Potentials.f90

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is atom3 affected by the potential,
       ! is it a 3-body potential
       if( is_active .and. n_targets == 3 .and. interaction%cutoff > distances(1) )then
          call potential_affects_atom(interaction,atom_triplet(test_index2),is_active,3)  ! in Potentials.f90
          if( is_active )then

             if( interaction%cutoff > distances(2) )then


                ! differentiate between energy, force, and electronegativity evaluation
                select case(calculation_type)
                case(energy_evaluation_index)

                   !***********************************!
                   ! 3-body energy (atom1-atom2-atom3) !
                   !***********************************!


                   ! If there is a bond order factor associated with the potential,
                   ! we add the contribution is brings:
                   !
                   ! V = \sum_ijk b_ijk v_ijk
                   ! b_ijk = (b_i + b_j + b_k) / 3
                   if(interaction%pot_index > -1)then
                      ! get b_i (for all i, they have been precalculated)
                      call core_get_bond_order_factors(interaction%pot_index,&
                           bo_factors)
                   else
                      bo_factors = 1.d0
                   end if

                   ! If a smooth cutoff is present, we add the
                   ! contribution it brings:
                   ! 
                   ! V = \sum_ijk v_ijk f(r_ij) f(r_ik)
                   if(interaction%smoothened)then
                      ! get f(r_ij)
                      call smoothening_factor(distances(1),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(1))
                      ! get f(r_ik)
                      call smoothening_factor(distances(2),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(2))
                   else
                      cut_factors(1:2) = 1.d0
                   end if

                   ! evaluate the 3-body energy involving atom2-atom1-atom3 interaction
                   call evaluate_energy(3,interaction%n_product,separations(1:3,1:2),distances(1:2),&
                        interaction,tmp_energy,atom_triplet)

                   ! add the term: b_ijk v_ijk f(r_ij) f(r_ik)
                   energy = energy + tmp_energy*cut_factors(1)*cut_factors(2)*&
                        (bo_factors(index1)+bo_factors(index2)+bo_factors(index3))/3.d0

                case(force_evaluation_index)

                   !**************!
                   ! 3-body force !
                   !**************!


                   ! We will need the energy contribution 
                   ! if smooth cutoffs or bond factors are used,
                   ! since we are mulplying the potential.
                   if((interaction%pot_index > -1) .or. &
                        interaction%smoothened)then
                      call evaluate_energy(3,interaction%n_product,separations(1:3,1:2),distances(1:2),&
                           interaction,tmp_energy,atom_triplet) ! in Potentials.f90
                   else
                      tmp_energy = 0.d0
                   end if

                   ! If a smooth cutoff is present, we add the
                   ! contribution it brings:
                   ! 
                   ! V = \sum_ijk v_ijk f(r_ij) f(r_ik)
                   ! F_a = - \nabla_a V 
                   !     = - \sum_ij ( v_ij f'(r_ij) f(r_ik) (\nabla_a r_ij) + 
                   !                   v_ij f(r_ij) f'(r_ik) (\nabla_a r_ik) +
                   !                   (\nabla_a v_ij) f(r_ij) f(r_ik) )
                   !
                   if(interaction%smoothened)then
                      ! get f(r_ij)
                      call smoothening_factor(distances(1),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(1)) ! in Potentials.f90
                      ! get f'(r_ij) (\nabla_a r_ij)
                      call smoothening_gradient(directions(1:3,1),distances(1),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_gradients(1:3,1)) ! in Potentials.f90 
                      ! get f(r_jk)
                      call smoothening_factor(distances(2),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(2)) ! in Potentials.f90
                      ! get f'(r_jk) (\nabla_a r_jk)
                      call smoothening_gradient(directions(1:3,2),distances(2),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_gradients(1:3,2)) ! in Potentials.f90
                   else
                      cut_factors(1:2) = 1.d0
                      cut_gradients(1:3,1:2) = 0.d0
                   end if

                   ! If there is a bond order factor associated with the potential,
                   ! we add the contribution is brings:
                   !
                   ! V = \sum_ijk b_ijk v_ijk
                   ! b_ijk = (b_i + b_j + b_k) / 3
                   ! F_a = - \nabla_a V 
                   !     = - \sum_ijk (\nabla_a b_ijk) v_ijk + b_ijk (\nabla_a v_ijk)
                   !     = - \sum_ijk (\nabla_a b_ijk) v_ijk + b_ij f_a,ijk
                   !
                   if(interaction%pot_index > -1)then
                      ! get b_i (for all i, they have been precalculated)
                      call core_get_bond_order_factors(interaction%pot_index,&
                           bo_factors)
                      ! get (\nabla_a b_i) (for all a)
                      call core_get_bond_order_gradients(interaction%pot_index,&
                           index1,& ! atom index
                           1, & ! slot_index
                           bo_gradients(1:3,1:n_atoms,1), &
                           bo_virial(1:6,1))
                      ! get (\nabla_a b_j) (for all a)
                      call core_get_bond_order_gradients(interaction%pot_index,&
                           index2,& ! atom index
                           2, & ! slot_index
                           bo_gradients(1:3,1:n_atoms,2), &
                           bo_virial(1:6,2))
                      ! get (\nabla_a b_k) (for all a)
                      call core_get_bond_order_gradients(interaction%pot_index,&
                           index3,& ! atom index
                           3, & ! slot_index
                           bo_gradients(1:3,1:n_atoms,3), &
                           bo_virial(1:6,3))

                      ! Add the bond order gradient terms involving the 
                      ! atom2-atom1-atom3 energy for all atoms.
                      ! That is, add the (\nabla_a b_ijk) v_ijk term with 
                      ! the given ijk (atom2,atom1,atom3) for all a.
                      forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) &
                           - tmp_energy*cut_factors(1)*cut_factors(2)*&
                           ( bo_gradients(1:3,1:n_atoms,1) &
                           + bo_gradients(1:3,1:n_atoms,2) &
                           + bo_gradients(1:3,1:n_atoms,3) )/3.d0

                      stress(1:6) = stress(1:6) - tmp_energy*cut_factors(1)*cut_factors(2)*&
                           ( bo_virial(1:6,1) &
                           + bo_virial(1:6,2) &
                           + bo_virial(1:6,3) )/3.d0

                   else
                      !bo_factors = 1.d0
                      !bo_sums = 0.d0
                      !bo_gradients = 0.d0
                   end if

                   ! evaluate the 3-body force
                   call evaluate_forces(3,interaction%n_product,separations(1:3,1:2),distances(1:2),interaction,&
                        tmp_forces(1:3,1:3),atom_triplet) ! in Potentials.f90


                   if(interaction%pot_index > -1)then
                    ! force on atom 1:
                    triplet_forces(1:3,1) = ( tmp_forces(1:3,1)*cut_factors(1)*cut_factors(2) + &
                        cut_gradients(1:3,1)*cut_factors(2)*tmp_energy ) * &
                        ( bo_factors(index1) &
                        + bo_factors(index2) &
                        + bo_factors(index3) )/3.d0

                    ! force on atom 2:
                    triplet_forces(1:3,2) = ( tmp_forces(1:3,2)*cut_factors(1)*cut_factors(2) + &
                        (-cut_gradients(1:3,1)*cut_factors(2) + &
                        cut_gradients(1:3,2)*cut_factors(1)) * tmp_energy ) * &
                        ( bo_factors(index1) &
                        + bo_factors(index2) &
                        + bo_factors(index3) )/3.d0

                    ! force on atom 3:
                    triplet_forces(1:3,3) = ( tmp_forces(1:3,3)*cut_factors(1)*cut_factors(2) - &
                        cut_gradients(1:3,2)*cut_factors(1)*tmp_energy ) * &
                        ( bo_factors(index1) &
                        + bo_factors(index2) &
                        + bo_factors(index3) )/3.d0
                   else
                    ! force on atom 1:
                    triplet_forces(1:3,1) = ( tmp_forces(1:3,1)*cut_factors(1)*cut_factors(2) + &
                        cut_gradients(1:3,1)*cut_factors(2)*tmp_energy ) 

                    ! force on atom 2:
                    triplet_forces(1:3,2) = ( tmp_forces(1:3,2)*cut_factors(1)*cut_factors(2) + &
                        (-cut_gradients(1:3,1)*cut_factors(2) + &
                        cut_gradients(1:3,2)*cut_factors(1)) * tmp_energy ) 

                    ! force on atom 3:
                    triplet_forces(1:3,3) = ( tmp_forces(1:3,3)*cut_factors(1)*cut_factors(2) - &
                        cut_gradients(1:3,2)*cut_factors(1)*tmp_energy ) 
                   
                   end if

                   forces(1:3,index1) = forces(1:3,index1) + triplet_forces(1:3,1)
                   forces(1:3,index2) = forces(1:3,index2) + triplet_forces(1:3,2)
                   forces(1:3,index3) = forces(1:3,index3) + triplet_forces(1:3,3)
                        
                   !***************!
                   ! stress tensor !
                   !***************!
                   
                   ! s_xx, s_yy, s_zz, s_yz, s_xz, s_xy:
                   stress(1) = stress(1) + separations(1,1) * triplet_forces(1,2) &
                        + (separations(1,1) + separations(1,2)) * triplet_forces(1,3)
                   stress(2) = stress(2) + separations(2,1) * triplet_forces(2,2) &
                        + (separations(2,1) + separations(2,2)) * triplet_forces(2,3)
                   stress(3) = stress(3) + separations(3,1) * triplet_forces(3,2) &
                        + (separations(3,1) + separations(3,2)) * triplet_forces(3,3)
                   stress(4) = stress(4) + separations(2,1) * triplet_forces(3,2) &
                        + (separations(2,1) + separations(2,2)) * triplet_forces(3,3)
                   stress(5) = stress(5) + separations(1,1) * triplet_forces(3,2) &
                        + (separations(1,1) + separations(1,2)) * triplet_forces(3,3)
                   stress(6) = stress(6) + separations(1,1) * triplet_forces(2,2) &
                        + (separations(1,1) + separations(1,2)) * triplet_forces(2,3)

                case(electronegativity_evaluation_index)

                   !**********************************************!
                   ! 3-body electronegativity (atom1-atom2-atom3) !
                   !**********************************************!

                   ! If a smooth cutoff is present, we add the
                   ! contribution it brings:
                   if(interaction%smoothened)then
                      ! get f(r_ij)
                      call smoothening_factor(distances(1),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(1)) ! in Potentials.f90
                      ! get f(r_ik)
                      call smoothening_factor(distances(2),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(2)) ! in Potentials.f90
                   else
                      cut_factors(1:2) = 1.d0
                   end if

                   ! If there is a bond order factor associated with the potential,
                   ! we add the contribution is brings:
                   if(interaction%pot_index > -1)then
                      ! get b_i (for all i, they have been precalculated)
                      call core_get_bond_order_factors(interaction%pot_index,&
                           bo_factors) ! in Potentials.f90
                   else
                      bo_factors = 1.d0
                   end if

                   ! evaluate the 3-body e-neg involving atom2-atom1-atom3 interaction
                   call evaluate_electronegativity(3,interaction%n_product,separations(1:3,1:2),&
                        distances(1:2),interaction,&
                        tmp_enegs(1:3),atom_triplet) ! in Potentials.f90

                   ! e-neg on atom 1:
                   enegs(index1) = enegs(index1) + &
                        ( tmp_enegs(1)*cut_factors(1)*cut_factors(2) ) * &
                        ( bo_factors(index1) &
                        + bo_factors(index2) &
                        + bo_factors(index3) )/3.d0

                   ! e-neg on atom 2:
                   enegs(index2) = enegs(index2) + &
                        ( tmp_enegs(2)*cut_factors(1)*cut_factors(2) ) * &
                        ( bo_factors(index1) &
                        + bo_factors(index2) &
                        + bo_factors(index3) )/3.d0

                   ! e-neg on atom 3:
                   enegs(index3) = enegs(index3) + &
                        ( tmp_enegs(3)*cut_factors(1)*cut_factors(2) ) * &
                        ( bo_factors(index1) &
                        + bo_factors(index2) &
                        + bo_factors(index3) )/3.d0


                end select


             end if ! interaction%cutoff > distances(2)

          end if ! is_active
       else if(n_targets > 3)then

          many_bodies_found = .true.

       end if ! is_active .and. n_targets == 3

    end do ! k = 1, size(interaction_indices)


  end subroutine core_evaluate_local_triplet





  subroutine core_evaluate_local_triplet_B(atom_triplet, &
       index1, index2, index3, &
       test_index1, test_index2, &
       interaction_indices, &
       separations, directions, distances, &
       calculation_type,energy,forces,enegs,stress, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: calculation_type, index1, index2, index3, test_index1, test_index2
    integer, pointer :: interaction_indices(:)
    double precision, intent(inout) :: energy, forces(:,:), enegs(:), stress(6)
    double precision, intent(in) :: separations(3,2), directions(3,2), distances(2)
    logical, intent(out) :: many_bodies_found
    type(atom), intent(in) :: atom_triplet(3)

    type(atom) :: atom1, atom2, atom3
    integer :: k, n_targets, n_atoms
    type(potential) :: interaction
    double precision :: tmp_energy, tmp_forces(3,3), tmp_enegs(3), &
         cut_factors(2), cut_gradients(3,2), &
         triplet_forces(3,3), bo_virial(6,3)
    logical :: is_active

    many_bodies_found = .false.
    atom1 = atom_triplet(1)
    atom2 = atom_triplet(2)
    atom3 = atom_triplet(3)

    ! loop over the potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))
       call get_number_of_targets_of_potential_index(interaction%type_index,&
            n_targets) ! in Potentials.f90
       call potential_affects_atom(interaction,atom_triplet(test_index1),is_active,2) ! in Potentials.f90

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is atom3 affected by the potential,
       ! is it a 3-body potential
       if( is_active .and. n_targets == 3 .and. interaction%cutoff > distances(1) )then
          call potential_affects_atom(interaction,atom_triplet(test_index2),is_active,3)  ! in Potentials.f90
          if( is_active )then

             if( interaction%cutoff > distances(2) )then


                ! differentiate between energy, force, and electronegativity evaluation
                select case(calculation_type)
                case(energy_evaluation_index)

                   !***********************************!
                   ! 3-body energy (atom1-atom2-atom3) !
                   !***********************************!


                   ! If there is a bond order factor associated with the potential,
                   ! we add the contribution is brings:
                   !
                   ! V = \sum_ijk b_ijk v_ijk
                   ! b_ijk = (b_i + b_j + b_k) / 3
                   if(interaction%pot_index > -1)then
                      ! get b_i (for all i, they have been precalculated)
                      call core_get_bond_order_factors(interaction%pot_index,&
                           bo_factors)
                   else

                   end if

                   ! If a smooth cutoff is present, we add the
                   ! contribution it brings:
                   ! 
                   ! V = \sum_ijk v_ijk f(r_ij) f(r_ik)
                   if(interaction%smoothened)then
                      ! get f(r_ij)
                      call smoothening_factor(distances(1),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(1))
                      ! get f(r_ik)
                      call smoothening_factor(distances(2),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(2))
                   else
                      cut_factors(1:2) = 1.d0
                   end if

                   ! evaluate the 3-body energy involving atom2-atom1-atom3 interaction
                   call evaluate_energy(3,interaction%n_product,separations(1:3,1:2),distances(1:2),&
                        interaction,tmp_energy,atom_triplet)

                   if(interaction%pot_index > -1)then
                      ! add the term: b_ijk v_ijk f(r_ij) f(r_ik)
                      energy = energy + tmp_energy*cut_factors(1)*cut_factors(2)*&
                           (bo_factors(index1)+bo_factors(index2)+bo_factors(index3))/3.d0
                   else
                      energy = energy + tmp_energy*cut_factors(1)*cut_factors(2)
                   end if

                case(force_evaluation_index)

                   !**************!
                   ! 3-body force !
                   !**************!


                   ! We will need the energy contribution 
                   ! if smooth cutoffs or bond factors are used,
                   ! since we are mulplying the potential.
                   if((interaction%pot_index > -1) .or. &
                        interaction%smoothened)then
                      call evaluate_energy(3,interaction%n_product,separations(1:3,1:2),distances(1:2),&
                           interaction,tmp_energy,atom_triplet) ! in Potentials.f90
                   else
                      tmp_energy = 0.d0
                   end if

                   ! If a smooth cutoff is present, we add the
                   ! contribution it brings:
                   ! 
                   ! V = \sum_ijk v_ijk f(r_ij) f(r_ik)
                   ! F_a = - \nabla_a V 
                   !     = - \sum_ij ( v_ij f'(r_ij) f(r_ik) (\nabla_a r_ij) + 
                   !                   v_ij f(r_ij) f'(r_ik) (\nabla_a r_ik) +
                   !                   (\nabla_a v_ij) f(r_ij) f(r_ik) )
                   !
                   if(interaction%smoothened)then
                      ! get f(r_ij)
                      call smoothening_factor(distances(1),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(1)) ! in Potentials.f90
                      ! get f'(r_ij) (\nabla_a r_ij)
                      call smoothening_gradient(directions(1:3,1),distances(1),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_gradients(1:3,1)) ! in Potentials.f90 
                      ! get f(r_jk)
                      call smoothening_factor(distances(2),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(2)) ! in Potentials.f90
                      ! get f'(r_jk) (\nabla_a r_jk)
                      call smoothening_gradient(directions(1:3,2),distances(2),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_gradients(1:3,2)) ! in Potentials.f90
                   else
                      cut_factors(1:2) = 1.d0
                      cut_gradients(1:3,1:2) = 0.d0
                   end if

                   ! If there is a bond order factor associated with the potential,
                   ! we add the contribution is brings:
                   !
                   ! V = \sum_ijk b_ijk v_ijk
                   ! b_ijk = (b_i + b_j + b_k) / 3
                   ! F_a = - \nabla_a V 
                   !     = - \sum_ijk (\nabla_a b_ijk) v_ijk + b_ijk (\nabla_a v_ijk)
                   !     = - \sum_ijk (\nabla_a b_ijk) v_ijk + b_ij f_a,ijk
                   !
                   if(interaction%pot_index > -1)then
                      ! get b_i (for all i, they have been precalculated)
                      call core_get_bond_order_factors(interaction%pot_index,&
                           bo_factors)
                      ! get (\nabla_a b_i) (for all a)
                      call core_get_bond_order_gradients(interaction%pot_index,&
                           index1,& ! atom index
                           1, & ! slot_index
                           bo_gradients(1:3,1:n_atoms,1), &
                           bo_virial(1:6,1))
                      ! get (\nabla_a b_j) (for all a)
                      call core_get_bond_order_gradients(interaction%pot_index,&
                           index2,& ! atom index
                           2, & ! slot_index
                           bo_gradients(1:3,1:n_atoms,2), &
                           bo_virial(1:6,2))
                      ! get (\nabla_a b_k) (for all a)
                      call core_get_bond_order_gradients(interaction%pot_index,&
                           index3,& ! atom index
                           3, & ! slot_index
                           bo_gradients(1:3,1:n_atoms,3), &
                           bo_virial(1:6,3))

                      ! Add the bond order gradient terms involving the 
                      ! atom2-atom1-atom3 energy for all atoms.
                      ! That is, add the (\nabla_a b_ijk) v_ijk term with 
                      ! the given ijk (atom2,atom1,atom3) for all a.
                      forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) &
                           - tmp_energy*cut_factors(1)*cut_factors(2)*&
                           ( bo_gradients(1:3,1:n_atoms,1) &
                           + bo_gradients(1:3,1:n_atoms,2) &
                           + bo_gradients(1:3,1:n_atoms,3) )/3.d0

                      stress(1:6) = stress(1:6) - tmp_energy*cut_factors(1)*cut_factors(2)*&
                           ( bo_virial(1:6,1) &
                           + bo_virial(1:6,2) &
                           + bo_virial(1:6,3) )/3.d0

                   else
                      !bo_factors = 1.d0
                      !bo_sums = 0.d0
                      !bo_gradients = 0.d0
                   end if

                   ! evaluate the 3-body force
                   call evaluate_forces(3,interaction%n_product,separations(1:3,1:2),distances(1:2),interaction,&
                        tmp_forces(1:3,1:3),atom_triplet) ! in Potentials.f90


                   if(interaction%pot_index > -1)then
                      ! force on atom 1:
                      triplet_forces(1:3,1) = ( tmp_forces(1:3,1)*cut_factors(1)*cut_factors(2) + &
                           cut_gradients(1:3,1)*cut_factors(2)*tmp_energy ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) )/3.d0
                      
                      ! force on atom 2:
                      triplet_forces(1:3,2) = ( tmp_forces(1:3,2)*cut_factors(1)*cut_factors(2) + &
                           (-cut_gradients(1:3,1)*cut_factors(2) + &
                           cut_gradients(1:3,2)*cut_factors(1)) * tmp_energy ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) )/3.d0
                      
                      ! force on atom 3:
                      triplet_forces(1:3,3) = ( tmp_forces(1:3,3)*cut_factors(1)*cut_factors(2) - &
                           cut_gradients(1:3,2)*cut_factors(1)*tmp_energy ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) )/3.d0
                   else
                      ! force on atom 1:
                      triplet_forces(1:3,1) = ( tmp_forces(1:3,1)*cut_factors(1)*cut_factors(2) + &
                           cut_gradients(1:3,1)*cut_factors(2)*tmp_energy ) 
                      
                      ! force on atom 2:
                      triplet_forces(1:3,2) = ( tmp_forces(1:3,2)*cut_factors(1)*cut_factors(2) + &
                           (-cut_gradients(1:3,1)*cut_factors(2) + &
                           cut_gradients(1:3,2)*cut_factors(1)) * tmp_energy ) 
                      
                      ! force on atom 3:
                      triplet_forces(1:3,3) = ( tmp_forces(1:3,3)*cut_factors(1)*cut_factors(2) - &
                           cut_gradients(1:3,2)*cut_factors(1)*tmp_energy ) 
                      
                   end if

                   forces(1:3,index1) = forces(1:3,index1) + triplet_forces(1:3,1)
                   forces(1:3,index2) = forces(1:3,index2) + triplet_forces(1:3,2)
                   forces(1:3,index3) = forces(1:3,index3) + triplet_forces(1:3,3)
                        
                   !***************!
                   ! stress tensor !
                   !***************!
                   
                   ! s_xx, s_yy, s_zz, s_yz, s_xz, s_xy:
                   stress(1) = stress(1) + separations(1,1) * triplet_forces(1,2) &
                        + (separations(1,1) + separations(1,2)) * triplet_forces(1,3)
                   stress(2) = stress(2) + separations(2,1) * triplet_forces(2,2) &
                        + (separations(2,1) + separations(2,2)) * triplet_forces(2,3)
                   stress(3) = stress(3) + separations(3,1) * triplet_forces(3,2) &
                        + (separations(3,1) + separations(3,2)) * triplet_forces(3,3)
                   stress(4) = stress(4) + separations(2,1) * triplet_forces(3,2) &
                        + (separations(2,1) + separations(2,2)) * triplet_forces(3,3)
                   stress(5) = stress(5) + separations(1,1) * triplet_forces(3,2) &
                        + (separations(1,1) + separations(1,2)) * triplet_forces(3,3)
                   stress(6) = stress(6) + separations(1,1) * triplet_forces(2,2) &
                        + (separations(1,1) + separations(1,2)) * triplet_forces(2,3)

                case(electronegativity_evaluation_index)

                   !**********************************************!
                   ! 3-body electronegativity (atom1-atom2-atom3) !
                   !**********************************************!

                   ! If a smooth cutoff is present, we add the
                   ! contribution it brings:
                   if(interaction%smoothened)then
                      ! get f(r_ij)
                      call smoothening_factor(distances(1),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(1)) ! in Potentials.f90
                      ! get f(r_ik)
                      call smoothening_factor(distances(2),&
                           interaction%cutoff,interaction%soft_cutoff,&
                           cut_factors(2)) ! in Potentials.f90
                   else
                      cut_factors(1:2) = 1.d0
                   end if

                   ! If there is a bond order factor associated with the potential,
                   ! we add the contribution is brings:
                   if(interaction%pot_index > -1)then
                      ! get b_i (for all i, they have been precalculated)
                      call core_get_bond_order_factors(interaction%pot_index,&
                           bo_factors) ! in Potentials.f90
                   else
                      !bo_factors = 1.d0
                   end if

                   ! evaluate the 3-body e-neg involving atom2-atom1-atom3 interaction
                   call evaluate_electronegativity(3,interaction%n_product,separations(1:3,1:2),&
                        distances(1:2),interaction,&
                        tmp_enegs(1:3),atom_triplet) ! in Potentials.f90

                   if(interaction%pot_index > -1)then
                      ! e-neg on atom 1:
                      enegs(index1) = enegs(index1) + &
                           ( tmp_enegs(1)*cut_factors(1)*cut_factors(2) ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) )/3.d0
                      
                      ! e-neg on atom 2:
                      enegs(index2) = enegs(index2) + &
                           ( tmp_enegs(2)*cut_factors(1)*cut_factors(2) ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) )/3.d0
                      
                      ! e-neg on atom 3:
                      enegs(index3) = enegs(index3) + &
                           ( tmp_enegs(3)*cut_factors(1)*cut_factors(2) ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) )/3.d0

                   else
                      ! e-neg on atom 1:
                      enegs(index1) = enegs(index1) + &
                           ( tmp_enegs(1)*cut_factors(1)*cut_factors(2) ) 
                      
                      ! e-neg on atom 2:
                      enegs(index2) = enegs(index2) + &
                           ( tmp_enegs(2)*cut_factors(1)*cut_factors(2) ) 
                      
                      ! e-neg on atom 3:
                      enegs(index3) = enegs(index3) + &
                           ( tmp_enegs(3)*cut_factors(1)*cut_factors(2) ) 

                   end if

                end select


             end if ! interaction%cutoff > distances(2)

          end if ! is_active
       else if(n_targets > 3)then

          many_bodies_found = .true.

       end if ! is_active .and. n_targets == 3

    end do ! k = 1, size(interaction_indices)


  end subroutine core_evaluate_local_triplet_B







  subroutine core_evaluate_local_quadruplet(n_atoms, &
       atom_quadruplet, &
       index1, index2, index3, index4, &
       test_index1, test_index2, test_index3, &
       interaction_indices, &
       separations, directions, distances, &
       calculation_type,energy,forces,enegs,stress, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: calculation_type, n_atoms, index1, index2, index3, index4, &
         test_index1, test_index2, test_index3
    integer, pointer :: interaction_indices(:)
    double precision, intent(out) :: energy, forces(3,n_atoms), enegs(n_atoms), stress(6)
    double precision, intent(in) :: separations(3,3), directions(3,3), distances(3)
    logical, intent(out) :: many_bodies_found
    type(atom), intent(in) :: atom_quadruplet(4)

    type(atom) :: atom1, atom2, atom3, atom4
    integer :: k, n_targets
    type(potential) :: interaction
    double precision :: bo_factors(n_atoms), bo_sums(n_atoms), bo_gradients(3,n_atoms,4), &
         tmp_energy, tmp_forces(3,4), tmp_enegs(4), &
         cut_factors(3), cut_gradients(3,3), quad_forces(3,4), bo_virial(6,4)
    logical :: is_active


    many_bodies_found = .false.
    atom1 = atom_quadruplet(1)
    atom2 = atom_quadruplet(2)
    atom3 = atom_quadruplet(3)
    atom4 = atom_quadruplet(3)

    ! loop over the potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))
       call get_number_of_targets_of_potential_index(interaction%type_index,&
            n_targets) ! in Potentials.f90
       call potential_affects_atom(interaction,atom_quadruplet(test_index1),is_active,2) ! in Potentials.f90

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is atom3 affected by the potential,
       ! is atom4 affected by the potential,
       ! is it a 4-body potential
       if( is_active .and. n_targets == 4 .and. interaction%cutoff > distances(1) )then
          call potential_affects_atom(interaction,atom_quadruplet(test_index2),is_active,3)  ! in Potentials.f90
          if( is_active )then
             call potential_affects_atom(interaction,atom_quadruplet(test_index3),is_active,4)  ! in Potentials.f90
             if( is_active )then
                if( interaction%cutoff > distances(2) .and. interaction%cutoff > distances(3) )then

                   ! differentiate between energy, force, and electronegativity evaluation
                   select case(calculation_type)
                   case(energy_evaluation_index)

                      !***************!
                      ! 4-body energy !
                      !***************!

                      ! If there is a bond order factor associated with the potential,
                      ! we add the contribution is brings:
                      !
                      ! V = \sum_ijkl b_ijkl v_ijkl
                      ! b_ijkl = (b_i + b_j + b_k + b_l) / 4
                      if(interaction%pot_index > -1)then
                         ! get b_i (for all i, they have been precalculated)
                         call core_get_bond_order_factors(interaction%pot_index,&
                              bo_factors)
                      else
                         bo_factors = 1.d0
                      end if

                      ! If a smooth cutoff is present, we add the
                      ! contribution it brings:
                      ! 
                      ! V = \sum_ijkl v_ijkl f(r_ij) f(r_jk) f(r_kl)
                      if(interaction%smoothened)then
                         ! get f(r_ij)
                         call smoothening_factor(distances(1),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(1))
                         ! get f(r_jk)
                         call smoothening_factor(distances(2),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(2))
                         ! get f(r_kl)
                         call smoothening_factor(distances(3),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(3))
                      else
                         cut_factors(1:3) = 1.d0
                      end if

                      ! evaluate the 4-body energy
                      call evaluate_energy(4,interaction%n_product,separations(1:3,1:3),distances(1:3),&
                           interaction,tmp_energy,atom_quadruplet)

                      ! add the term: b_ijkl v_ijkl f(r_ij) f(r_jk) f(r_lk)
                      energy = energy + tmp_energy*cut_factors(1)*cut_factors(2)*cut_factors(3)*&
                           (bo_factors(index1)+bo_factors(index2)+bo_factors(index3)+bo_factors(index3))/4.d0

                   case(force_evaluation_index)

                      !**************!
                      ! 4-body force !
                      !**************!

                      ! We will need the energy contribution 
                      ! if smooth cutoffs or bond factors are used,
                      ! since we are mulplying the potential.
                      if((interaction%pot_index > -1) .or. &
                           interaction%smoothened)then
                         call evaluate_energy(4,interaction%n_product,separations(1:3,1:3),distances(1:3),&
                              interaction,tmp_energy,atom_quadruplet) ! in Potentials.f90
                      else
                         tmp_energy = 0.d0
                      end if

                      ! If a smooth cutoff is present, we add the
                      ! contribution it brings:
                      ! 
                      ! V = \sum_ijkl v_ijkl f(r_ij) f(r_jk) f(r_kl)
                      ! F_a = - \nabla_a V 
                      !     = - \sum_ijkl ( v_ijkl f'(r_ij) f(r_jk) f(r_kl) (\nabla_a r_ij) + 
                      !                   v_ijkl f(r_ij) f'(r_jk) f(r_kl) (\nabla_a r_jk) + 
                      !                   v_ijkl f(r_ij) f(r_jk) f'(r_kl) (\nabla_a r_kl) +
                      !                   (\nabla_a v_ijkl) f(r_ij) f(r_jk) f(r_kl) )
                      !
                      if(interaction%smoothened)then
                         ! get f(r_ij)
                         call smoothening_factor(distances(1),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(1)) ! in Potentials.f90
                         ! get f'(r_ij) (\nabla_a r_ij)
                         call smoothening_gradient(directions(1:3,1),distances(1),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_gradients(1:3,1)) ! in Potentials.f90 
                         ! get f(r_jk)
                         call smoothening_factor(distances(2),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(2)) ! in Potentials.f90
                         ! get f'(r_jk) (\nabla_a r_ik)
                         call smoothening_gradient(directions(1:3,2),distances(2),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_gradients(1:3,2)) ! in Potentials.f90
                         ! get f(r_kl)
                         call smoothening_factor(distances(3),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(3)) ! in Potentials.f90
                         ! get f'(r_kl) (\nabla_a r_kl)
                         call smoothening_gradient(directions(1:3,3),distances(3),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_gradients(1:3,3)) ! in Potentials.f90
                      else
                         cut_factors(1:3) = 1.d0
                         cut_gradients(1:3,1:3) = 0.d0
                      end if


                      ! If there is a bond order factor associated with the potential,
                      ! we add the contribution is brings:
                      !
                      ! V = \sum_ijkl b_ijkl v_ijkl
                      ! b_ijkl = (b_i + b_j + b_k + b_l) / 4
                      ! F_a = - \nabla_a V 
                      !     = - \sum_ijkl (\nabla_a b_ijkl) v_ijkl + b_ijkl (\nabla_a v_ijkl)
                      !     = - \sum_ijkl (\nabla_a b_ijkl) v_ijkl + b_ijkl f_a,ijkl
                      !
                      if(interaction%pot_index > -1)then
                         ! get b_i (for all i, they have been precalculated)
                         call core_get_bond_order_factors(interaction%pot_index,&
                              bo_factors) ! in Potentials.f90
                         ! get (\nabla_a b_i) (for all a)
                         call core_get_bond_order_gradients(interaction%pot_index,&
                              index1,& ! atom index
                              1, & ! slot_index
                              bo_gradients(1:3,1:n_atoms,1), &
                              bo_virial(1:6,1)) ! in Potentials.f90
                         ! get (\nabla_a b_j) (for all a)
                         call core_get_bond_order_gradients(interaction%pot_index,&
                              index2,& ! atom index
                              2, & ! slot_index
                              bo_gradients(1:3,1:n_atoms,2), &
                              bo_virial(1:6,2)) ! in Potentials.f90
                         ! get (\nabla_a b_k) (for all a)
                         call core_get_bond_order_gradients(interaction%pot_index,&
                              index3,& ! atom index
                              3, & ! slot_index
                              bo_gradients(1:3,1:n_atoms,3), &
                              bo_virial(1:6,3)) ! in Potentials.f90
                         ! get (\nabla_a b_l) (for all a)
                         call core_get_bond_order_gradients(interaction%pot_index,&
                              index4,& ! atom index
                              4, & ! slot_index
                              bo_gradients(1:3,1:n_atoms,4), &
                              bo_virial(1:6,4)) ! in Potentials.f90

                         ! Add the bond order gradient terms involving the 
                         ! atom2-atom1-atom3 energy for all atoms.
                         ! That is, add the (\nabla_a b_ijk) v_ijk term with 
                         ! the given ijk (atom2,atom1,atom3) for all a.
                         forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) &
                              - tmp_energy*cut_factors(1)*cut_factors(2)*cut_factors(3)* &
                              ( bo_gradients(1:3,1:n_atoms,1) &
                              + bo_gradients(1:3,1:n_atoms,2) &
                              + bo_gradients(1:3,1:n_atoms,3) &
                              + bo_gradients(1:3,1:n_atoms,4) )/4.d0

                         stress(1:6) = stress(1:6) - tmp_energy*cut_factors(1)*cut_factors(2)*cut_factors(3)* &
                              ( bo_virial(1:6,1) &
                              + bo_virial(1:6,2) &
                              + bo_virial(1:6,3) &
                              + bo_virial(1:6,4) )/4.d0

                      else
                         bo_factors = 1.d0
                         bo_sums = 0.d0
                         bo_gradients = 0.d0
                      end if


                      ! evaluate the 4-body force
                      call evaluate_forces(4,interaction%n_product,separations(1:3,1:3),distances(1:3),interaction,&
                           tmp_forces(1:3,1:4),atom_quadruplet) ! in Potentials.f90

                      ! force on atom 1:
                      quad_forces(1:3,1) = ( tmp_forces(1:3,1)*cut_factors(1)*cut_factors(2)*cut_factors(3) + &
                           cut_gradients(1:3,1)*cut_factors(2)*cut_factors(3)*tmp_energy ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) &
                           + bo_factors(index4) )/4.d0

                      ! force on atom 2:
                      quad_forces(1:3,2) = ( tmp_forces(1:3,2)*cut_factors(1)*cut_factors(2)*cut_factors(3) + &
                           (-cut_gradients(1:3,1)*cut_factors(2)*cut_factors(3) + &
                           cut_gradients(1:3,2)*cut_factors(1)*cut_factors(3)) * tmp_energy ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) &
                           + bo_factors(index4) )/4.d0

                      ! force on atom 3:
                      quad_forces(1:3,3) = ( tmp_forces(1:3,3)*cut_factors(1)*cut_factors(2)*cut_factors(3) + &
                           (cut_gradients(1:3,3)*cut_factors(1)*cut_factors(2) - &
                           cut_gradients(1:3,2)*cut_factors(1)*cut_factors(3)) * tmp_energy ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) &
                           + bo_factors(index4) )/4.d0

                      ! force on atom 4:
                      quad_forces(1:3,4) = ( tmp_forces(1:3,4)*cut_factors(1)*cut_factors(2)*cut_factors(3) - &
                           cut_gradients(1:3,3)*cut_factors(1)*cut_factors(2)*tmp_energy ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) &
                           + bo_factors(index4) )/4.d0

                      forces(1:3,index1) = forces(1:3,index1) + quad_forces(1:3,1)
                      forces(1:3,index2) = forces(1:3,index2) + quad_forces(1:3,2)
                      forces(1:3,index3) = forces(1:3,index3) + quad_forces(1:3,3)
                      forces(1:3,index4) = forces(1:3,index4) + quad_forces(1:3,4)
                        
                   !***************!
                   ! stress tensor !
                   !***************!
                   
                   ! s_xx, s_yy, s_zz, s_yz, s_xz, s_xy:
                   stress(1) = stress(1) + separations(1,1) * quad_forces(1,2) &
                        + (separations(1,1) + separations(1,2)) * quad_forces(1,3) &
                        + (separations(1,1) + separations(1,2) + separations(1,3)) * quad_forces(1,4)
                   stress(2) = stress(2) + separations(2,1) * quad_forces(2,2) &
                        + (separations(2,1) + separations(2,2)) * quad_forces(2,3) &
                        + (separations(2,1) + separations(2,2) + separations(2,3)) * quad_forces(2,4)
                   stress(3) = stress(3) + separations(3,1) * quad_forces(3,2) &
                        + (separations(3,1) + separations(3,2)) * quad_forces(3,3) &
                        + (separations(3,1) + separations(3,2) + separations(3,3)) * quad_forces(3,4)
                   stress(4) = stress(4) + separations(2,1) * quad_forces(3,2) &
                        + (separations(2,1) + separations(2,2)) * quad_forces(3,3) &
                        + (separations(2,1) + separations(2,2) + separations(2,3)) * quad_forces(3,4)
                   stress(5) = stress(5) + separations(1,1) * quad_forces(3,2) &
                        + (separations(1,1) + separations(1,2)) * quad_forces(3,3) &
                        + (separations(1,1) + separations(1,2) + separations(1,3)) * quad_forces(3,4)
                   stress(6) = stress(6) + separations(1,1) * quad_forces(2,2) &
                        + (separations(1,1) + separations(1,2)) * quad_forces(2,3) &
                        + (separations(1,1) + separations(1,2) + separations(1,3)) * quad_forces(2,4)


                   case(electronegativity_evaluation_index)

                      !**************************!
                      ! 4-body electronegativity !
                      !**************************!

                      ! If a smooth cutoff is present, we add the
                      ! contribution it brings:
                      if(interaction%smoothened)then
                         ! get f(r_ij)
                         call smoothening_factor(distances(1),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(1)) ! in Potentials.f90
                         ! get f(r_jk)
                         call smoothening_factor(distances(2),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(2)) ! in Potentials.f90
                         ! get f(r_kl)
                         call smoothening_factor(distances(3),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(3)) ! in Potentials.f90
                      else
                         cut_factors(1:2) = 1.d0
                      end if

                      ! If there is a bond order factor associated with the potential,
                      ! we add the contribution is brings:
                      if(interaction%pot_index > -1)then
                         ! get b_i (for all i, they have been precalculated)
                         call core_get_bond_order_factors(interaction%pot_index,&
                              bo_factors) ! in Potentials.f90
                      else
                         bo_factors = 1.d0
                      end if

                      ! evaluate the 4-body e-neg
                      call evaluate_electronegativity(4,interaction%n_product,separations(1:3,1:3),&
                           distances(1:3),interaction,&
                           tmp_enegs(1:4),atom_quadruplet) ! in Potentials.f90


                      ! e-neg on atom 1:
                      enegs(index1) = enegs(index1) + &
                           ( tmp_enegs(1)*cut_factors(1)*cut_factors(2)*cut_factors(3) ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) &
                           + bo_factors(index4) )/4.d0

                      ! e-neg on atom 2:
                      enegs(index2) = enegs(index2) + &
                           ( tmp_enegs(2)*cut_factors(1)*cut_factors(2)*cut_factors(3) ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) &
                           + bo_factors(index4) )/4.d0

                      ! e-neg on atom 3:
                      enegs(index3) = enegs(index3) + &
                           ( tmp_enegs(3)*cut_factors(1)*cut_factors(2)*cut_factors(3) ) * &
                           ( bo_factors(index1) &
                           + bo_factors(index2) &
                           + bo_factors(index3) &
                           + bo_factors(index4) )/4.d0

                   end select

                end if
             end if

          end if
       end if
    end do

  end subroutine core_evaluate_local_quadruplet





  subroutine core_evaluate_local_quadruplet_B(atom_quadruplet, &
       index1, index2, index3, index4, &
       test_index1, test_index2, test_index3, &
       interaction_indices, &
       separations, directions, distances, &
       calculation_type,energy,forces,enegs,stress, &
       many_bodies_found)
    implicit none
    integer, intent(in) :: calculation_type, index1, index2, index3, index4, &
         test_index1, test_index2, test_index3
    integer, pointer :: interaction_indices(:)
    double precision, intent(inout) :: energy, forces(:,:), enegs(:), stress(6)
    double precision, intent(in) :: separations(3,3), directions(3,3), distances(3)
    logical, intent(out) :: many_bodies_found
    type(atom), intent(in) :: atom_quadruplet(4)

    type(atom) :: atom1, atom2, atom3, atom4
    integer :: k, n_targets, n_atoms
    type(potential) :: interaction
    double precision :: tmp_energy, tmp_forces(3,4), tmp_enegs(4), &
         cut_factors(3), cut_gradients(3,3), quad_forces(3,4), bo_virial(6,4)
    logical :: is_active

    n_atoms = size(atoms)
    many_bodies_found = .false.
    atom1 = atom_quadruplet(1)
    atom2 = atom_quadruplet(2)
    atom3 = atom_quadruplet(3)
    atom4 = atom_quadruplet(3)

    ! loop over the potentials affecting atom1
    do k = 1, size(interaction_indices)

       interaction = interactions(interaction_indices(k))
       call get_number_of_targets_of_potential_index(interaction%type_index,&
            n_targets) ! in Potentials.f90
       call potential_affects_atom(interaction,atom_quadruplet(test_index1),is_active,2) ! in Potentials.f90

       ! filter the potentials by:
       ! is atom2 affected by the potential,
       ! is atom3 affected by the potential,
       ! is atom4 affected by the potential,
       ! is it a 4-body potential
       if( is_active .and. n_targets == 4 .and. interaction%cutoff > distances(1) )then
          call potential_affects_atom(interaction,atom_quadruplet(test_index2),is_active,3)  ! in Potentials.f90
          if( is_active )then
             call potential_affects_atom(interaction,atom_quadruplet(test_index3),is_active,4)  ! in Potentials.f90
             if( is_active )then
                if( interaction%cutoff > distances(2) .and. interaction%cutoff > distances(3) )then

                   ! differentiate between energy, force, and electronegativity evaluation
                   select case(calculation_type)
                   case(energy_evaluation_index)

                      !***************!
                      ! 4-body energy !
                      !***************!

                      ! If there is a bond order factor associated with the potential,
                      ! we add the contribution is brings:
                      !
                      ! V = \sum_ijkl b_ijkl v_ijkl
                      ! b_ijkl = (b_i + b_j + b_k + b_l) / 4
                      if(interaction%pot_index > -1)then
                         ! get b_i (for all i, they have been precalculated)
                         call core_get_bond_order_factors(interaction%pot_index,&
                              bo_factors)
                      else
                         !bo_factors = 1.d0
                      end if

                      ! If a smooth cutoff is present, we add the
                      ! contribution it brings:
                      ! 
                      ! V = \sum_ijkl v_ijkl f(r_ij) f(r_jk) f(r_kl)
                      if(interaction%smoothened)then
                         ! get f(r_ij)
                         call smoothening_factor(distances(1),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(1))
                         ! get f(r_jk)
                         call smoothening_factor(distances(2),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(2))
                         ! get f(r_kl)
                         call smoothening_factor(distances(3),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(3))
                      else
                         cut_factors(1:3) = 1.d0
                      end if

                      ! evaluate the 4-body energy
                      call evaluate_energy(4,interaction%n_product,separations(1:3,1:3),distances(1:3),&
                           interaction,tmp_energy,atom_quadruplet)

                      if(interaction%pot_index > -1)then
                         ! add the term: b_ijkl v_ijkl f(r_ij) f(r_jk) f(r_lk)
                         energy = energy + tmp_energy*cut_factors(1)*cut_factors(2)*cut_factors(3)*&
                              (bo_factors(index1)+bo_factors(index2)+bo_factors(index3)+bo_factors(index3))/4.d0
                      else
                         energy = energy + tmp_energy*cut_factors(1)*cut_factors(2)*cut_factors(3)
                      end if


                   case(force_evaluation_index)

                      !**************!
                      ! 4-body force !
                      !**************!

                      ! We will need the energy contribution 
                      ! if smooth cutoffs or bond factors are used,
                      ! since we are mulplying the potential.
                      if((interaction%pot_index > -1) .or. &
                           interaction%smoothened)then
                         call evaluate_energy(4,interaction%n_product,separations(1:3,1:3),distances(1:3),&
                              interaction,tmp_energy,atom_quadruplet) ! in Potentials.f90
                      else
                         tmp_energy = 0.d0
                      end if

                      ! If a smooth cutoff is present, we add the
                      ! contribution it brings:
                      ! 
                      ! V = \sum_ijkl v_ijkl f(r_ij) f(r_jk) f(r_kl)
                      ! F_a = - \nabla_a V 
                      !     = - \sum_ijkl ( v_ijkl f'(r_ij) f(r_jk) f(r_kl) (\nabla_a r_ij) + 
                      !                   v_ijkl f(r_ij) f'(r_jk) f(r_kl) (\nabla_a r_jk) + 
                      !                   v_ijkl f(r_ij) f(r_jk) f'(r_kl) (\nabla_a r_kl) +
                      !                   (\nabla_a v_ijkl) f(r_ij) f(r_jk) f(r_kl) )
                      !
                      if(interaction%smoothened)then
                         ! get f(r_ij)
                         call smoothening_factor(distances(1),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(1)) ! in Potentials.f90
                         ! get f'(r_ij) (\nabla_a r_ij)
                         call smoothening_gradient(directions(1:3,1),distances(1),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_gradients(1:3,1)) ! in Potentials.f90 
                         ! get f(r_jk)
                         call smoothening_factor(distances(2),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(2)) ! in Potentials.f90
                         ! get f'(r_jk) (\nabla_a r_ik)
                         call smoothening_gradient(directions(1:3,2),distances(2),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_gradients(1:3,2)) ! in Potentials.f90
                         ! get f(r_kl)
                         call smoothening_factor(distances(3),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(3)) ! in Potentials.f90
                         ! get f'(r_kl) (\nabla_a r_kl)
                         call smoothening_gradient(directions(1:3,3),distances(3),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_gradients(1:3,3)) ! in Potentials.f90
                      else
                         cut_factors(1:3) = 1.d0
                         cut_gradients(1:3,1:3) = 0.d0
                      end if


                      ! If there is a bond order factor associated with the potential,
                      ! we add the contribution is brings:
                      !
                      ! V = \sum_ijkl b_ijkl v_ijkl
                      ! b_ijkl = (b_i + b_j + b_k + b_l) / 4
                      ! F_a = - \nabla_a V 
                      !     = - \sum_ijkl (\nabla_a b_ijkl) v_ijkl + b_ijkl (\nabla_a v_ijkl)
                      !     = - \sum_ijkl (\nabla_a b_ijkl) v_ijkl + b_ijkl f_a,ijkl
                      !
                      if(interaction%pot_index > -1)then
                         ! get b_i (for all i, they have been precalculated)
                         call core_get_bond_order_factors(interaction%pot_index,&
                              bo_factors) ! in Potentials.f90
                         ! get (\nabla_a b_i) (for all a)
                         call core_get_bond_order_gradients(interaction%pot_index,&
                              index1,& ! atom index
                              1, & ! slot_index
                              bo_gradients(1:3,1:n_atoms,1), &
                              bo_virial(1:6,1)) ! in Potentials.f90
                         ! get (\nabla_a b_j) (for all a)
                         call core_get_bond_order_gradients(interaction%pot_index,&
                              index2,& ! atom index
                              2, & ! slot_index
                              bo_gradients(1:3,1:n_atoms,2), &
                              bo_virial(1:6,2)) ! in Potentials.f90
                         ! get (\nabla_a b_k) (for all a)
                         call core_get_bond_order_gradients(interaction%pot_index,&
                              index3,& ! atom index
                              3, & ! slot_index
                              bo_gradients(1:3,1:n_atoms,3), &
                              bo_virial(1:6,3)) ! in Potentials.f90
                         ! get (\nabla_a b_l) (for all a)
                         call core_get_bond_order_gradients(interaction%pot_index,&
                              index4,& ! atom index
                              4, & ! slot_index
                              bo_gradients(1:3,1:n_atoms,4), &
                              bo_virial(1:6,4)) ! in Potentials.f90

                         ! Add the bond order gradient terms involving the 
                         ! atom2-atom1-atom3 energy for all atoms.
                         ! That is, add the (\nabla_a b_ijk) v_ijk term with 
                         ! the given ijk (atom2,atom1,atom3) for all a.
                         forces(1:3,1:n_atoms) = forces(1:3,1:n_atoms) &
                              - tmp_energy*cut_factors(1)*cut_factors(2)*cut_factors(3)* &
                              ( bo_gradients(1:3,1:n_atoms,1) &
                              + bo_gradients(1:3,1:n_atoms,2) &
                              + bo_gradients(1:3,1:n_atoms,3) &
                              + bo_gradients(1:3,1:n_atoms,4) )/4.d0

                         stress(1:6) = stress(1:6) - tmp_energy*cut_factors(1)*cut_factors(2)*cut_factors(3)* &
                              ( bo_virial(1:6,1) &
                              + bo_virial(1:6,2) &
                              + bo_virial(1:6,3) &
                              + bo_virial(1:6,4) )/4.d0

                      else
                         !bo_factors = 1.d0
                         !bo_sums = 0.d0
                         !bo_gradients = 0.d0
                      end if


                      ! evaluate the 4-body force
                      call evaluate_forces(4,interaction%n_product,separations(1:3,1:3),distances(1:3),interaction,&
                           tmp_forces(1:3,1:4),atom_quadruplet) ! in Potentials.f90

                      if(interaction%pot_index > -1)then
                         ! force on atom 1:
                         quad_forces(1:3,1) = ( tmp_forces(1:3,1)*cut_factors(1)*cut_factors(2)*cut_factors(3) + &
                              cut_gradients(1:3,1)*cut_factors(2)*cut_factors(3)*tmp_energy ) * &
                              ( bo_factors(index1) &
                              + bo_factors(index2) &
                              + bo_factors(index3) &
                              + bo_factors(index4) )/4.d0
                         
                         ! force on atom 2:
                         quad_forces(1:3,2) = ( tmp_forces(1:3,2)*cut_factors(1)*cut_factors(2)*cut_factors(3) + &
                              (-cut_gradients(1:3,1)*cut_factors(2)*cut_factors(3) + &
                              cut_gradients(1:3,2)*cut_factors(1)*cut_factors(3)) * tmp_energy ) * &
                              ( bo_factors(index1) &
                              + bo_factors(index2) &
                              + bo_factors(index3) &
                              + bo_factors(index4) )/4.d0
                         
                         ! force on atom 3:
                         quad_forces(1:3,3) = ( tmp_forces(1:3,3)*cut_factors(1)*cut_factors(2)*cut_factors(3) + &
                              (cut_gradients(1:3,3)*cut_factors(1)*cut_factors(2) - &
                              cut_gradients(1:3,2)*cut_factors(1)*cut_factors(3)) * tmp_energy ) * &
                              ( bo_factors(index1) &
                              + bo_factors(index2) &
                              + bo_factors(index3) &
                              + bo_factors(index4) )/4.d0
                         
                         ! force on atom 4:
                         quad_forces(1:3,4) = ( tmp_forces(1:3,4)*cut_factors(1)*cut_factors(2)*cut_factors(3) - &
                              cut_gradients(1:3,3)*cut_factors(1)*cut_factors(2)*tmp_energy ) * &
                              ( bo_factors(index1) &
                              + bo_factors(index2) &
                              + bo_factors(index3) &
                              + bo_factors(index4) )/4.d0

                      else
                         ! force on atom 1:
                         quad_forces(1:3,1) = ( tmp_forces(1:3,1)*cut_factors(1)*cut_factors(2)*cut_factors(3) + &
                              cut_gradients(1:3,1)*cut_factors(2)*cut_factors(3)*tmp_energy )
                         
                         ! force on atom 2:
                         quad_forces(1:3,2) = ( tmp_forces(1:3,2)*cut_factors(1)*cut_factors(2)*cut_factors(3) + &
                              (-cut_gradients(1:3,1)*cut_factors(2)*cut_factors(3) + &
                              cut_gradients(1:3,2)*cut_factors(1)*cut_factors(3)) * tmp_energy ) 
                         
                         ! force on atom 3:
                         quad_forces(1:3,3) = ( tmp_forces(1:3,3)*cut_factors(1)*cut_factors(2)*cut_factors(3) + &
                              (cut_gradients(1:3,3)*cut_factors(1)*cut_factors(2) - &
                              cut_gradients(1:3,2)*cut_factors(1)*cut_factors(3)) * tmp_energy ) 
                         
                         ! force on atom 4:
                         quad_forces(1:3,4) = ( tmp_forces(1:3,4)*cut_factors(1)*cut_factors(2)*cut_factors(3) - &
                              cut_gradients(1:3,3)*cut_factors(1)*cut_factors(2)*tmp_energy ) 

                      end if

                      forces(1:3,index1) = forces(1:3,index1) + quad_forces(1:3,1)
                      forces(1:3,index2) = forces(1:3,index2) + quad_forces(1:3,2)
                      forces(1:3,index3) = forces(1:3,index3) + quad_forces(1:3,3)
                      forces(1:3,index4) = forces(1:3,index4) + quad_forces(1:3,4)
                        
                   !***************!
                   ! stress tensor !
                   !***************!
                   
                   ! s_xx, s_yy, s_zz, s_yz, s_xz, s_xy:
                   stress(1) = stress(1) + separations(1,1) * quad_forces(1,2) &
                        + (separations(1,1) + separations(1,2)) * quad_forces(1,3) &
                        + (separations(1,1) + separations(1,2) + separations(1,3)) * quad_forces(1,4)
                   stress(2) = stress(2) + separations(2,1) * quad_forces(2,2) &
                        + (separations(2,1) + separations(2,2)) * quad_forces(2,3) &
                        + (separations(2,1) + separations(2,2) + separations(2,3)) * quad_forces(2,4)
                   stress(3) = stress(3) + separations(3,1) * quad_forces(3,2) &
                        + (separations(3,1) + separations(3,2)) * quad_forces(3,3) &
                        + (separations(3,1) + separations(3,2) + separations(3,3)) * quad_forces(3,4)
                   stress(4) = stress(4) + separations(2,1) * quad_forces(3,2) &
                        + (separations(2,1) + separations(2,2)) * quad_forces(3,3) &
                        + (separations(2,1) + separations(2,2) + separations(2,3)) * quad_forces(3,4)
                   stress(5) = stress(5) + separations(1,1) * quad_forces(3,2) &
                        + (separations(1,1) + separations(1,2)) * quad_forces(3,3) &
                        + (separations(1,1) + separations(1,2) + separations(1,3)) * quad_forces(3,4)
                   stress(6) = stress(6) + separations(1,1) * quad_forces(2,2) &
                        + (separations(1,1) + separations(1,2)) * quad_forces(2,3) &
                        + (separations(1,1) + separations(1,2) + separations(1,3)) * quad_forces(2,4)


                   case(electronegativity_evaluation_index)

                      !**************************!
                      ! 4-body electronegativity !
                      !**************************!

                      ! If a smooth cutoff is present, we add the
                      ! contribution it brings:
                      if(interaction%smoothened)then
                         ! get f(r_ij)
                         call smoothening_factor(distances(1),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(1)) ! in Potentials.f90
                         ! get f(r_jk)
                         call smoothening_factor(distances(2),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(2)) ! in Potentials.f90
                         ! get f(r_kl)
                         call smoothening_factor(distances(3),&
                              interaction%cutoff,interaction%soft_cutoff,&
                              cut_factors(3)) ! in Potentials.f90
                      else
                         cut_factors(1:2) = 1.d0
                      end if

                      ! If there is a bond order factor associated with the potential,
                      ! we add the contribution is brings:
                      if(interaction%pot_index > -1)then
                         ! get b_i (for all i, they have been precalculated)
                         call core_get_bond_order_factors(interaction%pot_index,&
                              bo_factors) ! in Potentials.f90
                      else
                         !bo_factors = 1.d0
                      end if

                      ! evaluate the 4-body e-neg
                      call evaluate_electronegativity(4,interaction%n_product,separations(1:3,1:3),&
                           distances(1:3),interaction,&
                           tmp_enegs(1:4),atom_quadruplet) ! in Potentials.f90


                      if(interaction%pot_index > -1)then
                         ! e-neg on atom 1:
                         enegs(index1) = enegs(index1) + &
                              ( tmp_enegs(1)*cut_factors(1)*cut_factors(2)*cut_factors(3) ) * &
                              ( bo_factors(index1) &
                              + bo_factors(index2) &
                              + bo_factors(index3) &
                              + bo_factors(index4) )/4.d0
                         
                         ! e-neg on atom 2:
                         enegs(index2) = enegs(index2) + &
                              ( tmp_enegs(2)*cut_factors(1)*cut_factors(2)*cut_factors(3) ) * &
                              ( bo_factors(index1) &
                              + bo_factors(index2) &
                              + bo_factors(index3) &
                              + bo_factors(index4) )/4.d0
                         
                         ! e-neg on atom 3:
                         enegs(index3) = enegs(index3) + &
                              ( tmp_enegs(3)*cut_factors(1)*cut_factors(2)*cut_factors(3) ) * &
                              ( bo_factors(index1) &
                              + bo_factors(index2) &
                              + bo_factors(index3) &
                              + bo_factors(index4) )/4.d0

                         ! e-neg on atom 4:
                         enegs(index4) = enegs(index4) + &
                              ( tmp_enegs(4)*cut_factors(1)*cut_factors(2)*cut_factors(3) ) * &
                              ( bo_factors(index1) &
                              + bo_factors(index2) &
                              + bo_factors(index3) &
                              + bo_factors(index4) )/4.d0
                      else
                         ! e-neg on atom 1:
                         enegs(index1) = enegs(index1) + &
                              ( tmp_enegs(1)*cut_factors(1)*cut_factors(2)*cut_factors(3) )
                         
                         ! e-neg on atom 2:
                         enegs(index2) = enegs(index2) + &
                              ( tmp_enegs(2)*cut_factors(1)*cut_factors(2)*cut_factors(3) ) 
                         
                         ! e-neg on atom 3:
                         enegs(index3) = enegs(index3) + &
                              ( tmp_enegs(3)*cut_factors(1)*cut_factors(2)*cut_factors(3) ) 

                         ! e-neg on atom 4:
                         enegs(index4) = enegs(index4) + &
                              ( tmp_enegs(4)*cut_factors(1)*cut_factors(2)*cut_factors(3) ) 
                      end if

                   end select

                end if
             end if

          end if
       end if
    end do

  end subroutine core_evaluate_local_quadruplet_B






! !!!: core_calculate_forces

  ! Calculates forces acting on all atoms of the system.
  !
  ! The routine calculates the potential gradient
  !
  ! .. math::
  !
  !    \mathbf{F}_\alpha = - \nabla_\alpha V
  !
  ! for all atoms :math:`\alpha`. This is done according to the
  ! the structure and potentials allocated in the core, so the
  ! routine does not accept arguments. Instead, the core modifying
  ! routines such as :func:`core_generate_atoms` must be called
  ! first to set up the calculation.
  !
  ! called from PyInterface: :func:`calculate_forces`
  !
  ! *total_forces an array containing the calculated forces for all atoms
  ! *total_stress as array containing the calculated stress tensor
  subroutine core_calculate_forces(total_forces,total_stress)
    implicit none
    double precision, intent(inout) :: total_forces(:,:), total_stress(6)
    double precision :: dummy_energy

    call core_loop_over_local_interactions(force_evaluation_index,&
         dummy_energy,total_forces,temp_enegs,total_stress)

  end subroutine core_calculate_forces





! !!!: core_calculate_electronegativities

  ! Calculates electronegativity forces acting on all atomic charges of the system.
  !
  ! The routine calculates the electronegativities
  !
  ! .. math::
  !
  !    \chi_{\alpha} = -\frac{\partial V}{\partial q_\alpha}
  !
  ! for all atoms :math:`\alpha`. This is done according to the
  ! the structure and potentials allocated in the core, so the
  ! routine does not accept arguments. Instead, the core modifying
  ! routines such as :func:`core_generate_atoms` must be called
  ! first to set up the calculation.
  !
  ! called from PyInterface: :func:`calculate_electronegativities`
  !
  ! *total_enegs an array containing the calculated charge forces for all atoms
  subroutine core_calculate_electronegativities(total_enegs)
    implicit none
    double precision, intent(inout) :: total_enegs(:)
    double precision :: dummy_energy, dummy_stress(6)

    call core_loop_over_local_interactions(electronegativity_evaluation_index,&
         dummy_energy,temp_forces,total_enegs,dummy_stress)
    return

  end subroutine core_calculate_electronegativities




! !!!: core_calculate_energy

  ! Calculates the total potential energy of the system.
  !
  ! This is done according to the
  ! the structure and potentials allocated in the core, so the
  ! routine does not accept arguments. Instead, the core modifying
  ! routines such as :func:`core_generate_atoms` must be called
  ! first to set up the calculation.
  !
  ! called from PyInterface: :func:`calculate_energy`
  !
  ! *total_energy calculated total potential energy
  subroutine core_calculate_energy(total_energy)
    implicit none
    double precision, intent(out) :: total_energy
    double precision :: dummy_stress(6)

    call core_loop_over_local_interactions(energy_evaluation_index,&
         total_energy,temp_forces,temp_enegs,dummy_stress)
    return

  end subroutine core_calculate_energy


  ! !!!: list_atoms

  ! Prints some information on the atoms stored in the core in stdout.
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


! !!!: list_cell

  ! Prints some information on the supercell stored in the core in stdout.
  subroutine list_cell()
    implicit none
    integer :: i

    write(*,*) "cell"
    do i = 1, 3
        write(*,'(F10.4,F10.4,F10.4)') cell%vectors(1:3,i)
    end do
    write(*,*) "pbc: ", cell%periodic(1:3)

  end subroutine list_cell


! !!!: list_interactions

  ! Prints some information on the potentials stored in the core in stdout.
  subroutine list_interactions()
    implicit none
    integer :: i, j
    
    write(*,*) "interactions"
    do i = 1, n_interactions
       write(*,'(A,I5,F10.4)') "type, cutoff ", interactions(i)%type_index, interactions(i)%cutoff
       write(*,*) "params ", interactions(i)%parameters
       if(interactions(i)%filter_elements)then
          write(*,*) "         symbols ", interactions(i)%apply_elements
          write(*,*) "original symbols ", interactions(i)%original_elements
       end if
       if(interactions(i)%filter_tags)then
          write(*,*) "         tags ", interactions(i)%apply_tags
          write(*,*) "original tags ", interactions(i)%original_tags
       end if
       if(interactions(i)%filter_indices)then
          write(*,*) "         indices ", interactions(i)%apply_indices
          write(*,*) "original indices ", interactions(i)%original_indices
       end if
       if(interactions(i)%pot_index >= 0)then
          write(*,*) "bond order index ", interactions(i)%pot_index
       end if
       if(size(interactions(i)%multipliers) > 0)then
          do j = 1, size(interactions(i)%multipliers)
             write(*,*) "      multiplier ", interactions(i)%multipliers(j)%type_index
          end do
       end if
       write(*,*) ""
    end do

  end subroutine list_interactions


! !!!: list_bonds

  ! Prints some information on the bond order factors stored in the core in stdout.
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

  ! Debug routine for Ewald
  subroutine core_get_ewald_energy(real_cut, k_cut, reciprocal_cut, sigma, epsilon, energy)
    implicit none
    double precision, intent(in) :: real_cut, k_cut, sigma, epsilon
    integer, intent(in) :: reciprocal_cut(3)
    double precision, intent(out) :: energy
    double precision :: scaler(size(atoms))

    scaler = 1.d0

    call calculate_ewald_energy(atoms,cell,real_cut,k_cut,reciprocal_cut,sigma,&
         epsilon,scaler,.true.,energy)

  end subroutine core_get_ewald_energy



  ! Sets the parameters for Ewald summation in the core. 
  !
  ! *real_cut the real-space cutoff
  ! *reciprocal_cut the k-space cutoffs
  ! *sigma the split parameter
  ! *epsilon electric constant
  ! *scaler scaling factors for the individual charges
  subroutine core_set_ewald_parameters(real_cut, k_radius, reciprocal_cut, sigma, epsilon, scaler)
    implicit none
    double precision, intent(in) :: real_cut, k_radius, sigma, epsilon, scaler(:)
    integer, intent(in) :: reciprocal_cut(3)

    evaluate_ewald = .true.
    ewald_k_cutoffs = reciprocal_cut
    ewald_k_radius = k_radius
    ewald_cutoff = real_cut
    ewald_sigma = sigma
    ewald_epsilon = epsilon
    if(ewald_allocated)then
       deallocate(ewald_scaler)
    else
       nullify(ewald_scaler)
    end if
    allocate(ewald_scaler(size(atoms)))
    ewald_scaler = scaler
    ewald_allocated = .true.
    
    call deallocate_ewald_arrays()
    call allocate_ewald_arrays(size(atoms))

  end subroutine core_set_ewald_parameters


  ! Partitions the simulation volume in subvolumes for fast neighbor searching
  !
  ! *max_cutoff the maximum cutoff radius for neighbor search
  subroutine core_create_space_partitioning(max_cutoff)
    implicit none
    double precision, intent(in) :: max_cutoff
    integer :: splits(3), i,j,k,in,jn,kn

    call get_optimal_splitting(cell,max_cutoff,splits)
    call divide_cell(cell,splits)
    do i = 1, size(atoms)
       call find_subcell_for_atom(cell,atoms(i))
    end do

  end subroutine core_create_space_partitioning

  ! Builds the neighbor lists in the core.
  ! The simulation cell must be partitioned with :func:`core_create_space_partitioning` 
  ! before this routine can be called.
  !
  ! *cutoffs list of cutoffs, atom by atom
  subroutine core_build_neighbor_lists(cutoffs)
    implicit none
    double precision, intent(in) :: cutoffs(:)
    integer :: cell_indices(3), i,j, i_n, j_n, k_n, neighbor_offset(3), atom1_index, atom2_index, &
         max_n_nbors=100, global_max, tmp_max, &
         atom1_wrap_offset(3), atom2_wrap_offset(3), n_atoms
    integer, pointer, save ::  nbors_and_offsets(:,:,:)
    type(subcell) :: atom_cell, neighbor_cell
    logical :: neighbor_include, first_run = .true.
    double precision :: separation(3), distance, dummy1(3)
    
    n_atoms = size(atoms) 
    if(first_run)then
       nullify(nbors_and_offsets)
       allocate(nbors_and_offsets(4,max_n_nbors,n_atoms))
       first_run = .false.
    else if(size(nbors_and_offsets(1,1,:)) /= n_atoms)then
       deallocate(nbors_and_offsets)
       allocate(nbors_and_offsets(4,max_n_nbors,n_atoms))
    end if

    n_nbs = 0
    nbors_and_offsets = 0

    do atom1_index = 1, size(atoms)

       if(is_my_atom(atom1_index))then
          
          cell_indices = atoms(atom1_index)%subcell_indices
          atom_cell = cell%subcells(cell_indices(1),cell_indices(2),cell_indices(3))
          call wrapped_coordinates(atoms(atom1_index)%position,cell,dummy1,atom1_wrap_offset)

          do k_n = -1,1
             do j_n = -1,1
                do i_n = -1,1
                   
                   cell_indices(1:3) = atom_cell%neighbors(1:3,i_n,j_n,k_n)

                   neighbor_cell = cell%subcells(cell_indices(1),cell_indices(2),cell_indices(3))
                   neighbor_offset = atom_cell%offsets(1:3,i_n,j_n,k_n)
                   neighbor_include = atom_cell%include(i_n,j_n,k_n)               

                   if(neighbor_include)then
                      
                      do j = 1, neighbor_cell%n_atoms
                         atom2_index = neighbor_cell%atoms(j)
                         call wrapped_coordinates(atoms(atom2_index)%position,cell,dummy1,atom2_wrap_offset)

                         ! prevent double counting
                         if(pick(atom1_index,atom2_index,neighbor_offset))then
                            call separation_vector(atoms(atom1_index)%position, &
                                 atoms(atom2_index)%position, &
                                 neighbor_offset - atom1_wrap_offset + atom2_wrap_offset, &
                                 cell, &
                                 separation)
                            ! distance squared
                            distance = separation.o.separation

                            ! atom2 is neighbor of atom1
                            if(distance < cutoffs(atom1_index)*cutoffs(atom1_index))then
                               n_nbs(atom1_index) = n_nbs(atom1_index)+1
                               nbors_and_offsets(1,n_nbs(atom1_index),atom1_index) = atom2_index
                               nbors_and_offsets(2:4,n_nbs(atom1_index),atom1_index) = neighbor_offset(1:3) &
                                     - atom1_wrap_offset(1:3) + atom2_wrap_offset(1:3)
                               if(n_nbs(atom1_index) == max_n_nbors)then
                                  call expand_neighbor_storage(nbors_and_offsets,max_n_nbors,max_n_nbors+50,n_atoms)
                                  max_n_nbors = max_n_nbors + 50
                               end if
                            end if

                            ! atom1 is neighbor of atom2
                            if(distance < cutoffs(atom2_index)*cutoffs(atom2_index))then
                               n_nbs(atom2_index) = n_nbs(atom2_index)+1
                               nbors_and_offsets(1,n_nbs(atom2_index),atom2_index) = atom1_index
                               nbors_and_offsets(2:4,n_nbs(atom2_index),atom2_index) = -neighbor_offset &
                                     + atom1_wrap_offset(1:3) - atom2_wrap_offset(1:3)
                               if(n_nbs(atom2_index) == max_n_nbors)then
                                  call expand_neighbor_storage(nbors_and_offsets,max_n_nbors,max_n_nbors+50,n_atoms)
                                  max_n_nbors = max_n_nbors + 50
                               end if
                            end if

                         end if ! pick

                      end do ! j = 1, neighbor_cell%n_atoms
                      
                   end if
                   
                end do
             end do
          end do
          
          
       end if ! is_my_atom
    end do ! atom1_index
    
#ifdef MPI

    ! get the total array size by summing the neighbors found by all cpus
    call mpi_allreduce(n_nbs,total_n_nbs,n_atoms,mpi_integer,mpi_sum,mpi_comm_world,mpistat)
    tmp_max = max(max_n_nbors, maxval(total_n_nbs))
    ! get the maximum array size (in case some array expanded during search)
    call mpi_allreduce(tmp_max,global_max,1,mpi_integer,mpi_max,mpi_comm_world,mpistat)
    if(global_max > max_n_nbors)then
       call expand_neighbor_storage(nbors_and_offsets,max_n_nbors,global_max,n_atoms)
       max_n_nbors = global_max
    end if

    ! stack the lists on cpu 0 and broadcast them to all cpus
    call mpi_stack(nbors_and_offsets,n_nbs,4,n_atoms,max_n_nbors)
    call mpi_bcast(nbors_and_offsets,4*n_atoms*max_n_nbors,mpi_integer,0,mpi_comm_world,mpistat)
    call mpi_bcast(n_nbs,n_atoms,mpi_integer,0,mpi_comm_world,mpistat)


#endif

    do atom1_index = 1, n_atoms
       call core_create_neighbor_list(n_nbs(atom1_index),atom1_index,&
            nbors_and_offsets(1,1:n_nbs(atom1_index),atom1_index),&
            nbors_and_offsets(2:4,1:n_nbs(atom1_index),atom1_index))
    end do

  end subroutine core_build_neighbor_lists


  ! Expands the allocated memory for storing neighbor lists
  subroutine expand_neighbor_storage(nbors_and_offsets,length,new_length,n_atoms)
    implicit none
    integer, intent(in) :: n_atoms
    integer, intent(in) :: length, new_length
    integer, pointer ::  nbors_and_offsets(:,:,:), tmp_nbors(:,:,:)

    nullify(tmp_nbors)
    allocate(tmp_nbors(4,length,n_atoms))
    tmp_nbors(1:4,1:length,1:n_atoms) = nbors_and_offsets(1:4,1:length,1:n_atoms)
    deallocate(nbors_and_offsets)
    allocate(nbors_and_offsets(4,new_length,n_atoms))
    nbors_and_offsets = 0
    nbors_and_offsets(1:4,1:length,1:n_atoms) = tmp_nbors(1:4,1:length,1:n_atoms)
    deallocate(tmp_nbors)

  end subroutine expand_neighbor_storage


  ! Returns the number of neighbors for an atom
  !
  ! *atom_index the index of the atoms
  ! *n_neighbors the number of neighbors
  subroutine core_get_number_of_neighbors(atom_index,n_neighbors)
    implicit none
    integer, intent(in) :: atom_index
    integer, intent(out) :: n_neighbors

    n_neighbors = atoms(atom_index)%neighbor_list%n_neighbors

  end subroutine core_get_number_of_neighbors


  ! Returns the list of neighbros for an atom
  !
  ! *atom_index the index of the atom whose neighbors are returned
  ! *n_neighbors the number of neighbors
  ! *neighbors the indices of the neighboring atoms
  ! *offsets the offsets for periodic boundaries
  subroutine core_get_neighbor_list_of_atom(atom_index, n_neighbors, neighbors, offsets)
    implicit none
    integer, intent(in) :: atom_index, n_neighbors
    integer, intent(out) :: neighbors(n_neighbors), offsets(3,n_neighbors)

    neighbors(1:n_neighbors) = atoms(atom_index)%neighbor_list%neighbors(1:n_neighbors)
    offsets(1:3,1:n_neighbors) = atoms(atom_index)%neighbor_list%pbc_offsets(1:3,1:n_neighbors)

  end subroutine core_get_neighbor_list_of_atom


  ! Write atomic coordinates and other info in a file.
  ! This is only for debugging.
  subroutine core_debug_dump(forces)
    implicit none
    double precision, intent(in) :: forces(:,:)
    integer :: step = 10
    character(len=13) :: filename
    integer :: i, j, k, channel, lowest, lowest_index
    double precision :: separation(3)
    type(neighbor_list) :: nbors
    logical, allocatable :: taken(:)
    

    step = step+1
    write(*,*) "atom 1 at step ",step,"on cpu ", cpu_id, atoms(1)%position(1:3)
    write(filename,'(A5,I1,A1,I2,A4)') "dump_",cpu_id,"_",step,".txt"
    channel = 1234+cpu_id*100+step*10
    open(channel,FILE=filename)

    do i = 1, size(atoms)
       write(channel,'(I10,F20.10,F20.10,F20.10)') i, atoms(i)%position(1:3)
    end do

    write(channel,*) ""

    do i = 1, size(atoms)
       write(channel,'(I10,F20.10,F20.10,F20.10)') i, forces(1:3,i)
    end do

    write(channel,*) ""

    do i = 1, size(atoms)
       nbors = atoms(i)%neighbor_list
       write(channel,'(I4,I5)',advance='no') i, nbors%n_neighbors
       allocate(taken(nbors%n_neighbors))
       taken = .false.
       do j = 1, nbors%n_neighbors      
          lowest = 999999999
          do k = 1, nbors%n_neighbors      
             if( nbors%neighbors(k) < lowest .and. .not.taken(k) )then
                lowest = nbors%neighbors(k)
                lowest_index = k
             end if
          end do
          taken(lowest_index) = .true.
          call separation_vector(atoms(i)%position, &
               atoms(nbors%neighbors(lowest_index))%position, &
               nbors%pbc_offsets(1:3,lowest_index), &
               cell, &
               separation)
          write(channel,'(I4,F9.3)',advance='no') nbors%neighbors(lowest_index), .norm.(separation)
       end do
       write(channel,'(A)') ""
       deallocate(taken)
    end do

    close(channel)    

  end subroutine core_debug_dump

end module pysic_core
