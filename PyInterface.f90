module pysic_interface
  use pysic_core
  use geometry
  use utility
  use potentials
  use mpi
  use mt95
  implicit none

  integer :: owner_id = -1

contains

  subroutine get_owner_id(id)
    implicit none
    integer, intent(out) :: id
    
    id = owner_id
    
  end subroutine get_owner_id


  subroutine set_owner_id(id)
    implicit none
    integer, intent(in) :: id

    owner_id = id

  end subroutine set_owner_id


  subroutine can_be_accessed(id, free)
    implicit none
    integer, intent(in) :: id
    logical, intent(out) :: free

    if(owner_id < 0)then
       free = .true.
       return
    else if(id == owner_id)then
       free = .true.
       return
    else
       free = .false.
       return
    end if

  end subroutine can_be_accessed


  ! Initialize Mersenne Twister rng
  subroutine start_rng(seed)
    implicit none
    integer, intent(in) :: seed
    integer :: synced_seed
    
    synced_seed = seed
    ! if we have many cpus it's important they all share the random seed
    call mpi_master_bcast_int(synced_seed)
    call genrand_init(synced_seed)

  end subroutine start_rng


  ! Initializes MPI for parallel calculations
  subroutine start_mpi()
    implicit none

    call mpi_initialize()

  end subroutine start_mpi


  ! Finishes MPI for parallel calculations
  subroutine finish_mpi()
    implicit none

    call mpi_finish()

  end subroutine finish_mpi


  subroutine sync_mpi()
    implicit none

    call mpi_sync()

  end subroutine sync_mpi

  ! Distributes atoms among the processors
  subroutine distribute_mpi(n_atoms)
    implicit none
    integer, intent(in) :: n_atoms

    call mpi_distribute(n_atoms)

  end subroutine distribute_mpi


  ! Returns the MPI cpu id number
  subroutine get_cpu_id(id)
    implicit none
    integer, intent(out) :: id

    id = cpu_id

  end subroutine get_cpu_id


  ! Returns the MPI cpu count
  subroutine get_number_of_cpus(ncpu)
    implicit none
    integer, intent(out) :: ncpu

    ncpu = n_cpus

  end subroutine get_number_of_cpus


  subroutine get_mpi_list_of_atoms(n_atoms,cpu_atoms)
    implicit none
    integer, intent(in) :: n_atoms
    logical, intent(out) :: cpu_atoms(n_atoms)
    
    cpu_atoms(1:n_atoms) = is_my_atom(1:n_atoms)

  end subroutine get_mpi_list_of_atoms


  ! Initializes the potentials module
  subroutine start_potentials()
    implicit none

    call initialize_potential_characterizers()

  end subroutine start_potentials


  ! Initializes the potentials module
  subroutine start_bond_order_factors()
    implicit none

    call initialize_bond_order_factor_characterizers()

  end subroutine start_bond_order_factors


  ! Creates a supercell for containing the calculation geometry
  subroutine create_cell(vectors,inverse,periodicity)
    implicit none
    double precision, intent(in) :: vectors(3,3), inverse(3,3)
    logical, intent(in) :: periodicity(3)

    call core_create_cell(vectors,inverse,periodicity)

  end subroutine create_cell


  subroutine get_cell_vectors(vectors)
    implicit none
    double precision, intent(out) :: vectors(3,3)

    call core_get_cell_vectors(vectors)
    
  end subroutine get_cell_vectors


  ! Creates atomic particles
  subroutine create_atoms(n_atoms,masses,charges,positions,momenta,tags,elements)
    implicit none
    integer, intent(in) :: n_atoms, elements(2,n_atoms), tags(n_atoms)
    double precision, intent(in) :: masses(n_atoms), charges(n_atoms), positions(3,n_atoms), &
         momenta(3,n_atoms)
    character(len=2) :: elements_str(n_atoms)
    integer :: i

    ! translates the integer-formatted labels to characters
    do i = 1, n_atoms
       call int2str(2,elements(:,i),elements_str(i))
    end do
    ! tells the Fortran core to create the atoms
    call core_generate_atoms(n_atoms,masses,charges,positions,momenta,tags,elements_str)

  end subroutine create_atoms


  ! Updates the positions and velocities of existing atoms
  subroutine update_atom_coordinates(n_atoms,positions,momenta)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(in) :: positions(3,n_atoms), momenta(3,n_atoms)

    call core_update_atom_coordinates(n_atoms,positions,momenta)

  end subroutine update_atom_coordinates


  ! Counts the number of atoms in the current core
  subroutine get_number_of_atoms(n_atoms)
    implicit none
    integer, intent(out) :: n_atoms

    call core_get_number_of_atoms(n_atoms)

  end subroutine get_number_of_atoms


  ! Returns the wanted properties of atomic particles
  subroutine examine_atoms()
    implicit none
    call list_atoms()
  end subroutine examine_atoms

  subroutine examine_cell()
    implicit none
    call list_cell()
  end subroutine examine_cell

  subroutine examine_potentials()
    implicit none
    call list_interactions()
  end subroutine examine_potentials


  ! Creates potentials for describing the atomic interactions
  subroutine allocate_potentials(n_pots)
    implicit none
    integer, intent(in) :: n_pots

    call core_allocate_potentials(n_pots)

  end subroutine allocate_potentials


  subroutine add_potential(n_targets,n_params,pot_name,parameters,cutoff,smooth_cut,&
       elements,tags,indices,orig_elements,orig_tags,orig_indices)
    implicit none
    integer, intent(in) :: n_targets, n_params
    character(len=*), intent(in) :: pot_name
    double precision, intent(in) :: parameters(n_params)
    double precision, intent(in) :: cutoff, smooth_cut
    integer, intent(in) :: elements(2,n_targets) ! label_length
    integer, intent(in) :: tags(n_targets), indices(n_targets)
    integer, intent(in) :: orig_elements(2,n_targets) ! label_length
    integer, intent(in) :: orig_tags(n_targets), orig_indices(n_targets)
    character(len=2) :: elements_str(n_targets), orig_elements_str(n_targets) ! label_length
    integer :: i

    ! translates the integer-formatted labels to characters
    do i = 1, n_targets
       call int2str(2,elements(:,i),elements_str(i))
       call int2str(2,orig_elements(:,i),orig_elements_str(i))
    end do
    ! indices +1 because fortran starts indexing from 1
    call core_add_potential(n_targets,n_params,pot_name,parameters,cutoff,smooth_cut,&
         elements_str,tags,indices+1,orig_elements_str,orig_tags,orig_indices+1)

  end subroutine add_potential


  ! Creates neighbor lists
  subroutine create_neighbor_list(n_nbs,atom_index,neighbors,offsets)
    implicit none
    integer, intent(in) :: n_nbs
    integer, intent(in) :: neighbors(n_nbs), offsets(3,n_nbs)
    integer, intent(in) :: atom_index

    ! add +1 to neighbors because python indexing begins from 0 and fortran from 1
    call core_create_neighbor_list(n_nbs,atom_index,neighbors+1,offsets)

  end subroutine create_neighbor_list


  subroutine create_potential_list()
    implicit none

    call core_assign_potential_indices()

  end subroutine create_potential_list


  subroutine calculate_coordinations(n_atoms,cutoffs,coordinations)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(in) :: cutoffs(2)
    double precision, intent(out) :: coordinations(n_atoms)

    call core_calculate_coordinations(n_atoms,cutoffs,coordinations)

  end subroutine calculate_coordinations



  subroutine calculate_bond_orders(n_atoms,group_index,bond_orders)
    implicit none
    integer, intent(in) :: n_atoms, group_index
    double precision, intent(out) :: bond_orders(n_atoms)

    call core_calculate_bond_orders(n_atoms,group_index,bond_orders)

  end subroutine calculate_bond_orders


  ! Calculates the total energy of the system
  subroutine calculate_energy(energy)
    implicit none
    double precision, intent(out) :: energy

    call core_calculate_energy(energy)

  end subroutine calculate_energy


  ! Calculates forces acting on the particles
  subroutine calculate_forces(n_atoms,forces)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(out) :: forces(3,n_atoms)

    call core_calculate_forces(n_atoms,forces)

  end subroutine calculate_forces


  ! Calculates the stress tensor of the cell
  subroutine calculate_stress()
    implicit none

  end subroutine calculate_stress


  ! Tells the number of differently named potentials the core knows
  subroutine number_of_potentials(n_pots)
    implicit none
    integer, intent(out) :: n_pots

    call get_number_of_potentials(n_pots)

  end subroutine number_of_potentials

  ! Tells the number of differently named bond order factors the core knows
  subroutine number_of_bond_order_factors(n_bonds)
    implicit none
    integer, intent(out) :: n_bonds

    call get_number_of_bond_order_factors(n_bonds)

  end subroutine number_of_bond_order_factors

  ! Tells whether a given keyword defines a potential or not
  subroutine is_potential(string,is_ok)
    implicit none
    character(len=*), intent(in) :: string
    logical, intent(out) :: is_ok

    call is_valid_potential(string,is_ok)

  end subroutine is_potential

  ! Tells whether a given keyword defines a bond order factor or not
  subroutine is_bond_order_factor(string,is_ok)
    implicit none
    character(len=*), intent(in) :: string
    logical, intent(out) :: is_ok

    call is_valid_bond_order_factor(string,is_ok)

  end subroutine is_bond_order_factor

  ! Lists all the keywords which define a potential
  subroutine list_valid_potentials(n_pots,potentials)
    implicit none
    integer, intent(in) :: n_pots
    integer, intent(out) :: potentials(11,n_pots) ! pot_name_length
    character(len=11), dimension(n_pots) :: pots ! pot_name_length
    integer :: i

    call list_potentials(n_pots,pots)
    do i = 1, n_pots
       call str2int(pot_name_length,pots(i),potentials(1:pot_name_length,i))
    end do

  end subroutine list_valid_potentials

  ! Lists all the keywords which define a bond order factor
  subroutine list_valid_bond_order_factors(n_bonds,bond_factors)
    implicit none
    integer, intent(in) :: n_bonds
    integer, intent(out) :: bond_factors(11,n_bonds) ! pot_name_length
    character(len=11), dimension(n_bonds) :: bonds ! pot_name_length
    integer :: i

    call list_bond_order_factors(n_bonds,bonds)
    do i = 1, n_bonds
       call str2int(pot_name_length,bonds(i),bond_factors(1:pot_name_length,i))
    end do

  end subroutine list_valid_bond_order_factors

  ! Tells how many targets a potential has, i.e., is it a many-body potential
  subroutine number_of_targets_of_potential(pot_name, n_target)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: n_target

    call get_number_of_targets_of_potential(pot_name,n_target)

  end subroutine number_of_targets_of_potential

  ! Tells how many targets a bond order factor has, i.e., is it many-body
  subroutine number_of_targets_of_bond_order_factor(bond_name, n_target)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(out) :: n_target

    call get_number_of_targets_of_bond_order_factor(bond_name,n_target)

  end subroutine number_of_targets_of_bond_order_factor

  ! Tells how many numeric parameters a potential incorporates
  subroutine number_of_parameters_of_potential(pot_name, n_params)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: n_params

    call get_number_of_parameters_of_potential(pot_name,n_params)

  end subroutine number_of_parameters_of_potential


  ! Tells how many numeric parameters a bond order factor incorporates
  subroutine number_of_parameters_of_bond_order_factor(bond_name, n_targets, n_params)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(in) :: n_targets
    integer, intent(out) :: n_params

    call get_number_of_parameters_of_bond_order_factor(bond_name,n_targets,n_params)

  end subroutine number_of_parameters_of_bond_order_factor


  subroutine names_of_parameters_of_potential(pot_name,param_names)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: param_names(10,4) ! param_name_length, n_max_param
    character(len=10), pointer :: param_name_str(:) ! param_name_length
    integer :: i
    
    call get_names_of_parameters_of_potential(pot_name,param_name_str)
    param_names = 0
    do i = 1, size(param_name_str(:))
       call str2int(param_name_length,param_name_str(i),param_names(1:param_name_length,i))
    end do

  end subroutine names_of_parameters_of_potential


  subroutine names_of_parameters_of_bond_order_factor(bond_name,n_targets,param_names)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(in) :: n_targets
    integer, intent(out) :: param_names(10,4) ! param_name_length, n_max_param
    character(len=10), pointer :: param_name_str(:) ! param_name_length
    integer :: i
    
    call get_names_of_parameters_of_bond_order_factor(bond_name, n_targets, param_name_str)
    param_names = 0
    do i = 1, size(param_name_str(:))
       call str2int(param_name_length,param_name_str(i),param_names(1:param_name_length,i))
    end do

  end subroutine names_of_parameters_of_bond_order_factor



  subroutine descriptions_of_parameters_of_potential(pot_name,param_notes)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: param_notes(100,3) ! param_note_length, n_max_param
    character(len=100), pointer :: param_note_str(:) ! param_note_length
    integer :: i
    
    call get_descriptions_of_parameters_of_potential(pot_name,param_note_str)
    param_notes = 0
    do i = 1, size(param_note_str(:))
       call str2int(param_note_length,param_note_str(i),param_notes(1:param_note_length,i))
    end do

  end subroutine descriptions_of_parameters_of_potential



  subroutine descriptions_of_parameters_of_bond_order_factor(bond_name,n_targets,param_notes)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(in) :: n_targets
    integer, intent(out) :: param_notes(100,4) ! param_note_length, n_max_param
    character(len=100), pointer :: param_note_str(:) ! param_note_length
    integer :: i
    
    call get_descriptions_of_parameters_of_bond_order_factor(bond_name,n_targets,param_note_str)
    param_notes = 0
    do i = 1, size(param_note_str(:))
       call str2int(param_note_length,param_note_str(i),param_notes(1:param_note_length,i))
    end do

  end subroutine descriptions_of_parameters_of_bond_order_factor



  subroutine description_of_potential(pot_name,description)
    implicit none
    character(len=*), intent(in) :: pot_name
    character(len=500), intent(out) :: description ! pot_note_length
    
    call get_description_of_potential(pot_name,description)

  end subroutine description_of_potential



  subroutine description_of_bond_order_factor(bond_name,description)
    implicit none
    character(len=*), intent(in) :: bond_name
    character(len=500), intent(out) :: description ! pot_note_length
    
    call get_description_of_bond_order_factor(bond_name,description)

  end subroutine description_of_bond_order_factor


  subroutine release()
    implicit none
    
    call core_release_all_memory()
    owner_id = -1

  end subroutine release


end module pysic_interface
