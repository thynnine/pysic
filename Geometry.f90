module geometry
  use quaternions
  implicit none
  integer, parameter :: label_length = 2

  ! Defines a list of neighbors for a single atom
  type neighbor_list
     integer, pointer :: neighbors(:), pbc_offsets(:,:)
     integer :: max_length, n_neighbors
  end type neighbor_list

  ! Defines an atomic particle
  type atom
     double precision :: mass, charge, position(3), momentum(3)
     integer :: index, tags, n_pots
     logical :: potentials_listed
     character(len=label_length) :: element
     type(neighbor_list) :: neighbor_list
     integer, pointer :: potential_indices(:)
  end type atom

  ! Defines the supercell
  type supercell
     double precision :: vectors(3,3), inverse_cell(3,3), vector_lengths(3)
     logical :: periodic(3)
  end type supercell

contains

  ! Creates the supercell
  subroutine generate_supercell(vectors,inverse,periodicity,cell)
    implicit none
    double precision, intent(in) :: vectors(3,3), inverse(3,3)
    logical, intent(in) :: periodicity(3)
    type(supercell), intent(out) :: cell
    integer :: i

    cell%vectors = vectors
    cell%inverse_cell = inverse
    cell%periodic = periodicity
    do i = 1, 3
       cell%vector_lengths(i) = (.norm.vectors(1:3,i))
    end do

  end subroutine generate_supercell



  ! Creates atoms to construct the system to be simulated
  subroutine generate_atoms(n_atoms,masses,charges,positions,momenta,tags,elements,atoms)
    implicit none
    integer, intent(in) :: n_atoms, tags(n_atoms)
    double precision, intent(in) :: masses(n_atoms), charges(n_atoms), positions(3,n_atoms), &
         momenta(3,n_atoms)
    character(len=label_length), intent(in) :: elements(n_atoms)
    type(atom), pointer :: atoms(:)
    integer :: i

    nullify(atoms)
    allocate(atoms(n_atoms))
    
    do i = 1, n_atoms
       atoms(i)%mass = masses(i)
       atoms(i)%charge = charges(i)
       atoms(i)%position(1:3) = positions(1:3,i)
       atoms(i)%momentum(1:3) = momenta(1:3,i)
       atoms(i)%tags = tags(i)
       atoms(i)%element = elements(i)
       atoms(i)%index = i
       atoms(i)%neighbor_list%max_length = 0
       atoms(i)%neighbor_list%n_neighbors = 0
       nullify(atoms(i)%neighbor_list%neighbors)
       nullify(atoms(i)%neighbor_list%pbc_offsets)
       atoms(i)%n_pots = 0
       atoms(i)%potentials_listed = .false.
       nullify(atoms(i)%potential_indices)
    end do

  end subroutine generate_atoms




  ! Updates the positions and momenta of atoms
  subroutine update_atomic_positions(n_atoms,positions,momenta,atoms)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(in) :: positions(3,n_atoms), momenta(3,n_atoms)
    type(atom), pointer :: atoms(:)
    integer :: i

    if(size(atoms) /= n_atoms)then
       write(*,*) "the number of atoms has changed, you should reinitialize the structure"
    else
       do i = 1, n_atoms
          atoms(i)%position(1:3) = positions(1:3,i)
          atoms(i)%momentum(1:3) = momenta(1:3,i)
       end do
    end if

  end subroutine update_atomic_positions



  subroutine assign_neighbor_list(n_nbs,nbor_list,neighbors,offsets)
    implicit none
    integer, intent(in) :: n_nbs
    integer, intent(in) :: neighbors(n_nbs), offsets(3,n_nbs)
    type(neighbor_list), intent(inout) :: nbor_list

    if(nbor_list%max_length < n_nbs)then
       if(nbor_list%max_length > 0)then
          deallocate(nbor_list%neighbors)
          deallocate(nbor_list%pbc_offsets)
       else
          nullify(nbor_list%neighbors)
          nullify(nbor_list%pbc_offsets)
       end if
       allocate(nbor_list%neighbors(2*n_nbs+10))
       allocate(nbor_list%pbc_offsets(3,2*n_nbs+10))
       nbor_list%max_length = 2*n_nbs+10
    end if
    nbor_list%neighbors = -1
    nbor_list%pbc_offsets = 0
    nbor_list%n_neighbors = n_nbs
    nbor_list%neighbors(1:n_nbs) = neighbors(1:n_nbs)
    nbor_list%pbc_offsets(1:3,1:n_nbs) = offsets(1:3,1:n_nbs)

  end subroutine assign_neighbor_list


  subroutine assign_potential_indices(n_pots,atom_in,indices)
    implicit none    
    integer, intent(in) :: n_pots
    type(atom), intent(inout) :: atom_in
    integer, intent(in) :: indices(n_pots)

    if(atom_in%potentials_listed)then
       if(atom_in%n_pots /= n_pots)then
          deallocate(atom_in%potential_indices)
          allocate(atom_in%potential_indices(size(indices)))
       end if
    else
       nullify(atom_in%potential_indices)
       allocate(atom_in%potential_indices(size(indices)))
    end if
    atom_in%potential_indices = indices
    atom_in%n_pots = size(indices)
    atom_in%potentials_listed = .true.

  end subroutine assign_potential_indices

  ! Calculates the minimum separation vector between two atoms, r2-r1, including possible periodicity
  subroutine separation_vector(r1,r2,offset,cell,separation)
    implicit none
    double precision, intent(in) :: r1(3), r2(3)
    integer, intent(in) :: offset(3)
    type(supercell), intent(in) :: cell
    double precision, intent(out) :: separation(3)
    integer :: i

    separation = r2 - r1
    do i = 1, 3
       separation = separation + offset(i)*cell%vectors(1:3,i)
    end do

  end subroutine separation_vector


  subroutine relative_coordinates(position,cell,relative)
    implicit none
    double precision, intent(in) :: position(3)
    type(supercell), intent(in) :: cell
    double precision, intent(out) :: relative(3)

    relative = matmul(cell%inverse_cell,position)

  end subroutine relative_coordinates


  subroutine absolute_coordinates(relative,cell,position)
    implicit none
    double precision, intent(out) :: position(3)
    type(supercell), intent(in) :: cell
    double precision, intent(in) :: relative(3)

    position = matmul(cell%vectors,relative)

  end subroutine absolute_coordinates


  ! Wraps a general coordinate inside the supercell if the system is periodic
  subroutine wrapped_coordinates(position,cell,wrapped)
    implicit none
    double precision, intent(in) :: position(3)
    type(supercell), intent(in) :: cell
    double precision, intent(out) :: wrapped(3)
    double precision :: relative(3)
    integer :: i
    
    call relative_coordinates(position,cell,relative)
    do i = 1, 3
       relative(i) = relative(i) - floor(relative(i))
    end do
    call absolute_coordinates(relative,cell,wrapped)

  end subroutine wrapped_coordinates



end module geometry
