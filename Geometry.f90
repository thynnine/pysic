!
! A module for handling the geometric structure of the system.
!
module geometry
  use quaternions
  implicit none

  ! *label_length the number of characters available for denoting chemical symbols
  integer, parameter :: label_length = 2

  ! Defines a list of neighbors for a single atom.
  ! The list contains the indices of the neighboring atoms
  ! as well as the periodic boundary condition (PBC) offsets.
  !
  ! The offsets are integer
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
  ! *neighbors indices of the neighboring atoms
  ! *pbc_offsets offsets for periodic boundaries for each neighbor
  ! *max_length The allocated length of the neighbor lists. To avoid deallocating and reallocating memory, extra space is reserved for the neighbors in case the number of neighbors increases during simulation (due to atoms moving).
  ! *n_neighbors the number of neighbors in the lists
  type neighbor_list
     integer, pointer :: neighbors(:), pbc_offsets(:,:)
     integer :: max_length, n_neighbors
  end type neighbor_list

  ! Defines an atomic particle.
  !
  ! *mass mass of th atom
  ! *charge charge of the atom
  ! *position coordinates of the atom
  ! *momentum momentum of the atom
  ! *index index of the atom
  ! *tags integer tag
  ! *n_pots number of potentials that may affect the atom
  ! *n_bonds number of bond order factors that may affect the atom
  ! *potentials_listed logical tag for checking if the potentials affecting the atom have been listed in potential_indices
  ! *bond_order_factors_listed logical tag for checking if the bond order factors affecting the atom have been listed in bond_indices
  ! *element the chemical symbol of the atom
  ! *neighbor_list the list of neighbors for the atom
  ! *potential_indices the indices of the potentials for which this atom is a valid target at first position (see :func:`potential_affects_atom`)
  ! *bond_indices the indices of the bond order factors for which this atom is a valid target at first position (see :func:`bond_order_factor_affects_atom`)
  type atom
     double precision :: mass, charge, position(3), momentum(3)
     integer :: index, tags, n_pots, n_bonds
     logical :: potentials_listed, bond_order_factors_listed
     character(len=label_length) :: element
     type(neighbor_list) :: neighbor_list
     integer, pointer :: potential_indices(:), bond_indices(:)
  end type atom

  ! Supercell containing the simulation.
  !
  ! The supercell is spanned by three vectors :math:`\mathbf{v}_1,\mathbf{v}_2,\mathbf{v}_3` stored as a 
  ! :math:`3 \times 3` matrix in format 
  ! 
  ! .. math::
  !
  !   \mathbf{M} = \left[ 
  !   \begin{array}{ccc}
  !   v_{1,x} & v_{1,y} & v_{1,z} \\
  !   v_{2,x} & v_{2,y} & v_{2,z} \\
  !   v_{3,x} & v_{3,y} & v_{3,z} 
  !   \end{array}
  !   \right].
  !
  ! Also the inverse cell matrix is kept for transformations between the absolute and fractional coordinates.
  !
  ! *vectors vectors spanning the supercell containing the system as a matrix :math:`\mathbf{M}`
  ! *inverse_cell the inverse of the cell matrix :math:`\mathbf{M}^{-1}`
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
       atoms(i)%bond_order_factors_listed = .false.
       nullify(atoms(i)%bond_indices)
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

    if(nbor_list%max_length <= n_nbs)then
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
    if(n_nbs > 0)then
       nbor_list%neighbors(1:n_nbs) = neighbors(1:n_nbs)
       nbor_list%pbc_offsets(1:3,1:n_nbs) = offsets(1:3,1:n_nbs)
    end if

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


  subroutine assign_bond_order_factor_indices(n_bonds,atom_in,indices)
    implicit none    
    integer, intent(in) :: n_bonds
    type(atom), intent(inout) :: atom_in
    integer, intent(in) :: indices(n_bonds)

    if(atom_in%bond_order_factors_listed)then
       if(atom_in%n_bonds /= n_bonds)then
          deallocate(atom_in%bond_indices)
          allocate(atom_in%bond_indices(size(indices)))
       end if
    else
       nullify(atom_in%bond_indices)
       allocate(atom_in%bond_indices(size(indices)))
    end if
    atom_in%bond_indices = indices
    atom_in%n_bonds = size(indices)
    atom_in%bond_order_factors_listed = .true.

  end subroutine assign_bond_order_factor_indices




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
