!
! Pysic_interface is an interface module between the Python and Fortran
! sides of pysic. When pysic is compiled, only this file is interfaced to
! Python via `f2py`_ while the rest of the Fortran source files are 
! directly compiled to .mod Fortran modules. The main reason for this
! is that F90 derived types are used extensively in the core and these
! are not yet (2011) supported by `f2py`_ (although the support for
! derived types is planned in third generation f2py). Because of this,
! most data is passed from :mod:`pysic` to 
! :ref:`pysic_interface` as NumPy arrays, and conversions from objects
! is needed on both sides. This is cumbersome and adds overhead, but
! it is not an efficiency issue since most of the information is only
! passed from Python to Fortran once and saved. Even during a molecular
! dynamics simulation only the forces, coordinates and momenta
! of atoms need to be communicated through the interface, which is 
! naturally and efficiently handled using just numeric arrays anyway.
!
! Another limitation in current `f2py`_ is handling of string arrays.
! To overcome this, string arrays are converted to integer arrays
! and back using simple mapping functions in :mod:`~pysic.pysic_utility`
! and :ref:`utility`.
!
! Due to the current limitations of `f2py`_, no derived types can appear
! in the module. This severly limits what the module can do, and therefore
! the module has been by design made to be very light in terms of 
! functionality: No data is stored in the module and almost all routines
! simply redirect the call to a submodule, most often :ref:`pysic_core`.
! In the descriptions of the routines in this documentation, 
! links are provided to the submodule routines that the call is directed
! to, if the routine is just a redirect of the call.
!
! .. _f2py: http://www.scipy.org/F2py
!
module pysic_interface
  use pysic_core
  use geometry
  use utility
  use potentials
  use mpi
  use mt95
  implicit none

contains

  ! Initialize Mersenne Twister random number generator.
  !
  ! A seed number has to be given. In case we run in MPI
  ! mode, the master cpu will broadcast its seed to all other
  ! cpus to ensure that the random number sequences match
  ! in all the cpus.
  !
  ! *seed a seed for the random number generator
  subroutine start_rng(seed)
    implicit none
    integer, intent(in) :: seed
    integer :: synced_seed
    
    synced_seed = seed
    ! if we have many cpus it's important they all share the random seed
    call mpi_master_bcast_int(synced_seed) ! in MPI.f90
    call genrand_init(synced_seed) ! in Mersenne.f90

  end subroutine start_rng

  !
  ! ToDo: routines for drawing random numbers from MT
  !

  ! Initializes MPI for parallel calculations.
  ! 
  ! Calls :func:`mpi_initialize`
  subroutine start_mpi()
    implicit none

    call mpi_initialize() ! in MPI.f90

  end subroutine start_mpi


  ! Finishes MPI for parallel calculations.
  ! 
  ! Calls :func:`mpi_finish`
  subroutine finish_mpi()
    implicit none

    call mpi_finish() ! in MPI.f90

  end subroutine finish_mpi


  ! Syncs MPI.
  ! This just calls mpi_barrier, so it makes all cpus
  ! wait until everyone is at this particular point in
  ! execution.
  ! 
  ! Calls :func:`mpi_sync`
  subroutine sync_mpi()
    implicit none

    call mpi_sync() ! in MPI.f90

  end subroutine sync_mpi

  ! Distributes atoms among the processors.
  ! In the MPI scheme, atoms are distributed among
  ! the cpus for force and energy calculations.
  ! This routine initializes the arrays that
  ! tell each cpu which atoms it has to calculate
  ! interactions for. It can be called before
  ! the atoms are created in the core but one has to
  ! make sure the number of atoms specified in the last call
  ! matches the number of atoms in the core when a
  ! calculation is invoked.
  ! 
  ! Calls :func:`mpi_distribute`
  !
  ! *n_atoms number of atoms
  subroutine distribute_mpi(n_atoms)
    implicit none
    integer, intent(in) :: n_atoms

    call mpi_distribute(n_atoms) ! in MPI.f90

  end subroutine distribute_mpi


  ! Returns the MPI cpu id number, which is an
  ! integer between 0 and :math:`n_\mathrm{cpus}-1`,
  ! where :math:`n_\mathrm{cpus}` is the total
  ! number of cpus.
  !
  ! *id cpu id number in MPI - 0 in serial mode
  subroutine get_cpu_id(id)
    implicit none
    integer, intent(out) :: id

    id = cpu_id ! in MPI.f90

  end subroutine get_cpu_id


  ! Returns the MPI cpu count
  !
  ! *ncpu the total number of cpus available
  subroutine get_number_of_cpus(ncpu)
    implicit none
    integer, intent(out) :: ncpu

    ncpu = n_cpus ! in MPI.f90

  end subroutine get_number_of_cpus


  ! Returns a logical array containing true for every
  ! atom that is allocated to this cpu, and false
  ! for all other atoms.
  !
  ! *n_atoms number of atoms
  ! *cpu_atoms array of logical values showing which atoms are marked to be handled by this cpu
  subroutine get_mpi_list_of_atoms(n_atoms,cpu_atoms)
    implicit none
    integer, intent(in) :: n_atoms
    logical, intent(out) :: cpu_atoms(n_atoms)
    
    cpu_atoms(1:n_atoms) = is_my_atom(1:n_atoms) ! in MPI.f90

  end subroutine get_mpi_list_of_atoms


  ! Initializes the potentials.
  ! A routine is called to generate descriptors for
  ! potentials. These descriptors are needed by the
  ! python interface in order to directly inquire
  ! the core on the types of potentials available.
  ! 
  ! Calls :func:`initialize_potential_characterizers`
  subroutine start_potentials()
    implicit none

    call initialize_potential_characterizers() ! in Potentials.f90

  end subroutine start_potentials


  ! Initializes the bond order factors.
  ! A routine is called to generate descriptors for
  ! potentials. These descriptors are needed by the
  ! python interface in order to directly inquire
  ! the core on the types of factors available.
  ! 
  ! Calls :func:`initialize_bond_order_factor_characterizers`
  subroutine start_bond_order_factors()
    implicit none

    call initialize_bond_order_factor_characterizers() ! in Potentials.f90

  end subroutine start_bond_order_factors


  ! Creates a supercell for containing the calculation geometry
  ! Also the inverse cell matrix (reciprocal cell) must be given,
  ! although it is not checked that the given inverse actually
  ! is the true inverse.
  ! 
  ! Calls :func:`core_create_cell`
  !
  ! *vectors A 3x3 matrix containing the vectors spanning the supercell. The first index runs over xyz and the second index runs over the three vectors.
  ! *inverse A 3x3 matrix containing the inverse matrix of the one given in vectors, i.e. :math:`A*B = I` for the two matrices. Since the latter represents a cell of non-zero volume, this inverse must exist. It is not tested that the given matrix actually is the inverse, the user must make sure it is.
  ! *periodicity A 3-element vector containing logical tags specifying if the system is periodic in the directions of the three vectors spanning the supercell.
  subroutine create_cell(vectors,inverse,periodicity)
    implicit none
    double precision, intent(in) :: vectors(3,3), inverse(3,3)
    logical, intent(in) :: periodicity(3)

    call core_create_cell(vectors,inverse,periodicity) ! in Core.f90

  end subroutine create_cell

  ! Returns the vectors defining the simulation supercell.
  ! 
  ! Calls :func:`core_get_cell_vectors`
  !
  ! *vectors A 3x3 matrix containing the vectors spanning the supercell. The first index runs over xyz and the second index runs over the three vectors.
  subroutine get_cell_vectors(vectors)
    implicit none
    double precision, intent(out) :: vectors(3,3)

    call core_get_cell_vectors(vectors) ! in Core.f90
    
  end subroutine get_cell_vectors


  ! Creates atomic particles.
  ! Atoms are handled as custom fortran types in the core. Currently
  ! f2py does not support direct creation of types from python, so instead
  ! all the necessary data is passed from python as arrays and reassembled
  ! as types in fortran. This is not much of an added overhead - the
  ! memory allocation itself already makes this a routine one does not
  ! wish to call repeatedly. Instead, one should call the routines
  ! for updating atoms whenever the actual atoms do not change
  ! (e.g., between MD timesteps).
  ! 
  ! Calls :func:`core_generate_atoms`
  !
  ! *n_atoms number of atoms
  ! *masses masses of atoms
  ! *charges electric charges of atoms
  ! *positions coordinates of atoms
  ! *momenta momenta of atoms
  ! *tags numeric tags for the atoms
  ! *elements atomic symbols of the atoms
  subroutine create_atoms(n_atoms,masses,charges,positions,momenta,tags,elements)
    implicit none
    integer, intent(in) :: n_atoms, elements(2,n_atoms), tags(n_atoms)
    double precision, intent(in) :: masses(n_atoms), charges(n_atoms), positions(3,n_atoms), &
         momenta(3,n_atoms)
    character(len=2) :: elements_str(n_atoms)
    integer :: i

    ! translates the integer-formatted labels to characters
    do i = 1, n_atoms
       call int2str(2,elements(:,i),elements_str(i)) ! in Utility.f90
    end do
    ! tells the Fortran core to create the atoms
    call core_generate_atoms(n_atoms,masses,charges,positions,momenta,tags,elements_str) ! in Core.f90

  end subroutine create_atoms


  ! Updates the positions and velocities of existing atoms.
  ! This method does not allocate memory and so the atoms
  ! must already exist in the core.
  ! 
  ! Calls :func:`core_update_atom_coordinates`
  !
  ! *n_atoms number of atoms
  ! *positions new coordinates for the atoms
  ! *momenta new momenta for the atoms
  subroutine update_atom_coordinates(n_atoms,positions,momenta)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(in) :: positions(3,n_atoms), momenta(3,n_atoms)

    call core_update_atom_coordinates(n_atoms,positions,momenta) ! in Core.f90

  end subroutine update_atom_coordinates


  ! Counts the number of atoms in the current core
  ! 
  ! Calls :func:`core_get_number_of_atoms`
  !
  ! *n_atoms number of atoms
  subroutine get_number_of_atoms(n_atoms)
    implicit none
    integer, intent(out) :: n_atoms

    call core_get_number_of_atoms(n_atoms) ! in Core.f90

  end subroutine get_number_of_atoms


  ! Prints some information about the atoms allocated in the core.
  ! This is mainly for debugging, as the python side should always
  ! dictate what is in the core.
  ! 
  ! Calls :func:`list_atoms`
  subroutine examine_atoms()
    implicit none
    call list_atoms() ! in Core.f90
  end subroutine examine_atoms

  ! Prints some information about the supercell allocated in the core.
  ! This is mainly for debugging, as the python side should always
  ! dictate what is in the core.
  ! 
  ! Calls :func:`list_cell`
  subroutine examine_cell()
    implicit none
    call list_cell() ! in Core.f90
  end subroutine examine_cell

  ! Prints some information about the potential allocated in the core.
  ! This is mainly for debugging, as the python side should always
  ! dictate what is in the core.
  ! 
  ! Calls :func:`list_interactions`
  subroutine examine_potentials()
    implicit none
    call list_interactions() ! in Core.f90
  end subroutine examine_potentials

  ! Prints some information about the bond order factors allocated in the core.
  ! This is mainly for debugging, as the python side should always
  ! dictate what is in the core.
  ! 
  ! Calls :func:`list_bonds`
  subroutine examine_bond_order_factors()
    implicit none
    call list_bonds() ! in Core.f90
  end subroutine examine_bond_order_factors


  ! Allocates memory for storing potentials for describing the atomic interactions.
  ! It is more convenient to loop through the potentials and format them in a
  ! suitable way in python than in fortran. Therefore the core is first called
  ! through this routine in order to allocate memory for the potentials.
  ! Then, each potential is created individually.
  ! 
  ! Calls :func:`core_allocate_potentials`
  !
  ! *n_pots number of potentials
  subroutine allocate_potentials(n_pots)
    implicit none
    integer, intent(in) :: n_pots

    call core_allocate_potentials(n_pots) ! in Core.f90

  end subroutine allocate_potentials

  ! Allocates memory for storing bond order parameters for describing the atomic interactions.
  ! Similar to the allocate_potentials routine.
  ! 
  ! Calls :func:`core_allocate_bond_order_factors`
  !
  ! *n_bonds number of bond order factors
  subroutine allocate_bond_order_factors(n_bonds)
    implicit none
    integer, intent(in) :: n_bonds

    call core_allocate_bond_order_factors(n_bonds) ! in Core.f90

  end subroutine allocate_bond_order_factors

  ! Allocates memory for storing bond order factors for describing the atomic interactions.
  ! The difference to allocate_bond_order_factors is that this method allocates
  ! space for arrays used in storing actual calculated bond order factors. The other
  ! routine allocates space for storing the parameters used in the calculations.
  ! 
  ! Calls :func:`core_allocate_bond_order_storage`
  !
  ! *n_atoms number of atoms
  ! *n_groups number of bond order groups
  ! *n_factors number of bond order parameters 
  subroutine allocate_bond_order_storage(n_atoms,n_groups,n_factors)
    implicit none
    integer, intent(in) :: n_atoms, n_groups, n_factors

    call core_allocate_bond_order_storage(n_atoms,n_groups,n_factors) ! in Core.f90

  end subroutine allocate_bond_order_storage

  ! Creates a potential in the core.
  ! The memory must have been allocated first using allocate_potentials.
  ! 
  ! Calls :func:`core_add_potential`
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
  !
  subroutine add_potential(n_targets,n_params,pot_name,parameters,cutoff,smooth_cut,&
       elements,tags,indices,orig_elements,orig_tags,orig_indices,pot_index)
    implicit none
    integer, intent(in) :: n_targets, n_params,pot_index
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
       call int2str(2,elements(:,i),elements_str(i)) ! in Utility.f90
       call int2str(2,orig_elements(:,i),orig_elements_str(i)) ! in Utility.f90
    end do
    ! indices +1 because fortran starts indexing from 1
    call core_add_potential(n_targets,n_params,pot_name,parameters,cutoff,smooth_cut,&
         elements_str,tags,indices+1,orig_elements_str,orig_tags,orig_indices+1,pot_index) ! in Core.f90

  end subroutine add_potential



  ! Creates a bond order factor in the core.
  ! The memory must have been allocated first using allocate_potentials.
  ! 
  ! Calls :func:`core_add_bond_order_factor`
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
  !
  subroutine add_bond_order_factor(n_targets,n_params,n_split,bond_name,parameters,param_split,&
       cutoff,smooth_cut,elements,orig_elements,group_index)
    implicit none
    integer, intent(in) :: n_targets, n_params, n_split
    integer, intent(in) :: param_split(n_split), group_index
    character(len=*), intent(in) :: bond_name
    double precision, intent(in) :: parameters(n_params)
    double precision, intent(in) :: cutoff, smooth_cut
    integer, intent(in) :: elements(2,n_targets) ! label_length
    integer, intent(in) :: orig_elements(2,n_targets) ! label_length
    character(len=2) :: elements_str(n_targets), orig_elements_str(n_targets) ! label_length
    integer :: i

    ! translates the integer-formatted labels to characters
    do i = 1, n_targets
       call int2str(2,elements(:,i),elements_str(i)) ! in Utility.f90
       call int2str(2,orig_elements(:,i),orig_elements_str(i)) ! in Utility.f90
    end do
    call core_add_bond_order_factor(n_targets,n_params,n_split,bond_name,parameters,param_split,&
         cutoff,smooth_cut,elements_str,orig_elements_str,group_index) ! in Core.f90

  end subroutine add_bond_order_factor




  ! Creates neighbor lists for a single atom 
  ! telling it which other atoms are in its
  ! immediate neighborhood.
  ! The neighbor list must be precalculated, this method only
  ! stores them in the core. The list must contain 
  ! an array storing the indices of the neighboring atoms
  ! as well as the supercell offsets. The offsets are integer
  ! triplets showing how many times must the supercell vectors
  ! be added to the position of the neighbor to find the
  ! neighboring image in a periodic system.
  ! Note that if the system is small, one atom can in
  ! principle appear several times in the neighbor list.
  ! 
  ! Calls :func:`core_create_neighbor_list`
  !
  ! *n_nbs number of neighbors
  ! *atom_index index of the atom for which the neighbor list is created
  ! *neighbors An array containing the indices of the neighboring atoms
  ! *offsets An array containing vectors specifying the offsets of the neighbors in periodic systems.
  subroutine create_neighbor_list(n_nbs,atom_index,neighbors,offsets)
    implicit none
    integer, intent(in) :: n_nbs
    integer, intent(in) :: neighbors(n_nbs), offsets(3,n_nbs)
    integer, intent(in) :: atom_index

    ! add +1 to neighbors because python indexing begins from 0 and fortran from 1
    call core_create_neighbor_list(n_nbs,atom_index,neighbors+1,offsets) ! in Core.f90

  end subroutine create_neighbor_list

  ! Creates a list of indices for all atoms showing which potentials
  ! act on them.
  ! The user may define many potentials to sum up the potential energy of the
  ! system. However, if some potentials only act on certain atoms, they will
  ! be redundant for the other atoms. The potential lists are lists
  ! given to each atom containing the potentials which can act on the
  ! atom.
  ! 
  ! Calls :func:`core_assign_potential_indices`
  subroutine create_potential_list()
    implicit none

    call core_assign_potential_indices() ! in Core.f90

  end subroutine create_potential_list

  ! Similarly to the potential lists, also list containing all the
  ! bond order factors that may affect an atom are stored in a list.
  ! 
  ! Calls :func:`core_assign_bond_order_factor_indices`
  subroutine create_bond_order_factor_list()
    implicit none

    call core_assign_bond_order_factor_indices() ! in Core.f90

  end subroutine create_bond_order_factor_list


  ! Returns bond order factors of the given group for all atoms.
  ! The group index is an identifier for the bond order parameters
  ! which are used for calculating one and the same factors.
  ! In practice, the Coordinators in pysic are indexed and this
  ! indexing is copied in the core. Thus the group index specifies
  ! the coordinator / potential.
  ! 
  ! Calls :func:`core_get_bond_order_factors`
  !
  ! *n_atoms number of atoms
  ! *group_index index for the bond order factor group
  ! *bond_orders the calculated bond order factors
  subroutine calculate_bond_order_factors(n_atoms,group_index,bond_orders)
    implicit none
    integer, intent(in) :: n_atoms, group_index
    double precision, intent(out) :: bond_orders(n_atoms)
    double precision :: raw_sums(n_atoms)

    call core_get_bond_order_factors(n_atoms,group_index,bond_orders) ! in Core.f90

  end subroutine calculate_bond_order_factors


  ! Returns bond order factors gradients of the given group.
  ! The gradients of all factors are given with respect to moving the given atom.
  ! The group index is an identifier for the bond order parameters
  ! which are used for calculating one and the same factors.
  ! In practice, the Coordinators in pysic are indexed and this
  ! indexing is copied in the core. Thus the group index specifies
  ! the coordinator / potential.
  !
  ! Calls :func:`core_get_bond_order_sums`
  ! Calls :func:`core_calculate_bond_order_gradients`
  !
  ! *n_atoms number of atoms
  ! *group_index an index denoting the potential to which the factor is connected
  ! *atom_index index of the atom with respect to which the factors are differentiated
  ! *gradients the calculated bond order gradients
  subroutine calculate_bond_order_gradients(n_atoms,group_index,atom_index,gradients)
    implicit none
    integer, intent(in) :: n_atoms, group_index, atom_index
    double precision, intent(out) :: gradients(3,n_atoms)
    double precision :: bond_orders(n_atoms)
 
    call core_get_bond_order_sums(n_atoms,group_index,bond_orders) ! in Core.f90
    call core_calculate_bond_order_gradients(n_atoms,group_index,atom_index,bond_orders,gradients) ! in Core.f90

  end subroutine calculate_bond_order_gradients



  ! Returns bond order factors gradients of the given group.
  ! The gradients of the given factors is given with respect to moving all atoms.
  ! The group index is an identifier for the bond order parameters
  ! which are used for calculating one and the same factors.
  ! In practice, the Coordinators in pysic are indexed and this
  ! indexing is copied in the core. Thus the group index specifies
  ! the coordinator / potential.
  !
  ! Calls :func:`core_get_bond_order_sums`
  ! Calls :func:`core_calculate_bond_order_gradients_of_factor`
  !
  ! *n_atoms number of atoms
  ! *group_index an index denoting the potential to which the factor is connected
  ! *atom_index index of the atom whose factor is differentiated
  ! *gradients the calculated bond order gradients
  subroutine calculate_bond_order_gradients_of_factor(n_atoms,group_index,atom_index,gradients)
    implicit none
    integer, intent(in) :: n_atoms, group_index, atom_index
    double precision, intent(out) :: gradients(3,n_atoms)
    double precision :: bond_orders(n_atoms)
 
    call core_get_bond_order_sums(n_atoms,group_index,bond_orders) ! in Core.f90
    call core_calculate_bond_order_gradients_of_factor(n_atoms,group_index,atom_index,bond_orders,gradients) ! in Core.f90

  end subroutine calculate_bond_order_gradients_of_factor


  ! Returns the total potential energy of the system
  ! 
  ! Calls :func:`core_calculate_energy`
  !
  ! *n_atoms number of atoms
  ! *energy total potential energy
  subroutine calculate_energy(n_atoms,energy)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(out) :: energy

    call core_calculate_energy(n_atoms,energy) ! in Core.f90

  end subroutine calculate_energy


  ! Returns forces acting on the particles
  ! 
  ! Calls :func:`core_calculate_forces`
  ! 
  ! *n_atoms number of atoms
  ! *forces array of forces on all atoms
  subroutine calculate_forces(n_atoms,forces)
    implicit none
    integer, intent(in) :: n_atoms
    double precision, intent(out) :: forces(3,n_atoms)

    call core_calculate_forces(n_atoms,forces) ! in Core.f90

  end subroutine calculate_forces


  ! Calculates the stress tensor of the cell
  !
  ! ToDo: implement this through force calculation and coordinates
  !
  subroutine calculate_stress()
    implicit none

  end subroutine calculate_stress


  ! Below are a group of routines used for inquiring the core
  ! about the potentials and bond order factors it knows.
  ! This framework was built so that all edits in the core
  ! would be automatically reflected in the python interface
  ! without the need to edit the python code.

  ! Tells the number of differently named potentials the core knows
  ! 
  ! Calls :func:`get_number_of_potentials`
  !
  ! *n_pots number of potentials
  subroutine number_of_potentials(n_pots)
    implicit none
    integer, intent(out) :: n_pots

    call get_number_of_potentials(n_pots) ! in Potentials.f90

  end subroutine number_of_potentials

  ! Tells the number of differently named bond order factors the core knows
  ! 
  ! Calls :func:`get_number_of_bond_order_factors`
  !
  ! *n_bonds number of bond order factors
  subroutine number_of_bond_order_factors(n_bonds)
    implicit none
    integer, intent(out) :: n_bonds

    call get_number_of_bond_order_factors(n_bonds) ! in Potentials.f90

  end subroutine number_of_bond_order_factors

  ! Tells whether a given keyword defines a potential or not
  ! 
  ! Calls :func:`is_valid_potential`
  !
  ! *string name of a potential
  ! *is_ok true if string is a name of a potential
  subroutine is_potential(string,is_ok)
    implicit none
    character(len=*), intent(in) :: string
    logical, intent(out) :: is_ok

    call is_valid_potential(string,is_ok) ! in Potentials.f90

  end subroutine is_potential

  ! Tells whether a given keyword defines a bond order factor or not
  ! 
  ! Calls :func:`is_valid_bond_order_factor`
  !
  ! *string name of a bond order factor
  ! *is_ok true if string is a name of a bond order factor
  subroutine is_bond_order_factor(string,is_ok)
    implicit none
    character(len=*), intent(in) :: string
    logical, intent(out) :: is_ok

    call is_valid_bond_order_factor(string,is_ok) ! in Potentials.f90

  end subroutine is_bond_order_factor

  ! Lists all the keywords which define a potential
  ! 
  ! Calls :func:`list_potentials`
  !
  ! *n_pots number of potential types
  ! *potentials names of the potential types
  subroutine list_valid_potentials(n_pots,potentials)
    implicit none
    integer, intent(in) :: n_pots
    integer, intent(out) :: potentials(11,n_pots) ! pot_name_length
    character(len=11), dimension(n_pots) :: pots ! pot_name_length
    integer :: i

    call list_potentials(n_pots,pots) ! in Potentials.f90
    do i = 1, n_pots
       call str2int(pot_name_length,pots(i),potentials(1:pot_name_length,i)) ! in Utility.f90
    end do

  end subroutine list_valid_potentials

  ! Lists all the keywords which define a bond order factor
  ! 
  ! Calls :func:`list_bond_order_factors`
  !
  ! *n_bonds number of bond order factor types
  ! *bond_factors names of the bond order factor types
  subroutine list_valid_bond_order_factors(n_bonds,bond_factors)
    implicit none
    integer, intent(in) :: n_bonds
    integer, intent(out) :: bond_factors(11,n_bonds) ! pot_name_length
    character(len=11), dimension(n_bonds) :: bonds ! pot_name_length
    integer :: i

    call list_bond_order_factors(n_bonds,bonds) ! in Potentials.f90
    do i = 1, n_bonds
       call str2int(pot_name_length,bonds(i),bond_factors(1:pot_name_length,i)) ! in Utility.f90
    end do

  end subroutine list_valid_bond_order_factors

  ! Tells how many targets a potential has, i.e., is it a many-body potential
  ! 
  ! Calls :func:`get_number_of_targets_of_potential`
  !
  ! *pot_name name of the potential
  ! *n_target number of targets
  subroutine number_of_targets_of_potential(pot_name, n_target)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: n_target

    call get_number_of_targets_of_potential(pot_name,n_target) ! in Potentials.f90

  end subroutine number_of_targets_of_potential

  ! Tells how many targets a bond order factor has, i.e., is it many-body
  ! 
  ! Calls :func:`get_number_of_targets_of_bond_order_factor`
  !
  ! *bond_name name of the bond order factor
  ! *n_target number of targets
  subroutine number_of_targets_of_bond_order_factor(bond_name, n_target)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(out) :: n_target

    call get_number_of_targets_of_bond_order_factor(bond_name,n_target) ! in Potentials.f90

  end subroutine number_of_targets_of_bond_order_factor

  ! Tells how many numeric parameters a potential incorporates
  ! 
  ! Calls :func:`get_number_of_parameters_of_potential`
  !
  ! *pot_name name of the potential
  ! *n_params number of parameters
  subroutine number_of_parameters_of_potential(pot_name, n_params)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: n_params

    call get_number_of_parameters_of_potential(pot_name,n_params) ! in Potentials.f90

  end subroutine number_of_parameters_of_potential

  ! Tells how many numeric parameters a bond order factor incorporates
  ! 
  ! Calls :func:`get_number_of_parameters_of_bond_order_factor`
  !
  ! *bond_name name of the bond order factor
  ! *n_targets number of targets
  ! *n_params number of parameters
  subroutine number_of_parameters_of_bond_order_factor(bond_name, n_targets, n_params)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(in) :: n_targets
    integer, intent(out) :: n_params

    call get_number_of_parameters_of_bond_order_factor(bond_name,n_targets,n_params) ! in Potentials.f90

  end subroutine number_of_parameters_of_bond_order_factor

  ! Lists the names of parameters the given potential knows.
  ! Output is an array of integers. This is because f2py doesn't
  ! currently support string arrays. So, the characters are translated to
  ! integers and back in fortran and python.
  ! This adds a bit of overhead, but the routine is only invoked
  ! on user command so it doesn't matter.
  ! 
  ! Calls :func:`get_names_of_parameters_of_potential`
  !
  ! *pot_name name of the potential
  ! *param_names names of the parameters
  subroutine names_of_parameters_of_potential(pot_name,param_names)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: param_names(10,4) ! param_name_length, n_max_param
    character(len=10), pointer :: param_name_str(:) ! param_name_length
    integer :: i
    
    call get_names_of_parameters_of_potential(pot_name,param_name_str) ! in Potentials.f90
    param_names = 0
    do i = 1, size(param_name_str(:))
       call str2int(param_name_length,param_name_str(i),param_names(1:param_name_length,i)) ! in Utility.f90
    end do

  end subroutine names_of_parameters_of_potential


  ! Lists the names of parameters the given bond order factor knows.
  ! Output is an array of integers. This is because f2py doesn't
  ! currently support string arrays. So, the characters are translated to
  ! integers and back in fortran and python.
  ! This adds a bit of overhead, but the routine is only invoked
  ! on user command so it doesn't matter.
  ! 
  ! Calls :func:`get_names_of_parameters_of_bond_order_factor`
  !
  ! *bond_name name of the bond order factor
  ! *n_targets number of targets
  ! *param_names names of the parameters
  subroutine names_of_parameters_of_bond_order_factor(bond_name,n_targets,param_names)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(in) :: n_targets
    integer, intent(out) :: param_names(10,4) ! param_name_length, n_max_param
    character(len=10), pointer :: param_name_str(:) ! param_name_length
    integer :: i
    
    call get_names_of_parameters_of_bond_order_factor(bond_name, n_targets, param_name_str) ! in Potentials.f90
    param_names = 0
    do i = 1, size(param_name_str(:))
       call str2int(param_name_length,param_name_str(i),param_names(1:param_name_length,i)) ! in Utility.f90
    end do

  end subroutine names_of_parameters_of_bond_order_factor


  
  ! Lists descriptions for parameters the given potential.
  ! Output is an array of integers. This is because f2py doesn't
  ! currently support string arrays. So, the characters are translated to
  ! integers and back in fortran and python.
  ! This adds a bit of overhead, but the routine is only invoked
  ! on user command so it doesn't matter.
  ! 
  ! Calls :func:`get_descriptions_of_parameters_of_potential`
  !
  ! *pot_name name of the potential
  ! *param_notes descriptions of the parameters
  subroutine descriptions_of_parameters_of_potential(pot_name,param_notes)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: param_notes(100,3) ! param_note_length, n_max_param
    character(len=100), pointer :: param_note_str(:) ! param_note_length
    integer :: i
    
    call get_descriptions_of_parameters_of_potential(pot_name,param_note_str) ! in Potentials.f90
    param_notes = 0
    do i = 1, size(param_note_str(:))
       call str2int(param_note_length,param_note_str(i),param_notes(1:param_note_length,i)) ! in Utility.f90
    end do

  end subroutine descriptions_of_parameters_of_potential


  ! Lists descriptions for parameters the given bond order factor.
  ! Output is an array of integers. This is because f2py doesn't
  ! currently support string arrays. So, the characters are translated to
  ! integers and back in fortran and python.
  ! This adds a bit of overhead, but the routine is only invoked
  ! on user command so it doesn't matter.
  ! 
  ! Calls :func:`get_descriptions_of_parameters_of_bond_order_factor`
  ! 
  ! *bond_name name of the bond order factor
  ! *n_targets number of targets
  ! *param_notes descriptions of the parameters
  subroutine descriptions_of_parameters_of_bond_order_factor(bond_name,n_targets,param_notes)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(in) :: n_targets
    integer, intent(out) :: param_notes(100,4) ! param_note_length, n_max_param
    character(len=100), pointer :: param_note_str(:) ! param_note_length
    integer :: i
    
    call get_descriptions_of_parameters_of_bond_order_factor(bond_name,n_targets,param_note_str) ! in Potentials.f90
    param_notes = 0
    do i = 1, size(param_note_str(:))
       call str2int(param_note_length,param_note_str(i),param_notes(1:param_note_length,i)) ! in Utility.f90
    end do

  end subroutine descriptions_of_parameters_of_bond_order_factor


  ! Returns a description of the given potential
  ! 
  ! Calls :func:`get_description_of_potential`
  !
  ! *pot_name name of the potential
  ! *description description of the potential
  subroutine description_of_potential(pot_name,description)
    implicit none
    character(len=*), intent(in) :: pot_name
    character(len=500), intent(out) :: description ! pot_note_length
    
    call get_description_of_potential(pot_name,description) ! in Potentials.f90

  end subroutine description_of_potential



  ! Returns a description of the given bond order factor
  ! 
  ! Calls :func:`get_description_of_bond_order_factor`
  !
  ! *bond_name name of the bond order factor
  ! *description description of the bond order actor
  subroutine description_of_bond_order_factor(bond_name,description)
    implicit none
    character(len=*), intent(in) :: bond_name
    character(len=500), intent(out) :: description ! pot_note_length
    
    call get_description_of_bond_order_factor(bond_name,description) ! in Potentials.f90

  end subroutine description_of_bond_order_factor

  ! Deallocates all the arrays in the core
  ! 
  ! Calls :func:`core_release_all_memory`
  subroutine release()
    implicit none
    
    call core_release_all_memory() ! in Core.f90

  end subroutine release


end module pysic_interface
