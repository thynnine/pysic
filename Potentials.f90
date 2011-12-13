module potentials
  use utility
  use geometry
  implicit none

  ! NB: When adding new potentials remember to update these parameters if necessary.
  ! At least n_potential_types must always be updated.
  ! Also note that in PyInterface, these constants are not used because f2py doesn't read them from here.
  ! search for their names to find them (in comments)
  integer, parameter :: pot_name_length = 11, &
       param_name_length = 10, &
       n_potential_types = 4, &
       n_max_params = 3, &
       pot_note_length = 500, &
       param_note_length = 100

  integer, parameter :: pair_lj_index = 1, &
       pair_spring_index = 2, &
       mono_const_index = 3, &
       tri_bend_index = 4

  character(len=label_length), parameter :: no_name = "xx"
  logical :: descriptors_created = .false.

  ! Contains a description of a type of a potential.
  ! The type contains the name and description of the potential
  ! and the parameters it contains.
  type potential_descriptor
     character(len=pot_name_length) :: name
     character(len=param_name_length), pointer :: parameter_names(:)
     integer :: n_parameters, n_targets, type_index
     character(len=pot_note_length) :: description
     character(len=param_note_length), pointer :: parameter_notes(:)
  end type potential_descriptor

  type(potential_descriptor), pointer :: potential_descriptors(:)

  ! Defines a particular potential.
  ! The potential should correspond to the description of some
  ! built-in type and hold actual numeric values for parameters. 
  ! In addition a real potential must have information on the 
  ! particles it acts on and the range it operates in.
  type potential
     integer :: type_index
     double precision, pointer :: parameters(:), derived_parameters(:)
     double precision :: cutoff, soft_cutoff
     character(len=2), pointer :: apply_elements(:) ! label_length
     integer, pointer :: apply_tags(:), apply_indices(:)
     character(len=2), pointer :: original_elements(:) ! label_length
     integer, pointer :: original_tags(:), original_indices(:)
     logical :: filter_elements, filter_tags, filter_indices, smoothened
  end type potential
  
contains

  subroutine create_potential(n_targets,n_params,pot_name,parameters,cutoff,soft_cutoff,&
       elements,tags,indices,orig_elements,orig_tags,orig_indices,new_potential)
    implicit none
    integer, intent(in) :: n_targets, n_params
    character(len=*), intent(in) :: pot_name
    double precision, intent(in) :: parameters(n_params)
    double precision, intent(in) :: cutoff, soft_cutoff
    character(len=2), intent(in) :: elements(n_targets) ! label_length
    integer, intent(in) :: tags(n_targets), indices(n_targets)
    character(len=2), intent(in) :: orig_elements(n_targets) ! label_length
    integer, intent(in) :: orig_tags(n_targets), orig_indices(n_targets)
    type(potential), intent(out) :: new_potential
    type(potential_descriptor) :: descriptor
    logical :: smoothen
    
    call get_descriptor(pot_name, descriptor)

    new_potential%type_index = descriptor%type_index
    nullify(new_potential%parameters)
    allocate(new_potential%parameters(n_params))
    new_potential%parameters = parameters
    new_potential%cutoff = cutoff
    if(soft_cutoff > 0.d0 .and. soft_cutoff < cutoff)then
       new_potential%soft_cutoff = soft_cutoff
       new_potential%smoothened = .true.
    else
       new_potential%soft_cutoff = cutoff
       new_potential%smoothened = .false.
    end if
    nullify(new_potential%apply_elements)
    allocate(new_potential%apply_elements(n_targets))
    new_potential%apply_elements = elements
    nullify(new_potential%apply_tags)
    nullify(new_potential%apply_indices)
    allocate(new_potential%apply_tags(n_targets))
    allocate(new_potential%apply_indices(n_targets))
    new_potential%apply_tags = tags
    new_potential%apply_indices = indices

    nullify(new_potential%original_elements)
    allocate(new_potential%original_elements(n_targets))
    new_potential%original_elements = orig_elements

    nullify(new_potential%original_tags)
    nullify(new_potential%original_indices)
    allocate(new_potential%original_tags(n_targets))
    allocate(new_potential%original_indices(n_targets))
    new_potential%original_tags = orig_tags
    new_potential%original_indices = orig_indices

    ! calculate derived parameters if necessary
    select case (new_potential%type_index)
    case(tri_bend_index) ! bond bending
       nullify(new_potential%derived_parameters)
       allocate(new_potential%derived_parameters(1))
       new_potential%derived_parameters(1) = cos(parameters(2))
    case default
       nullify(new_potential%derived_parameters)
       allocate(new_potential%derived_parameters(1))
       new_potential%derived_parameters(1) = 0.0
    end select
       

    if(n_targets < 1)then
       return
    end if

    if(elements(1) == no_name)then
       new_potential%filter_elements = .false.
    else
       new_potential%filter_elements = .true.       
    end if

    if(tags(1) < 0)then
       new_potential%filter_tags = .false.
    else
       new_potential%filter_tags = .true.       
    end if

    if(indices(1) < 1)then
       new_potential%filter_indices = .false.
    else
       new_potential%filter_indices = .true.       
    end if

  end subroutine create_potential


  subroutine evaluate_forces(n_targets,separations,distances,interaction,force,atoms)
    implicit none
    integer, intent(in) :: n_targets
    double precision, intent(in) :: separations(3,n_targets-1), distances(n_targets-1)
    type(potential), intent(in) :: interaction
    double precision, intent(out) :: force(3,n_targets)
    type(atom), optional, intent(in) :: atoms(n_targets)
    double precision :: r1, r2, r3, r4, r6, ratio
    double precision :: tmp1(3), tmp2(3), tmp3(3), tmp4(3)

    force = 0.d0

    select case (interaction%type_index)
    case (pair_lj_index)! lennard-jones

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          ratio = interaction%parameters(2) / r1
          r6 = ratio*ratio*ratio*ratio*ratio*ratio
          force(1:3,1) = interaction%parameters(1) * ( 6.d0*r6 - 12.d0*r6*r6 ) * separations(1:3,1) / (r1*r1)
          force(1:3,2) = -force(1:3,1)
       end if

    case (pair_spring_index)! spring-potential

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          r2 = (r1 - interaction%parameters(2))
          force(1:3,1) = interaction%parameters(1) * r2 * separations(1:3,1) / r1
          force(1:3,2) = -force(1:3,1)
       end if

    case (mono_const_index)! constant force

       force(1:3,1) = interaction%parameters(1:3)

    case (tri_bend_index)! bond bending

       ! the bond bending is applied to all ordered triplets a1--a2--a3
       ! for which the central atom (a2) is of the correct type
       ! note that the core passes all triplets here for filtering

       if(.not.present(atoms))then
          return
       end if
       if(interaction%filter_elements)then
          if(interaction%original_elements(2) /= atoms(2)%element)then
             return
          end if
       end if
       if(interaction%filter_tags)then
          if(interaction%original_tags(2) /= atoms(2)%tags)then
             return
          end if
       end if
       if(interaction%filter_indices)then
          if(interaction%original_indices(2) /= atoms(2)%index)then
             return
          end if
       end if

       r1 = distances(1)
       r2 = distances(2)

       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          if(r2 < interaction%cutoff .and. r2 > 0.d0)then

             tmp1 = separations(1:3,1)
             tmp2 = separations(1:3,2)
             r3 = tmp1 .o. tmp2
             tmp3 = tmp2 - r3 / ( r1 * r1 ) * tmp1 ! = r_23 - (r_21 . r_23)/|r_21|^2 * r_21
             tmp4 = tmp1 - r3 / ( r2 * r2 ) * tmp2 ! = r_21 - (r_21 . r_23)/|r_23|^2 * r_23
             ratio = 1.d0 / ( r1 * r2 )

             ! cos theta = (r_21 . r_23) / ( |r_21| |r_23| ) = r3*ratio
             ! k ( cos theta - cos theta_0) =
             r6 = interaction%parameters(1) * (r3 * ratio - interaction%derived_parameters(1))
             
             force(1:3,1) = -tmp3*ratio * r6 
             force(1:3,2) = (tmp3+tmp4)*ratio * r6
             force(1:3,3) = -tmp4*ratio * r6

          end if
       end if

    end select

  end subroutine evaluate_forces


  subroutine evaluate_energy(n_targets,separations,distances,interaction,energy,atoms)
    implicit none
    integer, intent(in) :: n_targets
    double precision, intent(in) :: separations(3,n_targets-1), distances(n_targets-1)
    type(potential), intent(in) :: interaction
    double precision, intent(out) :: energy
    type(atom), optional, intent(in) :: atoms(n_targets)
    double precision :: r1, r2, r3, r6, ratio

    energy = 0.d0

    select case (interaction%type_index)
    case (pair_lj_index) ! lennard-jones

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then          
          ratio = interaction%parameters(2) / r1
          r6 = ratio*ratio*ratio*ratio*ratio*ratio
          energy = interaction%parameters(1) * (r6*r6 - r6)
       end if

    case (pair_spring_index)! spring-potential

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          r2 = (r1 - interaction%parameters(2))
          r6 = (interaction%cutoff - interaction%parameters(2))
          energy = interaction%parameters(1) * 0.5d0 * (r2*r2 - r6*r6)
       end if

    case (mono_const_index)! constant force

       energy = - interaction%parameters(1:3) .o. atoms(1)%position(1:3)

    case (tri_bend_index)! bond-bending
       
       ! the bond bending is applied to all ordered triplets a1--a2--a3
       ! for which the central atom (a2) is of the correct type
       ! note that the core passes all triplets here for filtering
       if(.not.present(atoms))then
          return
       end if
       if(interaction%filter_elements)then
          if(interaction%original_elements(2) /= atoms(2)%element)then
             return
          end if
       end if
       if(interaction%filter_tags)then
          if(interaction%original_tags(2) /= atoms(2)%tags)then
             return
          end if
       end if
       if(interaction%filter_indices)then
          if(interaction%original_indices(2) /= atoms(2)%index)then
             return
          end if
       end if

       r1 = distances(1)
       r2 = distances(2)
       
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          if(r2 < interaction%cutoff .and. r2 > 0.d0)then

             r3 = separations(1:3,1) .o. separations(1:3,2)
             ! cos theta = (r_21 . r_23) / ( |r_21| |r_23| )
             r6 = (r3  / ( r1 * r2 )  - interaction%derived_parameters(1))
             ! k/2 ( cos theta - cos theta_0 )^2
             energy = 0.5d0 * interaction%parameters(1) * r6*r6

          end if
       end if

    end select

  end subroutine evaluate_energy



  subroutine smoothening_factor(r,hard_cut,soft_cut,factor)
    implicit none
    double precision, intent(in) :: r, hard_cut, soft_cut
    double precision, intent(out) :: factor

    if(r > hard_cut)then
       factor = 0.d0
       return
    else if(r < soft_cut)then
       factor = 1.d0
       return
    else
       factor = 0.5d0*( 1.d0 + cos( pi * (r-soft_cut) / (hard_cut-soft_cut) ) )
       return
    end if

  end subroutine smoothening_factor



  subroutine smoothening_derivative(r,hard_cut,soft_cut,factor)
    implicit none
    double precision, intent(in) :: r, hard_cut, soft_cut
    double precision, intent(out) :: factor

    if(r > hard_cut)then
       factor = 0.d0
       return
    else if(r < soft_cut)then
       factor = 0.d0
       return
    else
       factor = - 0.5d0 * pi / (hard_cut-soft_cut) * sin( pi * (r-soft_cut) / (hard_cut-soft_cut) ) 
       return
    end if

  end subroutine smoothening_derivative



  subroutine smoothening_gradient(unit_vector,r,hard_cut,soft_cut,gradient)
    implicit none
    double precision, intent(in) :: unit_vector(3), r, hard_cut, soft_cut
    double precision, intent(out) :: gradient(3)
    double precision :: factor
    
    call smoothening_derivative(r,hard_cut,soft_cut,factor)
    gradient = unit_vector * factor

  end subroutine smoothening_gradient



  subroutine clear_potential_characterizers()
    implicit none

    if(descriptors_created)then
       deallocate(potential_descriptors)
    else
       nullify(potential_descriptors)
    end if

  end subroutine clear_potential_characterizers



  subroutine initialize_potential_characterizers()
    implicit none
    integer :: index

    call clear_potential_characterizers()
    allocate(potential_descriptors(0:n_potential_types))
    index = 0

    ! Lennard-Jones potential
    index = index+1
    if(pair_lj_index /= index)then
       write(*,*) "Potential indices in the core no not match!"
    end if
    potential_descriptors(index)%type_index = index
    call pad_string('pair_LJ', pot_name_length,potential_descriptors(index)%name)
    potential_descriptors(index)%n_parameters = 2
    potential_descriptors(index)%n_targets = 2
    allocate(potential_descriptors(index)%parameter_names(potential_descriptors(index)%n_parameters))
    allocate(potential_descriptors(index)%parameter_notes(potential_descriptors(index)%n_parameters))
    call pad_string('epsilon', param_name_length,potential_descriptors(index)%parameter_names(1))
    call pad_string('energy scale constant', param_note_length,potential_descriptors(index)%parameter_notes(1))
    call pad_string('sigma', param_name_length,potential_descriptors(index)%parameter_names(2))
    call pad_string('length scale constant', param_note_length,potential_descriptors(index)%parameter_notes(2))
    call pad_string('A standard Lennard-Jones potential: V(r) = epsilon * ( (sigma/r)^12 - (sigma/r)^6 )', &
         pot_note_length,potential_descriptors(index)%description)

    ! spring potential
    index = index+1
    if(pair_spring_index /= index)then
       write(*,*) "Potential indices in the core no not match!"
    end if
    potential_descriptors(index)%type_index = index
    call pad_string('pair_spring', pot_name_length,potential_descriptors(index)%name)
    potential_descriptors(index)%n_parameters = 2
    potential_descriptors(index)%n_targets = 2
    allocate(potential_descriptors(index)%parameter_names(potential_descriptors(index)%n_parameters))
    allocate(potential_descriptors(index)%parameter_notes(potential_descriptors(index)%n_parameters))
    call pad_string('k', param_name_length,potential_descriptors(index)%parameter_names(1))
    call pad_string('spring constant', param_note_length,potential_descriptors(index)%parameter_notes(1))
    call pad_string('R_0', param_name_length,potential_descriptors(index)%parameter_names(2))
    call pad_string('equilibrium separation', param_note_length,potential_descriptors(index)%parameter_notes(2))
    call pad_string('A standard spring potential: V(r) = k/2 (r-R_0)^2', &
         pot_note_length,potential_descriptors(index)%description)

    ! constant force potential
    index = index+1
    if(mono_const_index /= index)then
       write(*,*) "Potential indices in the core no not match!"
    end if
    potential_descriptors(index)%type_index = index
    call pad_string('mono_const', pot_name_length,potential_descriptors(index)%name)
    potential_descriptors(index)%n_parameters = 3
    potential_descriptors(index)%n_targets = 1
    allocate(potential_descriptors(index)%parameter_names(potential_descriptors(index)%n_parameters))
    allocate(potential_descriptors(index)%parameter_notes(potential_descriptors(index)%n_parameters))
    call pad_string('Fx', param_name_length,potential_descriptors(index)%parameter_names(1))
    call pad_string('x-component', param_note_length,potential_descriptors(index)%parameter_notes(1))
    call pad_string('Fy', param_name_length,potential_descriptors(index)%parameter_names(2))
    call pad_string('y-component', param_note_length,potential_descriptors(index)%parameter_notes(2))
    call pad_string('Fz', param_name_length,potential_descriptors(index)%parameter_names(3))
    call pad_string('z-component', param_note_length,potential_descriptors(index)%parameter_notes(3))
    call pad_string('A constant force: F = [Fx, Fy, Fz]', &
         pot_note_length,potential_descriptors(index)%description)
    
    ! bond-bending potential
    index = index+1
    if(tri_bend_index /= index)then
       write(*,*) "Potential indices in the core no not match!"
    end if
    potential_descriptors(index)%type_index = index
    call pad_string('tri_bend', pot_name_length,potential_descriptors(index)%name)
    potential_descriptors(index)%n_parameters = 2
    potential_descriptors(index)%n_targets = 3
    allocate(potential_descriptors(index)%parameter_names(potential_descriptors(index)%n_parameters))
    allocate(potential_descriptors(index)%parameter_notes(potential_descriptors(index)%n_parameters))
    call pad_string('k', param_name_length,potential_descriptors(index)%parameter_names(1))
    call pad_string('bond angle spring constant', param_note_length,potential_descriptors(index)%parameter_notes(1))
    call pad_string('theta_0', param_name_length,potential_descriptors(index)%parameter_names(2))
    call pad_string('equilibrium bond angle', param_note_length,potential_descriptors(index)%parameter_notes(2))
    call pad_string('Bond bending potential: V(theta) = k/2 (cos theta - cos theta_0)^2', &
         pot_note_length,potential_descriptors(index)%description)

    ! null potential, for errorneous inquiries
    index = 0
    potential_descriptors(index)%type_index = index
    call pad_string('null',pot_name_length,potential_descriptors(index)%name)
    potential_descriptors(index)%n_parameters = 1
    potential_descriptors(index)%n_targets = 0
    allocate(potential_descriptors(index)%parameter_names(potential_descriptors(index)%n_parameters))
    allocate(potential_descriptors(index)%parameter_notes(potential_descriptors(index)%n_parameters))
    call pad_string('null',param_name_length,potential_descriptors(index)%parameter_names(1))
    call pad_string('a dummy parameter',param_note_length,potential_descriptors(index)%parameter_notes(1))
    call pad_string('No such potential is available.'&
         ,pot_note_length,potential_descriptors(index)%description)

    descriptors_created = .true.

  end subroutine initialize_potential_characterizers
  


  subroutine get_number_of_potentials(n_pots)
    implicit none
    integer, intent(out) :: n_pots

    n_pots = n_potential_types

  end subroutine get_number_of_potentials



  subroutine list_potentials(pots)
    implicit none
    character(len=pot_name_length), dimension(n_potential_types), intent(out) :: pots
    integer :: i

    do i = 1, n_potential_types
       pots(i) = potential_descriptors(i)%name
    end do

  end subroutine list_potentials



  subroutine get_index_of_potential(pot_name, index)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: index
    integer :: i

    index = 0
    do i = 1, n_potential_types
       if(potential_descriptors(i)%name == pot_name)then
          index = i
          return
       end if
    end do

  end subroutine get_index_of_potential



  subroutine get_descriptor(pot_name,descriptor)
    implicit none
    character(len=*), intent(in) :: pot_name
    type(potential_descriptor), intent(out) :: descriptor
    integer :: i

    do i = 1, n_potential_types
       if(potential_descriptors(i)%name == pot_name)then
          descriptor = potential_descriptors(i)
          return
       end if
    end do
    descriptor = potential_descriptors(0)

  end subroutine get_descriptor



  subroutine get_names_of_parameters_of_potential(pot_name, param_names)
    implicit none
    character(len=*), intent(in) :: pot_name
    character(len=param_name_length), pointer :: param_names(:)
    type(potential_descriptor) :: descriptor
    integer :: i

    call get_descriptor(pot_name,descriptor)
    nullify(param_names)
    allocate(param_names(descriptor%n_parameters))
    do i = 1, descriptor%n_parameters
       param_names(i) = descriptor%parameter_names(i)
    end do

  end subroutine get_names_of_parameters_of_potential



  subroutine get_descriptions_of_parameters_of_potential(pot_name, param_notes)
    implicit none
    character(len=*), intent(in) :: pot_name
    character(len=param_note_length), pointer :: param_notes(:)
    type(potential_descriptor) :: descriptor
    integer :: i

    call get_descriptor(pot_name,descriptor)
    nullify(param_notes)
    allocate(param_notes(descriptor%n_parameters))
    do i = 1, descriptor%n_parameters
       param_notes(i) = descriptor%parameter_notes(i)
    end do

  end subroutine get_descriptions_of_parameters_of_potential



  subroutine get_number_of_parameters_of_potential(pot_name,n_params)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: n_params
    type(potential_descriptor) :: descriptor

    call get_descriptor(pot_name,descriptor)
    n_params = descriptor%n_parameters

  end subroutine get_number_of_parameters_of_potential



  subroutine get_number_of_targets_of_potential(pot_name,n_target)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: n_target
    type(potential_descriptor) :: descriptor

    call get_descriptor(pot_name,descriptor)
    n_target = descriptor%n_targets

  end subroutine get_number_of_targets_of_potential



  subroutine get_number_of_targets_of_potential_index(pot_index,n_target)
    implicit none
    integer, intent(in) :: pot_index
    integer, intent(out) :: n_target

    n_target = potential_descriptors(pot_index)%n_targets

  end subroutine get_number_of_targets_of_potential_index


  subroutine is_valid_potential(string,is_valid)
    implicit none
    character(len=*), intent(in) :: string
    logical, intent(out) :: is_valid
    integer :: i

    is_valid = .false.
    do i = 1, n_potential_types
       if(potential_descriptors(i)%name == string)then
          is_valid = .true.
          return
       end if
    end do

  end subroutine is_valid_potential



  subroutine get_index_of_parameter_of_potential(pot_name,param_name,index)
    implicit none
    character(len=*), intent(in) :: pot_name
    character(len=*), intent(in) :: param_name
    integer, intent(out) :: index
    type(potential_descriptor) :: descriptor
    integer :: i

    call get_descriptor(pot_name,descriptor)
    
    index = 0
    do i = 1, descriptor%n_parameters
       if (descriptor%parameter_names(i) == param_name) then
          index = i
          return
       end if
    end do

  end subroutine get_index_of_parameter_of_potential



  subroutine get_description_of_potential(pot_name,description)
    implicit none
    character(len=*), intent(in) :: pot_name
    character(len=pot_note_length), intent(out) :: description
    type(potential_descriptor) :: descriptor

    call get_descriptor(pot_name,descriptor)
    description = descriptor%description

  end subroutine get_description_of_potential


  subroutine potential_affects_atom(interaction,atom_in,affects,position)
    implicit none
    type(potential), intent(in) :: interaction
    type(atom), intent(in) :: atom_in
    logical, intent(out) :: affects
    integer, optional, intent(in) :: position
    integer :: i

    affects = .false.

    if(present(position) .and. position <= size(interaction%apply_elements) .and. position > 0)then

       if(interaction%filter_elements)then
          if(interaction%apply_elements(position) == atom_in%element)then
             affects = .true.
             return
          end if
       end if

       if(interaction%filter_tags)then
          if(interaction%apply_tags(position) == atom_in%tags)then
             affects = .true.
             return
          end if
       end if

       if(interaction%filter_indices)then
          if(interaction%apply_indices(position) == atom_in%index)then
             affects = .true.
             return
          end if
       end if

    else

       if(interaction%filter_elements)then
          do i = 1, size(interaction%apply_elements)
             if(interaction%apply_elements(i) == atom_in%element)then
                affects = .true.
                return
             end if
          end do
       end if

       if(interaction%filter_tags)then
          do i = 1, size(interaction%apply_tags)
             if(interaction%apply_tags(i) == atom_in%tags)then
                affects = .true.
                return
             end if
          end do
       end if

       if(interaction%filter_indices)then
          do i = 1, size(interaction%apply_indices)
             if(interaction%apply_indices(i) == atom_in%index)then
                affects = .true.
                return
             end if
          end do
       end if

    end if

  end subroutine potential_affects_atom


end module potentials
