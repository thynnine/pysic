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
       n_potential_types = 5, &
       n_bond_order_types = 3, &
       n_max_params = 4, &
       n_max_targets = 2, &
       pot_note_length = 500, &
       param_note_length = 100

  integer, parameter :: pair_lj_index = 1, &
       pair_spring_index = 2, &
       mono_const_index = 3, &
       tri_bend_index = 4, &
       mono_none_index = 5

  integer, parameter :: coordination_index = 1, &
       tersoff_index = 2, &
       c_scale_index = 3

  character(len=label_length), parameter :: no_name = "xx"
  logical :: descriptors_created = .false., bond_descriptors_created = .false.

  ! Contains a description of a type of a potential.
  ! The type contains the name and description of the potential
  ! and the parameters it contains.
  type potential_descriptor
     character(len=pot_name_length) :: name
     character(len=param_name_length), pointer :: parameter_names(:)
     integer :: n_parameters, n_targets, type_index
     character(len=pot_note_length) :: description
     character(len=param_note_length), pointer :: parameter_notes(:)
     logical :: bond_order_argument
  end type potential_descriptor

  type bond_order_descriptor
     character(len=pot_name_length) :: name
     character(len=param_name_length), pointer :: parameter_names(:,:)
     integer :: n_targets, type_index
     integer, pointer :: n_parameters(:)
     character(len=pot_note_length) :: description
     character(len=param_note_length), pointer :: parameter_notes(:,:)
     logical :: includes_post_processing
  end type bond_order_descriptor

  type(potential_descriptor), pointer :: potential_descriptors(:)
  type(bond_order_descriptor), pointer :: bond_order_descriptors(:)

  ! Defines a particular potential.
  ! The potential should correspond to the description of some
  ! built-in type and hold actual numeric values for parameters. 
  ! In addition a real potential must have information on the 
  ! particles it acts on and the range it operates in.
  type potential
     integer :: type_index, pot_index
     double precision, pointer :: parameters(:), derived_parameters(:)
     double precision :: cutoff, soft_cutoff
     character(len=2), pointer :: apply_elements(:) ! label_length
     integer, pointer :: apply_tags(:), apply_indices(:)
     character(len=2), pointer :: original_elements(:) ! label_length
     integer, pointer :: original_tags(:), original_indices(:)
     logical :: filter_elements, filter_tags, filter_indices, smoothened, bond_order_argument
  end type potential

  ! Defines parameters for bond order factor calculation.
  type bond_order_parameters
     integer :: type_index, group_index ! group index connects the parameters to a coordinator entity in the python side of the simulator
     double precision, pointer :: parameters(:,:), derived_parameters(:,:)
     double precision :: cutoff, soft_cutoff
     integer, pointer :: n_params(:)
     character(len=2), pointer :: apply_elements(:) ! label_length
     character(len=2), pointer :: original_elements(:) ! label_length
     logical :: includes_post_processing
  end type bond_order_parameters
  
contains

  subroutine create_potential(n_targets,n_params,pot_name,parameters,cutoff,soft_cutoff,&
       elements,tags,indices,orig_elements,orig_tags,orig_indices,pot_index,new_potential)
    implicit none
    integer, intent(in) :: n_targets, n_params, pot_index
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

    new_potential%pot_index = pot_index
    new_potential%bond_order_argument = descriptor%bond_order_argument


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



  subroutine create_bond_order_factor(n_targets,n_params,n_split,bond_name,parameters,param_split,&
       cutoff,soft_cutoff,elements,orig_elements,group_index,new_bond)
    implicit none
    integer, intent(in) :: n_targets, n_params, n_split, group_index
    integer, intent(in) :: param_split(n_split)
    character(len=*), intent(in) :: bond_name
    double precision, intent(in) :: parameters(n_params)
    double precision, intent(in) :: cutoff, soft_cutoff
    character(len=2), intent(in) :: elements(n_targets) ! label_length
    character(len=2), intent(in) :: orig_elements(n_targets) ! label_length
    type(bond_order_parameters), intent(out) :: new_bond
    type(bond_order_descriptor) :: descriptor
    integer :: max_split, accumulated_split, i

    call get_bond_descriptor(bond_name, descriptor)
    
    new_bond%group_index = group_index
    new_bond%type_index = descriptor%type_index
    nullify(new_bond%parameters)
    nullify(new_bond%n_params)
    max_split = MAXVAL(param_split)
    allocate(new_bond%parameters(max_split,n_split))
    allocate(new_bond%n_params(n_split))
    new_bond%parameters = 0.d0
    accumulated_split = 0
    do i = 1, n_split
       new_bond%n_params(i) = param_split(i)
       new_bond%parameters(1:param_split(i),i) = &
            parameters(accumulated_split+1:accumulated_split+param_split(i))
       accumulated_split = accumulated_split + param_split(i)
    end do
    new_bond%cutoff = cutoff
    new_bond%soft_cutoff = soft_cutoff
    nullify(new_bond%apply_elements)
    allocate(new_bond%apply_elements(n_targets))
    new_bond%apply_elements = elements
    nullify(new_bond%original_elements)
    allocate(new_bond%original_elements(n_targets))
    new_bond%original_elements = orig_elements
    new_bond%includes_post_processing = descriptor%includes_post_processing

    ! calculate derived parameters if necessary
    select case (new_bond%type_index)
    case(tersoff_index) ! tersoff
       nullify(new_bond%derived_parameters)
       allocate(new_bond%derived_parameters(0,n_split))
    case default
       nullify(new_bond%derived_parameters)
       allocate(new_bond%derived_parameters(0,n_split))
    end select

  end subroutine create_bond_order_factor


  subroutine post_process_bond_order_factor(raw_sum, bond_params, factor_out)
    implicit none
    double precision, intent(in) :: raw_sum
    type(bond_order_parameters), intent(in) :: bond_params
    double precision, intent(out) :: factor_out
    double precision :: beta, eta, dN
    
    select case(bond_params%type_index)
    case(tersoff_index) ! tersoff factor
       beta = bond_params%parameters(1,1)
       eta = bond_params%parameters(2,1)
       factor_out = ( 1 + ( beta*raw_sum )**eta )**(-1.d0/(2.d0*eta))

    case (c_scale_index) ! coordination correction scaling function

       dN = ( raw_sum - bond_params%parameters(2,1) ) * bond_params%parameters(3,1)
       factor_out = bond_params%parameters(1,1) * dN / (1.d0 + exp(bond_params%parameters(4,1)*dN))

    case default
       factor_out = raw_sum
    end select

  end subroutine post_process_bond_order_factor

  subroutine post_process_bond_order_gradient(raw_sum, raw_gradient, bond_params, factor_out)
    implicit none
    double precision, intent(in) :: raw_sum, raw_gradient(3)
    type(bond_order_parameters), intent(in) :: bond_params
    double precision, intent(out) :: factor_out(3)
    double precision :: beta, eta, inv_eta, dN, expo, inv_exp
    
    select case(bond_params%type_index)
    case(tersoff_index)
       beta = bond_params%parameters(1,1)
       eta = bond_params%parameters(2,1)
       inv_eta = -1.d0/(2.d0*eta)
       factor_out = inv_eta * ( 1.d0 + ( beta*raw_sum )**eta )**(inv_eta-1.d0) * &
            eta * (beta*raw_sum)**(eta-1.d0) * beta * raw_gradient

    case (c_scale_index) ! coordination correction scaling function

       dN = ( raw_sum - bond_params%parameters(2,1) ) * bond_params%parameters(3,1)
       expo = exp( bond_params%parameters(4,1)*dN )
       inv_exp = 1.d0 / (1.d0 + expo)
       factor_out = bond_params%parameters(1,1) * &
            bond_params%parameters(3,1) * ( inv_exp - &
            dN * inv_exp*inv_exp * bond_params%parameters(4,1) * expo ) * raw_gradient

    case default
       factor_out = raw_gradient
    end select

  end subroutine post_process_bond_order_gradient

  subroutine evaluate_bond_order_factor(n_targets,separations,distances,bond_params,factor,atoms)
    implicit none
    integer, intent(in) :: n_targets
    double precision, intent(in) :: separations(3,n_targets-1), distances(n_targets-1)
    type(bond_order_parameters), intent(in) :: bond_params(n_targets-1)
    double precision, intent(out) :: factor(n_targets)
    type(atom), optional, intent(in) :: atoms(n_targets)
    double precision :: r1, r2, cosine, decay1, decay2, xi1, xi2, gee1, gee2, &
         mu, a1, a2, cc1, dd1, h1, cc2, dd2, h2
    double precision :: tmp1(3), tmp2(3)

    factor = 0.d0

    select case (bond_params(1)%type_index)
    case(coordination_index) ! number of neighbors

       r1 = distances(1)
       if(r1 < bond_params(1)%cutoff .and. r1 > 0.d0)then
          call smoothening_factor(r1,bond_params(1)%cutoff,bond_params(1)%soft_cutoff,decay1)
          factor = decay1 ! symmetric, so factor(1) = factor(2)
       end if

    case(tersoff_index) ! tersoff bond-order factor

       ! note that the given distances and separation vectors must be for index pairs
       ! ij and ik (in the notation described in the documentation) since these are needed.
       !
       ! bond_params(1) should contain the ij parameters and bond_params(2) the ik ones,
       ! but it is only checked that bond_params(1) is actually of tersoff type since only
       ! the cutoffs of bond_params(2) are used

       if(.not.present(atoms))then
          return
       end if
       ! check that bond_params(1) is for indices (ij)
       if(bond_params(1)%original_elements(1) /= atoms(2)%element)then
          return
       end if
       if(bond_params(1)%original_elements(2) /= atoms(1)%element)then
          return
       end if
       ! check that bond_params(2) is for indices (ik)
       if(bond_params(2)%original_elements(1) /= atoms(2)%element)then
          return
       end if
       if(bond_params(2)%original_elements(2) /= atoms(3)%element)then
          return
       end if

       r1 = distances(1)
       r2 = distances(2)

       if( ( r2 < bond_params(2)%cutoff .and. r2 > 0.d0 ) .and. &
           ( r1 < bond_params(1)%cutoff .and. r1 > 0.d0 ) )then

          tmp1 = separations(1:3,1)
          tmp2 = separations(1:3,2)
          cosine = (tmp1 .o. tmp2) / ( r1 * r2 ) ! angle between ij and ik
          
          ! bond_params: beta_i, eta_i, mu_i, a_ij, c_ij, d_ij, h_ij
          mu = bond_params(1)%parameters(3,1)
          a1 = bond_params(1)%parameters(1,2)
          a2 = bond_params(2)%parameters(1,2)
          cc1 = bond_params(1)%parameters(2,2)*bond_params(1)%parameters(2,2)
          dd1 = bond_params(1)%parameters(3,2)*bond_params(1)%parameters(3,2)
          h1 = bond_params(1)%parameters(4,2)
          cc2 = bond_params(2)%parameters(2,2)*bond_params(2)%parameters(2,2)
          dd2 = bond_params(2)%parameters(3,2)*bond_params(2)%parameters(3,2)
          h2 = bond_params(2)%parameters(4,2)
          
          call smoothening_factor(r2,bond_params(2)%cutoff,bond_params(2)%soft_cutoff,decay2)
          call smoothening_factor(r1,bond_params(1)%cutoff,bond_params(1)%soft_cutoff,decay1)
          
          xi1 = decay1 * decay2 * exp( (a1 * (r1-r2))**mu ) 
          xi2 = decay1 * decay2 * exp( (a2 * (r2-r1))**mu ) 
          gee1 = 1 + cc1/dd1 - cc1/(dd1+(h1-cosine)*(h1-cosine))
          gee2 = 1 + cc2/dd2 - cc2/(dd2+(h2-cosine)*(h2-cosine))
          
          ! only the middle atom gets a contribution, so factor(1) = factor(3) = 0.0
          factor(2) = xi1*gee1 + xi2*gee2
          
       end if
    case (c_scale_index) ! coordination correction scaling function
       ! there is no contribution to raw sums
    case default
       ! if we have an invalid case, do nothing
    end select

  end subroutine evaluate_bond_order_factor


  ! Returns the gradients of bond order terms with respect to moving an atom.
  ! The returned array has three layers:
  ! gradient( coordinates, atom with respect to which we differentiate, atom whose factor is differentiated )
  ! So for example, for a three body term atom1-atom2-atom3, gradient(1,2,3) contains
  ! the x-coordinate (1), of the factor for atom2 (2), with respect to moving atom3 (3).
  subroutine evaluate_bond_order_gradient(n_targets,separations,distances,bond_params,gradient,atoms)
    implicit none
    integer, intent(in) :: n_targets
    double precision, intent(in) :: separations(3,n_targets-1), distances(n_targets-1)
    type(bond_order_parameters), intent(in) :: bond_params(n_targets-1)
    double precision, intent(out) :: gradient(3,n_targets,n_targets)
    type(atom), optional, intent(in) :: atoms(n_targets)
    double precision :: r1, r2, nablar1(3,3), nablar2(3,3), cosine, nablacosine(3,3), decay1, decay2,&
         nabladecay(3,3), unitvector(3,2), xi1, gee1, nablaxi1(3,3), nablagee1(3,3), &
         xi2, gee2, nablaxi2(3,3), nablagee2(3,3), &
         mu, a1, a2, cc1, dd1, h1, cc2, dd2, h2, dot, &
         ratio, exponent1a, exponent1b, exponent2a, exponent2b
    double precision :: tmp1(3), tmp2(3), tmp3(3), tmp4(3), tmp5(3), tmp6(3), tmpmat1(3,3), tmpmat2(3,3)

    gradient = 0.d0

    select case (bond_params(1)%type_index)
    case(coordination_index) ! number of neighbors

       r1 = distances(1)
       if(r1 < bond_params(1)%cutoff .and. r1 > 0.d0)then
          call smoothening_gradient(separations(1:3,1) / r1,r1,&
               bond_params(1)%cutoff,&
               bond_params(1)%soft_cutoff,&
               tmp1)
          gradient(1:3,1,1) = -tmp1(1:3)
          gradient(1:3,2,1) = -tmp1(1:3)
          gradient(1:3,1,2) = tmp1(1:3)
          gradient(1:3,2,2) = tmp1(1:3)
       end if

    case(tersoff_index) ! tersoff bond-order factor

       if(.not.present(atoms))then
          return
       end if
       ! check that bond_params(1) is for indices (ij)
       if(bond_params(1)%original_elements(1) /= atoms(2)%element)then
          return
       end if
       if(bond_params(1)%original_elements(2) /= atoms(1)%element)then
          return
       end if
       ! check that bond_params(2) is for indices (ik)
       if(bond_params(2)%original_elements(1) /= atoms(2)%element)then
          return
       end if
       if(bond_params(2)%original_elements(2) /= atoms(3)%element)then
          return
       end if

       r1 = distances(1)
       r2 = distances(2)

       if( ( r2 < bond_params(2)%cutoff .and. r2 > 0.d0 ) .and. &
           ( r1 < bond_params(1)%cutoff .and. r1 > 0.d0 ) )then

             tmp1 = separations(1:3,1)
             tmp2 = separations(1:3,2)
             unitvector(1:3,1) = tmp1 / r1
             unitvector(1:3,2) = tmp2 / r2
             dot = (tmp1 .o. tmp2)
             ratio = 1.d0 / ( r1 * r2 )
             cosine = dot * ratio ! angle between ij and ik

             ! gradients of the r_ij and r_ik vectors with respect to the positions 
             ! of the three particles i, j, and k

             ! gradient 1 affects atoms ij, so atoms 2 and 1
             nablar1 = 0.d0
             nablar1(1:3,2) = -unitvector(1:3,1)
             nablar1(1:3,1) = unitvector(1:3,1)
             ! gradient 2 affects atoms ik, so atoms 2 and 3
             nablar2 = 0.d0
             nablar2(1:3,2) = -unitvector(1:3,2)
             nablar2(1:3,3) = unitvector(1:3,2)

             ! gradient of the cos theta_ijk factor
             nablacosine = 0.d0
             tmp3 = tmp2 - dot / (r1*r1) * tmp1
             tmp4 = tmp1 - dot / (r2*r2) * tmp2
             nablacosine(1:3,1) = tmp3*ratio
             nablacosine(1:3,2) = -(tmp3+tmp4)*ratio
             nablacosine(1:3,3) = tmp4*ratio

             ! bond_params: beta_i, eta_i, mu_i, a_ij, c_ij, d_ij, h_ij
             mu = bond_params(1)%parameters(3,1)
             a1 = bond_params(1)%parameters(1,2)
             a2 = bond_params(2)%parameters(1,2)
             cc1 = bond_params(1)%parameters(2,2)*bond_params(1)%parameters(2,2)
             dd1 = bond_params(1)%parameters(3,2)*bond_params(1)%parameters(3,2)
             h1 = bond_params(1)%parameters(4,2)
             cc2 = bond_params(2)%parameters(2,2)*bond_params(2)%parameters(2,2)
             dd2 = bond_params(2)%parameters(3,2)*bond_params(2)%parameters(3,2)
             h2 = bond_params(2)%parameters(4,2)

             call smoothening_factor(r2,bond_params(2)%cutoff,bond_params(2)%soft_cutoff,decay2)
             call smoothening_factor(r1,bond_params(1)%cutoff,bond_params(2)%soft_cutoff,decay1)
             ! store the gradients of decay1 and decay2 in tmp5 and tmp6
             call smoothening_gradient(unitvector(1:3,2),r2,&
                  bond_params(2)%cutoff,bond_params(2)%soft_cutoff,tmp6)
             call smoothening_gradient(unitvector(1:3,1),r1,&
                  bond_params(1)%cutoff,bond_params(1)%soft_cutoff,tmp5)

             ! gradient 1 (tmp5) affects atoms ij, so atoms 2 and 1
             tmpmat1 = 0.d0
             tmpmat1(1:3,1) = tmp5
             tmpmat1(1:3,2) = -tmp5
             ! gradient 2 (tmp6) affects atoms ik, so atoms 2 and 3
             tmpmat2 = 0.d0
             tmpmat2(1:3,2) = -tmp6
             tmpmat2(1:3,3) = tmp6

             exponent1a = (a1 * (r1-r2))**(mu-1)
             exponent2a = (a2 * (r2-r1))**(mu-1)
             exponent1b = (a1 * (r1-r2))*exponent1a
             exponent2b = (a2 * (r2-r1))*exponent2a
             xi1 = exp( exponent1b ) 
             xi2 = exp( exponent2b )
             gee1 = 1 + cc1/dd1 - cc1/(dd1+(h1-cosine)*(h1-cosine))
             gee2 = 1 + cc2/dd2 - cc2/(dd2+(h2-cosine)*(h2-cosine))
             nablaxi1 = (tmpmat1 * decay2 + tmpmat2 * decay1) * xi1 + &
                  decay1*decay2*( exp( exponent1b ) * a1 * mu * exponent1a ) * (nablar1 - nablar2)
             nablaxi2 = (tmpmat2 * decay1 + tmpmat1 * decay2) * xi2 + &
                  decay1*decay2*( exp( exponent2b ) * a2 * mu * exponent2a ) * (nablar2 - nablar1)
             xi1 = xi1 * decay1 * decay2
             xi2 = xi2 * decay1 * decay2
             nablagee1 = - cc1 / ( (dd1+(h1-cosine)*(h1-cosine))*(dd1+(h1-cosine)*(h1-cosine)) ) * &
                  2.d0 * (h1-cosine) * nablacosine
             nablagee2 = - cc2 / ( (dd2+(h2-cosine)*(h2-cosine))*(dd2+(h2-cosine)*(h2-cosine)) ) * &
                  2.d0 * (h2-cosine) * nablacosine

             ! there is only contribution to the factor of the middle atom, so
             ! also the gradients of atoms 1 and 3 are 0.
             gradient(1:3,2,1) = nablaxi1(1:3,1) * gee1 + xi1 * nablagee1(1:3,1) + &
                  nablaxi2(1:3,1) * gee2 + xi2 * nablagee2(1:3,1)
             gradient(1:3,2,2) = nablaxi1(1:3,2) * gee1 + xi1 * nablagee1(1:3,2) + &
                  nablaxi2(1:3,2) * gee2 + xi2 * nablagee2(1:3,2)
             gradient(1:3,2,3) = nablaxi1(1:3,3) * gee1 + xi1 * nablagee1(1:3,3) + &
                  nablaxi2(1:3,3) * gee2 + xi2 * nablagee2(1:3,3)

       end if
       
    case (c_scale_index) ! coordination correction scaling function
       ! there is no contribution to raw sums
    end select

  end subroutine evaluate_bond_order_gradient


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

    case (mono_none_index) ! constant potential

       force(1:3,1) = 0.d0

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

    case (mono_none_index) ! constant potential

       energy = interaction%parameters(1)

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
    descriptors_created = .false.

  end subroutine clear_potential_characterizers


  subroutine clear_bond_order_factor_characterizers()
    implicit none

    if(bond_descriptors_created)then
       deallocate(bond_order_descriptors)
    else
       nullify(bond_order_descriptors)
    end if
    bond_descriptors_created = .false.

  end subroutine clear_bond_order_factor_characterizers


  subroutine initialize_potential_characterizers()
    implicit none
    integer :: index

    call clear_potential_characterizers()
    allocate(potential_descriptors(0:n_potential_types))
    index = 0

    ! Lennard-Jones potential
    index = index+1
    if(pair_lj_index /= index)then
       write(*,*) "Potential indices in the core do not match!"
    end if
    potential_descriptors(index)%type_index = index
    potential_descriptors(index)%bond_order_argument = .false.
    call pad_string('LJ', pot_name_length,potential_descriptors(index)%name)
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
       write(*,*) "Potential indices in the core do not match!"
    end if
    potential_descriptors(index)%type_index = index
    potential_descriptors(index)%bond_order_argument = .false.
    call pad_string('spring', pot_name_length,potential_descriptors(index)%name)
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
       write(*,*) "Potential indices in the core do not match!"
    end if
    potential_descriptors(index)%type_index = index
    potential_descriptors(index)%bond_order_argument = .false.
    call pad_string('force', pot_name_length,potential_descriptors(index)%name)
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
       write(*,*) "Potential indices in the core do not match!"
    end if
    potential_descriptors(index)%type_index = index
    potential_descriptors(index)%bond_order_argument = .false.
    call pad_string('bond_bend', pot_name_length,potential_descriptors(index)%name)
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

    ! constant potential
    index = index+1
    if(mono_none_index /= index)then
       write(*,*) "Potential indices in the core do not match!"
    end if
    potential_descriptors(index)%type_index = index
    potential_descriptors(index)%bond_order_argument = .false.
    call pad_string('constant', pot_name_length,potential_descriptors(index)%name)
    potential_descriptors(index)%n_parameters = 1
    potential_descriptors(index)%n_targets = 1
    allocate(potential_descriptors(index)%parameter_names(potential_descriptors(index)%n_parameters))
    allocate(potential_descriptors(index)%parameter_notes(potential_descriptors(index)%n_parameters))
    call pad_string('V', param_name_length,potential_descriptors(index)%parameter_names(1))
    call pad_string('potential value', param_note_length,potential_descriptors(index)%parameter_notes(1))
    call pad_string('Constant potential: V(r) = V', &
         pot_note_length,potential_descriptors(index)%description)

    ! null potential, for errorneous inquiries
    index = 0
    potential_descriptors(index)%type_index = index
    potential_descriptors(index)%bond_order_argument = .false.
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
  

  subroutine initialize_bond_order_factor_characterizers()
    implicit none
    integer :: index, i, max_params

    call clear_bond_order_factor_characterizers()
    allocate(bond_order_descriptors(0:n_bond_order_types))

    index = 0

    ! Coordination
    index = index+1
    if(coordination_index /= index)then
       write(*,*) "Bond-order indices in the core do not match!"
    end if
    bond_order_descriptors(index)%type_index = index
    call pad_string('neighbors',pot_name_length,bond_order_descriptors(index)%name)
    allocate(bond_order_descriptors(index)%n_parameters(2))
    bond_order_descriptors(index)%n_parameters(1) = 0
    bond_order_descriptors(index)%n_parameters(2) = 0
    bond_order_descriptors(index)%n_targets = 2
    bond_order_descriptors(index)%includes_post_processing = .false.
!    max_params = 1
!    allocate(bond_order_descriptors(index)%parameter_names(max_params,bond_order_descriptors(index)%n_targets))
!    allocate(bond_order_descriptors(index)%parameter_notes(max_params,bond_order_descriptors(index)%n_targets))
!    call pad_string('none',param_name_length,bond_order_descriptors(index)%parameter_names(1,1))
!    call pad_string('no parameters',param_note_length,bond_order_descriptors(index)%parameter_notes(1,1))
    call pad_string('Counter for the number of neighbors: b_i = sum f(r_ij)', &
         pot_note_length,bond_order_descriptors(index)%description)

    ! Tersoff bond-order factor
    index = index+1
    if(tersoff_index /= index)then
       write(*,*) "Bond-order indices in the core do not match!"
    end if
    bond_order_descriptors(index)%type_index = index
    call pad_string('tersoff',pot_name_length,bond_order_descriptors(index)%name)
    allocate(bond_order_descriptors(index)%n_parameters(3))
    bond_order_descriptors(index)%n_parameters(1) = 3
    bond_order_descriptors(index)%n_parameters(2) = 4
    bond_order_descriptors(index)%n_parameters(3) = 0
    max_params = 4
    bond_order_descriptors(index)%n_targets = 3
    bond_order_descriptors(index)%includes_post_processing = .true.
    allocate(bond_order_descriptors(index)%parameter_names(max_params,bond_order_descriptors(index)%n_targets))
    allocate(bond_order_descriptors(index)%parameter_notes(max_params,bond_order_descriptors(index)%n_targets))
    call pad_string('beta',param_name_length,bond_order_descriptors(index)%parameter_names(1,1))
    call pad_string('prefactor',param_note_length,bond_order_descriptors(index)%parameter_notes(1,1))
    call pad_string('eta',param_name_length,bond_order_descriptors(index)%parameter_names(2,1))
    call pad_string('overall exponent',param_note_length,bond_order_descriptors(index)%parameter_notes(2,1))
    call pad_string('mu',param_name_length,bond_order_descriptors(index)%parameter_names(3,1))
    call pad_string('decay exponent',param_note_length,bond_order_descriptors(index)%parameter_notes(3,1))
    call pad_string('a',param_name_length,bond_order_descriptors(index)%parameter_names(1,2))
    call pad_string('inverse decay factor',param_note_length,bond_order_descriptors(index)%parameter_notes(1,2))
    call pad_string('c',param_name_length,bond_order_descriptors(index)%parameter_names(2,2))
    call pad_string('angle term nominator',param_note_length,bond_order_descriptors(index)%parameter_notes(2,2))
    call pad_string('d',param_name_length,bond_order_descriptors(index)%parameter_names(3,2))
    call pad_string('angle term denominator 1',param_note_length,bond_order_descriptors(index)%parameter_notes(3,2))
    call pad_string('h',param_name_length,bond_order_descriptors(index)%parameter_names(4,2))
    call pad_string('angle term denominator 2',param_note_length,bond_order_descriptors(index)%parameter_notes(4,2))
!    call pad_string('none',param_name_length,bond_order_descriptors(index)%parameter_names(1,3))
!    call pad_string('no parameters',param_note_length,bond_order_descriptors(index)%parameter_notes(1,3))
    call pad_string('Tersoff bond-order: b_i = [ 1+( beta_i sum xi_ijk*g_ijk)^eta_i ]^(-1/eta_i), '//&
         'xi_ijk = f(r_ik)*exp(alpha_ij^m_i (r_ij-r_ik)^m_i), '// &
         'g_ijk = 1+c_ij^2/d_ij^2-c_ij^2/(d_ij^2+(h_ij^2-cos theta_ijk))', &
         pot_note_length,bond_order_descriptors(index)%description)
    

    ! Coordination correction
    index = index+1
    if(c_scale_index /= index)then
       write(*,*) "Bond-order indices in the core do not match!"
    end if
    bond_order_descriptors(index)%type_index = index
    call pad_string('c_scale',pot_name_length,bond_order_descriptors(index)%name)
    allocate(bond_order_descriptors(index)%n_parameters(1))
    bond_order_descriptors(index)%n_parameters(1) = 4
    bond_order_descriptors(index)%n_targets = 1
    bond_order_descriptors(index)%includes_post_processing = .true.
    max_params = 4
    allocate(bond_order_descriptors(index)%parameter_names(max_params,bond_order_descriptors(index)%n_targets))
    allocate(bond_order_descriptors(index)%parameter_notes(max_params,bond_order_descriptors(index)%n_targets))
    call pad_string('epsilon',param_name_length,bond_order_descriptors(index)%parameter_names(1,1))
    call pad_string('energy scale constant',param_note_length,bond_order_descriptors(index)%parameter_notes(1,1))
    call pad_string('N',param_name_length,bond_order_descriptors(index)%parameter_names(2,1))
    call pad_string('normal coordination number',param_note_length,bond_order_descriptors(index)%parameter_notes(2,1))
    call pad_string('C',param_name_length,bond_order_descriptors(index)%parameter_names(3,1))
    call pad_string('coordination difference scale constant',param_note_length,bond_order_descriptors(index)%parameter_notes(3,1))
    call pad_string('gamma',param_name_length,bond_order_descriptors(index)%parameter_names(4,1))
    call pad_string('exponential decay constant',param_note_length,bond_order_descriptors(index)%parameter_notes(4,1))
    call pad_string('Coordination difference scaling function: '//&
         'b_i(n_i) = epsilon_i dN_i / (1 + exp(gamma_i dN_i)); dN_i = C_i(n_i - N_i)', &
         pot_note_length,bond_order_descriptors(index)%description)


    ! null, for error inquiries
    index = 0
    bond_order_descriptors(index)%type_index = index
    call pad_string('null',pot_name_length,bond_order_descriptors(index)%name)
    allocate(bond_order_descriptors(index)%n_parameters(1))
    bond_order_descriptors(index)%n_parameters(1) = 1
    bond_order_descriptors(index)%n_targets = 1
    allocate(bond_order_descriptors(index)%parameter_names(bond_order_descriptors(index)%n_parameters(1),1))
    allocate(bond_order_descriptors(index)%parameter_notes(bond_order_descriptors(index)%n_parameters(1),1))
    call pad_string('null',param_name_length,bond_order_descriptors(index)%parameter_names(1,1))
    call pad_string('a dummy parameter',param_note_length,bond_order_descriptors(index)%parameter_notes(1,1))
    call pad_string('No such bond-order factor is available', &
         pot_note_length,bond_order_descriptors(index)%description)
    
    bond_descriptors_created = .true.

  end subroutine initialize_bond_order_factor_characterizers


  subroutine get_number_of_potentials(n_pots)
    implicit none
    integer, intent(out) :: n_pots

    n_pots = n_potential_types

  end subroutine get_number_of_potentials


  subroutine get_number_of_bond_order_factors(n_bond)
    implicit none
    integer, intent(out) :: n_bond

    n_bond = n_bond_order_types

  end subroutine get_number_of_bond_order_factors


  subroutine list_potentials(n_pots,pots)
    implicit none
    integer, intent(in) :: n_pots
    character(len=pot_name_length), dimension(n_pots), intent(out) :: pots
    integer :: i

    do i = 1, n_pots
       if (i > n_potential_types)then
          pots(i) = potential_descriptors(0)%name
       else
          pots(i) = potential_descriptors(i)%name
       end if
    end do

  end subroutine list_potentials

  subroutine list_bond_order_factors(n_bonds,bonds)
    implicit none
    integer, intent(in) :: n_bonds
    character(len=pot_name_length), dimension(n_bonds), intent(out) :: bonds
    integer :: i

    do i = 1, n_bonds
       if (i > n_bond_order_types)then
          bonds(i) = bond_order_descriptors(0)%name
       else
          bonds(i) = bond_order_descriptors(i)%name
       end if
    end do

  end subroutine list_bond_order_factors


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

  subroutine get_index_of_bond_order_factor(bond_name, index)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(out) :: index
    integer :: i

    index = 0
    do i = 1, n_bond_order_types
       if(bond_order_descriptors(i)%name == bond_name)then
          index = i
          return
       end if
    end do

  end subroutine get_index_of_bond_order_factor


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

  subroutine get_bond_descriptor(bond_name,descriptor)
    implicit none
    character(len=*), intent(in) :: bond_name
    type(bond_order_descriptor), intent(out) :: descriptor
    integer :: i

    do i = 1, n_bond_order_types
       if(bond_order_descriptors(i)%name == bond_name)then
          descriptor = bond_order_descriptors(i)
          return
       end if
    end do
    descriptor = bond_order_descriptors(0)

  end subroutine get_bond_descriptor


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

  subroutine get_names_of_parameters_of_bond_order_factor(bond_name, n_targets, param_names)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(in) :: n_targets
    character(len=param_name_length), pointer :: param_names(:)
    type(bond_order_descriptor) :: descriptor
    integer :: i

    call get_bond_descriptor(bond_name,descriptor)
    if(descriptor%n_targets < n_targets)then
       return
    end if
    nullify(param_names)
    allocate(param_names(descriptor%n_parameters(n_targets)))
    do i = 1, descriptor%n_parameters(n_targets)
       param_names(i) = descriptor%parameter_names(i,n_targets)
    end do

  end subroutine get_names_of_parameters_of_bond_order_factor


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

  subroutine get_descriptions_of_parameters_of_bond_order_factor(bond_name, &
       n_targets, param_notes)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(in) :: n_targets
    character(len=param_note_length), pointer :: param_notes(:)
    type(bond_order_descriptor) :: descriptor
    integer :: i

    call get_bond_descriptor(bond_name,descriptor)
    if(descriptor%n_targets < n_targets)then
       return
    end if
    nullify(param_notes)
    allocate(param_notes(descriptor%n_parameters(n_targets)))
    do i = 1, descriptor%n_parameters(n_targets)
       param_notes(i) = descriptor%parameter_notes(i,n_targets)
    end do

  end subroutine get_descriptions_of_parameters_of_bond_order_factor


  subroutine get_number_of_parameters_of_potential(pot_name,n_params)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: n_params
    type(potential_descriptor) :: descriptor

    call get_descriptor(pot_name,descriptor)
    n_params = descriptor%n_parameters

  end subroutine get_number_of_parameters_of_potential


  subroutine get_number_of_parameters_of_bond_order_factor(bond_name,n_targets,n_params)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(in) :: n_targets
    integer, intent(out) :: n_params
    type(bond_order_descriptor) :: descriptor

    call get_bond_descriptor(bond_name,descriptor)
    n_params = descriptor%n_parameters(n_targets)

  end subroutine get_number_of_parameters_of_bond_order_factor


  subroutine get_number_of_targets_of_potential(pot_name,n_target)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: n_target
    type(potential_descriptor) :: descriptor

    call get_descriptor(pot_name,descriptor)
    n_target = descriptor%n_targets

  end subroutine get_number_of_targets_of_potential

  subroutine get_number_of_targets_of_bond_order_factor(bond_name,n_target)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(out) :: n_target
    type(bond_order_descriptor) :: descriptor

    call get_bond_descriptor(bond_name,descriptor)
    n_target = descriptor%n_targets

  end subroutine get_number_of_targets_of_bond_order_factor

  subroutine get_number_of_targets_of_potential_index(pot_index,n_target)
    implicit none
    integer, intent(in) :: pot_index
    integer, intent(out) :: n_target

    n_target = potential_descriptors(pot_index)%n_targets

  end subroutine get_number_of_targets_of_potential_index

  subroutine get_number_of_targets_of_bond_order_factor_index(bond_index,n_target)
    implicit none
    integer, intent(in) :: bond_index
    integer, intent(out) :: n_target

    n_target = bond_order_descriptors(bond_index)%n_targets

  end subroutine get_number_of_targets_of_bond_order_factor_index


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


  subroutine is_valid_bond_order_factor(string,is_valid)
    implicit none
    character(len=*), intent(in) :: string
    logical, intent(out) :: is_valid
    integer :: i

    is_valid = .false.
    do i = 1, n_bond_order_types
       if(bond_order_descriptors(i)%name == string)then
          is_valid = .true.
          return
       end if
    end do

  end subroutine is_valid_bond_order_factor


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

  subroutine get_index_of_parameter_of_bond_order_factor(bond_name,param_name,index,n_targets)
    implicit none
    character(len=*), intent(in) :: bond_name
    character(len=*), intent(in) :: param_name
    integer, intent(out) :: index
    integer, intent(out) :: n_targets
    type(bond_order_descriptor) :: descriptor
    integer :: i,j

    call get_bond_descriptor(bond_name,descriptor)
    
    index = 0
    n_targets = 0
    do j = 1, descriptor%n_targets
       do i = 1, descriptor%n_parameters(j)
          if (descriptor%parameter_names(i,j) == param_name) then
             index = i
             n_targets = j
             return
          end if
       end do
    end do

  end subroutine get_index_of_parameter_of_bond_order_factor



  subroutine get_description_of_potential(pot_name,description)
    implicit none
    character(len=*), intent(in) :: pot_name
    character(len=pot_note_length), intent(out) :: description
    type(potential_descriptor) :: descriptor

    call get_descriptor(pot_name,descriptor)
    description = descriptor%description

  end subroutine get_description_of_potential

  subroutine get_description_of_bond_order_factor(bond_name,description)
    implicit none
    character(len=*), intent(in) :: bond_name
    character(len=pot_note_length), intent(out) :: description
    type(bond_order_descriptor) :: descriptor

    call get_bond_descriptor(bond_name,descriptor)
    description = descriptor%description

  end subroutine get_description_of_bond_order_factor

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

  subroutine bond_order_factor_is_in_group(factor,group_index,in_group)
    implicit none
    type(bond_order_parameters), intent(in) :: factor
    integer, intent(in) :: group_index
    logical, intent(out) :: in_group

    in_group = ( factor%group_index == group_index )

  end subroutine bond_order_factor_is_in_group

  subroutine bond_order_factor_affects_atom(factor,atom_in,affects,position)
    implicit none
    type(bond_order_parameters), intent(in) :: factor
    type(atom), intent(in) :: atom_in
    logical, intent(out) :: affects
    integer, optional, intent(in) :: position
    integer :: i

    affects = .false.

    if(present(position) .and. position <= size(factor%apply_elements) .and. position > 0)then

       if(factor%apply_elements(position) == atom_in%element)then
          affects = .true.
          return
       end if

    else

       do i = 1, size(factor%apply_elements)
          if(factor%apply_elements(i) == atom_in%element)then
             affects = .true.
             return
          end if
       end do
       
    end if

  end subroutine bond_order_factor_affects_atom




end module potentials
