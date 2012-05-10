! 
! Potentials contains the low-level routines for handling interactions.
! The module defines custom types for both describing the types of
! potentials and bond order factors (:data:`potential_descriptor`, :data:`bond_order_descriptor`) 
! as well as for storing the parameters of actual interactions in use
! for the Fortran calculations (:data:`potential`, :data:`bond_order_parameters`).
! Tools for creating the custom datatypes (:func:`create_potential`, :func:`create_bond_order_factor`)
! are provided.
!
! The types of potentials and bond order factors are defined using the types
! :data:`potential_descriptor` and :data:`bond_order_descriptor`. 
! These should be created at start-up and remain untouched during simulation.
! They are used by the Fortran core for checking the types of parameters a potential
! needs, for instance, but they are also accessible from the Python interface.
! Especially, upon creation of :class:`~pysic.Potential` and :class:`~pysic.BondOrderParameters` 
! instances, one needs to specify the type as a keyword. This keyword is then compared to the list of 
! characterizers in the core to determine the type of the interaction.
!
! The basic routines for calculating the actual forces and energies are also defined in
! this module (:func:`evaluate_energy`, :func:`evaluate_forces`, :func:`evaluate_bond_order_factor`,
! :func:`evaluate_bond_order_gradient`). However, these routines do not calculate the total potential
! energy of the system, :math:`V`, or the total forces acting on the particles,
! :math:`\mathbf{F} = -\nabla_\alpha V`. Instead, the routines evaluate the contributions from individual
! atoms, atom pairs, atom triplets, etc. For instance, let the total energy of the system be
!
! .. math::
!
!   V = \sum_p \left( \sum_i v^p_i + \sum_{(i,j)} v^p_{ij} + \sum_{(i,j,k)} v^p_{ijk} \right),
!
! where :math:`p` sums over the different potentials acting on the system and :math:`i`, :math:`(i,j)` and
! :math:`(i,j,k)` sum over all atoms, pairs and triplet, respectively. Then the energy terms :math:`v`
! are obtained from :func:`evaluate_energy`. In pseudo-code,
!
! .. math::
!
!   v^p_{S} = \mathtt{evaluate\_energy}(S,p),
!
! where :math:`S` is a set of atoms. The summation over potentials and atoms is done in :ref:`pysic_core`
! in :func:`calculate_energy`. Similarly for forces, the summation is carried out in :func:`calculate_forces`.
!
! The reason for separating the calculation of individual interaction terms to :ref:`potentials`
! and the overall summation to :ref:`pysic_core` is that only the core knows the current structure and
! interactions of the system. 
! It is the task of this module to tell the core how all the potentials behave given
! any local structure, but the overall system information is kept in the core. So during energy
! evaluation, :ref:`pysic_core` finds all local structures that possibly contribute with an interaction
! and asks :ref:`potentials` to calculate this contribution.
!
! Bond order factors are potential modifiers, not direct interactions themselves.
! In general, the factors are scalar functions defined per atom, for instance,
!
! .. math::
!
!    b^p_i = s^p_i\left( \sum_{(i,j)} c^p_{ij} + \sum_{(i,j,k)} c^p_{ijk} \right)
! 
! for a three-body factor, where :math:`c^p` are local contributions 
! (usually representing chemical bonds) and :math:`s^p_i` is a per atom scaling function. 
! The bond factors multiply the potentials :math:`p`
! leading to the total energy
!
! .. math::
!
!    V = \sum_p \left( \sum_i b^p_i v^p_i + \sum_{(i,j)} \frac{1}{2}(b^p_i+b^p_j) v^p_{ij} + \sum_{(i,j,k)} \frac{1}{2}(b^p_i+b^p_j+b^p_k) v^p_{ijk} \right).
!
! The corresponding total force on atom :math:`\alpha` is then
!
! .. math::
!
!    \mathbf{F}_{\alpha} = - \nabla_\alpha V = - \sum_p \left( \sum_i ((\nabla_\alpha b^p_i) v^p_i + b^p_i (\nabla_\alpha v^p_i) ) + \ldots \right).
!
! The contributions :math:`\mathbf{f}^p_\alpha = -\nabla_\alpha v^p`, :math:`c^p`, 
! and :math:`\nabla_\alpha c^p` are
! calculated in :func:`evaluate_forces`, :func:`evaluate_bond_order_factor`, 
! and :func:`evaluate_bond_order_gradient`.
! Application of the scaling functions :math:`s_i` and :math:`s_i'` on the sums 
! :math:`\sum_{(i,j)} c^p_{ij} + \sum_{(i,j,k)} c^p_{ijk}` is done in the routines
! :func:`post_process_bond_order_factor` and :func:`post_process_bond_order_gradient` to
! produce the actual bond order factors :math:`b^p_i` and gradients :math:`\nabla_\alpha b^p_i`.
! These sums, similarly to the energy and force summations, are evaluated with 
! :func:`core_calculate_bond_order_factors` in :ref:`pysic_core`.
!
! Note when adding potentials or bond order factors in the source code:
!
! The parameters defined in Potentials.f90 are used for determining the maximum sizes of arrays, 
! numbers of potentials and bond factors, and the internally used indices
! for them. When adding new potentials of bond factors, make sure to update
! the relevant numbers. Especially the number of potentials (:data:`n_potential_types`)
! or number of bond order factors (:data:`n_bond_order_types`) must be increased
! when more types are defined.
! 
! Also note that in :ref:`pysic_core`, some of these parameters are used for
! determining array sizes. However, the actual parameters are not used
! because f2py does not read the values from here. Therefore if you change
! a parameter here, search for its name in :ref:`pysic_core` to see if the
! name appears in a comment. That is an indicator that a numeric value
! must be updated accordingly.
!
module potentials
  use geometry
  use quaternions
  use utility
  use mpi
  implicit none

  !*********************************!
  ! EDIT WHEN ADDING NEW POTENTIALS !
  !*********************************!

  !***********************************!
  ! EDIT WHEN ADDING NEW BOND FACTORS !
  !***********************************!

  ! Parameters for potential descriptors
  ! *pot_name_length maximum length allowed for the names of potentials
  ! *parameter_name_length maximum length allowed for the names of parameters
  ! *n_potential_types number of different types of potentials known
  ! *n_bond_order_types number of different types of bond order factors known
  ! *n_max_param maximum number of parameters for a potential or for a bond factor, per number of targets
  ! *pot_note_length maximum lenght allowed for the description of the potential
  ! *param_note_length maximum length allowed for the descriptions of parameters
  integer, parameter :: pot_name_length = 11, &
       param_name_length = 10, &
       n_potential_types = 6, &
       n_bond_order_types = 3, &
       n_max_params = 12, &
       pot_note_length = 500, &
       param_note_length = 100


  !*********************************!
  ! EDIT WHEN ADDING NEW POTENTIALS !
  !*********************************!

  ! Indices for potential types.
  ! These are used internally so that when
  ! a potential is evaluated, its type index
  ! is compared against these values.
  ! It makes the code easier to read when
  ! it says 
  !    if(type == pair_lj_index)then
  ! instead  of
  !    if(type == 1)then
  !
  ! *pair_lj_index internal index for the Lennard-Jones potential
  ! *pair_spring_index internal index for the spring potential
  ! *mono_const_index internal index for the constant force potential
  ! *tri_bend_index internal index for the bond bending potential
  ! *mono_none_index internal index for the constant potential
  integer, parameter :: pair_lj_index = 1, &
       pair_spring_index = 2, &
       mono_const_index = 3, &
       tri_bend_index = 4, &
       pair_exp_index = 5, &
       mono_none_index = 6


  !***********************************!
  ! EDIT WHEN ADDING NEW BOND FACTORS !
  !***********************************!

  ! Indices for bond order factor types.
  ! Similar to potential types above.
  !
  ! *coordiantion_index internal index for the coordination bondorder factor
  ! *tersoff_index internal index for the Tersoff bond order factor
  ! *c_scale_index internal index for the coordination scaling potential
  integer, parameter :: coordination_index = 1, &
       tersoff_index = 2, &
       c_scale_index = 3

  ! *no_name The label for unlabeled atoms. In other words, there are routines that expect atomic symbols as arguments, but if there are no symbols to pass, this should be given to mark an empty entry.
  character(len=label_length), parameter :: no_name = "xx"

  ! *descriptors_created logical tag used for managing pointer allocations for potential descriptors
  ! *bond_descriptors_created logical tag used for managing pointer allocations for bond order factor descriptors
  logical :: descriptors_created = .false., bond_descriptors_created = .false.

  ! Description of a type of a potential.
  ! The type contains the name and description of the potential
  ! and the parameters it contains.
  ! The descriptors contain the information that the inquiry methods in
  ! the python interface fetch.
  !
  ! *name The name of the potential: this is a keyword according to which the potentials may be recognized.
  ! *parameter_names The names of the parameters of the potential: these are keywords according to which the parameters may be recognized.
  ! *n_parameters number of parameters
  ! *n_targets number of targets, i.e., interacting bodies
  ! *type_index The internal index of the potential. This can also be used for recognizing the potential and must therefore match the name. For instance, if name = 'LJ', type_index = :data:`pair_lj_index`.
  ! *description A description of the potential. This should contain the mathematical formulation as well as a short verbal explanation.
  ! *parameter_notes Descriptions of the parameters. The descriptions should be very short indicators such as 'spring constant' or 'energy coefficient'. For more detailed explanations, the proper documentation should be used.
  type potential_descriptor
     character(len=pot_name_length) :: name
     character(len=param_name_length), pointer :: parameter_names(:)
     integer :: n_parameters, n_targets, type_index
     character(len=pot_note_length) :: description
     character(len=param_note_length), pointer :: parameter_notes(:)
  end type potential_descriptor

  ! Description of a type of a potential.
  ! The type contains the name and description of the potential
  ! and the parameters it contains.
  ! The descriptors contain the information that the inquiry methods in
  ! the python interface fetch.
  !
  ! *name The name of the bond order factor: this is a keyword according to which the factor may be recognized.
  ! *parameter_names The names of the parameters of the bond order factor: these are keywords according to which the parameters may be recognized.
  ! *n_parameters number of parameters for each number of bodies (1-body parameters, 2-body parameters etc.)
  ! *n_targets number of targets, i.e., interacting bodies
  ! *type_index The internal index of the bond order factor. This can also be used for recognizing the factor and must therefore match the name. For instance, if name = 'neighbors', type_index = :data:`coordination_index`.
  ! *description A description of the bond order factor. This should contain the mathematical formulation as well as a short verbal explanation.
  ! *parameter_notes Descriptions of the parameters. The descriptions should be very short indicators such as 'spring constant' or 'energy coefficient'. For more detailed explanations, the proper documentation should be used.
  ! *includes_post_processing a logical tag specifying if there is a scaling function :math:`s_i` attached to the factor.
  type bond_order_descriptor
     character(len=pot_name_length) :: name
     character(len=param_name_length), pointer :: parameter_names(:,:)
     integer :: n_targets, type_index
     integer, pointer :: n_parameters(:)
     character(len=pot_note_length) :: description
     character(len=param_note_length), pointer :: parameter_notes(:,:)
     logical :: includes_post_processing
  end type bond_order_descriptor

  ! *potential_descriptors an array for storing descriptors for the different *types* of potentials
  type(potential_descriptor), pointer :: potential_descriptors(:)
  ! *bond_order_descriptors an array for storing descriptors for the different *types* of bond order factors
  type(bond_order_descriptor), pointer :: bond_order_descriptors(:)

  ! Defines a particular potential.
  ! The potential should correspond to the description of some
  ! built-in type and hold actual numeric values for parameters. 
  ! In addition, a real potential must have information on the 
  ! particles it acts on and the range it operates in.
  ! These are to be created based on the :class:`~pysic.Potential` objects in the Python
  ! interface when calculations are invoked.
  !
  ! *type_index The internal index of the potential *type*. This is used for recognizing the potential. Note that the potential instance does not have a name. If the name is needed, it can be obtained from the :data:`potential_descriptor` of the correct index.
  ! *pot_index The internal index of the *actual potential*. This is needed when bond order factors are included so that the factors may be joint with the correct potentials.
  ! *parameters numerical values for parameters
  ! *derived_parameters numerical values for parameters calculated based on the parameters specified by the user
  ! *cutoff The hard cutoff for the potential. If the atoms are farther away from each other than this, the potential does not affect them.
  ! *soft_cutoff The soft cutoff for the potential. If this is smaller than the hard cutoff, the potential is scaled to zero continuously when the interatomic distances grow from the soft to the hard cutoff.
  ! *apply_elements A list of elements (atomic symbols) the potential affects. E.g., for a Si-O potential, it would be ('Si','O'). Note that unlike in the Python interface, a single :data:`potential` only has one set of targets, and for multiple target options several :data:`potential` instances are created.
  ! *apply_tags A list of atom tags the potential affects. Note that unlike in the Python interface, a single :data:`potential` only has one set of targets, and for multiple target options several :data:`potential` instances are created.
  ! *apply_indices A list of atom indices the potential affects. Note that unlike in the Python interface, a single :data:`potential` only has one set of targets, and for multiple target options several :data:`potential` instances are created.
  ! *original_elements The list of elements (atomic symbols) of the original :class:`~pysic.Potential` in the Python interface from which this potential was created. Whereas the apply_* lists are used for finding all pairs and triplets of atoms for which the potential could act on, the original_* lists specify the roles of atoms in the interaction.
  ! *original_tags The list of atom tags of the original :class:`~pysic.Potential` in the Python interface from which this potential was created. Whereas the apply_* lists are used for finding all pairs and triplets of atoms for which the potential could act on, the original_* lists specify the roles of atoms in the interaction.
  ! *original_indices The list of atom indices of the original :class:`~pysic.Potential` in the Python interface from which this potential was created. Whereas the apply_* lists are used for finding all pairs and triplets of atoms for which the potential could act on, the original_* lists specify the roles of atoms in the interaction.
  ! *filter_elements a logical switch specifying whether the potential targets atoms based on the atomic symbols
  ! *filter_tags a logical switch specifying whether the potential targets atoms based on the atom tags
  ! *filter_indices a logical switch specifying whether the potential targets atoms based on the atom indices
  ! *smoothened logical switch specifying if a smooth cutoff is applied to the potential
  type potential
     integer :: type_index, pot_index
     double precision, pointer :: parameters(:), derived_parameters(:)
     double precision :: cutoff, soft_cutoff
     character(len=2), pointer :: apply_elements(:) ! label_length
     integer, pointer :: apply_tags(:), apply_indices(:)
     character(len=2), pointer :: original_elements(:) ! label_length
     integer, pointer :: original_tags(:), original_indices(:)
     logical :: filter_elements, filter_tags, filter_indices, smoothened
  end type potential

  ! Defines a particular bond order factor.
  ! The factor should correspond to the description of some
  ! built-in type and hold actual numeric values for parameters. 
  ! In addition a real bond order factor must have information on the 
  ! particles it acts on and the range it operates in.
  ! These are created based on the :class:`~pysic.BondOrderParameters` objects in the Python
  ! interface when calculations are invoked.
  !
  ! *type_index The internal index of the bond order factor *type*. This is used for recognizing the factor. Note that the bond order parameters instance does not have a name. If the name is needed, it can be obtained from the :data:`bond_order_descriptor` of the correct index.
  ! *group_index The internal index of the *potential* the bond order factor is modifying.
  ! *parameters numerical values for parameters
  ! *n_params array containing the numbers of parameters for different number of targets (1-body parameters, 2-body parameters, etc.)
  ! *derived_parameters numerical values for parameters calculated based on the parameters specified by the user
  ! *cutoff The hard cutoff for the bond order factor. If the atoms are farther away from each other than this, they do not contribute to the total bond order factor does not affect them.
  ! *soft_cutoff The soft cutoff for the bond order factor. If this is smaller than the hard cutoff, the bond contribution is scaled to zero continuously when the interatomic distances grow from the soft to the hard cutoff.
  ! *apply_elements A list of elements (atomic symbols) the factor affects. E.g., for Si-O bonds, it would be ('Si','O'). Note that unlike in the Python interface, a single :data:`bond_order_parameters` only has one set of targets, and for multiple target options several :data:`bond_order_parameters` instances are created.
  ! *original_elements The list of elements (atomic symbols) of the original :class:`~pysic.BondOrderParameters` in the Python interface from which this factor was created. Whereas the apply_elements lists are used for finding all pairs and triplets of atoms which could contribute to the bond order factor, the original_elements lists specify the roles of atoms in the factor.
  ! *includes_post_processing a logical switch specifying if there is a scaling function :math:`s_i` attached to the factor
  type bond_order_parameters
     integer :: type_index, group_index
     double precision, pointer :: parameters(:,:), derived_parameters(:,:)
     double precision :: cutoff, soft_cutoff
     integer, pointer :: n_params(:)
     character(len=2), pointer :: apply_elements(:) ! label_length
     character(len=2), pointer :: original_elements(:) ! label_length
     logical :: includes_post_processing
  end type bond_order_parameters
  
contains


  ! !!!: create_potential

  ! Returns a :data:`potential`.
  !
  ! The routine takes as arguments all the necessary parameters
  ! and returns a potential type wrapping them in one package.
  !
  ! *n_targets number of targets, i.e., interacting bodies
  ! *n_params number of parameters
  ! *pot_name name of the potential - a keyword that must match a name of one of the :data:`potential_descriptors`
  ! *parameters array of numerical values for the parameters
  ! *cutoff the hard cutoff for the potential
  ! *soft_cutoff the soft cutoff for the potential
  ! *elements the elements (atomic symbols) the potential acts on
  ! *tags the atom tags the potential acts on
  ! *indices the atom indices the potential acts on
  ! *orig_elements The elements (atomic symbols) in the :class:`~pysic.Potential` used for generating the potential. This is needed to specify the roles of the atoms in the interaction.
  ! *orig_tags The atom tags in the :class:`~pysic.Potential` used for generating the potential. This is needed to specify the roles of the atoms in the interaction.
  ! *orig_indices The atom indices in the :class:`~pysic.Potential` used for generating the potential. This is needed to specify the roles of the atoms in the interaction.
  ! *pot_index the internal index of the potential
  ! *new_potential the created :data:`potential`
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

    !*********************************!
    ! EDIT WHEN ADDING NEW POTENTIALS !
    !*********************************!

    ! calculate derived parameters if necessary
    select case (new_potential%type_index)
    case(tri_bend_index) ! bond bending
       nullify(new_potential%derived_parameters)
       allocate(new_potential%derived_parameters(1))
       new_potential%derived_parameters(1) = cos(parameters(2))
    case(pair_exp_index) ! charge-dep. exp.
       nullify(new_potential%derived_parameters)
       allocate(new_potential%derived_parameters(4))
       ! eta_i = [ln R_i,max/(R_i,max - R_i,min) ] / [ln Q_i,max - ln(Q_i,max - Q_i,min) ] 
       new_potential%derived_parameters(1) = ( log(parameters(3)/(parameters(3)-parameters(4))) ) / &
            ( log(parameters(5)/(parameters(5)-parameters(6))) ) 
       new_potential%derived_parameters(2) = ( log(parameters(7)/(parameters(7)-parameters(8))) ) / &
            ( log(parameters(9)/(parameters(9)-parameters(10))) )
       ! beta_i = (R_i,min - R_i,max)^(1/eta_i) / (Q_i,max - Q_i,min)
       new_potential%derived_parameters(3) = (parameters(4) - parameters(3))**(1.d0/new_potential%derived_parameters(1)) / &
            (parameters(5) - parameters(6))
       new_potential%derived_parameters(4) = (parameters(8) - parameters(7))**(1.d0/new_potential%derived_parameters(2)) / &
            (parameters(9) - parameters(10))

    case default
       nullify(new_potential%derived_parameters)
       allocate(new_potential%derived_parameters(1))
       new_potential%derived_parameters(1) = 0.0
    end select

    new_potential%pot_index = pot_index


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


  ! !!!: create_bond_order_factor

  ! Returns a :data:`bond_order_parameters`.
  !
  ! The routine takes as arguments all the necessary parameters
  ! and returns a bond order parameters type wrapping them in one package.
  !
  ! *group_index The internal index of the *potential* the bond order factor is modifying.
  ! *parameters numerical values for parameters as a one-dimensional array
  ! *n_targets number of targets, i.e., interacting bodies
  ! *n_params array containing the numbers of parameters for different number of targets (1-body parameters, 2-body parameters, etc.)
  ! *n_split number of groupings in the list of parameters, per number of bodies - should equal n_targets
  ! *param_split Array containing the numbers of 1-body, 2-body, etc. parameters. The parameters are given as a list, but a bond order factor may have parameters separately for different numbers of targets. This list specifies the number of parameters for each.
  ! *bond_name name of the bond order factor - a keyword that must match a name of one of the :data:`bond_order_descriptors`
  ! *cutoff The hard cutoff for the bond order factor. If the atoms are farther away from each other than this, they do not contribute to the total bond order factor does not affect them.
  ! *soft_cutoff The soft cutoff for the bond order factor. If this is smaller than the hard cutoff, the bond contribution is scaled to zero continuously when the interatomic distances grow from the soft to the hard cutoff.
  ! *elements a list of elements (atomic symbols) the factor affects
  ! *orig_elements the list of elements (atomic symbols) of the original :class:`~pysic.BondOrderParameters` in the Python interface from which this factor was created
  ! *new_bond the created :data:`bond_order_parameters`
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


    !***********************************!
    ! EDIT WHEN ADDING NEW BOND FACTORS !
    !***********************************!

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


  ! !!!: post_process_bond_order_factor

  ! Bond-order post processing, i.e., 
  ! application of per-atom scaling functions.
  !
  ! By post processing, we mean any operations done after calculating the
  ! sum of pair- and many-body terms. That is, if a factor is, say,
  !
  ! .. math::
  !
  !      b_i = f(\sum_j c_{ij}) = 1 + \sum_j c_{ij},
  !
  ! the :math:`\sum_j c_ij` would have been calculated already and the 
  ! operation :math:`f(x) = 1 + x` remains to be carried out.
  ! The post processing is done per atom regardless of if the
  ! bond factor is of a pair or many body type.
  ! 
  ! This routine applies the scaling function on the given
  ! bond order sum accoding to the given parameters.
  !
  ! *raw_sum the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
  ! *bond_params a :data:`bond_order_parameters` specifying the parameters
  ! *factor_out the calculated bond order factor :math:`b_i`
  subroutine post_process_bond_order_factor(raw_sum, bond_params, factor_out)
    implicit none
    double precision, intent(in) :: raw_sum
    type(bond_order_parameters), intent(in) :: bond_params
    double precision, intent(out) :: factor_out
    double precision :: beta, eta, dN


    !***********************************!
    ! EDIT WHEN ADDING NEW BOND FACTORS !
    !***********************************!
    
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


  ! !!!: post_process_bond_order_gradient

  ! Bond-order post processing, i.e., 
  ! application of per-atom scaling functions.
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
  ! This routine applies the scaling function on the given
  ! bond order sum and gradient accoding to the given parameters.
  !
  ! *raw_sum the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
  ! *raw_gradient the precalculated bond order gradient sum, :math:`\nabla_\alpha \sum_j c_ij` in the above example
  ! *bond_params a :data:`bond_order_parameters` specifying the parameters
  ! *factor_out the calculated bond order factor :math:`\nabla_\alpha b_i`
  subroutine post_process_bond_order_gradient(raw_sum, raw_gradient, bond_params, factor_out)
    implicit none
    double precision, intent(in) :: raw_sum, raw_gradient(3)
    type(bond_order_parameters), intent(in) :: bond_params
    double precision, intent(out) :: factor_out(3)
    double precision :: beta, eta, inv_eta, dN, expo, inv_exp

    !***********************************!
    ! EDIT WHEN ADDING NEW BOND FACTORS !
    !***********************************!
    
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


  ! !!!: evaluate_bond_order_factor

  ! Returns a bond order factor term.
  ! 
  ! By a bond order factor term, we mean the contribution from
  ! specific atoms, :math:`c_{ijk}`, appearing in the factor
  !
  ! .. math::
  !
  !       b_i = f(\sum_{jk} c_{ijk})
  ! 
  ! This routine evaluates the term :math:`c_{ij}` or :math:`c_{ijk}` for the given
  ! atoms :math:`ij` or :math:`ijk` according to the given parameters.
  !
  ! *n_targets number of targets
  ! *separations atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
  ! *distances atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
  ! *bond_params a :data:`bond_order_parameters` containing the parameters
  ! *atoms a list of the actual :data:`atom` objects for which the term is calculated
  ! *factor the calculated bond order term :math:`c`
  subroutine evaluate_bond_order_factor(n_targets,separations,distances,bond_params,factor,atoms)
    implicit none
    integer, intent(in) :: n_targets
    double precision, intent(in) :: separations(3,n_targets-1), &
         distances(n_targets-1)
    type(bond_order_parameters), intent(in) :: bond_params(n_targets-1)
    double precision, intent(out) :: factor(n_targets)
    type(atom), optional, intent(in) :: atoms(n_targets)
    double precision :: r1, r2, cosine, decay1, decay2, &
         xi1, xi2, gee1, gee2, &
         mu, a1, a2, cc1, dd1, h1, cc2, dd2, h2
    double precision :: tmp1(3), tmp2(3)

    factor = 0.d0

    !***********************************!
    ! EDIT WHEN ADDING NEW BOND FACTORS !
    !***********************************!

    select case (bond_params(1)%type_index)
    case(coordination_index) ! number of neighbors

       r1 = distances(1)
       if(r1 < bond_params(1)%cutoff .and. r1 > 0.d0)then
          ! coordination is simply the number of neighbors
          ! calculated as a sum of cutoff functions
          call smoothening_factor(r1,bond_params(1)%cutoff,&
               bond_params(1)%soft_cutoff,decay1)
          factor = decay1 ! symmetric, so factor(1) = factor(2)
       end if

    case(tersoff_index) ! tersoff bond-order factor

       ! note that the given distances and separation vectors 
       ! must be for index pairs ij and ik (in the notation 
       ! described in the documentation) since these are needed.
       !
       ! bond_params(1) should contain the ij parameters and 
       ! bond_params(2) the ik ones

       if(bond_params(1)%type_index /= bond_params(2)%type_index)then
          return
       end if

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

          ! tmp1 and tmp2 are the vectors r_ij, r_ik
          tmp1 = separations(1:3,1)
          tmp2 = separations(1:3,2)
          ! cosine of the angle between ij and ik
          cosine = (tmp1 .o. tmp2) / ( r1 * r2 )
          
          ! bond_params: mu_i, a_ij, a_ik, &
          !              c_ij^2, d_ij^2, h_ij, 
          !              c_ik^2, d_ik^2, h_ik
          mu = bond_params(1)%parameters(3,1)
          a1 = bond_params(1)%parameters(1,2)
          a2 = bond_params(2)%parameters(1,2)
          cc1 = bond_params(1)%parameters(2,2)*bond_params(1)%parameters(2,2)
          dd1 = bond_params(1)%parameters(3,2)*bond_params(1)%parameters(3,2)
          h1 = bond_params(1)%parameters(4,2)
          cc2 = bond_params(2)%parameters(2,2)*bond_params(2)%parameters(2,2)
          dd2 = bond_params(2)%parameters(3,2)*bond_params(2)%parameters(3,2)
          h2 = bond_params(2)%parameters(4,2)
          
          call smoothening_factor(r2,bond_params(2)%cutoff,&
               bond_params(2)%soft_cutoff,decay2)
          call smoothening_factor(r1,bond_params(1)%cutoff,&
               bond_params(1)%soft_cutoff,decay1)
          
          xi1 = decay1 * decay2 * exp( (a1 * (r1-r2))**mu ) 
          xi2 = decay1 * decay2 * exp( (a2 * (r2-r1))**mu ) 
          gee1 = 1 + cc1/dd1 - cc1/(dd1+(h1-cosine)*(h1-cosine))
          gee2 = 1 + cc2/dd2 - cc2/(dd2+(h2-cosine)*(h2-cosine))
          
          ! only the middle atom gets a contribution, 
          ! so factor(1) = factor(3) = 0.0
          factor(2) = xi1*gee1 + xi2*gee2
          
       end if
    case (c_scale_index) ! coordination correction scaling function
       ! there is no contribution to raw sums
    case default
       ! if we have an invalid case, do nothing
    end select

  end subroutine evaluate_bond_order_factor


  ! !!!: evaluate_bond_order_gradient

  ! Returns the gradients of bond order terms with respect to moving an atom.
  ! 
  ! By a bond order factor term, we mean the contribution from
  ! specific atoms, c_ijk, appearing in the factor
  !
  ! .. math::
  !
  !       b_i = f(\sum_{jk} c_{ijk})
  ! 
  ! This routine evaluates the gradient term :math:`\nabla_\alpha c_{ij}` or 
  ! :math:`\nabla_\alpha c_{ijk}` for the given atoms :math:`ij` or :math:`ijk` according to the given parameters.
  !
  ! The returned array has three dimensions:
  ! gradient( coordinates, atom with respect to which we differentiate, atom whose factor is differentiated )
  ! So for example, for a three body term atom1-atom2-atom3, gradient(1,2,3) contains
  ! the x-coordinate (1), of the factor for atom2 (2), with respect to moving atom3 (3).
  !
  ! *n_targets number of targets
  ! *separations atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
  ! *distances atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
  ! *bond_params a :data:`bond_order_parameters` containing the parameters
  ! *atoms a list of the actual :data:`atom` objects for which the term is calculated
  ! *gradient the calculated bond order term :math:`\nabla_\alpha c`
  subroutine evaluate_bond_order_gradient(n_targets,separations,distances,bond_params,gradient,atoms)
    implicit none
    integer, intent(in) :: n_targets
    double precision, intent(in) :: separations(3,n_targets-1), &
         distances(n_targets-1)
    type(bond_order_parameters), intent(in) :: bond_params(n_targets-1)
    double precision, intent(out) :: gradient(3,n_targets,n_targets)
    type(atom), optional, intent(in) :: atoms(n_targets)
    double precision :: r1, r2, nablar1(3,3), nablar2(3,3), &
         cosine, nablacosine(3,3), decay1, decay2,&
         nabladecay(3,3), unitvector(3,2), xi1, gee1, &
         nablaxi1(3,3), nablagee1(3,3), &
         xi2, gee2, nablaxi2(3,3), nablagee2(3,3), &
         mu, a1, a2, cc1, dd1, h1, cc2, dd2, h2, dot, &
         ratio, exponent1a, exponent1b, exponent2a, exponent2b
    double precision :: tmp1(3), tmp2(3), tmp3(3), tmp4(3), &
         tmp5(3), tmp6(3), tmpmat1(3,3), tmpmat2(3,3)

    gradient = 0.d0

    !***********************************!
    ! EDIT WHEN ADDING NEW BOND FACTORS !
    !***********************************!

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

       if(bond_params(1)%type_index /= bond_params(2)%type_index)then
          return
       end if

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
             ! cosine of the angle between ij and ik, theta_ijk
             cosine = dot * ratio 

             ! gradients of the r_ij and r_ik vectors with 
             ! respect to the positions 
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

          
             ! bond_params: mu_i, a_ij, a_ik, &
             !              c_ij^2, d_ij^2, h_ij, 
             !              c_ik^2, d_ik^2, h_ik
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

             ! tmpmat1 ad tmpmat2 store the gradients of decay1
             ! and decay2 with respect to moving any atom i, j, or k

             ! gradient 1 (tmp5) affects atoms ij, so it affects atoms 2 and 1
             tmpmat1 = 0.d0
             tmpmat1(1:3,1) = tmp5
             tmpmat1(1:3,2) = -tmp5
             ! gradient 2 (tmp6) affects atoms ik, so it affects atoms 2 and 3
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


  ! !!!: evaluate_electronegativity

  ! If a potential, say, :math:`U_{ijk}` depends on the charges of atoms :math:`q_i` 
  ! it will not only create a force,
  ! but also a difference in chemical potential :math:`\mu_i` for the atomic partial charges.
  ! Similarly to :func:`evaluate_forces`, this function evaluates the chemical
  ! 'force' on the atomic charges
  !
  ! .. math::
  !
  !    \chi_{\alpha,ijk} = -\mu_{\alpha,ijk} = -\frac{\partial U_{ijk}}{\partial q_\alpha}
  !
  ! To be consist the forces returned by :func:`evaluate_electronegativity` must be
  ! derivatives of the energies returned by :func:`evaluate_energy`.
  !
  ! *n_targets number of targets
  ! *separations atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
  ! *distances atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
  ! *interaction a :data:`potential` containing the parameters
  ! *atoms a list of the actual :data:`atom` objects for which the term is calculated
  ! *eneg the calculated electronegativity component :math:`\chi_{\alpha,ijk}`
  subroutine evaluate_electronegativity(n_targets,separations,distances,interaction,eneg,atoms)
    implicit none
    integer, intent(in) :: n_targets
    double precision, intent(in) :: separations(3,n_targets-1), distances(n_targets-1)
    type(potential), intent(in) :: interaction
    double precision, intent(out) :: eneg(n_targets)
    type(atom), intent(in) :: atoms(n_targets)
    double precision :: r1, r2, r3, r4, r5, r6, r7, r8, r9, ratio, d1, d2, d3, d4, d5, d6
    logical :: inverse_params
    
    eneg = 0.d0


    !*********************************!
    ! EDIT WHEN ADDING NEW POTENTIALS !
    !*********************************!

    ! The interaction type is decided based on the type index.
    ! This decides what kind of a function is applied and what 
    ! the parameters provided actually mean.
    select case (interaction%type_index)
    case (pair_exp_index) ! charge-dependent exponential potential
       
       ! If we have parameters for a pair A-B and we get a pair B-A, we
       ! must flip the per atom parameters
       inverse_params = .true.
       if(interaction%filter_elements)then
          if(interaction%original_elements(1) == atoms(1)%element)then
             inverse_params = .false.
          end if
       end if
       if(interaction%filter_tags)then
          if(interaction%original_tags(1) == atoms(1)%tags)then
             inverse_params = .false.
          end if
       end if
       if(interaction%filter_indices)then
          if(interaction%original_indices(1) == atoms(1)%index)then
             inverse_params = .false.
          end if
       end if

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          if(inverse_params)then
             r2 = interaction%parameters(7) ! R_i,max
             r3 = interaction%parameters(3) ! R_j,max
             r4 = interaction%parameters(9) ! Q_i,max
             r5 = interaction%parameters(5) ! Q_j,max
             r6 = interaction%derived_parameters(2) ! eta_i
             r7 = interaction%derived_parameters(1) ! eta_j
             r8 = interaction%derived_parameters(4) ! beta_i
             r9 = interaction%derived_parameters(3) ! beta_j
             d3 = interaction%parameters(12) ! xi_i
             d4 = interaction%parameters(11) ! xi_j
          else
             r2 = interaction%parameters(3) ! R_i,max
             r3 = interaction%parameters(7) ! R_j,max
             r4 = interaction%parameters(5) ! Q_i,max
             r5 = interaction%parameters(9) ! Q_j,max
             r6 = interaction%derived_parameters(1) ! eta_i
             r7 = interaction%derived_parameters(2) ! eta_j
             r8 = interaction%derived_parameters(3) ! beta_i
             r9 = interaction%derived_parameters(4) ! beta_j
             d3 = interaction%parameters(11) ! xi_i
             d4 = interaction%parameters(12) ! xi_j
          end if

          d5 = r8 * (r4 - atoms(1)%charge)
          d6 = r9 * (r5 - atoms(2)%charge)
          d1 = r2 + abs(d5)**r6
          d2 = r3 + abs(d6)**r7
          eneg(1) = interaction%parameters(1)* &
               exp(-interaction%parameters(2)*r1 + &
               0.5d0*(d3*d1 + d4*d2) ) ! this is just the energy
          eneg(2) = eneg(1) * d4 * 0.5d0 * r7 * abs(d6)**(r7-1.d0) * sign(1.d0,d6) * r9
          eneg(1) = eneg(1) * d3 * 0.5d0 * r6 * abs(d5)**(r6-1.d0) * sign(1.d0,d5) * r8
       end if

    end select

  end subroutine evaluate_electronegativity



  ! !!!: evaluate_forces

  ! Evaluates the forces due to an interaction between the given
  ! atoms. In other words, if the total force on atom :math:`\alpha` is
  !
  ! .. math::
  !
  !    \mathbf{F}_\alpha = \sum_{ijk} -\nabla_\alpha v_{ijk} = \sum \mathbf{f}_{\alpha,ijk},
  !
  ! this routine evaluates :math:`\mathbf{f}_{\alpha,ijk}` for :math:`\alpha = (i,j,k)` for the given
  ! atoms i, j, and k.
  !
  ! To be consist the forces returned by :func:`evaluate_forces` must be
  ! gradients of the energies returned by :func:`evaluate_energy`.
  !
  ! *n_targets number of targets
  ! *separations atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
  ! *distances atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
  ! *interaction a :data:`potential` containing the parameters
  ! *atoms a list of the actual :data:`atom` objects for which the term is calculated
  ! *force the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
  subroutine evaluate_forces(n_targets,separations,distances,interaction,force,atoms)
    implicit none
    integer, intent(in) :: n_targets
    double precision, intent(in) :: separations(3,n_targets-1), distances(n_targets-1)
    type(potential), intent(in) :: interaction
    double precision, intent(out) :: force(3,n_targets)
    type(atom), optional, intent(in) :: atoms(n_targets)
    double precision :: r1, r2, r3, r4, r5, r6, r7, r8, r9, ratio, d1, d2, d3, d4
    double precision :: tmp1(3), tmp2(3), tmp3(3), tmp4(3)
    logical :: inverse_params

    force = 0.d0

    !*********************************!
    ! EDIT WHEN ADDING NEW POTENTIALS !
    !*********************************!

    ! The interaction type is decided based on the type index.
    ! This decides what kind of a function is applied and what 
    ! the parameters provided actually mean.
    select case (interaction%type_index)
    case (pair_lj_index) ! lennard-jones

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          ratio = interaction%parameters(2) / r1
          r6 = ratio*ratio*ratio*ratio*ratio*ratio
          force(1:3,1) = interaction%parameters(1) * ( 6.d0*r6 - 12.d0*r6*r6 ) * separations(1:3,1) / (r1*r1)
          force(1:3,2) = -force(1:3,1)
       end if

    case (pair_spring_index) ! spring-potential

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          r2 = (r1 - interaction%parameters(2))
          force(1:3,1) = interaction%parameters(1) * r2 * separations(1:3,1) / r1
          force(1:3,2) = -force(1:3,1)
       end if

    case (mono_const_index) ! constant force

       force(1:3,1) = interaction%parameters(1:3)

    case (tri_bend_index) ! bond bending

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

    case (pair_exp_index) ! charge-dependent exponential potential
       
       ! If we have parameters for a pair A-B and we get a pair B-A, we
       ! must flip the per atom parameters
       inverse_params = .true.
       if(.not.present(atoms))then
          return
       end if
       if(interaction%filter_elements)then
          if(interaction%original_elements(1) == atoms(1)%element)then
             inverse_params = .false.
          end if
       end if
       if(interaction%filter_tags)then
          if(interaction%original_tags(1) == atoms(1)%tags)then
             inverse_params = .false.
          end if
       end if
       if(interaction%filter_indices)then
          if(interaction%original_indices(1) == atoms(1)%index)then
             inverse_params = .false.
          end if
       end if

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          if(inverse_params)then
             r2 = interaction%parameters(7) ! R_i,max
             r3 = interaction%parameters(3) ! R_j,max
             r4 = interaction%parameters(9) ! Q_i,max
             r5 = interaction%parameters(5) ! Q_j,max
             r6 = interaction%derived_parameters(2) ! eta_i
             r7 = interaction%derived_parameters(1) ! eta_j
             r8 = interaction%derived_parameters(4) ! beta_i
             r9 = interaction%derived_parameters(3) ! beta_j
             d3 = interaction%parameters(12) ! xi_i
             d4 = interaction%parameters(11) ! xi_j
          else
             r2 = interaction%parameters(3) ! R_i,max
             r3 = interaction%parameters(7) ! R_j,max
             r4 = interaction%parameters(5) ! Q_i,max
             r5 = interaction%parameters(9) ! Q_j,max
             r6 = interaction%derived_parameters(1) ! eta_i
             r7 = interaction%derived_parameters(2) ! eta_j
             r8 = interaction%derived_parameters(3) ! beta_i
             r9 = interaction%derived_parameters(4) ! beta_j
             d3 = interaction%parameters(11) ! xi_i
             d4 = interaction%parameters(12) ! xi_j
          end if

          d1 = r2 + abs(r8 * (r4 - atoms(1)%charge))**r6
          d2 = r3 + abs(r9 * (r5 - atoms(2)%charge))**r7
          force(1:3,1) = -interaction%parameters(2)*interaction%parameters(1)* &
               exp(-interaction%parameters(2)*r1 + &
               0.5d0*(d3*d1 + d4*d2) ) * &
               separations(1:3,1) / r1
          force(1:3,2) = -force(1:3,1)

       end if

    case (mono_none_index) ! constant potential

       force(1:3,1) = 0.d0

    end select

  end subroutine evaluate_forces

  ! !!!: evaluate_energy

  ! Evaluates the potential energy due to an interaction between the given
  ! atoms. In other words, if the total potential energy is
  !
  ! .. math::
  !
  !    E = \sum_{ijk} v_{ijk}
  !
  ! this routine evaluates :math:`v_{ijk}` for the given
  ! atoms i, j, and k.
  ! 
  ! To be consist the forces returned by :func:`evaluate_forces` must be
  ! gradients of the energies returned by :func:`evaluate_energy`.
  !
  ! *n_targets number of targets
  ! *separations atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
  ! *distances atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
  ! *interaction a :data:`bond_order_parameters` containing the parameters
  ! *atoms a list of the actual :data:`atom` objects for which the term is calculated
  ! *energy the calculated energy :math:`v_{ijk}`
  subroutine evaluate_energy(n_targets,separations,distances,interaction,energy,atoms)
    implicit none
    integer, intent(in) :: n_targets
    double precision, intent(in) :: separations(3,n_targets-1), distances(n_targets-1)
    type(potential), intent(in) :: interaction
    double precision, intent(out) :: energy
    type(atom), optional, intent(in) :: atoms(n_targets)
    double precision :: r1, r2, r3, r4, r5, r6, r7, r8, r9, ratio, d1, d2, d3, d4
    logical :: inverse_params
    integer :: i

    energy = 0.d0

    !*********************************!
    ! EDIT WHEN ADDING NEW POTENTIALS !
    !*********************************!

    ! The interaction type is decided based on the type index.
    ! This decides what kind of a function is applied and what 
    ! the parameters provided actually mean.
    select case (interaction%type_index)
    case (pair_lj_index) ! lennard-jones

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then          
          ratio = interaction%parameters(2) / r1
          r6 = ratio*ratio*ratio*ratio*ratio*ratio
          energy = interaction%parameters(1) * (r6*r6 - r6)
       end if

    case (pair_spring_index) ! spring-potential

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          r2 = (r1 - interaction%parameters(2))
          r6 = (interaction%cutoff - interaction%parameters(2))
          energy = interaction%parameters(1) * 0.5d0 * (r2*r2 - r6*r6)
       end if

    case (mono_const_index) ! constant force

       energy = - interaction%parameters(1:3) .o. atoms(1)%position(1:3)

    case (tri_bend_index) ! bond-bending
       
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

    case (pair_exp_index) ! charge-dependent exponential potential
       
       ! If we have parameters for a pair A-B and we get a pair B-A, we
       ! must flip the per atom parameters
       inverse_params = .true.
       if(.not.present(atoms))then
          return
       end if
       if(interaction%filter_elements)then
          if(interaction%original_elements(1) == atoms(1)%element)then
             inverse_params = .false.
          end if
       end if
       if(interaction%filter_tags)then
          if(interaction%original_tags(1) == atoms(1)%tags)then
             inverse_params = .false.
          end if
       end if
       if(interaction%filter_indices)then
          if(interaction%original_indices(1) == atoms(1)%index)then
             inverse_params = .false.
          end if
       end if

       r1 = distances(1)
       if(r1 < interaction%cutoff .and. r1 > 0.d0)then
          if(inverse_params)then
             r2 = interaction%parameters(7) ! R_i,max
             r3 = interaction%parameters(3) ! R_j,max
             r4 = interaction%parameters(9) ! Q_i,max
             r5 = interaction%parameters(5) ! Q_j,max
             r6 = interaction%derived_parameters(2) ! eta_i
             r7 = interaction%derived_parameters(1) ! eta_j
             r8 = interaction%derived_parameters(4) ! beta_i
             r9 = interaction%derived_parameters(3) ! beta_j
             d3 = interaction%parameters(12) ! xi_i
             d4 = interaction%parameters(11) ! xi_j
          else
             r2 = interaction%parameters(3) ! R_i,max
             r3 = interaction%parameters(7) ! R_j,max
             r4 = interaction%parameters(5) ! Q_i,max
             r5 = interaction%parameters(9) ! Q_j,max
             r6 = interaction%derived_parameters(1) ! eta_i
             r7 = interaction%derived_parameters(2) ! eta_j
             r8 = interaction%derived_parameters(3) ! beta_i
             r9 = interaction%derived_parameters(4) ! beta_j
             d3 = interaction%parameters(11) ! xi_i
             d4 = interaction%parameters(12) ! xi_j
          end if

          d1 = r2 + abs(r8 * (r4 - atoms(1)%charge))**r6
          d2 = r3 + abs(r9 * (r5 - atoms(2)%charge))**r7
          energy = interaction%parameters(1)* &
               exp(-interaction%parameters(2)*r1 + &
               0.5d0*(d3*d1 + d4*d2) )
       end if

    case (mono_none_index) ! constant potential

       energy = interaction%parameters(1)

    end select

  end subroutine evaluate_energy


  ! !!!: smoothening_factor

  ! Function for smoothening potential and bond order cutoffs.
  ! In principle any "nice" function which goes from 1 to 0
  ! in a finite interval could be used. Here, we choose
  !
  ! .. math::
  !
  !   f(r) = \frac{1}{2} ( 1 + \cos \pi \frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}} ) 
  !
  ! for :math:`r \in [r_\mathrm{soft},r_\mathrm{hard}]`.
  ! 
  ! This routine takes as arguments :math:`r`, :math:`r_\mathrm{soft}`, and :math:`r_\mathrm{hard}`, and
  ! returns the value of the smoothening function.
  !
  ! *r distance :math:`r`
  ! *hard_cut the hard cutoff :math:`r_\mathrm{hard}`
  ! *soft_cut the soft cutoff :math:`r_\mathrm{soft}`
  ! *factor the calculated smoothening factor
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


  ! !!!: smoothening_derivative

  ! Derivative of the function for smoothening potential 
  ! and bond order cutoffs.
  ! In principle any "nice" function which goes from 1 to 0
  ! in a finite interval could be used. Here, we choose
  !
  ! .. math::
  !
  !   f(r) = \frac{1}{2} ( 1 + \cos \pi \frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}} ) 
  !
  ! for :math:`r \in [r_\mathrm{soft},r_\mathrm{hard}]`.
  ! The derivative is then
  !
  ! .. math::
  !
  !   f'(r) = \frac{\pi}{2 (r_\mathrm{soft}-r_\mathrm{hard})} \sin \pi \frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}}.
  !
  ! This routine takes as arguments :math:`r`, :math:`r_\mathrm{soft}`, and :math:`r_\mathrm{hard}`, and
  ! returns the value of the derivative of the smoothening function.
  !
  ! *r distance :math:`r`
  ! *hard_cut the hard cutoff :math:`r_\mathrm{hard}`
  ! *soft_cut the soft cutoff :math:`r_\mathrm{soft}`
  ! *factor the calculated derivative of the smoothening factor
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


  ! !!!: smoothening_gradient

  ! Gradient of the function for smoothening potential 
  ! and bond order cutoffs.
  ! In principle any "nice" function which goes from 1 to 0
  ! in a finite interval could be used. Here, we choose
  !
  ! .. math::
  !
  !   f(r) = \frac{1}{2} ( 1 + \cos \pi \frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}} ) 
  !
  ! for :math:`r \in [r_\mathrm{soft},r_\mathrm{hard}]`.
  ! The derivative is then
  !
  ! .. math::
  !
  !   f'(r) = \frac{\pi}{2 (r_\mathrm{soft}-r_\mathrm{hard})} \sin \pi \frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}}.
  !
  ! and the gradient with respect to :math:`r`
  !
  ! .. math::
  !
  !   \nabla f(r) = f'(r) \nabla r = f'(r) \hat{r}
  !
  ! where :math:`\hat{r}` is the unit vector in the direction of :math:`\mathbf{r}`.
  !
  ! This routine takes as arguments :math:`\hat{r}`, :math:`r`, :math:`r_\mathrm{soft}`, and :math:`r_\mathrm{hard}`, and
  ! returns the value of the gradient of the smoothening function.
  ! 
  ! *unit_vector the vector :math:`\hat{r}`
  ! *r distance :math:`r`
  ! *hard_cut the hard cutoff :math:`r_\mathrm{hard}`
  ! *soft_cut the soft cutoff :math:`r_\mathrm{soft}`
  ! *gradient the calculated derivative of the smoothening factor
  subroutine smoothening_gradient(unit_vector,r,hard_cut,soft_cut,gradient)
    implicit none
    double precision, intent(in) :: unit_vector(3), r, hard_cut, soft_cut
    double precision, intent(out) :: gradient(3)
    double precision :: factor
    
    call smoothening_derivative(r,hard_cut,soft_cut,factor)
    gradient = unit_vector * factor

  end subroutine smoothening_gradient


  ! !!!: clear_potential_characterizers

  ! Deallocates all memory associated with potential characterizes.
  subroutine clear_potential_characterizers()
    implicit none

    if(descriptors_created)then
       deallocate(potential_descriptors)
    else
       nullify(potential_descriptors)
    end if
    descriptors_created = .false.

  end subroutine clear_potential_characterizers


  ! !!!: clear_bond_order_factor_characterizers

  ! Deallocates all memory associated with bond order factor characterizes.
  subroutine clear_bond_order_factor_characterizers()
    implicit none

    if(bond_descriptors_created)then
       deallocate(bond_order_descriptors)
    else
       nullify(bond_order_descriptors)
    end if
    bond_descriptors_created = .false.

  end subroutine clear_bond_order_factor_characterizers


  ! !!!: initialize_potential_characterizers

  ! Creates potential characterizers.
  !
  ! This routine is meant to be run once, as pysic is
  ! imported, to create the characterizers for
  ! potentials. Once created, they are accessible
  ! by both the python and fortran sides of pysic
  ! as a tool for describing the general structure
  ! of potential objects.
  subroutine initialize_potential_characterizers()
    implicit none
    integer :: index

    call clear_potential_characterizers()
    allocate(potential_descriptors(0:n_potential_types))
    index = 0



    !*********************************!
    ! EDIT WHEN ADDING NEW POTENTIALS !
    !*********************************!



    ! **** Lennard-Jones potential ****

    ! Index is a helper variable that is increased by one
    ! after each descriptor is created. It is stored as
    ! the type index of the potential and thus must be
    ! unique for each descriptor. The reason for
    ! using it here is that it allows one to just copy the 
    ! structure of the code creating one descriptor
    ! to make another one without having to change all the
    ! indices by hand.
    index = index+1

    ! This is a check-up that the indexing is consistent.
    ! Should this test fail, it means that the named indices
    ! pair_lj_index etc. are wrong, which will likely lead
    ! to errors in interpreting potentials elsewhere.
    ! This error cannot result from a user error, it is
    ! a bug in programming.
    if(pair_lj_index /= index)then
       write(*,*) "Potential indices in the core do not match!"
    end if

    ! Record type index
    potential_descriptors(index)%type_index = index

    ! Record the name of the potential.
    ! This is a keyword used for accessing the type of potential
    ! in pysic, also in the python interface.
    call pad_string('LJ', pot_name_length,potential_descriptors(index)%name)

    ! Record the number of parameters
    potential_descriptors(index)%n_parameters = 2

    ! Record the number of targets (i.e., is the potential 1-body, 2-body etc.)
    potential_descriptors(index)%n_targets = 2

    ! Allocate space for storing the parameter names and descriptions.
    allocate(potential_descriptors(index)%parameter_names(potential_descriptors(index)%n_parameters))
    allocate(potential_descriptors(index)%parameter_notes(potential_descriptors(index)%n_parameters))

    ! Record parameter names and descriptions.
    ! Names are keywords with which one can intuitively 
    ! and easily access the parameters in the python
    ! interface.
    ! Descriptions are short descriptions of the
    ! physical or mathematical meaning of the parameters.
    ! They can be viewed from the python interface to
    ! remind the user how to parameterize the potential.
    call pad_string('epsilon', param_name_length,potential_descriptors(index)%parameter_names(1))
    call pad_string('energy scale constant', param_note_length,potential_descriptors(index)%parameter_notes(1))
    call pad_string('sigma', param_name_length,potential_descriptors(index)%parameter_names(2))
    call pad_string('length scale constant', param_note_length,potential_descriptors(index)%parameter_notes(2))

    ! Record a description of the entire potential.
    ! This description can also be viewed in the python
    ! interface as a reminder of the properties of the
    ! potential.
    ! The description should contain the mathematical
    ! formulation of the potential as well as a short
    ! verbal description.
    call pad_string('A standard Lennard-Jones potential: V(r) = epsilon * ( (sigma/r)^12 - (sigma/r)^6 )', &
         pot_note_length,potential_descriptors(index)%description)

    ! **** spring potential ****

    index = index+1
    if(pair_spring_index /= index)then
       write(*,*) "Potential indices in the core do not match!"
    end if
    potential_descriptors(index)%type_index = index
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

    ! **** constant force potential ****

    index = index+1
    if(mono_const_index /= index)then
       write(*,*) "Potential indices in the core do not match!"
    end if
    potential_descriptors(index)%type_index = index
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
    
    ! **** bond-bending potential ****

    index = index+1
    if(tri_bend_index /= index)then
       write(*,*) "Potential indices in the core do not match!"
    end if
    potential_descriptors(index)%type_index = index
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


    ! **** charge dependent exp potential ****

    index = index+1
    if(pair_exp_index /= index)then
       write(*,*) "Potential indices in the core do not match!"
    end if
    potential_descriptors(index)%type_index = index
    call pad_string('exponential', pot_name_length,potential_descriptors(index)%name)
    potential_descriptors(index)%n_parameters = 12
    potential_descriptors(index)%n_targets = 2
    allocate(potential_descriptors(index)%parameter_names(potential_descriptors(index)%n_parameters))
    allocate(potential_descriptors(index)%parameter_notes(potential_descriptors(index)%n_parameters))
    call pad_string('epsilon', param_name_length,potential_descriptors(index)%parameter_names(1))
    call pad_string('energy scale constant', param_note_length,potential_descriptors(index)%parameter_notes(1))
    call pad_string('zeta', param_name_length,potential_descriptors(index)%parameter_names(2))
    call pad_string('inverse decay length', param_note_length,potential_descriptors(index)%parameter_notes(2))
    call pad_string('Rmax1', param_name_length,potential_descriptors(index)%parameter_names(3))
    call pad_string('maximum covalent radius for atom i', &
         param_note_length,potential_descriptors(index)%parameter_notes(3))
    call pad_string('Rmin1', param_name_length,potential_descriptors(index)%parameter_names(4))
    call pad_string('minimum covalent radius for atom i', &
         param_note_length,potential_descriptors(index)%parameter_notes(4))
    call pad_string('Qmax1', param_name_length,potential_descriptors(index)%parameter_names(5))
    call pad_string('maximum covalent charge for atom i', &
         param_note_length,potential_descriptors(index)%parameter_notes(5))
    call pad_string('Qmin1', param_name_length,potential_descriptors(index)%parameter_names(6))
    call pad_string('minimum covalent charge for atom i', &
         param_note_length,potential_descriptors(index)%parameter_notes(6))
    call pad_string('Rmax2', param_name_length,potential_descriptors(index)%parameter_names(7))
    call pad_string('maximum covalent radius for atom j', &
         param_note_length,potential_descriptors(index)%parameter_notes(7))
    call pad_string('Rmin2', param_name_length,potential_descriptors(index)%parameter_names(8))
    call pad_string('minimum covalent radius for atom j', &
         param_note_length,potential_descriptors(index)%parameter_notes(8))
    call pad_string('Qmax2', param_name_length,potential_descriptors(index)%parameter_names(9))
    call pad_string('maximum covalent charge for atom j', &
         param_note_length,potential_descriptors(index)%parameter_notes(9))
    call pad_string('Qmin2', param_name_length,potential_descriptors(index)%parameter_names(10))
    call pad_string('minimum covalent charge for atom j', &
         param_note_length,potential_descriptors(index)%parameter_notes(10))
    call pad_string('xi1', param_name_length,potential_descriptors(index)%parameter_names(11))
    call pad_string('inverse charge decay for atom i', &
         param_note_length,potential_descriptors(index)%parameter_notes(11))
    call pad_string('xi2', param_name_length,potential_descriptors(index)%parameter_names(12))
    call pad_string('inverse charge decay for atom j', &
         param_note_length,potential_descriptors(index)%parameter_notes(12))
    call pad_string('Charge-dependent exponential potential: '//&
         'V(r,q) = epsilon exp(-zeta r + (xi_i D_i(q_i) + xi_j D_j(q_j))/2 ) '//&
         'D_i(q) = R_i,max + |beta_i (Q_i,max - q)|^eta_i '//&
         'beta_i = (R_i,min - R_i,max)^(1/eta_i) / (Q_i,max - Q_i,min) '//&
         'eta_i = [ ln R_i,max - ln (R_i,max - R_i,min) ] / [ ln Q_i,max - ln (Q_i,max - Q_i,min) ]',&
         pot_note_length,potential_descriptors(index)%description)

    ! **** constant potential ****

    index = index+1
    if(mono_none_index /= index)then
       write(*,*) "Potential indices in the core do not match!"
    end if
    potential_descriptors(index)%type_index = index
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
  

  ! !!!: initialize_bond_order_factor_characterizers

  ! Creates bond order factor characterizers.
  !
  ! This routine is meant to be run once, as pysic is
  ! imported, to create the characterizers for
  ! bond order factors. Once created, they are accessible
  ! by both the python and fortran sides of pysic
  ! as a tool for describing the general structure
  ! of bond order factor objects.
  subroutine initialize_bond_order_factor_characterizers()
    implicit none
    integer :: index, i, max_params

    call clear_bond_order_factor_characterizers()
    allocate(bond_order_descriptors(0:n_bond_order_types))

    index = 0


    !***********************************!
    ! EDIT WHEN ADDING NEW BOND FACTORS !
    !***********************************!


    ! **** Coordination ****

    ! Index is a helper variable that is increased by one
    ! after each descriptor is created. It is stored as
    ! the type index of the potential and thus must be
    ! unique for each descriptor. The reason for
    ! using it here is that it allows one to just copy the 
    ! structure of the code creating one descriptor
    ! to make another one without having to change all the
    ! indices by hand.
    index = index+1

    ! This is a check-up that the indexing is consistent.
    ! Should this test fail, it means that the named indices
    ! coordination_index etc. are wrong, which will likely lead
    ! to errors in interpreting potentials elsewhere.
    ! This error cannot result from a user error, it is
    ! a bug in programming.
    if(coordination_index /= index)then
       write(*,*) "Bond-order indices in the core do not match!"
    end if

    ! Record type index
    bond_order_descriptors(index)%type_index = index

    ! Record the name of the bond order factor.
    ! This is a keyword used for accessing the type of factor
    ! in pysic, also in the python interface.
    call pad_string('neighbors',pot_name_length,bond_order_descriptors(index)%name)

    ! Record the number of targets
    bond_order_descriptors(index)%n_targets = 2 

    ! Allocate space for storing the numbers of parameters (for 1-body, 2-body etc.).
    allocate(bond_order_descriptors(index)%n_parameters(bond_order_descriptors(index)%n_targets))

    ! Record the number of parameters
    bond_order_descriptors(index)%n_parameters(1) = 0 ! no 1-body params
    bond_order_descriptors(index)%n_parameters(2) = 0 ! no 2-body params

    ! Record if the bond factor includes post processing (per-atom scaling)
    bond_order_descriptors(index)%includes_post_processing = .false.

    ! Record a description of the entire bond order factor.
    ! This description can also be viewed in the python
    ! interface as a reminder of the properties of the
    ! bond order factor.
    ! The description should contain the mathematical
    ! formulation of the factor as well as a short
    ! verbal description.
    call pad_string('Counter for the number of neighbors: b_i = sum f(r_ij)', &
         pot_note_length,bond_order_descriptors(index)%description)


    ! **** Tersoff bond-order factor ****

    index = index+1
    if(tersoff_index /= index)then
       write(*,*) "Bond-order indices in the core do not match!"
    end if
    bond_order_descriptors(index)%type_index = index
    call pad_string('tersoff',pot_name_length,bond_order_descriptors(index)%name)
    bond_order_descriptors(index)%n_targets = 3
    allocate(bond_order_descriptors(index)%n_parameters(bond_order_descriptors(index)%n_targets))
    bond_order_descriptors(index)%n_parameters(1) = 3 ! 3 1-body parameters
    bond_order_descriptors(index)%n_parameters(2) = 4 ! 4 2-body parameters
    bond_order_descriptors(index)%n_parameters(3) = 0 ! no 3-body parameters
    max_params = 4 ! maxval(bond_order_descriptors(index)%n_parameters)
    bond_order_descriptors(index)%includes_post_processing = .true.

    ! Allocate space for storing the parameter names and descriptions.
    allocate(bond_order_descriptors(index)%parameter_names(max_params,bond_order_descriptors(index)%n_targets))
    allocate(bond_order_descriptors(index)%parameter_notes(max_params,bond_order_descriptors(index)%n_targets))

    ! Record parameter names and descriptions.
    ! Names are keywords with which one can intuitively 
    ! and easily access the parameters in the python
    ! interface.
    ! Descriptions are short descriptions of the
    ! physical or mathematical meaning of the parameters.
    ! They can be viewed from the python interface to
    ! remind the user how to parameterize the potential.
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

    call pad_string('Tersoff bond-order: b_i = [ 1+( beta_i sum xi_ijk*g_ijk)^eta_i ]^(-1/eta_i), \n'//&
         'xi_ijk = f(r_ik)*exp(alpha_ij^m_i (r_ij-r_ik)^m_i), \n'// &
         'g_ijk = 1+c_ij^2/d_ij^2-c_ij^2/(d_ij^2+(h_ij^2-cos theta_ijk))', &
         pot_note_length,bond_order_descriptors(index)%description)
    

    ! **** Coordination correction ****
    index = index+1
    if(c_scale_index /= index)then
       write(*,*) "Bond-order indices in the core do not match!"
    end if

    bond_order_descriptors(index)%type_index = index
    call pad_string('c_scale',pot_name_length,bond_order_descriptors(index)%name)

    bond_order_descriptors(index)%n_targets = 1
    allocate(bond_order_descriptors(index)%n_parameters(bond_order_descriptors(index)%n_targets))
    bond_order_descriptors(index)%n_parameters(1) = 4 ! 4 1-body parameters

    bond_order_descriptors(index)%includes_post_processing = .true.
    max_params = 4 ! maxval(bond_order_descriptors(index)%n_parameters)

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

  
  !******************************************************!
  !                                                      !
  ! Below are the inquiry routines with which one can    !
  ! access the information stored in the potential and   !
  ! bond order factor descriptors from the python        !
  ! interface.                                           !
  !                                                      !
  ! These routines all follow the information path of    !
  ! pysic.py method                                      !
  !   -> PyInterface.f90 routine                         !
  !     -> Potentials.f90 routine                        !
  !       -> Potentials.f90 descriptors                  !
  !                                                      !
  !*******************************************************


  ! !!!: get_number_of_potentials

  ! Return the number of :data:`potential_descriptor`  known.
  !
  ! *n_pots number of potential types
  subroutine get_number_of_potentials(n_pots)
    implicit none
    integer, intent(out) :: n_pots

    n_pots = n_potential_types

  end subroutine get_number_of_potentials


  ! !!!: get_number_of_bond_order_factors

  ! Returns the number of :data:`bond_order_descriptor` known.
  !
  ! *n_bond number of bond order factor types
  subroutine get_number_of_bond_order_factors(n_bond)
    implicit none
    integer, intent(out) :: n_bond

    n_bond = n_bond_order_types

  end subroutine get_number_of_bond_order_factors


  ! !!!: list_potentials

  ! Returns the names of :data:`potential_descriptor` objects.
  !
  ! *n_pots number of potential types
  ! *pots names of the potential types
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


  ! !!!: list_bond_order_factors

  ! Returns the names of :data:`bond_order_descriptor` objects.
  !
  ! *n_bonds number of bond order factor types
  ! *bonds names of the bond order factor types
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


  ! !!!: get_index_of_potential

  ! Returns the index of a :data:`potential_descriptor` in the internal list of potential types :data:`potential_descriptors`.
  !
  ! *pot_name name of the potential - a keyword
  ! *index index of the potential in the internal array
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


  ! !!!: get_index_of_bond_order_factor

  ! Returns the index of a :data:`bond_order_descriptor` in the internal list of bond order factor types :data:`bond_order_descriptors`.
  !
  ! *bond_name name of the bond order factor - a keyword
  ! *index index of the potential in the internal array
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


  ! !!!: get_descriptor

  ! Returns the :data:`potential_descriptor` of a given name.
  !
  ! *pot_name name of the potential
  ! *descriptor the matching :data:`potential_descriptor` 
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



  ! !!!: get_bond_descriptor

  ! Returns the :data:`bond_order_descriptor` of a given name.
  !
  ! *bond_name name of the bond order factor
  ! *descriptor the matching :data:`bond_order_descriptor`
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



  ! !!!: get_names_of_parameters_of_potential

  ! Returns the names of parameters of a potential as a list of strings.
  !
  ! *pot_name name of the potential
  ! *param_names names of the parameters
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



  ! !!!: get_names_of_parameters_of_bond_order_factor

  ! Returns the names of parameters of a bond order factor as a list of strings.
  !
  ! *bond_name name of the bond order factor
  ! *n_targets number of targets
  ! *param_names names of the parameters
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



  ! !!!: get_descriptions_of_parameters_of_potentials

  ! Returns the descriptions of the parameters of a potential
  ! as a list of strings.
  !
  ! *pot_name name of the potential
  ! *param_notes descriptions of the parameters
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


  ! !!!: get_descriptions_of_parameters_of_bond_order_factor

  ! Returns the descriptions of the parameters of a bond order factor
  ! as a list of strings.
  ! 
  ! *bond_name name of the bond order factor
  ! *n_targets number of targets
  ! *param_notes descriptions of the parameters
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



  ! !!!: get_number_of_parameters_of_potential

  ! Returns the number of parameters of a potential.
  !
  ! *pot_name name of the potential
  ! *n_params number of parameters
  subroutine get_number_of_parameters_of_potential(pot_name,n_params)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: n_params
    type(potential_descriptor) :: descriptor

    call get_descriptor(pot_name,descriptor)
    n_params = descriptor%n_parameters

  end subroutine get_number_of_parameters_of_potential

  
  ! !!!: get_number_of_parameters_of_bond_order_factor

  ! Returns the number of parameters of a bond order factor as a list of strings,
  ! each element showing the number of parameters for that number of bodies.
  !
  ! *bond_name name of the bond order factor
  ! *n_targets number of targets
  ! *n_params number of parameters
  subroutine get_number_of_parameters_of_bond_order_factor(bond_name,n_targets,n_params)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(in) :: n_targets
    integer, intent(out) :: n_params
    type(bond_order_descriptor) :: descriptor

    call get_bond_descriptor(bond_name,descriptor)
    n_params = descriptor%n_parameters(n_targets)

  end subroutine get_number_of_parameters_of_bond_order_factor


  ! !!!: get_number_of_targets_of_potential

  ! Returns the number of targets (i.e., bodies) of a potential.
  !
  ! *pot_name name of the potential
  ! *n_target number of targets
  subroutine get_number_of_targets_of_potential(pot_name,n_target)
    implicit none
    character(len=*), intent(in) :: pot_name
    integer, intent(out) :: n_target
    type(potential_descriptor) :: descriptor

    call get_descriptor(pot_name,descriptor)
    n_target = descriptor%n_targets

  end subroutine get_number_of_targets_of_potential


  ! !!!: get_number_of_targets_of_bond_order_factor

  ! Returns the number of targets (i.e., bodies) of a bond order factor.
  !
  ! *bond_name name of the bond order factor
  ! *n_target number of targets
  subroutine get_number_of_targets_of_bond_order_factor(bond_name,n_target)
    implicit none
    character(len=*), intent(in) :: bond_name
    integer, intent(out) :: n_target
    type(bond_order_descriptor) :: descriptor

    call get_bond_descriptor(bond_name,descriptor)
    n_target = descriptor%n_targets

  end subroutine get_number_of_targets_of_bond_order_factor



  ! !!!: get_number_of_targets_of_potential_index

  ! Returns the number of targets (i.e., bodies) of a potential 
  ! specified by its index.
  ! 
  ! *pot_index index of the potential
  ! *n_target numner of targets
  subroutine get_number_of_targets_of_potential_index(pot_index,n_target)
    implicit none
    integer, intent(in) :: pot_index
    integer, intent(out) :: n_target

    n_target = potential_descriptors(pot_index)%n_targets

  end subroutine get_number_of_targets_of_potential_index


  ! !!!: get_number_of_targets_of_bond_order_factor_index

  ! Returns the number of targets (i.e., bodies) of a bond order factor
  ! specified by its index.
  !
  ! *bond_index index of the bond order factor
  ! *n_target number of targets
  subroutine get_number_of_targets_of_bond_order_factor_index(bond_index,n_target)
    implicit none
    integer, intent(in) :: bond_index
    integer, intent(out) :: n_target

    n_target = bond_order_descriptors(bond_index)%n_targets

  end subroutine get_number_of_targets_of_bond_order_factor_index


  ! !!!: is_valid_potential

  ! Returns true if the given keyword is the name of a potential
  ! and false otherwise.
  !
  ! *string name of a potential
  ! *is_valid true if string is a name of a potential
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


  ! !!!: is_valid_bond_order_factor

  ! Returns true if the given keyword is the name of a bond order factor
  ! and false otherwise.
  !
  ! *string name of a bond order factor
  ! *is_valid true if string is a name of a bond order factor
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


  ! !!!: get_index_of_parameter_of_potential

  ! Returns the index of a parameter of a potential in the
  ! internal list of parameters.
  !
  ! *pot_name name of the potential
  ! *param_name name of the parameter
  ! *index the index of the parameter
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



  ! !!!: get_index_of_parameter_of_bond_order_factor

  ! Returns the index of a parameter of a bond order factor in the
  ! internal list of parameters. Since bond order factors can have
  ! parameters for different number of targets, also the number of
  ! targets of this parameter is returned.
  !
  ! *bond_name name of the bond order factor
  ! *param_name name of the parameter
  ! *index index of the parameter
  ! *n_targets number of targets of the parameter
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



  ! !!!: get_description_of_potential

  ! Returns the description of a potential.
  !
  ! *pot_name name of the potential
  ! *description description of the potential
  subroutine get_description_of_potential(pot_name,description)
    implicit none
    character(len=*), intent(in) :: pot_name
    character(len=pot_note_length), intent(out) :: description
    type(potential_descriptor) :: descriptor

    call get_descriptor(pot_name,descriptor)
    description = descriptor%description

  end subroutine get_description_of_potential



  ! !!!: get_description_of_bond_order_factors

  ! Returns the description of a bond order factor.
  !
  ! *bond_name name of the bond order factor
  ! *description description of the bond order factor
  subroutine get_description_of_bond_order_factor(bond_name,description)
    implicit none
    character(len=*), intent(in) :: bond_name
    character(len=pot_note_length), intent(out) :: description
    type(bond_order_descriptor) :: descriptor

    call get_bond_descriptor(bond_name,descriptor)
    description = descriptor%description

  end subroutine get_description_of_bond_order_factor


  ! !!!: potential_affects_atom

  ! Tests whether the given potential affects the specific atom.
  !
  ! For potentials, the atoms are specified as valid targets by 
  ! the atomic symbol, index, or tag.
  !
  ! If position is not given, then the routine returns true if
  ! the atom can appear in the potential in any role.
  ! If position is given, then true is returned only if the atom
  ! is valid for that particular position.
  !
  ! For instance, in a 3-body potential A-B-C, the potential
  ! May be specified so that only certain elements are valid for
  ! positions A and C while some other elements are valid for B.
  ! In a water molecule, for instance, we could have an H-O-H
  ! bond bending potential, but no H-H-O potentials.
  !
  ! *interaction the :data:`potential`
  ! *atom_in the :data:`atom`
  ! *affects true if the potential affects the atom
  ! *position specifies the particular role of the atom in the interaction
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



  ! !!!: bond_order_factor_is_in_group

  ! Tests whether the given bond order factor is a member of a specific group,
  ! i.e., if it affects the potential specifiesd by the group index.
  !
  ! *factor the :data:`bond_order_parameters`
  ! *group_index the index for the potential
  ! *in_group true if the factor is a member of the group
  subroutine bond_order_factor_is_in_group(factor,group_index,in_group)
    implicit none
    type(bond_order_parameters), intent(in) :: factor
    integer, intent(in) :: group_index
    logical, intent(out) :: in_group

    in_group = ( factor%group_index == group_index )

  end subroutine bond_order_factor_is_in_group



  ! !!!: bond_order_factor_affects_atom

  ! Tests whether the given bond order factor affects the specific atom.
  !
  ! For bond order factors, the atoms are specified as valid targets by 
  ! the atomic symbol only.
  !
  ! If position is not given, then the routine returns true if
  ! the atom can appear in the bond order factor in any role.
  ! If position is given, then true is returned only if the atom
  ! is valid for that particular position.
  !
  ! For instance, we may want to calculate the coordination of
  ! Cu-O bonds for Cu but not for O.
  !
  ! *factor the :data:`bond_order_parameters`
  ! *atom_in the :data:`atom`
  ! *affects true if the bond order factor is affected by the atom
  ! *position specifies the particular role of the atom in the bond order factor
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

  ! Calculates the energy of :math:`\frac{1}{r}` potentials through Ewald summation.
  !
  ! If a periodic system contains charges interacting via the :math:`\frac{1}{r}` Coulomb potential,
  ! direct summation of the interactions
  !
  ! .. math::
  !    :label: direct_sum
  !
  !    E = \sum_{(i,j)} \frac{1}{4\pi\epsilon_0}\frac{q_i q_j}{r_{ij}},
  !
  ! where the sum is over pairs of charges :math:`q_i, q_j` 
  ! (charges of the entire system, not just the simulation cell) and the distance between the charges is 
  ! :math:`r_{ij} = |\mathbf{r}_j - \mathbf{r}_i|`,
  ! does not work in general because the sum :eq:`direct_sum` converges very slowly. [#]_ Therefore truncating the
  ! sum may lead to severe errors.
  !
  ! The standard technique for overcoming this problem is the so called Ewald summation method.
  ! The idea is to split the long ranged and singular Coulomb potential to a short ranged singular and
  ! long ranged smooth parts, and calculate the long ranged part in reciprocal space via Fourier transformations.
  ! This is possible since the system is periodic and the same supercell repeats infinitely in all directions.
  ! In practice the calculation can be done by adding (and subtracting) Gaussian charge densities over the
  ! point charges to screen the
  ! potential in real space. That is, the original charge density 
  ! :math:`\rho(\mathbf{r}) = \sum_i q_i \delta(\mathbf{r} - \mathbf{r}_i)` is split by
  !
  ! .. math::
  !   :nowrap:
  !
  !   \begin{eqnarray}
  !   \rho(\mathbf{r}) & = & \rho_s(\mathbf{r}) + \rho_l(\mathbf{r}) \\
  !   \rho_s(\mathbf{r}) & = & \sum_i \left[ q_i \delta(\mathbf{r} - \mathbf{r}_i) - q_i G_\sigma(\mathbf{r} - \mathbf{r}_i) \right] \\
  !   \rho_l(\mathbf{r}) & = & \sum_i q_i G_\sigma(\mathbf{r} - \mathbf{r}_i) \\
  !   G_\sigma(\mathbf{r}) & = & \frac{1}{(2 \pi \sigma^2)^{3/2}} \exp\left( -\frac{|\mathbf{r}|^2}{2 \sigma^2} \right)
  !   \end{eqnarray}
  !
  ! Here :math:`\rho_l` generates a long range interaction since at large distances the Gaussian densities
  ! :math:`G_\sigma` appear the same as point charges 
  ! (:math:`\lim_{\sigma/r \to 0} G_\sigma(\mathbf{r}) = \delta(\mathbf{r})`). 
  ! Since the charge density is smooth, so will be the potential it creates.
  ! The density :math:`\rho_s` exhibits short ranged interactions for the same reason: 
  ! At distances longer than the width of the
  ! Gaussians the point charges are screened by the Gaussians which exactly cancel them
  ! (:math:`\lim_{\sigma/r \to 0} \delta(\mathbf{r}) - G_\sigma(\mathbf{r}) = 0`).
  !
  ! The short ranged interactions are directly calculated in real space
  !
  ! .. math::
  !    :nowrap:
  !
  !    \begin{eqnarray}
  !    E_s & = & \frac{1}{4 \pi \varepsilon_0} \int \frac{\rho_s(\mathbf{r}) \rho_s(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|} \mathrm{d}^3 r \mathrm{d}^3 r' \\
  !        & = & \frac{1}{4 \pi \varepsilon_0} \sum_{(i,j)} \frac{q_i q_j}{r_{ij}} \mathrm{erfc} \left( \frac{r_{ij}}{\sigma \sqrt{2}} \right).
  !    \end{eqnarray}
  !
  ! The complementary error function :math:`\mathrm{erfc}(r) = 1 - \mathrm{erf}(r) = 1 - \frac{2}{\sqrt{\pi}} \int_0^r e^{-t^2/2} \mathrm{d}t` makes the sum converge rapidly as :math:`\frac{r_{ij}}{\sigma} \to \infty`.
  !
  ! The long ranged interaction can be calculated in reciprocal space by Fourier transformation. The result is
  !
  ! .. math::
  !    :nowrap:
  !
  !    \begin{eqnarray}
  !    E_l & = & \frac{1}{2 V \varepsilon_0} \sum_{\mathbf{k} \ne 0} \frac{e^{-\sigma^2 k^2 / 2}}{k^2} |S(\mathbf{k})|^2 - \frac{1}{4 \pi \varepsilon_0} \frac{1}{\sqrt{2 \pi} \sigma} \sum_i^N q_i^2\\
  !    S(\mathbf{k}) & = & \sum_i^N q_i e^{\mathrm{i} \mathbf{k} \cdot \mathbf{r}_i} 
  !    \end{eqnarray}
  ! 
  ! The first sum in :math:`E_l` runs over the reciprocal lattice 
  ! :math:`\mathbf{k} = k_1 \mathbf{b}_1 + k_2 \mathbf{b}_2 + k_3 \mathbf{b}_3` where :math:`\mathbf{b}_i`
  ! are the vectors spanning the reciprocal cell (:math:`[\mathbf{b}_1 \mathbf{b}_2 \mathbf{b}_3] = ([\mathbf{v}_1 \mathbf{v}_2 \mathbf{v}_3]^{-1})^T` where :math:`\mathbf{v}_i` are the real space cell vectors).
  ! The latter sum is the self energy of each point charge in the potential of the particular Gaussian that 
  ! screens the charge, and the sum runs
  ! over all charges in the supercell spanning the periodic system. 
  ! (The self energy must be removed because it is present in the first sum even though when evaluating
  ! the potential at the position of a charge
  ! due to the other charges, no screening Gaussian function should be placed over the charge itself.)
  ! Likewise the sum in the structure factor :math:`S(\mathbf{k})` runs over all charges in the supercell.
  ! 
  ! The total energy is then the sum of the short and long range energies
  !
  ! .. math::
  !
  !    E = E_s + E_l.
  !
  ! .. [#] In fact, the sum converges only conditionally meaning the result depends on the order of summation. Physically this is not a problem, because one never has infinite lattices.
  !
  ! *n_atoms number of atoms
  ! *atoms list of atoms
  ! *cell the supercell containing the system
  ! *real_cutoff Cutoff radius of real-space interactions. Note that the neighbor lists stored in the atoms are used for neighbor finding so the cutoff cannot exceed the cutoff for the neighbor lists. (Or, it can, but the neighbors not in the lists will not be found.)
  ! *reciprocal_cutoff The number of cells to be included in the reciprocal sum in the directions of the reciprocal cell vectors. For example, if ``reciprocal_cutoff = [3,4,5]``, the reciprocal sum will be truncated as :math:`\sum_{\mathbf{k} \ne 0} = \sum_{k_1=-3}^3 \sum_{k_2=-4}^4 \sum_{k_3 = -5,(k_1,k_2,k_3) \ne (0,0,0)}^5`.
  ! *gaussian_width The :math:`\sigma` parameter, i.e., the distribution width of the screening Gaussians. This should not influence the actual value of the energy, but it does influence the convergence of the summation. If :math:`\sigma` is large, the real space sum :math:`E_s` converges slowly and a large real space cutoff is needed. If it is small, the reciprocal term :math:`E_l` converges slowly and the sum over the reciprocal lattice has to be evaluated over several cell lengths.
  ! *electric_constant The electic constant, i.e., vacuum permittivity :math:`\varepsilon_0`. In atomic units, it is :math:`\varepsilon_0 = 0.00552635 \frac{e^2}{\mathrm{Å\ eV}}`, but if one wishes to scale the results to some other unit system (such as reduced units with :math:`\varepsilon_0 = 1`), that is possible as well.
  ! *filter a list of logical values, one per atom, false for the atoms that should be ignored in the calculation
  ! *scaler a list of numerical values to scale the individual charges of the atoms
  ! *include_dipole_correction if true, a dipole correction term is included in the energy
  ! *total_energy the calculated energy
  subroutine calculate_ewald_energy(n_atoms,atoms,cell,real_cutoff,reciprocal_cutoff,gaussian_width,&
       electric_constant,filter,scaler,include_dipole_correction,total_energy)
    implicit none
    type(atom), intent(in) :: atoms(n_atoms)
    type(supercell), intent(in) :: cell
    double precision, intent(in) :: real_cutoff, gaussian_width, electric_constant, scaler(n_atoms)
    integer, intent(in) :: n_atoms, reciprocal_cutoff(3)
    double precision, intent(out) :: total_energy
    logical, intent(in) :: filter(n_atoms), include_dipole_correction
    double precision :: energy(7), tmp_energy(7), charge1, charge2, inv_eps_2v, inv_eps_4pi, &
         separation(3), distance, inv_sigma_sqrt_2pi, inv_sigma_sqrt_2, &
         s_factor(2, -reciprocal_cutoff(1):reciprocal_cutoff(1), &
         -reciprocal_cutoff(2):reciprocal_cutoff(2), &
         -reciprocal_cutoff(3):reciprocal_cutoff(3)), &
         tmp_factor(2, -reciprocal_cutoff(1):reciprocal_cutoff(1), &
         -reciprocal_cutoff(2):reciprocal_cutoff(2), &
         -reciprocal_cutoff(3):reciprocal_cutoff(3)), &
         k_vector(3), dot, coords2(3)
    integer :: index1, index2, j, k1, k2, k3, offset(3)
    type(atom) :: atom1, atom2
    type(neighbor_list) :: nbors1
    logical :: evaluate

    energy = 0.d0
    tmp_energy = 0.d0
    total_energy = 0.d0
    s_factor = 0.d0
    tmp_factor = 0.d0    

    inv_eps_4pi = 1.d0 / (4.d0 * pi * electric_constant)
    inv_eps_2v = 1.d0 / (2.d0 * cell%volume * electric_constant)
    inv_sigma_sqrt_2 = 1.d0 / (sqrt(2.d0) * gaussian_width)
    inv_sigma_sqrt_2pi = 1.d0 / (sqrt(2.d0 * pi) * gaussian_width)

    ! loop over atoms
    do index1 = 1, n_atoms
       if(is_my_atom(index1) .and. filter(index1))then
          atom1 = atoms(index1)
          nbors1 = atom1%neighbor_list
          charge1 = atom1%charge*scaler(index1)

          !
          ! calculate the real space sum
          !
          if(charge1 /= 0.d0)then
             ! loop over neighbors
             do j = 1, nbors1%n_neighbors
                
                ! neighboring atom
                index2 = nbors1%neighbors(j)
                
                ! prevent double counting (pick index1 < index2)
                if(pick(index1,index2,nbors1%pbc_offsets(1:3,j)))then

                   atom2 = atoms(index2)
                   charge2 = atom2%charge*scaler(index2)
                   ! calculate atom1-atom2 separation vector
                   ! and distance
                   call separation_vector(atom1%position, &
                        atom2%position, &
                        nbors1%pbc_offsets(1:3,j), &
                        cell, &
                        separation) ! in Geometry.f90
                   distance = .norm.separation
                   
                   if(distance < real_cutoff)then
                      !write(*,'(A,I8,I8,F8.2,F8.2,F8.2,F10.3)') "pair ", index1, index2, separation, distance

                      ! q_i q_j / r * erfc(r / (sqrt(2) sigma))
                      energy(1) = energy(1) + charge1*charge2/distance*(1.d0 - erf(distance*inv_sigma_sqrt_2))
                   end if

                end if

             end do

             !
             ! calculate the structure factors
             !
             do k1 = -reciprocal_cutoff(1), reciprocal_cutoff(1)
                do k2 = -reciprocal_cutoff(2), reciprocal_cutoff(2)
                   do k3 = -reciprocal_cutoff(3), reciprocal_cutoff(3)
                      if(k1 /= 0 .or. k2 /= 0 .or. k3 /= 0)then
                         
                         k_vector = k1*cell%reciprocal_cell(1:3,1) + &
                              k2*cell%reciprocal_cell(1:3,2) + &
                              k3*cell%reciprocal_cell(1:3,3)
                         dot = k_vector.o.atom1%position ! .o. in Quaternions.f90

                         ! S(k) = q exp(i k.r) = q [cos(k.r) + i sin(k.r)]
                         tmp_factor(1:2,k1,k2,k3) = tmp_factor(1:2,k1,k2,k3) + &
                              (/ charge1 * cos(dot), charge1 * sin(dot) /)
                         
                      end if
                   end do
                end do
             end do

             !
             ! calculate the self energy
             !
             energy(2) = energy(2) + charge1*charge1
             
             !
             ! calculate the charged background correction
             !
             energy(4) = energy(4) + charge1 
             
             !
             ! calculate the dipole correction
             !
             if(include_dipole_correction)then
                energy(5:7) = energy(5:7) + charge1 * atom1%position(1:3)
             end if
             
          end if

       end if
    end do

#ifdef MPI
    ! collect structure factors from all cpus in MPI
    call mpi_allreduce(tmp_factor,s_factor,size(s_factor),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
#else
    s_factor = tmp_factor
#endif
    
    !
    ! calculate the reciprocal space sum
    !
    do k1 = -reciprocal_cutoff(1), reciprocal_cutoff(1)
       do k2 = -reciprocal_cutoff(2), reciprocal_cutoff(2)
          do k3 = -reciprocal_cutoff(3), reciprocal_cutoff(3)
             if(k1 /= 0 .or. k2 /= 0 .or. k3 /= 0)then
                      
                k_vector = k1*cell%reciprocal_cell(1:3,1) + &
                     k2*cell%reciprocal_cell(1:3,2) + &
                     k3*cell%reciprocal_cell(1:3,3)
                distance = .norm.k_vector

                ! exp(- sigma^2 k^2 / 2) / k^2 |S(k)|^2
                ! |z|^2 = x^2 + y^2
                energy(3) = energy(3) + exp(-gaussian_width*gaussian_width*distance*distance*0.5d0) / &
                     (distance*distance) * &
                     (s_factor(1,k1,k2,k3)*s_factor(1,k1,k2,k3) + s_factor(2,k1,k2,k3)*s_factor(2,k1,k2,k3))
    
             end if
          end do
       end do
    end do

#ifdef MPI
    ! collect energies from all cpus in MPI (energy -> tmp_energy)
    call mpi_allreduce(energy,tmp_energy,size(energy),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
#else
    tmp_energy = energy
#endif

    ! multiply with leading coefficients
    energy(1) = tmp_energy(1) * inv_eps_4pi ! real space
    energy(2) = -tmp_energy(2) * inv_eps_4pi * inv_sigma_sqrt_2pi ! self energy
    energy(3) = tmp_energy(3) * inv_eps_2v ! reciprocal space
    energy(4) = - tmp_energy(4)*tmp_energy(4) * 0.5d0 * inv_eps_2v * gaussian_width*gaussian_width ! charged background
    energy(5) = (tmp_energy(5:7).o.tmp_energy(5:7)) * inv_eps_2v / 3.d0 ! dipole correction
    energy(6) = 0.d0
    energy(7) = 0.d0

    !write(*,'(A,F10.5,F10.5,F10.5,F10.5,F10.5)') "components", energy(1:5)

    total_energy = total_energy + sum(energy)


  end subroutine calculate_ewald_energy
 

  ! Calculates the forces due to long ranged :math:`\frac{1}{r}` potentials.
  ! These forces are the gradients of the energies :math:`U` given by :func:`calculate_ewald_energy`
  !
  ! .. math::
  !
  !    \mathbf{F}_\alpha = - \nabla_\alpha U
  !
  ! *n_atoms number of atoms
  ! *atoms list of atoms
  ! *cell the supercell containing the system
  ! *real_cutoff Cutoff radius of real-space interactions. Note that the neighbor lists stored in the atoms are used for neighbor finding so the cutoff cannot exceed the cutoff for the neighbor lists. (Or, it can, but the neighbors not in the lists will not be found.)
  ! *reciprocal_cutoff The number of cells to be included in the reciprocal sum in the directions of the reciprocal cell vectors. For example, if ``reciprocal_cutoff = [3,4,5]``, the reciprocal sum will be truncated as :math:`\sum_{\mathbf{k} \ne 0} = \sum_{k_1=-3}^3 \sum_{k_2=-4}^4 \sum_{k_3 = -5,(k_1,k_2,k_3) \ne (0,0,0)}^5`.
  ! *gaussian_width The :math:`\sigma` parameter, i.e., the distribution width of the screening Gaussians. This should not influence the actual value of the energy, but it does influence the convergence of the summation. If :math:`\sigma` is large, the real space sum :math:`E_s` converges slowly and a large real space cutoff is needed. If it is small, the reciprocal term :math:`E_l` converges slowly and the sum over the reciprocal lattice has to be evaluated over several cell lengths.
  ! *electric_constant The electic constant, i.e., vacuum permittivity :math:`\varepsilon_0`. In atomic units, it is :math:`\varepsilon_0 = 0.00552635 \frac{e^2}{\mathrm{Å\ eV}}`, but if one wishes to scale the results to some other unit system (such as reduced units with :math:`\varepsilon_0 = 1`), that is possible as well.
  ! *filter a list of logical values, one per atom, false for the atoms that should be ignored in the calculation
  ! *scaler a list of numerical values to scale the individual charges of the atoms
  ! *include_dipole_correction if true, a dipole correction term is included in the energy
  ! *total_forces the calculated forces
  subroutine calculate_ewald_forces(n_atoms,atoms,cell,real_cutoff,reciprocal_cutoff,gaussian_width,&
       electric_constant,filter,scaler,include_dipole_correction,total_forces)
    implicit none
    type(atom), intent(in) :: atoms(n_atoms)
    type(supercell), intent(in) :: cell
    double precision, intent(in) :: real_cutoff, gaussian_width, electric_constant, scaler(n_atoms)
    integer, intent(in) :: n_atoms, reciprocal_cutoff(3)
    double precision, intent(out) :: total_forces(3,n_atoms)
    logical, intent(in) :: filter(n_atoms), include_dipole_correction
    double precision :: forces(3,3,n_atoms), sum_forces(3,3,n_atoms), tmp_forces(3), charge1, charge2, &
         inv_eps_2v, inv_eps_4pi, sin_dot, cos_dot, &
         separation(3), distance, inv_dist, inv_sigma_sqrt_2pi, &
         inv_sigma_sqrt_2, inv_sigma_sqrt_2perpi, inv_sigma_sq_2, &
         s_factor(2, -reciprocal_cutoff(1):reciprocal_cutoff(1), &
         -reciprocal_cutoff(2):reciprocal_cutoff(2), &
         -reciprocal_cutoff(3):reciprocal_cutoff(3)), &
         tmp_factor(2, -reciprocal_cutoff(1):reciprocal_cutoff(1), &
         -reciprocal_cutoff(2):reciprocal_cutoff(2), &
         -reciprocal_cutoff(3):reciprocal_cutoff(3)), &
         nabla_factor(3, 2, n_atoms, -reciprocal_cutoff(1):reciprocal_cutoff(1), &
         -reciprocal_cutoff(2):reciprocal_cutoff(2), &
         -reciprocal_cutoff(3):reciprocal_cutoff(3)), &
         tmp_nabla_factor(3, 2, n_atoms, -reciprocal_cutoff(1):reciprocal_cutoff(1), &
         -reciprocal_cutoff(2):reciprocal_cutoff(2), &
         -reciprocal_cutoff(3):reciprocal_cutoff(3)), &
         k_vector(3), dot, &
         dipole(1:3), tmp_dipole(1:3)
    integer :: index1, index2, j, k1, k2, k3
    type(atom) :: atom1, atom2
    type(neighbor_list) :: nbors1

    forces = 0.d0
    tmp_forces = 0.d0
    total_forces = 0.d0
    s_factor = 0.d0
    tmp_factor = 0.d0
    nabla_factor = 0.d0
    tmp_nabla_factor = 0.d0

    inv_eps_4pi = 1.d0 / (4.d0 * pi * electric_constant)
    inv_eps_2v = 1.d0 / (2.d0 * cell%volume * electric_constant)
    inv_sigma_sqrt_2 = 1.d0 / (sqrt(2.d0) * gaussian_width)
    inv_sigma_sqrt_2pi = 1.d0 / (sqrt(2.d0 * pi) * gaussian_width)
    inv_sigma_sqrt_2 = 1.d0 / (sqrt(2.d0) * gaussian_width)
    inv_sigma_sq_2 = inv_sigma_sqrt_2*inv_sigma_sqrt_2
    inv_sigma_sqrt_2perpi = 2.d0 * inv_sigma_sqrt_2pi

    !
    ! calculate the real space sum
    !

    ! loop over atoms
    do index1 = 1, n_atoms
       if(is_my_atom(index1) .and. filter(index1))then
          atom1 = atoms(index1)
          nbors1 = atom1%neighbor_list
          charge1 = atom1%charge*scaler(index1)

          if(charge1 /= 0.d0)then
             ! loop over neighbors
             do j = 1, nbors1%n_neighbors
                
                ! neighboring atom
                index2 = nbors1%neighbors(j)
                
                ! prevent double counting
                if(pick(index1,index2,nbors1%pbc_offsets(1:3,j)))then
                   atom2 = atoms(index2)
                   charge2 = atom2%charge*scaler(index2)
                   
                   ! calculate atom1-atom2 separation vector
                   ! and distance
                   call separation_vector(atom1%position, &
                        atom2%position, &
                        nbors1%pbc_offsets(1:3,j), &
                        cell, &
                        separation) ! in Geometry.f90
                   distance = .norm.separation
                   inv_dist = 1.d0 / distance
                   
                   if(distance < real_cutoff)then
                      ! q_i q_j * 
                      ! ( erfc(r/(sigma*sqrt(2)))/r^2 + 
                      ! 1/sigma sqrt(2/pi) exp(-r^2/(2 sigma^2))/r ) \hat{r} 
                      tmp_forces =  -inv_eps_4pi * charge1*charge2 * ( &
                           ( 1 - erf(distance*inv_sigma_sqrt_2) ) * inv_dist*inv_dist + &
                           inv_sigma_sqrt_2perpi * exp( -distance*distance*inv_sigma_sq_2 ) * inv_dist  ) * &
                           separation * inv_dist
                      
                      forces(1:3,1,index1) = forces(1:3,1,index1) + tmp_forces
                      forces(1:3,1,index2) = forces(1:3,1,index2) - tmp_forces
                   end if

                end if
                
             end do

             !
             ! calculate the dipole correction terms
             !
             if(include_dipole_correction)then
                tmp_dipole(1:3) = tmp_dipole(1:3) + charge1 * atom1%position(1:3)
             end if

             !
             ! calculate the structure factors
             !             
             do k1 = -reciprocal_cutoff(1), reciprocal_cutoff(1)
                do k2 = -reciprocal_cutoff(2), reciprocal_cutoff(2)
                   do k3 = -reciprocal_cutoff(3), reciprocal_cutoff(3)
                      if(k1 /= 0 .or. k2 /= 0 .or. k3 /= 0)then
                         
                         k_vector = k1*cell%reciprocal_cell(1:3,1) + &
                              k2*cell%reciprocal_cell(1:3,2) + &
                              k3*cell%reciprocal_cell(1:3,3)
                         dot = k_vector.o.atom1%position ! .o. in Quaternions.f90
                         sin_dot = sin(dot)
                         cos_dot = cos(dot)
                         
                         ! this is the complex conjugate of S(k):
                         ! S*(k) = q exp(- i k.r) = q [cos(k.r) - i sin(k.r)]
                         tmp_factor(1:2,k1,k2,k3) = tmp_factor(1:2,k1,k2,k3) + &
                              (/ charge1 * cos_dot, - charge1 * sin_dot /)
                         
                         ! \nabla S(k) = i q k [cos(k.r) + i sin(k.r)]
                         tmp_nabla_factor(1:3,1,index1,k1,k2,k3) = &
                              -charge1 * sin_dot * k_vector ! real part
                         tmp_nabla_factor(1:3,2,index1,k1,k2,k3) = &
                              charge1 * cos_dot * k_vector ! imaginary part              
                         
                      end if
                   end do
                end do
             end do
          
          end if

          ! self energy or charged background do not contribute to forces
          
       end if
    end do

#ifdef MPI
    ! collect structure factors from all cpus in MPI (tmp_factor -> factor)
    call mpi_allreduce(tmp_factor,s_factor,size(s_factor),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
    call mpi_allreduce(tmp_nabla_factor,nabla_factor,size(nabla_factor),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
    if(include_dipole_correction)then
       call mpi_allreduce(tmp_dipole,dipole,1,mpi_double_precision,&
            mpi_sum,mpi_comm_world,mpistat)
    end if
#else
    s_factor = tmp_factor
    nabla_factor = tmp_nabla_factor
    dipole = tmp_dipole
#endif

    do index1 = 1, n_atoms
    
       !
       ! calculate the dipole correction
       ! 
       if(include_dipole_correction)then
          forces(1:3,3,index1) = atoms(index1)%charge*scaler(index1) * inv_eps_2v / (-1.5d0) * dipole(1:3)
       end if


       !
       ! calculate the reciprocal space sum
       !
       do k1 = -reciprocal_cutoff(1), reciprocal_cutoff(1)
          do k2 = -reciprocal_cutoff(2), reciprocal_cutoff(2)
             do k3 = -reciprocal_cutoff(3), reciprocal_cutoff(3)
                if(k1 /= 0 .or. k2 /= 0 .or. k3 /= 0)then
                   
                   k_vector = k1*cell%reciprocal_cell(1:3,1) + &
                        k2*cell%reciprocal_cell(1:3,2) + &
                        k3*cell%reciprocal_cell(1:3,3)
                   distance = .norm.k_vector
                   
                   ! - exp(- sigma^2 k^2 / 2) / k^2 2 Re[ S*(k) \nabla S(k) ]
                   ! Re[ z1 z2 ] = x1 x2 - y1 y2
                   forces(1:3,2,index1) = forces(1:3,2,index1) - inv_eps_2v * &
                        exp(-gaussian_width*gaussian_width*distance*distance*0.5d0) / &
                        (distance*distance) * 2.d0 * &
                        (s_factor(1,k1,k2,k3)*nabla_factor(1:3,1,index1,k1,k2,k3) - &
                        s_factor(2,k1,k2,k3)*nabla_factor(1:3,2,index1,k1,k2,k3))
                   
                end if
             end do
          end do
       end do
    end do

#ifdef MPI
    ! collect forces from all cpus in MPI
    call mpi_allreduce(forces,sum_forces,size(forces),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
    total_forces(1:3,1:n_atoms) = sum_forces(1:3,1,1:n_atoms) + &
         sum_forces(1:3,2,1:n_atoms) + sum_forces(1:3,3,1:n_atoms)
#else
    total_forces(1:3,1:n_atoms) = forces(1:3,1,1:n_atoms) + &
         forces(1:3,2,1:n_atoms) + forces(1:3,3,1:n_atoms)
#endif


  end subroutine calculate_ewald_forces


! ToDo: calculate_ewald_electronegativities


  ! Calculates the electronegativities due to long ranged :math:`\frac{1}{r}` potentials.
  ! These electronegativities are the derivatives of the energies :math:`U` given by :func:`calculate_ewald_energy`
  !
  ! .. math::
  !
  !    \chi_\alpha = - \frac{\partial U}{\partial q_\alpha}
  !
  ! *n_atoms number of atoms
  ! *atoms list of atoms
  ! *cell the supercell containing the system
  ! *real_cutoff Cutoff radius of real-space interactions. Note that the neighbor lists stored in the atoms are used for neighbor finding so the cutoff cannot exceed the cutoff for the neighbor lists. (Or, it can, but the neighbors not in the lists will not be found.)
  ! *reciprocal_cutoff The number of cells to be included in the reciprocal sum in the directions of the reciprocal cell vectors. For example, if ``reciprocal_cutoff = [3,4,5]``, the reciprocal sum will be truncated as :math:`\sum_{\mathbf{k} \ne 0} = \sum_{k_1=-3}^3 \sum_{k_2=-4}^4 \sum_{k_3 = -5,(k_1,k_2,k_3) \ne (0,0,0)}^5`.
  ! *gaussian_width The :math:`\sigma` parameter, i.e., the distribution width of the screening Gaussians. This should not influence the actual value of the energy, but it does influence the convergence of the summation. If :math:`\sigma` is large, the real space sum :math:`E_s` converges slowly and a large real space cutoff is needed. If it is small, the reciprocal term :math:`E_l` converges slowly and the sum over the reciprocal lattice has to be evaluated over several cell lengths.
  ! *electric_constant The electic constant, i.e., vacuum permittivity :math:`\varepsilon_0`. In atomic units, it is :math:`\varepsilon_0 = 0.00552635 \frac{e^2}{\mathrm{Å\ eV}}`, but if one wishes to scale the results to some other unit system (such as reduced units with :math:`\varepsilon_0 = 1`), that is possible as well.
  ! *filter a list of logical values, one per atom, false for the atoms that should be ignored in the calculation
  ! *scaler a list of numerical values to scale the individual charges of the atoms
  ! *include_dipole_correction if true, a dipole correction term is included
  ! *total_enegs the calculated electronegativities
  subroutine calculate_ewald_electronegativities(n_atoms,atoms,cell,real_cutoff,reciprocal_cutoff,gaussian_width,&
       electric_constant,filter,scaler,include_dipole_correction,total_enegs)
    implicit none
    type(atom), intent(in) :: atoms(n_atoms)
    type(supercell), intent(in) :: cell
    double precision, intent(in) :: real_cutoff, gaussian_width, electric_constant, scaler(n_atoms)
    integer, intent(in) :: n_atoms, reciprocal_cutoff(3)
    double precision, intent(out) :: total_enegs(n_atoms)
    logical, intent(in) :: filter(n_atoms), include_dipole_correction
    double precision :: tmp, qsum(4), tmp_qsum(4), &
         enegs(n_atoms), tmp_enegs(n_atoms), charge1, charge2, inv_eps_2v, inv_eps_4pi, &
         separation(3), distance, inv_sigma_sqrt_2pi, inv_sigma_sqrt_2, &
         s_factor(2, -reciprocal_cutoff(1):reciprocal_cutoff(1), &
         -reciprocal_cutoff(2):reciprocal_cutoff(2), &
         -reciprocal_cutoff(3):reciprocal_cutoff(3)), &
         tmp_factor(2, -reciprocal_cutoff(1):reciprocal_cutoff(1), &
         -reciprocal_cutoff(2):reciprocal_cutoff(2), &
         -reciprocal_cutoff(3):reciprocal_cutoff(3)), &
         diff_factor(2, n_atoms, -reciprocal_cutoff(1):reciprocal_cutoff(1), &
         -reciprocal_cutoff(2):reciprocal_cutoff(2), &
         -reciprocal_cutoff(3):reciprocal_cutoff(3)), &
         tmp_diff_factor(2, n_atoms, -reciprocal_cutoff(1):reciprocal_cutoff(1), &
         -reciprocal_cutoff(2):reciprocal_cutoff(2), &
         -reciprocal_cutoff(3):reciprocal_cutoff(3)), &
         k_vector(3), dot, sin_dot, cos_dot, coords2(3)
    integer :: index1, index2, j, k1, k2, k3, offset(3)
    type(atom) :: atom1, atom2
    type(neighbor_list) :: nbors1
    logical :: evaluate

    enegs = 0.d0
    tmp_enegs = 0.d0
    total_enegs = 0.d0
    s_factor = 0.d0
    tmp_factor = 0.d0
    diff_factor = 0.d0
    tmp_diff_factor = 0.d0
    qsum = 0.d0

    inv_eps_4pi = 1.d0 / (4.d0 * pi * electric_constant)
    inv_eps_2v = 1.d0 / (2.d0 * cell%volume * electric_constant)
    inv_sigma_sqrt_2 = 1.d0 / (sqrt(2.d0) * gaussian_width)
    inv_sigma_sqrt_2pi = 1.d0 / (sqrt(2.d0 * pi) * gaussian_width)

    ! loop over atoms
    do index1 = 1, n_atoms
       if(is_my_atom(index1) .and. filter(index1))then
          atom1 = atoms(index1)
          nbors1 = atom1%neighbor_list
          charge1 = atom1%charge*scaler(index1)

          !
          ! calculate the real space sum
          !
          if(charge1 /= 0.d0)then
             ! loop over neighbors
             do j = 1, nbors1%n_neighbors
                
                ! neighboring atom
                index2 = nbors1%neighbors(j)
                
                ! prevent double counting (pick index1 < index2)
                if(pick(index1,index2,nbors1%pbc_offsets(1:3,j)))then

                   atom2 = atoms(index2)
                   charge2 = atom2%charge*scaler(index2)
                   ! calculate atom1-atom2 separation vector
                   ! and distance
                   call separation_vector(atom1%position, &
                        atom2%position, &
                        nbors1%pbc_offsets(1:3,j), &
                        cell, &
                        separation) ! in Geometry.f90
                   distance = .norm.separation
                   
                   if(distance < real_cutoff)then

                      ! q_j / r * erfc(r / (sqrt(2) sigma))
                      tmp = -inv_eps_4pi/distance*(1.d0 - erf(distance*inv_sigma_sqrt_2))

                      enegs(index1) = enegs(index1) + charge2*tmp
                      enegs(index2) = enegs(index2) + charge1*tmp
                      
                   end if

                end if

             end do

             !
             ! calculate the structure factors
             !
             do k1 = -reciprocal_cutoff(1), reciprocal_cutoff(1)
                do k2 = -reciprocal_cutoff(2), reciprocal_cutoff(2)
                   do k3 = -reciprocal_cutoff(3), reciprocal_cutoff(3)
                      if(k1 /= 0 .or. k2 /= 0 .or. k3 /= 0)then
                         
                         k_vector = k1*cell%reciprocal_cell(1:3,1) + &
                              k2*cell%reciprocal_cell(1:3,2) + &
                              k3*cell%reciprocal_cell(1:3,3)
                         dot = k_vector.o.atom1%position ! .o. in Quaternions.f90
                         sin_dot = sin(dot)
                         cos_dot = cos(dot)

                         ! this is the complex conjugate of S(k):
                         ! S*(k) = q exp(- i k.r) = q [cos(k.r) - i sin(k.r)]
                         tmp_factor(1:2,k1,k2,k3) = tmp_factor(1:2,k1,k2,k3) + &
                              (/ charge1 * cos_dot, - charge1 * sin_dot /)

                         ! d S(k) = exp(i k.r) = [cos(k.r) + i sin(k.r)]
                         tmp_diff_factor(1:2,index1,k1,k2,k3) = tmp_diff_factor(1:2,index1,k1,k2,k3) + &
                              (/ cos_dot, sin_dot /)
                         
                      end if
                   end do
                end do
             end do

             !
             ! calculate the self energy term
             !
             enegs(index1) = enegs(index1) + 2.d0*charge1 * inv_eps_4pi * inv_sigma_sqrt_2pi 
             
             !
             ! calculate the total charge (charged background term)
             !
             tmp_qsum(1) = tmp_qsum(1) + charge1 
             
             !
             ! calculate the dipole correction
             !
             if(include_dipole_correction)then
                tmp_qsum(2:4) = tmp_qsum(2:4) + charge1 * atom1%position(1:3)
             end if
             
          end if

       end if
    end do

#ifdef MPI
    ! collect structure factors from all cpus in MPI
    call mpi_allreduce(tmp_factor,s_factor,size(s_factor),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
    call mpi_allreduce(tmp_diff_factor,diff_factor,size(diff_factor),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
    call mpi_allreduce(tmp_qsum,qsum,size(qsum),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
#else
    s_factor = tmp_factor
    diff_factor = tmp_diff_factor
    qsum = tmp_qsum
#endif
    
    do index1 = 1, n_atoms

       !
       ! calculate the total charged background and dipole correction terms
       !
       enegs(index1) = enegs(index1) + scaler(index1) * inv_eps_2v * &
            ( - gaussian_width*gaussian_width *qsum(1) )

    end do

    if(include_dipole_correction)then
       !
       ! dipole correction terms
       !
       enegs(index1) = enegs(index1) + scaler(index1) * inv_eps_2v * &
            ( -0.666666666666666667d0 ) * &
            ( atom1%position(1) * qsum(2) + &
            atom1%position(2) * qsum(3) + &
            atom1%position(3) * qsum(4) ) 
    end if

    do index1 = 1, n_atoms
       !
       ! calculate the reciprocal space sum
       !
       do k1 = -reciprocal_cutoff(1), reciprocal_cutoff(1)
          do k2 = -reciprocal_cutoff(2), reciprocal_cutoff(2)
             do k3 = -reciprocal_cutoff(3), reciprocal_cutoff(3)
                if(k1 /= 0 .or. k2 /= 0 .or. k3 /= 0)then
                   
                   k_vector = k1*cell%reciprocal_cell(1:3,1) + &
                        k2*cell%reciprocal_cell(1:3,2) + &
                        k3*cell%reciprocal_cell(1:3,3)
                   distance = .norm.k_vector
                   
                   ! - exp(- sigma^2 k^2 / 2) / k^2 2 Re[ S*(k) d S(k) ]
                   ! Re[ z1 z2 ] = x1 x2 - y1 y2
                   enegs(index1) = enegs(index1) - inv_eps_2v * &
                        exp(-gaussian_width*gaussian_width*distance*distance*0.5d0) / &
                        (distance*distance) * 2.d0 * &
                        (s_factor(1,k1,k2,k3)*diff_factor(1,index1,k1,k2,k3) - &
                        s_factor(2,k1,k2,k3)*diff_factor(2,index1,k1,k2,k3))
                   
                end if
             end do
          end do
       end do
       
    end do

#ifdef MPI
    ! collect energies from all cpus in MPI (energy -> tmp_energy)
    call mpi_allreduce(enegs,total_enegs,size(enegs),mpi_double_precision,&
         mpi_sum,mpi_comm_world,mpistat)
#else
    ! enegs: 1 - real space, 2 - self, 3 - reciprocal, 4 - dipole
    total_enegs = enegs
#endif

  end subroutine calculate_ewald_electronegativities
 


end module potentials
