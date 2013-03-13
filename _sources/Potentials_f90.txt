
.. _potentials:
        
===================================================
potentials (Potentials.f90)
===================================================



Potentials contains the low-level routines for handling interactions.
The module defines custom types for both describing the types of
potentials and bond order factors (:data:`potential_descriptor`, :data:`bond_order_descriptor`)
as well as for storing the parameters of actual interactions in use
for the Fortran calculations (:data:`potential`, :data:`bond_order_parameters`).
Tools for creating the custom datatypes (:func:`create_potential`, :func:`create_bond_order_factor`)
are provided.

The types of potentials and bond order factors are defined using the types
:data:`potential_descriptor` and :data:`bond_order_descriptor`.
These should be created at start-up and remain untouched during simulation.
They are used by the Fortran core for checking the types of parameters a potential
needs, for instance, but they are also accessible from the Python interface.
Especially, upon creation of :class:`~pysic.Potential` and :class:`~pysic.BondOrderParameters`
instances, one needs to specify the type as a keyword. This keyword is then compared to the list of
characterizers in the core to determine the type of the interaction.

The basic routines for calculating the actual forces and energies are also defined in
this module (:func:`evaluate_energy`, :func:`evaluate_forces`, :func:`evaluate_bond_order_factor`,
:func:`evaluate_bond_order_gradient`). However, these routines do not calculate the total potential
energy of the system, :math:`V`, or the total forces acting on the particles,
:math:`\mathbf{F} = -\nabla_\alpha V`. Instead, the routines evaluate the contributions from individual
atoms, atom pairs, atom triplets, etc. For instance, let the total energy of the system be

.. math::

  V = \sum_p \left( \sum_i v^p_i + \sum_{(i,j)} v^p_{ij} + \sum_{(i,j,k)} v^p_{ijk} \right),

where :math:`p` sums over the different potentials acting on the system and :math:`i`, :math:`(i,j)` and
:math:`(i,j,k)` sum over all atoms, pairs and triplet, respectively. Then the energy terms :math:`v`
are obtained from :func:`evaluate_energy`. In pseudo-code,

.. math::

  v^p_{S} = \mathtt{evaluate\_energy}(S,p),

where :math:`S` is a set of atoms. The summation over potentials and atoms is done in :ref:`pysic_core`
in :func:`calculate_energy`. Similarly for forces, the summation is carried out in :func:`calculate_forces`.

The reason for separating the calculation of individual interaction terms to :ref:`potentials`
and the overall summation to :ref:`pysic_core` is that only the core knows the current structure and
interactions of the system.
It is the task of this module to tell the core how all the potentials behave given
any local structure, but the overall system information is kept in the core. So during energy
evaluation, :ref:`pysic_core` finds all local structures that possibly contribute with an interaction
and asks :ref:`potentials` to calculate this contribution.

Bond order factors are potential modifiers, not direct interactions themselves.
In general, the factors are scalar functions defined per atom, for instance,

.. math::

   b^p_i = s^p_i\left( \sum_{(i,j)} c^p_{ij} + \sum_{(i,j,k)} c^p_{ijk} \right)

for a three-body factor, where :math:`c^p` are local contributions
(usually representing chemical bonds) and :math:`s^p_i` is a per atom scaling function.
The bond factors multiply the potentials :math:`p`
leading to the total energy

.. math::

   V = \sum_p \left( \sum_i b^p_i v^p_i + \sum_{(i,j)} \frac{1}{2}(b^p_i+b^p_j) v^p_{ij} + \sum_{(i,j,k)} \frac{1}{3}(b^p_i+b^p_j+b^p_k) v^p_{ijk} \right).

The corresponding total force on atom :math:`\alpha` is then

.. math::

   \mathbf{F}_{\alpha} = - \nabla_\alpha V = - \sum_p \left( \sum_i ((\nabla_\alpha b^p_i) v^p_i + b^p_i (\nabla_\alpha v^p_i) ) + \ldots \right).

The contributions :math:`\mathbf{f}^p_\alpha = -\nabla_\alpha v^p`, :math:`c^p`,
and :math:`\nabla_\alpha c^p` are
calculated in :func:`evaluate_forces`, :func:`evaluate_bond_order_factor`,
and :func:`evaluate_bond_order_gradient`.
Application of the scaling functions :math:`s_i` and :math:`s_i'` on the sums
:math:`\sum_{(i,j)} c^p_{ij} + \sum_{(i,j,k)} c^p_{ijk}` is done in the routines
:func:`post_process_bond_order_factor` and :func:`post_process_bond_order_gradient` to
produce the actual bond order factors :math:`b^p_i` and gradients :math:`\nabla_\alpha b^p_i`.
These sums, similarly to the energy and force summations, are evaluated with
:func:`core_calculate_bond_order_factors` in :ref:`pysic_core`.

Note when adding potentials or bond order factors in the source code:

The parameters defined in Potentials.f90 are used for determining the maximum sizes of arrays,
numbers of potentials and bond factors, and the internally used indices
for them. When adding new potentials of bond factors, make sure to update
the relevant numbers. Especially the number of potentials (:data:`n_potential_types`)
or number of bond order factors (:data:`n_bond_order_types`) must be increased
when more types are defined.

Also note that in :ref:`pysic_interface`, some of these parameters are used for
determining array sizes. However, the actual parameters are not used
because f2py does not read the values from here. Therefore if you change
a parameter here, search for its name in :ref:`pysic_interface` to see if the
name appears in a comment. That is an indicator that a numeric value
must be updated accordingly.


.. only:: html


    Modules used by potentials
    --------------------------
    - :ref:`geometry`
    - :ref:`mpi`
    - :ref:`quaternions`
    - :ref:`utility`

    List of global variables in potentials
    --------------------------------------
    - :data:`bond_descriptors_created`
    - :data:`bond_order_descriptors`
    - :data:`c_scale_index`
    - :data:`coordination_index`
    - :data:`descriptors_created`
    - :data:`mono_const_index`
    - :data:`mono_none_index`
    - :data:`mono_qself_index`
    - :data:`n_bond_order_types`
    - :data:`n_max_params`
    - :data:`n_potential_types`
    - :data:`no_name`
    - :data:`pair_buck_index`
    - :data:`pair_exp_index`
    - :data:`pair_lj_index`
    - :data:`pair_power_index`
    - :data:`pair_qabs_index`
    - :data:`pair_qexp_index`
    - :data:`pair_qpair_index`
    - :data:`pair_shift_index`
    - :data:`pair_spring_index`
    - :data:`pair_table_index`
    - :data:`param_name_length`
    - :data:`param_note_length`
    - :data:`pot_name_length`
    - :data:`pot_note_length`
    - :data:`potential_descriptors`
    - :data:`power_index`
    - :data:`quad_dihedral_index`
    - :data:`sqrt_scale_index`
    - :data:`table_bond_index`
    - :data:`table_prefix`
    - :data:`table_scale_index`
    - :data:`table_suffix`
    - :data:`tersoff_index`
    - :data:`tri_bend_index`
    - :data:`triplet_index`

    List of custom types in potentials
    ----------------------------------
    - :data:`bond_order_descriptor`
    - :data:`bond_order_parameters`
    - :data:`potential`
    - :data:`potential_descriptor`

    List of subroutines in potentials
    ---------------------------------
        
    - :func:`bond_order_factor_affects_atom`
    - :func:`bond_order_factor_is_in_group`
    - :func:`calculate_derived_parameters_bond_bending`
    - :func:`calculate_derived_parameters_charge_exp`
    - :func:`calculate_derived_parameters_dihedral`
    - :func:`calculate_ewald_electronegativities`
    - :func:`calculate_ewald_energy`
    - :func:`calculate_ewald_forces`
    - :func:`clear_bond_order_factor_characterizers`
    - :func:`clear_potential_characterizers`
    - :func:`create_bond_order_factor`
    - :func:`create_bond_order_factor_characterizer_coordination`
    - :func:`create_bond_order_factor_characterizer_power`
    - :func:`create_bond_order_factor_characterizer_scaler_1`
    - :func:`create_bond_order_factor_characterizer_scaler_sqrt`
    - :func:`create_bond_order_factor_characterizer_scaler_table`
    - :func:`create_bond_order_factor_characterizer_table`
    - :func:`create_bond_order_factor_characterizer_tersoff`
    - :func:`create_bond_order_factor_characterizer_triplet`
    - :func:`create_potential`
    - :func:`create_potential_characterizer_LJ`
    - :func:`create_potential_characterizer_bond_bending`
    - :func:`create_potential_characterizer_buckingham`
    - :func:`create_potential_characterizer_charge_exp`
    - :func:`create_potential_characterizer_charge_pair`
    - :func:`create_potential_characterizer_charge_pair_abs`
    - :func:`create_potential_characterizer_charge_self`
    - :func:`create_potential_characterizer_constant_force`
    - :func:`create_potential_characterizer_constant_potential`
    - :func:`create_potential_characterizer_dihedral`
    - :func:`create_potential_characterizer_exp`
    - :func:`create_potential_characterizer_power`
    - :func:`create_potential_characterizer_shift`
    - :func:`create_potential_characterizer_spring`
    - :func:`create_potential_characterizer_table`
    - :func:`evaluate_bond_order_factor`
    - :func:`evaluate_bond_order_factor_coordination`
    - :func:`evaluate_bond_order_factor_power`
    - :func:`evaluate_bond_order_factor_table`
    - :func:`evaluate_bond_order_factor_triplet`
    - :func:`evaluate_bond_order_gradient`
    - :func:`evaluate_bond_order_gradient_coordination`
    - :func:`evaluate_bond_order_gradient_power`
    - :func:`evaluate_bond_order_gradient_table`
    - :func:`evaluate_bond_order_gradient_triplet`
    - :func:`evaluate_electronegativity`
    - :func:`evaluate_electronegativity_charge_exp`
    - :func:`evaluate_electronegativity_charge_pair`
    - :func:`evaluate_electronegativity_charge_pair_abs`
    - :func:`evaluate_electronegativity_charge_self`
    - :func:`evaluate_electronegativity_component`
    - :func:`evaluate_energy`
    - :func:`evaluate_energy_LJ`
    - :func:`evaluate_energy_bond_bending`
    - :func:`evaluate_energy_buckingham`
    - :func:`evaluate_energy_charge_exp`
    - :func:`evaluate_energy_charge_pair`
    - :func:`evaluate_energy_charge_pair_abs`
    - :func:`evaluate_energy_charge_self`
    - :func:`evaluate_energy_component`
    - :func:`evaluate_energy_constant_force`
    - :func:`evaluate_energy_constant_potential`
    - :func:`evaluate_energy_dihedral`
    - :func:`evaluate_energy_exp`
    - :func:`evaluate_energy_power`
    - :func:`evaluate_energy_shift`
    - :func:`evaluate_energy_spring`
    - :func:`evaluate_energy_table`
    - :func:`evaluate_force_LJ`
    - :func:`evaluate_force_bond_bending`
    - :func:`evaluate_force_buckingham`
    - :func:`evaluate_force_component`
    - :func:`evaluate_force_constant_force`
    - :func:`evaluate_force_constant_potential`
    - :func:`evaluate_force_dihedral`
    - :func:`evaluate_force_exp`
    - :func:`evaluate_force_power`
    - :func:`evaluate_force_shift`
    - :func:`evaluate_force_spring`
    - :func:`evaluate_force_table`
    - :func:`evaluate_forces`
    - :func:`evaluate_pair_bond_order_factor`
    - :func:`evaluate_pair_bond_order_factor_tersoff`
    - :func:`evaluate_pair_bond_order_gradient`
    - :func:`evaluate_pair_bond_order_gradient_tersoff`
    - :func:`get_bond_descriptor`
    - :func:`get_description_of_bond_order_factor`
    - :func:`get_description_of_potential`
    - :func:`get_descriptions_of_parameters_of_bond_order_factor`
    - :func:`get_descriptions_of_parameters_of_potential`
    - :func:`get_descriptor`
    - :func:`get_index_of_bond_order_factor`
    - :func:`get_index_of_parameter_of_bond_order_factor`
    - :func:`get_index_of_parameter_of_potential`
    - :func:`get_index_of_potential`
    - :func:`get_names_of_parameters_of_bond_order_factor`
    - :func:`get_names_of_parameters_of_potential`
    - :func:`get_number_of_bond_order_factors`
    - :func:`get_number_of_parameters_of_bond_order_factor`
    - :func:`get_number_of_parameters_of_potential`
    - :func:`get_number_of_potentials`
    - :func:`get_number_of_targets_of_bond_order_factor`
    - :func:`get_number_of_targets_of_bond_order_factor_index`
    - :func:`get_number_of_targets_of_potential`
    - :func:`get_number_of_targets_of_potential_index`
    - :func:`initialize_bond_order_factor_characterizers`
    - :func:`initialize_potential_characterizers`
    - :func:`is_valid_bond_order_factor`
    - :func:`is_valid_potential`
    - :func:`list_bond_order_factors`
    - :func:`list_potentials`
    - :func:`post_process_bond_order_factor`
    - :func:`post_process_bond_order_factor_scaler_1`
    - :func:`post_process_bond_order_factor_scaler_sqrt`
    - :func:`post_process_bond_order_factor_scaler_table`
    - :func:`post_process_bond_order_factor_tersoff`
    - :func:`post_process_bond_order_gradient`
    - :func:`post_process_bond_order_gradient_scaler_1`
    - :func:`post_process_bond_order_gradient_scaler_sqrt`
    - :func:`post_process_bond_order_gradient_scaler_table`
    - :func:`post_process_bond_order_gradient_tersoff`
    - :func:`potential_affects_atom`
    - :func:`smoothening_derivative`
    - :func:`smoothening_factor`
    - :func:`smoothening_gradient`


Full documentation of global variables in potentials
----------------------------------------------------
        
        
  .. data:: bond_descriptors_created

    logical    *scalar*    

    *initial value* = .false.
    
    logical tag used for managing pointer allocations for bond order factor descriptors
    
  .. data:: bond_order_descriptors

    type(bond_order_descriptor)  *pointer*  *size(:)*    
    
    an array for storing descriptors for the different *types* of bond order factors
    
  .. data:: c_scale_index

    integer    *scalar*  *parameter*  

    *initial value* = 3
    
    internal index for the coordination scaling function
    
  .. data:: coordination_index

    integer    *scalar*  *parameter*  

    *initial value* = 1
    
    
    
  .. data:: descriptors_created

    logical    *scalar*    

    *initial value* = .false.
    
    logical tag used for managing pointer allocations for potential descriptors
    
  .. data:: mono_const_index

    integer    *scalar*  *parameter*  

    *initial value* = 3
    
    internal index for the constant force potential
    
  .. data:: mono_none_index

    integer    *scalar*  *parameter*  

    *initial value* = 6
    
    internal index for the constant potential
    
  .. data:: mono_qself_index

    integer    *scalar*  *parameter*  

    *initial value* = 11
    
    
    
  .. data:: n_bond_order_types

    integer    *scalar*  *parameter*  

    *initial value* = 8
    
    number of different types of bond order factors known
    
  .. data:: n_max_params

    integer    *scalar*  *parameter*  

    *initial value* = 12
    
    
    
  .. data:: n_potential_types

    integer    *scalar*  *parameter*  

    *initial value* = 15
    
    number of different types of potentials known
    
  .. data:: no_name

    character(len=label_length)    *scalar*  *parameter*  

    *initial value* = "xx"
    
    The label for unlabeled atoms. In other words, there are routines that expect atomic symbols as arguments, but if there are no symbols to pass, this should be given to mark an empty entry.
    
  .. data:: pair_buck_index

    integer    *scalar*  *parameter*  

    *initial value* = 7
    
    internal index for the Buckingham potential
    
  .. data:: pair_exp_index

    integer    *scalar*  *parameter*  

    *initial value* = 5
    
    
    
  .. data:: pair_lj_index

    integer    *scalar*  *parameter*  

    *initial value* = 1
    
    internal index for the Lennard-Jones potential
    
  .. data:: pair_power_index

    integer    *scalar*  *parameter*  

    *initial value* = 9
    
    internal index for the power law potential
    
  .. data:: pair_qabs_index

    integer    *scalar*  *parameter*  

    *initial value* = 14
    
    
    
  .. data:: pair_qexp_index

    integer    *scalar*  *parameter*  

    *initial value* = 13
    
    
    
  .. data:: pair_qpair_index

    integer    *scalar*  *parameter*  

    *initial value* = 12
    
    
    
  .. data:: pair_shift_index

    integer    *scalar*  *parameter*  

    *initial value* = 15
    
    
    
  .. data:: pair_spring_index

    integer    *scalar*  *parameter*  

    *initial value* = 2
    
    internal index for the spring potential
    
  .. data:: pair_table_index

    integer    *scalar*  *parameter*  

    *initial value* = 10
    
    internal index for the tabulated potential
    
  .. data:: param_name_length

    integer    *scalar*  *parameter*  

    *initial value* = 10
    
    
    
  .. data:: param_note_length

    integer    *scalar*  *parameter*  

    *initial value* = 100
    
    maximum length allowed for the descriptions of parameters
    
  .. data:: pot_name_length

    integer    *scalar*  *parameter*  

    *initial value* = 11
    
    maximum length allowed for the names of potentials
    
  .. data:: pot_note_length

    integer    *scalar*  *parameter*  

    *initial value* = 500
    
    maximum lenght allowed for the description of the potential
    
  .. data:: potential_descriptors

    type(potential_descriptor)  *pointer*  *size(:)*    
    
    an array for storing descriptors for the different *types* of potentials
    
  .. data:: power_index

    integer    *scalar*  *parameter*  

    *initial value* = 5
    
    internal index for the power law bond order factor
    
  .. data:: quad_dihedral_index

    integer    *scalar*  *parameter*  

    *initial value* = 8
    
    internal index for the dihedral angle potential
    
  .. data:: sqrt_scale_index

    integer    *scalar*  *parameter*  

    *initial value* = 6
    
    internal index for the square root scaling function
    
  .. data:: table_bond_index

    integer    *scalar*  *parameter*  

    *initial value* = 7
    
    internal index for the tabulated bond order factor
    
  .. data:: table_prefix

    character(len=6)    *scalar*  *parameter*  

    *initial value* = "table_"
    
    prefix for filenames for storing tables
    
  .. data:: table_scale_index

    integer    *scalar*  *parameter*  

    *initial value* = 8
    
    internal index for the tabulated scaling function
    
  .. data:: table_suffix

    character(len=4)    *scalar*  *parameter*  

    *initial value* = ".txt"
    
    
    
  .. data:: tersoff_index

    integer    *scalar*  *parameter*  

    *initial value* = 2
    
    internal index for the Tersoff bond order factor
    
  .. data:: tri_bend_index

    integer    *scalar*  *parameter*  

    *initial value* = 4
    
    internal index for the bond bending potential
    
  .. data:: triplet_index

    integer    *scalar*  *parameter*  

    *initial value* = 4
    
    internal index for the triplet bond order factor
    

Full documentation of custom types in potentials
------------------------------------------------
        
        
  .. data:: bond_order_descriptor

    Description of a type of a bond order factor.
    The type contains the name and description of the bond order factor
    and the parameters it contains.
    The descriptors contain the information that the inquiry methods in
    the python interface fetch.
    

    Contained data:

    parameter_notes: character(len=param_note_length)  *pointer*  *size(:, :)*
        Descriptions of the parameters. The descriptions should be very short indicators such as 'spring constant' or 'energy coefficient'. For more detailed explanations, the proper documentation should be used.
    n_parameters: integer  *pointer*  *size(:)*
        number of parameters for each number of bodies (1-body parameters, 2-body parameters etc.)
    n_level: integer    *scalar*
        
    name: character(len=pot_name_length)    *scalar*
        The name of the bond order factor: this is a keyword according to which the factor may be recognized.
    description: character(len=pot_note_length)    *scalar*
        A description of the bond order factor. This should contain the mathematical formulation as well as a short verbal explanation.
    includes_post_processing: logical    *scalar*
        a logical tag specifying if there is a scaling function :math:`s_i` attached to the factor.
    type_index: integer    *scalar*
        The internal index of the bond order factor. This can also be used for recognizing the factor and must therefore match the name. For instance, if name = 'neighbors', type_index = :data:`coordination_index`.
    n_targets: integer    *scalar*
        number of targets, i.e., interacting bodies
    parameter_names: character(len=param_name_length)  *pointer*  *size(:, :)*
        The names of the parameters of the bond order factor: these are keywords according to which the parameters may be recognized.
  .. data:: bond_order_parameters

    Defines a particular bond order factor.
    The factor should correspond to the description of some
    built-in type and hold actual numeric values for parameters.
    In addition a real bond order factor must have information on the
    particles it acts on and the range it operates in.
    These are created based on the :class:`~pysic.BondOrderParameters` objects in the Python
    interface when calculations are invoked.
    

    Contained data:

    cutoff: double precision    *scalar*
        The hard cutoff for the bond order factor. If the atoms are farther away from each other than this, they do not contribute to the total bond order factor does not affect them.
    n_level: integer    *scalar*
        
    soft_cutoff: double precision    *scalar*
        The soft cutoff for the bond order factor. If this is smaller than the hard cutoff, the bond contribution is scaled to zero continuously when the interatomic distances grow from the soft to the hard cutoff.
    parameters: double precision  *pointer*  *size(:, :)*
        numerical values for parameters
    group_index: integer    *scalar*
        The internal index of the *potential* the bond order factor is modifying.
    includes_post_processing: logical    *scalar*
        a logical switch specifying if there is a scaling function :math:`s_i` attached to the factor
    type_index: integer    *scalar*
        The internal index of the bond order factor *type*. This is used for recognizing the factor. Note that the bond order parameters instance does not have a name. If the name is needed, it can be obtained from the :data:`bond_order_descriptor` of the correct index.
    n_params: integer  *pointer*  *size(:)*
        array containing the numbers of parameters for different number of targets (1-body parameters, 2-body parameters, etc.)
    table: double precision  *pointer*  *size(:, :)*
        array for storing tabulated values
    original_elements: character(len=2)  *pointer*  *size(:)*
        The list of elements (atomic symbols) of the original :class:`~pysic.BondOrderParameters` in the Python interface from which this factor was created. Whereas the apply_elements lists are used for finding all pairs and triplets of atoms which could contribute to the bond order factor, the original_elements lists specify the roles of atoms in the factor.
    derived_parameters: double precision  *pointer*  *size(:, :)*
        numerical values for parameters calculated based on the parameters specified by the user
    apply_elements: character(len=2)  *pointer*  *size(:)*
        A list of elements (atomic symbols) the factor affects. E.g., for Si-O bonds, it would be ('Si','O'). Note that unlike in the Python interface, a single :data:`bond_order_parameters` only has one set of targets, and for multiple target options several :data:`bond_order_parameters` instances are created.
  .. data:: potential

    Defines a particular potential.
    The potential should correspond to the description of some
    built-in type and hold actual numeric values for parameters.
    In addition, a real potential must have information on the
    particles it acts on and the range it operates in.
    These are to be created based on the :class:`~pysic.Potential` objects in the Python
    interface when calculations are invoked.
    

    Contained data:

    pot_index: integer    *scalar*
        The internal index of the *actual potential*. This is needed when bond order factors are included so that the factors may be joint with the correct potentials.
    smoothened: logical    *scalar*
        logical switch specifying if a smooth cutoff is applied to the potential
    n_product: integer    *scalar*
        number of multipliers for a product potential
    filter_elements: logical    *scalar*
        a logical switch specifying whether the potential targets atoms based on the atomic symbols
    parameters: double precision  *pointer*  *size(:)*
        numerical values for parameters
    cutoff: double precision    *scalar*
        The hard cutoff for the potential. If the atoms are farther away from each other than this, the potential does not affect them.
    soft_cutoff: double precision    *scalar*
        The soft cutoff for the potential. If this is smaller than the hard cutoff, the potential is scaled to zero continuously when the interatomic distances grow from the soft to the hard cutoff.
    apply_tags: integer  *pointer*  *size(:)*
        A list of atom tags the potential affects. Note that unlike in the Python interface, a single :data:`potential` only has one set of targets, and for multiple target options several :data:`potential` instances are created.
    original_indices: integer  *pointer*  *size(:)*
        The list of atom indices of the original :class:`~pysic.Potential` in the Python interface from which this potential was created. Whereas the apply_* lists are used for finding all pairs and triplets of atoms for which the potential could act on, the original_* lists specify the roles of atoms in the interaction.
    apply_indices: integer  *pointer*  *size(:)*
        A list of atom indices the potential affects. Note that unlike in the Python interface, a single :data:`potential` only has one set of targets, and for multiple target options several :data:`potential` instances are created.
    filter_indices: logical    *scalar*
        a logical switch specifying whether the potential targets atoms based on the atom indices
    type_index: integer    *scalar*
        The internal index of the potential *type*. This is used for recognizing the potential. Note that the potential instance does not have a name. If the name is needed, it can be obtained from the :data:`potential_descriptor` of the correct index.
    original_tags: integer  *pointer*  *size(:)*
        The list of atom tags of the original :class:`~pysic.Potential` in the Python interface from which this potential was created. Whereas the apply_* lists are used for finding all pairs and triplets of atoms for which the potential could act on, the original_* lists specify the roles of atoms in the interaction.
    table: double precision  *pointer*  *size(:, :)*
        array for storing tabulated values
    original_elements: character(len=2)  *pointer*  *size(:)*
        The list of elements (atomic symbols) of the original :class:`~pysic.Potential` in the Python interface from which this potential was created. Whereas the apply_* lists are used for finding all pairs and triplets of atoms for which the potential could act on, the original_* lists specify the roles of atoms in the interaction.
    multipliers: type(potential)  *pointer*  *size(:)*
        additional potentials with the same targets and cutoff, for potential multiplication
    derived_parameters: double precision  *pointer*  *size(:)*
        numerical values for parameters calculated based on the parameters specified by the user
    apply_elements: character(len=2)  *pointer*  *size(:)*
        A list of elements (atomic symbols) the potential affects. E.g., for a Si-O potential, it would be ('Si','O'). Note that unlike in the Python interface, a single :data:`potential` only has one set of targets, and for multiple target options several :data:`potential` instances are created.
    filter_tags: logical    *scalar*
        a logical switch specifying whether the potential targets atoms based on the atom tags
  .. data:: potential_descriptor

    Description of a type of a potential.
    The type contains the name and description of the potential
    and the parameters it contains.
    The descriptors contain the information that the inquiry methods in
    the python interface fetch.
    

    Contained data:

    parameter_notes: character(len=param_note_length)  *pointer*  *size(:)*
        Descriptions of the parameters. The descriptions should be very short indicators such as 'spring constant' or 'energy coefficient'. For more detailed explanations, the proper documentation should be used.
    n_parameters: integer    *scalar*
        number of parameters
    description: character(len=pot_note_length)    *scalar*
        A description of the potential. This should contain the mathematical formulation as well as a short verbal explanation.
    type_index: integer    *scalar*
        The internal index of the potential. This can also be used for recognizing the potential and must therefore match the name. For instance, if name = 'LJ', type_index = :data:`pair_lj_index`.
    n_targets: integer    *scalar*
        number of targets, i.e., interacting bodies
    parameter_names: character(len=param_name_length)  *pointer*  *size(:)*
        The names of the parameters of the potential: these are keywords according to which the parameters may be recognized.
    name: character(len=pot_name_length)    *scalar*
        The name of the potential: this is a keyword according to which the potentials may be recognized.

Full documentation of subroutines in potentials
-----------------------------------------------
        
        
            
  .. function:: bond_order_factor_affects_atom(factor, atom_in, affects, position)

    Tests whether the given bond order factor affects the specific atom.
    
    For bond order factors, the atoms are specified as valid targets by
    the atomic symbol only.
    
    If position is not given, then the routine returns true if
    the atom can appear in the bond order factor in any role.
    If position is given, then true is returned only if the atom
    is valid for that particular position.
    
    For instance, we may want to calculate the coordination of
    Cu-O bonds for Cu but not for O.
    

    Parameters:

    factor: type(bond_order_parameters)  *intent(in)*    *scalar*  
        the :data:`bond_order_parameters`
    atom_in: type(atom)  *intent(in)*    *scalar*  
        the :data:`atom`
    **affects**: logical  **intent(out)**    *scalar*  
        true if the bond order factor is affected by the atom
    position: integer  *intent(in)*    *scalar*  *optional*
        specifies the particular role of the atom in the bond order factor
            
  .. function:: bond_order_factor_is_in_group(factor, group_index, in_group)

    Tests whether the given bond order factor is a member of a specific group,
    i.e., if it affects the potential specifiesd by the group index.
    

    Parameters:

    factor: type(bond_order_parameters)  *intent(in)*    *scalar*  
        the :data:`bond_order_parameters`
    group_index: integer  *intent(in)*    *scalar*  
        the index for the potential
    **in_group**: logical  **intent(out)**    *scalar*  
        true if the factor is a member of the group
            
  .. function:: calculate_derived_parameters_bond_bending(n_params, parameters, new_potential)

    Bond bending derived parameters
    

    Parameters:

    n_params: integer  *intent(in)*    *scalar*  
        
    parameters: double precision  *intent(in)*    *size(n_params)*  
        
    **new_potential**: type(potential)  **intent(inout)**    *scalar*  
        the potential object for which the parameters are calculated
            
  .. function:: calculate_derived_parameters_charge_exp(n_params, parameters, new_potential)

    Charge exponential derived parameters
    

    Parameters:

    n_params: integer  *intent(in)*    *scalar*  
        
    parameters: double precision  *intent(in)*    *size(n_params)*  
        
    **new_potential**: type(potential)  **intent(inout)**    *scalar*  
        the potential object for which the parameters are calculated
            
  .. function:: calculate_derived_parameters_dihedral(n_params, parameters, new_potential)

    Dihedral angle derived parameters
    

    Parameters:

    n_params: integer  *intent(in)*    *scalar*  
        
    parameters: double precision  *intent(in)*    *size(n_params)*  
        
    **new_potential**: type(potential)  **intent(inout)**    *scalar*  
        the potential object for which the parameters are calculated
            
  .. function:: calculate_ewald_electronegativities(n_atoms, atoms, cell, real_cutoff, reciprocal_cutoff, gaussian_width, electric_constant, filter, scaler, include_dipole_correction, total_enegs)

    Calculates the electronegativities due to long ranged :math:`\frac{1}{r}` potentials.
    These electronegativities are the derivatives of the energies :math:`U` given by :func:`calculate_ewald_energy`
    
    .. math::
    
       \chi_\alpha = - \frac{\partial U}{\partial q_\alpha}
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    atoms: type(atom)  *intent(in)*    *size(n_atoms)*  
        list of atoms
    cell: type(supercell)  *intent(in)*    *scalar*  
        the supercell containing the system
    real_cutoff: double precision  *intent(in)*    *scalar*  
        Cutoff radius of real-space interactions. Note that the neighbor lists stored in the atoms are used for neighbor finding so the cutoff cannot exceed the cutoff for the neighbor lists. (Or, it can, but the neighbors not in the lists will not be found.)
    reciprocal_cutoff: integer  *intent(in)*    *size(3)*  
        The number of cells to be included in the reciprocal sum in the directions of the reciprocal cell vectors. For example, if ``reciprocal_cutoff = [3,4,5]``, the reciprocal sum will be truncated as :math:`\sum_{\mathbf{k} \ne 0} = \sum_{k_1=-3}^3 \sum_{k_2=-4}^4 \sum_{k_3 = -5,(k_1,k_2,k_3) \ne (0,0,0)}^5`.
    gaussian_width: double precision  *intent(in)*    *scalar*  
        The :math:`\sigma` parameter, i.e., the distribution width of the screening Gaussians. This should not influence the actual value of the energy, but it does influence the convergence of the summation. If :math:`\sigma` is large, the real space sum :math:`E_s` converges slowly and a large real space cutoff is needed. If it is small, the reciprocal term :math:`E_l` converges slowly and the sum over the reciprocal lattice has to be evaluated over several cell lengths.
    electric_constant: double precision  *intent(in)*    *scalar*  
        The electic constant, i.e., vacuum permittivity :math:`\varepsilon_0`. In atomic units, it is :math:`\varepsilon_0 = 0.00552635 \frac{e^2}{\mathrm{Å\ eV}}`, but if one wishes to scale the results to some other unit system (such as reduced units with :math:`\varepsilon_0 = 1`), that is possible as well.
    filter: logical  *intent(in)*    *size(n_atoms)*  
        a list of logical values, one per atom, false for the atoms that should be ignored in the calculation
    scaler: double precision  *intent(in)*    *size(n_atoms)*  
        a list of numerical values to scale the individual charges of the atoms
    include_dipole_correction: logical  *intent(in)*    *scalar*  
        if true, a dipole correction term is included
    **total_enegs**: double precision  **intent(out)**    *size(n_atoms)*  
        the calculated electronegativities
            
  .. function:: calculate_ewald_energy(n_atoms, atoms, cell, real_cutoff, reciprocal_cutoff, gaussian_width, electric_constant, filter, scaler, include_dipole_correction, total_energy)

    Calculates the energy of :math:`\frac{1}{r}` potentials through Ewald summation.
    
    If a periodic system contains charges interacting via the :math:`\frac{1}{r}` Coulomb potential,
    direct summation of the interactions
    
    .. math::
       :label: direct_sum
    
       E = \sum_{(i,j)} \frac{1}{4\pi\epsilon_0}\frac{q_i q_j}{r_{ij}},
    
    where the sum is over pairs of charges :math:`q_i, q_j`
    (charges of the entire system, not just the simulation cell) and the distance between the charges is
    :math:`r_{ij} = |\mathbf{r}_j - \mathbf{r}_i|`,
    does not work in general because the sum :eq:`direct_sum` converges very slowly. [#]_ Therefore truncating the
    sum may lead to severe errors.
    
    The standard technique for overcoming this problem is the so called Ewald summation method.
    The idea is to split the long ranged and singular Coulomb potential to a short ranged singular and
    long ranged smooth parts, and calculate the long ranged part in reciprocal space via Fourier transformations.
    This is possible since the system is periodic and the same supercell repeats infinitely in all directions.
    In practice the calculation can be done by adding (and subtracting) Gaussian charge densities over the
    point charges to screen the
    potential in real space. That is, the original charge density
    :math:`\rho(\mathbf{r}) = \sum_i q_i \delta(\mathbf{r} - \mathbf{r}_i)` is split by
    
    .. math::
      :nowrap:
    
      \begin{eqnarray}
      \rho(\mathbf{r}) & = & \rho_s(\mathbf{r}) + \rho_l(\mathbf{r}) \\
      \rho_s(\mathbf{r}) & = & \sum_i \left[ q_i \delta(\mathbf{r} - \mathbf{r}_i) - q_i G_\sigma(\mathbf{r} - \mathbf{r}_i) \right] \\
      \rho_l(\mathbf{r}) & = & \sum_i q_i G_\sigma(\mathbf{r} - \mathbf{r}_i) \\
      G_\sigma(\mathbf{r}) & = & \frac{1}{(2 \pi \sigma^2)^{3/2}} \exp\left( -\frac{|\mathbf{r}|^2}{2 \sigma^2} \right)
      \end{eqnarray}
    
    Here :math:`\rho_l` generates a long range interaction since at large distances the Gaussian densities
    :math:`G_\sigma` appear the same as point charges
    (:math:`\lim_{\sigma/r \to 0} G_\sigma(\mathbf{r}) = \delta(\mathbf{r})`).
    Since the charge density is smooth, so will be the potential it creates.
    The density :math:`\rho_s` exhibits short ranged interactions for the same reason:
    At distances longer than the width of the
    Gaussians the point charges are screened by the Gaussians which exactly cancel them
    (:math:`\lim_{\sigma/r \to 0} \delta(\mathbf{r}) - G_\sigma(\mathbf{r}) = 0`).
    
    The short ranged interactions are directly calculated in real space
    
    .. math::
       :nowrap:
    
       \begin{eqnarray}
       E_s & = & \frac{1}{4 \pi \varepsilon_0} \int \frac{\rho_s(\mathbf{r}) \rho_s(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|} \mathrm{d}^3 r \mathrm{d}^3 r' \\
           & = & \frac{1}{4 \pi \varepsilon_0} \sum_{(i,j)} \frac{q_i q_j}{r_{ij}} \mathrm{erfc} \left( \frac{r_{ij}}{\sigma \sqrt{2}} \right).
       \end{eqnarray}
    
    The complementary error function :math:`\mathrm{erfc}(r) = 1 - \mathrm{erf}(r) = 1 - \frac{2}{\sqrt{\pi}} \int_0^r e^{-t^2/2} \mathrm{d}t` makes the sum converge rapidly as :math:`\frac{r_{ij}}{\sigma} \to \infty`.
    
    The long ranged interaction can be calculated in reciprocal space by Fourier transformation. The result is
    
    .. math::
       :nowrap:
    
       \begin{eqnarray}
       E_l & = & \frac{1}{2 V \varepsilon_0} \sum_{\mathbf{k} \ne 0} \frac{e^{-\sigma^2 k^2 / 2}}{k^2} |S(\mathbf{k})|^2 - \frac{1}{4 \pi \varepsilon_0} \frac{1}{\sqrt{2 \pi} \sigma} \sum_i^N q_i^2\\
       S(\mathbf{k}) & = & \sum_i^N q_i e^{\mathrm{i} \mathbf{k} \cdot \mathbf{r}_i}
       \end{eqnarray}
    
    The first sum in :math:`E_l` runs over the reciprocal lattice
    :math:`\mathbf{k} = k_1 \mathbf{b}_1 + k_2 \mathbf{b}_2 + k_3 \mathbf{b}_3` where :math:`\mathbf{b}_i`
    are the vectors spanning the reciprocal cell (:math:`[\mathbf{b}_1 \mathbf{b}_2 \mathbf{b}_3] = ([\mathbf{v}_1 \mathbf{v}_2 \mathbf{v}_3]^{-1})^T` where :math:`\mathbf{v}_i` are the real space cell vectors).
    The latter sum is the self energy of each point charge in the potential of the particular Gaussian that
    screens the charge, and the sum runs
    over all charges in the supercell spanning the periodic system.
    (The self energy must be removed because it is present in the first sum even though when evaluating
    the potential at the position of a charge
    due to the other charges, no screening Gaussian function should be placed over the charge itself.)
    Likewise the sum in the structure factor :math:`S(\mathbf{k})` runs over all charges in the supercell.
    
    The total energy is then the sum of the short and long range energies
    
    .. math::
    
       E = E_s + E_l.
    
    .. [#] In fact, the sum converges only conditionally meaning the result depends on the order of summation. Physically this is not a problem, because one never has infinite lattices.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    atoms: type(atom)  *intent(in)*    *size(n_atoms)*  
        list of atoms
    cell: type(supercell)  *intent(in)*    *scalar*  
        the supercell containing the system
    real_cutoff: double precision  *intent(in)*    *scalar*  
        Cutoff radius of real-space interactions. Note that the neighbor lists stored in the atoms are used for neighbor finding so the cutoff cannot exceed the cutoff for the neighbor lists. (Or, it can, but the neighbors not in the lists will not be found.)
    reciprocal_cutoff: integer  *intent(in)*    *size(3)*  
        The number of cells to be included in the reciprocal sum in the directions of the reciprocal cell vectors. For example, if ``reciprocal_cutoff = [3,4,5]``, the reciprocal sum will be truncated as :math:`\sum_{\mathbf{k} \ne 0} = \sum_{k_1=-3}^3 \sum_{k_2=-4}^4 \sum_{k_3 = -5,(k_1,k_2,k_3) \ne (0,0,0)}^5`.
    gaussian_width: double precision  *intent(in)*    *scalar*  
        The :math:`\sigma` parameter, i.e., the distribution width of the screening Gaussians. This should not influence the actual value of the energy, but it does influence the convergence of the summation. If :math:`\sigma` is large, the real space sum :math:`E_s` converges slowly and a large real space cutoff is needed. If it is small, the reciprocal term :math:`E_l` converges slowly and the sum over the reciprocal lattice has to be evaluated over several cell lengths.
    electric_constant: double precision  *intent(in)*    *scalar*  
        The electic constant, i.e., vacuum permittivity :math:`\varepsilon_0`. In atomic units, it is :math:`\varepsilon_0 = 0.00552635 \frac{e^2}{\mathrm{Å\ eV}}`, but if one wishes to scale the results to some other unit system (such as reduced units with :math:`\varepsilon_0 = 1`), that is possible as well.
    filter: logical  *intent(in)*    *size(n_atoms)*  
        a list of logical values, one per atom, false for the atoms that should be ignored in the calculation
    scaler: double precision  *intent(in)*    *size(n_atoms)*  
        a list of numerical values to scale the individual charges of the atoms
    include_dipole_correction: logical  *intent(in)*    *scalar*  
        if true, a dipole correction term is included in the energy
    **total_energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy
            
  .. function:: calculate_ewald_forces(n_atoms, atoms, cell, real_cutoff, reciprocal_cutoff, gaussian_width, electric_constant, filter, scaler, include_dipole_correction, total_forces, total_stress)

    Calculates the forces due to long ranged :math:`\frac{1}{r}` potentials.
    These forces are the gradients of the energies :math:`U` given by :func:`calculate_ewald_energy`
    
    .. math::
    
       \mathbf{F}_\alpha = - \nabla_\alpha U
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    atoms: type(atom)  *intent(in)*    *size(n_atoms)*  
        list of atoms
    cell: type(supercell)  *intent(in)*    *scalar*  
        the supercell containing the system
    real_cutoff: double precision  *intent(in)*    *scalar*  
        Cutoff radius of real-space interactions. Note that the neighbor lists stored in the atoms are used for neighbor finding so the cutoff cannot exceed the cutoff for the neighbor lists. (Or, it can, but the neighbors not in the lists will not be found.)
    reciprocal_cutoff: integer  *intent(in)*    *size(3)*  
        The number of cells to be included in the reciprocal sum in the directions of the reciprocal cell vectors. For example, if ``reciprocal_cutoff = [3,4,5]``, the reciprocal sum will be truncated as :math:`\sum_{\mathbf{k} \ne 0} = \sum_{k_1=-3}^3 \sum_{k_2=-4}^4 \sum_{k_3 = -5,(k_1,k_2,k_3) \ne (0,0,0)}^5`.
    gaussian_width: double precision  *intent(in)*    *scalar*  
        The :math:`\sigma` parameter, i.e., the distribution width of the screening Gaussians. This should not influence the actual value of the energy, but it does influence the convergence of the summation. If :math:`\sigma` is large, the real space sum :math:`E_s` converges slowly and a large real space cutoff is needed. If it is small, the reciprocal term :math:`E_l` converges slowly and the sum over the reciprocal lattice has to be evaluated over several cell lengths.
    electric_constant: double precision  *intent(in)*    *scalar*  
        The electic constant, i.e., vacuum permittivity :math:`\varepsilon_0`. In atomic units, it is :math:`\varepsilon_0 = 0.00552635 \frac{e^2}{\mathrm{Å\ eV}}`, but if one wishes to scale the results to some other unit system (such as reduced units with :math:`\varepsilon_0 = 1`), that is possible as well.
    filter: logical  *intent(in)*    *size(n_atoms)*  
        a list of logical values, one per atom, false for the atoms that should be ignored in the calculation
    scaler: double precision  *intent(in)*    *size(n_atoms)*  
        a list of numerical values to scale the individual charges of the atoms
    include_dipole_correction: logical  *intent(in)*    *scalar*  
        if true, a dipole correction term is included in the energy
    **total_forces**: double precision  **intent(out)**    *size(3, n_atoms)*  
        the calculated forces
    **total_stress**: double precision  **intent(out)**    *size(6)*  
        the calculated stress
            
  .. function:: clear_bond_order_factor_characterizers()

    Deallocates all memory associated with bond order factor characterizes.

            
  .. function:: clear_potential_characterizers()

    Deallocates all memory associated with potential characterizes.

            
  .. function:: create_bond_order_factor(n_targets, n_params, n_split, bond_name, parameters, param_split, cutoff, soft_cutoff, elements, orig_elements, group_index, new_bond, success)

    Returns a :data:`bond_order_parameters`.
    
    The routine takes as arguments all the necessary parameters
    and returns a bond order parameters type wrapping them in one package.
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets, i.e., interacting bodies
    n_params: integer  *intent(in)*    *scalar*  
        array containing the numbers of parameters for different number of targets (1-body parameters, 2-body parameters, etc.)
    n_split: integer  *intent(in)*    *scalar*  
        number of groupings in the list of parameters, per number of bodies - should equal n_targets
    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor - a keyword that must match a name of one of the :data:`bond_order_descriptors`
    parameters: double precision  *intent(in)*    *size(n_params)*  
        numerical values for parameters as a one-dimensional array
    param_split: integer  *intent(in)*    *size(n_split)*  
        Array containing the numbers of 1-body, 2-body, etc. parameters. The parameters are given as a list, but a bond order factor may have parameters separately for different numbers of targets. This list specifies the number of parameters for each.
    cutoff: double precision  *intent(in)*    *scalar*  
        The hard cutoff for the bond order factor. If the atoms are farther away from each other than this, they do not contribute to the total bond order factor does not affect them.
    soft_cutoff: double precision  *intent(in)*    *scalar*  
        The soft cutoff for the bond order factor. If this is smaller than the hard cutoff, the bond contribution is scaled to zero continuously when the interatomic distances grow from the soft to the hard cutoff.
    elements: character(len=2)  *intent(in)*    *size(n_targets)*  
        a list of elements (atomic symbols) the factor affects
    orig_elements: character(len=2)  *intent(in)*    *size(n_targets)*  
        the list of elements (atomic symbols) of the original :class:`~pysic.BondOrderParameters` in the Python interface from which this factor was created
    group_index: integer  *intent(in)*    *scalar*  
        The internal index of the *potential* the bond order factor is modifying.
    **new_bond**: type(bond_order_parameters)  **intent(out)**    *scalar*  
        the created :data:`bond_order_parameters`
    **success**: logical  **intent(out)**    *scalar*  
        logical tag specifying if creation of the factor succeeded
            
  .. function:: create_bond_order_factor_characterizer_coordination(index)

    Coordination characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the bond order factor
            
  .. function:: create_bond_order_factor_characterizer_power(index)

    Power decay characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the bond order factor
            
  .. function:: create_bond_order_factor_characterizer_scaler_1(index)

    Scaler characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the bond order factor
            
  .. function:: create_bond_order_factor_characterizer_scaler_sqrt(index)

    Square root scaler characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the bond order factor
            
  .. function:: create_bond_order_factor_characterizer_scaler_table(index)

    Square root scaler characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the bond order factor
            
  .. function:: create_bond_order_factor_characterizer_table(index)

    Tabulated characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the bond order factor
            
  .. function:: create_bond_order_factor_characterizer_tersoff(index)

    Tersoff characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the bond order factor
            
  .. function:: create_bond_order_factor_characterizer_triplet(index)

    Triplet characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the bond order factor
            
  .. function:: create_potential(n_targets, n_params, pot_name, parameters, cutoff, soft_cutoff, elements, tags, indices, orig_elements, orig_tags, orig_indices, pot_index, n_multi, multipliers, new_potential, success)

    Returns a :data:`potential`.
    
    The routine takes as arguments all the necessary parameters
    and returns a potential type wrapping them in one package.
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets, i.e., interacting bodies
    n_params: integer  *intent(in)*    *scalar*  
        number of parameters
    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential - a keyword that must match a name of one of the :data:`potential_descriptors`
    parameters: double precision  *intent(in)*    *size(n_params)*  
        array of numerical values for the parameters
    cutoff: double precision  *intent(in)*    *scalar*  
        the hard cutoff for the potential
    soft_cutoff: double precision  *intent(in)*    *scalar*  
        the soft cutoff for the potential
    elements: character(len=2)  *intent(in)*    *size(n_targets)*  
        the elements (atomic symbols) the potential acts on
    tags: integer  *intent(in)*    *size(n_targets)*  
        the atom tags the potential acts on
    indices: integer  *intent(in)*    *size(n_targets)*  
        the atom indices the potential acts on
    orig_elements: character(len=2)  *intent(in)*    *size(n_targets)*  
        The elements (atomic symbols) in the :class:`~pysic.Potential` used for generating the potential. This is needed to specify the roles of the atoms in the interaction.
    orig_tags: integer  *intent(in)*    *size(n_targets)*  
        The atom tags in the :class:`~pysic.Potential` used for generating the potential. This is needed to specify the roles of the atoms in the interaction.
    orig_indices: integer  *intent(in)*    *size(n_targets)*  
        The atom indices in the :class:`~pysic.Potential` used for generating the potential. This is needed to specify the roles of the atoms in the interaction.
    pot_index: integer  *intent(in)*    *scalar*  
        the internal index of the potential
    n_multi: integer  *intent(in)*    *scalar*  
        number of multiplying potentials
    multipliers: type(potential)  *intent(in)*    *size(n_multi)*  
        the multiplying potentials
    **new_potential**: type(potential)  **intent(out)**    *scalar*  
        the created :data:`potential`
    **success**: logical  **intent(out)**    *scalar*  
        logical tag specifying if creation of the potential succeeded
            
  .. function:: create_potential_characterizer_LJ(index)

    LJ characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_bond_bending(index)

    bond-bending characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_buckingham(index)

    Buckingham characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_charge_exp(index)

    charge exponential characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_charge_pair(index)

    charge self energy characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_charge_pair_abs(index)

    charge abs energy characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_charge_self(index)

    charge self energy characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_constant_force(index)

    constant F characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_constant_potential(index)

    constant potential characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_dihedral(index)

    dihedral angle characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_exp(index)

    exponential characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_power(index)

    Power law characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_shift(index)

    Shift characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_spring(index)

    spring characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: create_potential_characterizer_table(index)

    Tabulated characterizer initialization
    

    Parameters:

    index: integer  *intent(in)*    *scalar*  
        index of the potential
            
  .. function:: evaluate_bond_order_factor(n_targets, separations, distances, bond_params, factor, atoms)

    Returns a bond order factor term.
    
    By a bond order factor term, we mean the contribution from
    specific atoms, :math:`c_{ijk}`, appearing in the factor
    
    .. math::
    
          b_i = f(\sum_{jk} c_{ijk})
    
    This routine evaluates the term :math:`c_{ij}` or :math:`c_{ijk}` for the given
    atoms :math:`ij` or :math:`ijk` according to the given parameters.
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    separations: double precision  *intent(in)*    *size(3, n_targets-1)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(n_targets-1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *size(n_targets-1)*  
        a :data:`bond_order_parameters` containing the parameters
    **factor**: double precision  **intent(out)**    *size(n_targets)*  
        the calculated bond order term :math:`c`
    atoms: type(atom)  *intent(in)*    *size(n_targets)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_bond_order_factor_coordination(separations, distances, bond_params, factor)

    Coordination bond order factor
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *size(1)*  
        a :data:`bond_order_parameters` containing the parameters
    **factor**: double precision  **intent(out)**    *size(2)*  
        the calculated bond order term :math:`c`
            
  .. function:: evaluate_bond_order_factor_power(separations, distances, bond_params, factor)

    Power bond order factor
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *size(1)*  
        a :data:`bond_order_parameters` containing the parameters
    **factor**: double precision  **intent(out)**    *size(2)*  
        the calculated bond order term :math:`c`
            
  .. function:: evaluate_bond_order_factor_table(separations, distances, bond_params, factor)

    Tabulated bond order factor
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *size(1)*  
        a :data:`bond_order_parameters` containing the parameters
    **factor**: double precision  **intent(out)**    *size(2)*  
        the calculated bond order term :math:`c`
            
  .. function:: evaluate_bond_order_factor_triplet(separations, distances, bond_params, factor, atoms)

    Triplet bond factor
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 2)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(2)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *size(2)*  
        a :data:`bond_order_parameters` containing the parameters
    **factor**: double precision  **intent(out)**    *size(3)*  
        the calculated bond order term :math:`c`
    atoms: type(atom)  *intent(in)*    *size(3)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_bond_order_gradient(n_targets, separations, distances, bond_params, gradient, atoms)

    Returns the gradients of bond order terms with respect to moving an atom.
    
    By a bond order factor term, we mean the contribution from
    specific atoms, c_ijk, appearing in the factor
    
    .. math::
    
          b_i = f(\sum_{jk} c_{ijk})
    
    This routine evaluates the gradient term :math:`\nabla_\alpha c_{ij}` or
    :math:`\nabla_\alpha c_{ijk}` for the given atoms :math:`ij` or :math:`ijk` according to the given parameters.
    
    The returned array has three dimensions:
    gradient( coordinates, atom whose factor is differentiated, atom with respect to which we differentiate )
    So for example, for a three body term atom1-atom2-atom3, gradient(1,2,3) contains
    the x-coordinate (1), of the factor for atom2 (2), with respect to moving atom3 (3).
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    separations: double precision  *intent(in)*    *size(3, n_targets-1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(n_targets-1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *size(n_targets-1)*  
        a :data:`bond_order_parameters` containing the parameters
    **gradient**: double precision  **intent(out)**    *size(3, n_targets, n_targets)*  
        the calculated bond order term :math:`\nabla_\alpha c`
    atoms: type(atom)  *intent(in)*    *size(n_targets)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_bond_order_gradient_coordination(separations, distances, bond_params, gradient)

    Coordination bond order factor gradient
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *size(1)*  
        a :data:`bond_order_parameters` containing the parameters
    **gradient**: double precision  **intent(out)**    *size(3, 2, 2)*  
        the calculated bond order term :math:`c`
            
  .. function:: evaluate_bond_order_gradient_power(separations, distances, bond_params, gradient)

    Power bond order factor gradient
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *size(1)*  
        a :data:`bond_order_parameters` containing the parameters
    **gradient**: double precision  **intent(out)**    *size(3, 2, 2)*  
        the calculated bond order term :math:`c`
            
  .. function:: evaluate_bond_order_gradient_table(separations, distances, bond_params, gradient)

    Tabulated bond order factor gradient
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *size(1)*  
        a :data:`bond_order_parameters` containing the parameters
    **gradient**: double precision  **intent(out)**    *size(3, 2, 2)*  
        the calculated bond order term :math:`c`
            
  .. function:: evaluate_bond_order_gradient_triplet(separations, distances, bond_params, gradient, atoms)

    Coordination bond order factor gradient
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 2)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(2)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *size(2)*  
        a :data:`bond_order_parameters` containing the parameters
    **gradient**: double precision  **intent(out)**    *size(3, 3, 3)*  
        the calculated bond order term :math:`c`
    atoms: type(atom)  *intent(in)*    *size(3)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_electronegativity(n_targets, n_product, separations, distances, interaction, eneg, atoms)

    If a potential, say, :math:`U_{ijk}` depends on the charges of atoms :math:`q_i`
    it will not only create a force,
    but also a difference in chemical potential :math:`\mu_i` for the atomic partial charges.
    Similarly to :func:`evaluate_forces`, this function evaluates the chemical
    'force' on the atomic charges
    
    .. math::
    
       \chi_{\alpha,ijk} = -\mu_{\alpha,ijk} = -\frac{\partial U_{ijk}}{\partial q_\alpha}
    
    This routine can evaluate the contribution from a product potential.
    
    To be consist the forces returned by :func:`evaluate_electronegativity` must be
    derivatives of the energies returned by :func:`evaluate_energy`.
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    n_product: integer  *intent(in)*    *scalar*  
        the number of potentials for a product potential
    separations: double precision  *intent(in)*    *size(3, n_targets-1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(n_targets-1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **eneg**: double precision  **intent(out)**    *size(n_targets)*  
        the calculated electronegativity component :math:`\chi_{\alpha,ijk}`
    atoms: type(atom)  *intent(in)*    *size(n_targets)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_electronegativity_charge_exp(interaction, eneg, atoms)

    Charge exp electronegativity
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **eneg**: double precision  **intent(out)**    *size(2)*  
        the calculated electronegativity component :math:`\chi_{\alpha,ijk}`
    atoms: type(atom)  *intent(in)*    *size(2)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_electronegativity_charge_pair(interaction, eneg, atoms)

    Charge pair energy electronegativity
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **eneg**: double precision  **intent(out)**    *size(2)*  
        the calculated electronegativity component :math:`\chi_{\alpha,ijk}`
    atoms: type(atom)  *intent(in)*    *size(2)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_electronegativity_charge_pair_abs(interaction, eneg, atoms)

    Charge pair abs energy electronegativity
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **eneg**: double precision  **intent(out)**    *size(2)*  
        the calculated electronegativity component :math:`\chi_{\alpha,ijk}`
    atoms: type(atom)  *intent(in)*    *size(2)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_electronegativity_charge_self(interaction, eneg, atoms)

    Charge self energy electronegativity
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **eneg**: double precision  **intent(out)**    *size(1)*  
        the calculated electronegativity component :math:`\chi_{\alpha,ijk}`
    atoms: type(atom)  *intent(in)*    *size(1)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_electronegativity_component(n_targets, separations, distances, interaction, eneg, atoms)

    If a potential, say, :math:`U_{ijk}` depends on the charges of atoms :math:`q_i`
    it will not only create a force,
    but also a difference in chemical potential :math:`\mu_i` for the atomic partial charges.
    Similarly to :func:`evaluate_forces`, this function evaluates the chemical
    'force' on the atomic charges
    
    .. math::
    
       \chi_{\alpha,ijk} = -\mu_{\alpha,ijk} = -\frac{\partial U_{ijk}}{\partial q_\alpha}
    
    This routine evaluates an elemental :math:`\chi_{\alpha,ijk}`.
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    separations: double precision  *intent(in)*    *size(3, n_targets-1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(n_targets-1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **eneg**: double precision  **intent(out)**    *size(n_targets)*  
        the calculated electronegativity component :math:`\chi_{\alpha,ijk}`
    atoms: type(atom)  *intent(in)*    *size(n_targets)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_energy(n_targets, n_product, separations, distances, interaction, energy, atoms)

    Evaluates the potential energy due to an interaction between the given
    atoms. In other words, if the total potential energy is
    
    .. math::
    
       E = \sum_{ijk} v_{ijk}
    
    this routine evaluates :math:`v_{ijk}` for the given
    atoms i, j, and k.
    
    This routine can evaluate the contribution from a product potential.
    
    To be consist the forces returned by :func:`evaluate_forces` must be
    gradients of the energies returned by :func:`evaluate_energy`.
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    n_product: integer  *intent(in)*    *scalar*  
        the number of potentials for a product potential
    separations: double precision  *intent(in)*    *size(3, n_targets-1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(n_targets-1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
    atoms: type(atom)  *intent(in)*    *size(n_targets)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_energy_LJ(separations, distances, interaction, energy)

    LJ energy
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
            
  .. function:: evaluate_energy_bond_bending(separations, distances, interaction, energy, atoms)

    Bond bending energy
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 2)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(2)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
    atoms: type(atom)  *intent(in)*    *size(3)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_energy_buckingham(separations, distances, interaction, energy)

    Buckingham energy
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
            
  .. function:: evaluate_energy_charge_exp(interaction, energy, atoms)

    Charge exp energy
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
    atoms: type(atom)  *intent(in)*    *size(2)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_energy_charge_pair(interaction, energy, atoms)

    Charge pair energy
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
    atoms: type(atom)  *intent(in)*    *size(2)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_energy_charge_pair_abs(interaction, energy, atoms)

    Charge pair abs energy
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
    atoms: type(atom)  *intent(in)*    *size(2)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_energy_charge_self(interaction, energy, atoms)

    Charge self energy
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
    atoms: type(atom)  *intent(in)*    *size(1)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_energy_component(n_targets, separations, distances, interaction, energy, atoms)

    Evaluates the potential energy due to an interaction between the given
    atoms. In other words, if the total potential energy is
    
    .. math::
    
       E = \sum_{ijk} v_{ijk}
    
    this routine evaluates :math:`v_{ijk}` for the given
    atoms i, j, and k.
    
    This routine evaluates an elemental :math:`v_{\alpha,ijk}`.
    
    To be consist the forces returned by :func:`evaluate_forces` must be
    gradients of the energies returned by :func:`evaluate_energy`.
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    separations: double precision  *intent(in)*    *size(3, n_targets-1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(n_targets-1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
    atoms: type(atom)  *intent(in)*    *size(n_targets)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_energy_constant_force(interaction, energy, atoms)

    constant force energy
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
    atoms: type(atom)  *intent(in)*    *size(1)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_energy_constant_potential(interaction, energy)

    Constant potential energy
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
            
  .. function:: evaluate_energy_dihedral(separations, distances, interaction, energy, atoms)

    Dihedral angle energy
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 3)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(3)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
    atoms: type(atom)  *intent(in)*    *size(4)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_energy_exp(separations, distances, interaction, energy)

    Exp energy
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
            
  .. function:: evaluate_energy_power(separations, distances, interaction, energy)

    Power energy
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
            
  .. function:: evaluate_energy_shift(separations, distances, interaction, energy)

    Shift energy
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
            
  .. function:: evaluate_energy_spring(separations, distances, interaction, energy)

    spring energy
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
            
  .. function:: evaluate_energy_table(separations, distances, interaction, energy)

    Tabulated energy
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **energy**: double precision  **intent(out)**    *scalar*  
        the calculated energy :math:`v_{ijk}`
            
  .. function:: evaluate_force_LJ(separations, distances, interaction, force)

    LJ force
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 2)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
            
  .. function:: evaluate_force_bond_bending(separations, distances, interaction, force, atoms)

    Bond bending force
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 2)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(2)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 3)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
    atoms: type(atom)  *intent(in)*    *size(3)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_force_buckingham(separations, distances, interaction, force)

    Buckingham force
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 2)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
            
  .. function:: evaluate_force_component(n_targets, separations, distances, interaction, force, atoms)

    Evaluates the forces due to an interaction between the given
    atoms. In other words, if the total force on atom :math:`\alpha` is
    
    .. math::
    
       \mathbf{F}_\alpha = \sum_{ijk} -\nabla_\alpha v_{ijk} = \sum \mathbf{f}_{\alpha,ijk},
    
    this routine evaluates :math:`\mathbf{f}_{\alpha,ijk}` for :math:`\alpha = (i,j,k)` for the given
    atoms i, j, and k.
    
    This routine evaluates an elemental :math:`\mathbf{f}_{\alpha,ijk}`.
    
    To be consist the forces returned by :func:`evaluate_forces` must be
    gradients of the energies returned by :func:`evaluate_energy`.
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    separations: double precision  *intent(in)*    *size(3, n_targets-1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(n_targets-1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, n_targets)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
    atoms: type(atom)  *intent(in)*    *size(n_targets)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_force_constant_force(interaction, force)

    constant force
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 1)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
            
  .. function:: evaluate_force_constant_potential(interaction, force)

    constant force
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 1)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
            
  .. function:: evaluate_force_dihedral(separations, distances, interaction, force, atoms)

    Dihedral angle force
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 3)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(3)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 4)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
    atoms: type(atom)  *intent(in)*    *size(4)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_force_exp(separations, distances, interaction, force)

    Exp force
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 2)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
            
  .. function:: evaluate_force_power(separations, distances, interaction, force)

    Power force
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 2)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
            
  .. function:: evaluate_force_shift(separations, distances, interaction, force)

    Shift force
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 2)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
            
  .. function:: evaluate_force_spring(separations, distances, interaction, force)

    spring force
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 2)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
            
  .. function:: evaluate_force_table(separations, distances, interaction, force)

    Tabulated force
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, 2)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
            
  .. function:: evaluate_forces(n_targets, n_product, separations, distances, interaction, force, atoms)

    Evaluates the forces due to an interaction between the given
    atoms. In other words, if the total force on atom :math:`\alpha` is
    
    .. math::
    
       \mathbf{F}_\alpha = \sum_{ijk} -\nabla_\alpha v_{ijk} = \sum \mathbf{f}_{\alpha,ijk},
    
    this routine evaluates :math:`\mathbf{f}_{\alpha,ijk}` for :math:`\alpha = (i,j,k)` for the given
    atoms i, j, and k.
    
    This routine can evaluate the contribution from a product potential.
    
    To be consist the forces returned by :func:`evaluate_forces` must be
    gradients of the energies returned by :func:`evaluate_energy`.
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    n_product: integer  *intent(in)*    *scalar*  
        the number of potentials for a product potential
    separations: double precision  *intent(in)*    *size(3, n_targets-1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(n_targets-1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    interaction: type(potential)  *intent(in)*    *scalar*  
        a :data:`potential` containing the parameters
    **force**: double precision  **intent(out)**    *size(3, n_targets)*  
        the calculated force component :math:`\mathbf{f}_{\alpha,ijk}`
    atoms: type(atom)  *intent(in)*    *size(n_targets)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_pair_bond_order_factor(n_targets, separations, distances, bond_params, factor, atoms)

    Returns a bond order factor term.
    
    By a bond order factor term, we mean the contribution from
    specific atoms, :math:`c_{ijk}`, appearing in the factor
    
    .. math::
    
          b_ij = f(\sum_{k} c_{ijk})
    
    This routine evaluates the term :math:`c_{ijk}` for the given
    atoms :math:`ijk` according to the given parameters.
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    separations: double precision  *intent(in)*    *size(3, n_targets-1)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(n_targets-1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **factor**: double precision  **intent(out)**    *scalar*  
        the calculated bond order term :math:`c`
    atoms: type(atom)  *intent(in)*    *size(n_targets)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_pair_bond_order_factor_tersoff(separations, distances, bond_params, factor, atoms)


    Parameters:

    separations: double precision  *intent(in)*    *size(3, 2)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(2)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **factor**: double precision  **intent(out)**    *scalar*  
        the calculated bond order term :math:`c`
    atoms: type(atom)  *intent(in)*    *size(3)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_pair_bond_order_gradient(n_targets, separations, distances, bond_params, gradient, atoms)

    Returns the gradients of bond order terms with respect to moving an atom.
    
    By a bond order factor term, we mean the contribution from
    specific atoms, c_ijk, appearing in the factor
    
    .. math::
    
          b_ij = f(\sum_{k} c_{ijk})
    
    This routine evaluates the gradient term
    :math:`\nabla_\alpha c_{ijk}` for the given atoms :math:`ijk` according to the given parameters.
    
    The returned array has two dimensions:
    gradient( coordinates, atom with respect to which we differentiate ).
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    separations: double precision  *intent(in)*    *size(3, n_targets-1)*  
        atom-atom separation vectors :math:`\mathrm{r}_{12}`, :math:`\mathrm{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(n_targets-1)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **gradient**: double precision  **intent(out)**    *size(3, n_targets)*  
        the calculated bond order term :math:`\nabla_\alpha c`
    atoms: type(atom)  *intent(in)*    *size(n_targets)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: evaluate_pair_bond_order_gradient_tersoff(separations, distances, bond_params, gradient, atoms)

    Tersoff bond order factor gradient
    

    Parameters:

    separations: double precision  *intent(in)*    *size(3, 2)*  
        atom-atom separation vectors :math:`\mathbf{r}_{12}`, :math:`\mathbf{r}_{23}` etc. for the atoms 123...
    distances: double precision  *intent(in)*    *size(2)*  
        atom-atom distances :math:`r_{12}`, :math:`r_{23}` etc. for the atoms 123..., i.e., the norms of the separation vectors.
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` containing the parameters
    **gradient**: double precision  **intent(out)**    *size(3, 3)*  
        the calculated bond order term :math:`c`
    atoms: type(atom)  *intent(in)*    *size(3)*  
        a list of the actual :data:`atom` objects for which the term is calculated
            
  .. function:: get_bond_descriptor(bond_name, descriptor)

    Returns the :data:`bond_order_descriptor` of a given name.
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    **descriptor**: type(bond_order_descriptor)  **intent(out)**    *scalar*  
        the matching :data:`bond_order_descriptor`
            
  .. function:: get_description_of_bond_order_factor(bond_name, description)

    Returns the description of a bond order factor.
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    **description**: character(len=pot_note_length)  **intent(out)**    *scalar*  
        description of the bond order factor
            
  .. function:: get_description_of_potential(pot_name, description)

    Returns the description of a potential.
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    **description**: character(len=pot_note_length)  **intent(out)**    *scalar*  
        description of the potential
            
  .. function:: get_descriptions_of_parameters_of_bond_order_factor(bond_name, n_targets, param_notes)

    Returns the descriptions of the parameters of a bond order factor
    as a list of strings.
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    param_notes: character(len=param_note_length)  *intent()*  *pointer*  *size(:)*  
        descriptions of the parameters
            
  .. function:: get_descriptions_of_parameters_of_potential(pot_name, param_notes)

    Returns the descriptions of the parameters of a potential
    as a list of strings.
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    param_notes: character(len=param_note_length)  *intent()*  *pointer*  *size(:)*  
        descriptions of the parameters
            
  .. function:: get_descriptor(pot_name, descriptor)

    Returns the :data:`potential_descriptor` of a given name.
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    **descriptor**: type(potential_descriptor)  **intent(out)**    *scalar*  
        the matching :data:`potential_descriptor`
            
  .. function:: get_index_of_bond_order_factor(bond_name, index)

    Returns the index of a :data:`bond_order_descriptor` in the internal list of bond order factor types :data:`bond_order_descriptors`.
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor - a keyword
    **index**: integer  **intent(out)**    *scalar*  
        index of the potential in the internal array
            
  .. function:: get_index_of_parameter_of_bond_order_factor(bond_name, param_name, index, n_targets)

    Returns the index of a parameter of a bond order factor in the
    internal list of parameters. Since bond order factors can have
    parameters for different number of targets, also the number of
    targets of this parameter is returned.
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    param_name: character(len=*)  *intent(in)*    *scalar*  
        name of the parameter
    **index**: integer  **intent(out)**    *scalar*  
        index of the parameter
    **n_targets**: integer  **intent(out)**    *scalar*  
        number of targets of the parameter
            
  .. function:: get_index_of_parameter_of_potential(pot_name, param_name, index)

    Returns the index of a parameter of a potential in the
    internal list of parameters.
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    param_name: character(len=*)  *intent(in)*    *scalar*  
        name of the parameter
    **index**: integer  **intent(out)**    *scalar*  
        the index of the parameter
            
  .. function:: get_index_of_potential(pot_name, index)

    Returns the index of a :data:`potential_descriptor` in the internal list of potential types :data:`potential_descriptors`.
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential - a keyword
    **index**: integer  **intent(out)**    *scalar*  
        index of the potential in the internal array
            
  .. function:: get_names_of_parameters_of_bond_order_factor(bond_name, n_targets, param_names)

    Returns the names of parameters of a bond order factor as a list of strings.
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    param_names: character(len=param_name_length)  *intent()*  *pointer*  *size(:)*  
        names of the parameters
            
  .. function:: get_names_of_parameters_of_potential(pot_name, param_names)

    Returns the names of parameters of a potential as a list of strings.
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    param_names: character(len=param_name_length)  *intent()*  *pointer*  *size(:)*  
        names of the parameters
            
  .. function:: get_number_of_bond_order_factors(n_bond)

    Returns the number of :data:`bond_order_descriptor` known.
    

    Parameters:

    **n_bond**: integer  **intent(out)**    *scalar*  
        number of bond order factor types
            
  .. function:: get_number_of_parameters_of_bond_order_factor(bond_name, n_targets, n_params)

    Returns the number of parameters of a bond order factor as a list of strings,
    each element showing the number of parameters for that number of bodies.
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    **n_params**: integer  **intent(out)**    *scalar*  
        number of parameters
            
  .. function:: get_number_of_parameters_of_potential(pot_name, n_params)

    Returns the number of parameters of a potential.
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    **n_params**: integer  **intent(out)**    *scalar*  
        number of parameters
            
  .. function:: get_number_of_potentials(n_pots)

    Return the number of :data:`potential_descriptor`  known.
    

    Parameters:

    **n_pots**: integer  **intent(out)**    *scalar*  
        number of potential types
            
  .. function:: get_number_of_targets_of_bond_order_factor(bond_name, n_target)

    Returns the number of targets (i.e., bodies) of a bond order factor.
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    **n_target**: integer  **intent(out)**    *scalar*  
        number of targets
            
  .. function:: get_number_of_targets_of_bond_order_factor_index(bond_index, n_target)

    Returns the number of targets (i.e., bodies) of a bond order factor
    specified by its index.
    

    Parameters:

    bond_index: integer  *intent(in)*    *scalar*  
        index of the bond order factor
    **n_target**: integer  **intent(out)**    *scalar*  
        number of targets
            
  .. function:: get_number_of_targets_of_potential(pot_name, n_target)

    Returns the number of targets (i.e., bodies) of a potential.
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    **n_target**: integer  **intent(out)**    *scalar*  
        number of targets
            
  .. function:: get_number_of_targets_of_potential_index(pot_index, n_target)

    Returns the number of targets (i.e., bodies) of a potential
    specified by its index.
    

    Parameters:

    pot_index: integer  *intent(in)*    *scalar*  
        index of the potential
    **n_target**: integer  **intent(out)**    *scalar*  
        numner of targets
            
  .. function:: initialize_bond_order_factor_characterizers()

    Creates bond order factor characterizers.
    
    This routine is meant to be run once, as pysic is
    imported, to create the characterizers for
    bond order factors. Once created, they are accessible
    by both the python and fortran sides of pysic
    as a tool for describing the general structure
    of bond order factor objects.

            
  .. function:: initialize_potential_characterizers()

    Creates potential characterizers.
    
    This routine is meant to be run once, as pysic is
    imported, to create the characterizers for
    potentials. Once created, they are accessible
    by both the python and fortran sides of pysic
    as a tool for describing the general structure
    of potential objects.

            
  .. function:: is_valid_bond_order_factor(string, is_valid)

    Returns true if the given keyword is the name of a bond order factor
    and false otherwise.
    

    Parameters:

    string: character(len=*)  *intent(in)*    *scalar*  
        name of a bond order factor
    **is_valid**: logical  **intent(out)**    *scalar*  
        true if string is a name of a bond order factor
            
  .. function:: is_valid_potential(string, is_valid)

    Returns true if the given keyword is the name of a potential
    and false otherwise.
    

    Parameters:

    string: character(len=*)  *intent(in)*    *scalar*  
        name of a potential
    **is_valid**: logical  **intent(out)**    *scalar*  
        true if string is a name of a potential
            
  .. function:: list_bond_order_factors(n_bonds, bonds)

    Returns the names of :data:`bond_order_descriptor` objects.
    

    Parameters:

    n_bonds: integer  *intent(in)*    *scalar*  
        number of bond order factor types
    **bonds**: character(len=pot_name_length)  **intent(out)**    *size(n_bonds)*  
        names of the bond order factor types
            
  .. function:: list_potentials(n_pots, pots)

    Returns the names of :data:`potential_descriptor` objects.
    

    Parameters:

    n_pots: integer  *intent(in)*    *scalar*  
        number of potential types
    **pots**: character(len=pot_name_length)  **intent(out)**    *size(n_pots)*  
        names of the potential types
            
  .. function:: post_process_bond_order_factor(raw_sum, bond_params, factor_out)

    Bond-order post processing, i.e.,
    application of per-atom scaling functions.
    
    By post processing, we mean any operations done after calculating the
    sum of pair- and many-body terms. That is, if a factor is, say,
    
    .. math::
    
         b_i = f(\sum_j c_{ij}) = 1 + \sum_j c_{ij},
    
    the :math:`\sum_j c_ij` would have been calculated already and the
    operation :math:`f(x) = 1 + x` remains to be carried out.
    The post processing is done per atom regardless of if the
    bond factor is of a pair or many body type.
    
    This routine applies the scaling function on the given
    bond order sum accoding to the given parameters.
    

    Parameters:

    raw_sum: double precision  *intent(in)*    *scalar*  
        the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` specifying the parameters
    **factor_out**: double precision  **intent(out)**    *scalar*  
        the calculated bond order factor :math:`b_i`
            
  .. function:: post_process_bond_order_factor_scaler_1(raw_sum, bond_params, factor_out)

    Scaler post process factor
    

    Parameters:

    raw_sum: double precision  *intent(in)*    *scalar*  
        the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` specifying the parameters
    **factor_out**: double precision  **intent(out)**    *scalar*  
        the calculated bond order factor :math:`b_i`
            
  .. function:: post_process_bond_order_factor_scaler_sqrt(raw_sum, bond_params, factor_out)

    Square root scaler post process factor
    

    Parameters:

    raw_sum: double precision  *intent(in)*    *scalar*  
        the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` specifying the parameters
    **factor_out**: double precision  **intent(out)**    *scalar*  
        the calculated bond order factor :math:`b_i`
            
  .. function:: post_process_bond_order_factor_scaler_table(raw_sum, bond_params, factor_out)

    Tabulated scaler post process factor
    

    Parameters:

    raw_sum: double precision  *intent(in)*    *scalar*  
        the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` specifying the parameters
    **factor_out**: double precision  **intent(out)**    *scalar*  
        the calculated bond order factor :math:`b_i`
            
  .. function:: post_process_bond_order_factor_tersoff(raw_sum, bond_params, factor_out)

    Tersoff post process factor
    

    Parameters:

    raw_sum: double precision  *intent(in)*    *scalar*  
        the precalculated bond order sum, :math:`\sum_k c_{ijk}` in the above example
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` specifying the parameters
    **factor_out**: double precision  **intent(out)**    *scalar*  
        the calculated bond order factor :math:`b_{ij}`
            
  .. function:: post_process_bond_order_gradient(raw_sum, raw_gradient, bond_params, factor_out)

    Bond-order post processing, i.e.,
    application of per-atom scaling functions.
    
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
    
    This routine applies the scaling function on the given
    bond order sum and gradient accoding to the given parameters.
    

    Parameters:

    raw_sum: double precision  *intent(in)*    *scalar*  
        the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
    raw_gradient: double precision  *intent(in)*    *size(3)*  
        the precalculated bond order gradient sum, :math:`\nabla_\alpha \sum_j c_ij` in the above example
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` specifying the parameters
    **factor_out**: double precision  **intent(out)**    *size(3)*  
        the calculated bond order factor :math:`\nabla_\alpha b_i`
            
  .. function:: post_process_bond_order_gradient_scaler_1(raw_sum, raw_gradient, bond_params, factor_out)

    Scaler post process gradient
    

    Parameters:

    raw_sum: double precision  *intent(in)*    *scalar*  
        the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
    raw_gradient: double precision  *intent(in)*    *size(3)*  
        the precalculated bond order gradient sum, :math:`\nabla_\alpha \sum_j c_ij` in the above example
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` specifying the parameters
    **factor_out**: double precision  **intent(out)**    *size(3)*  
        the calculated bond order factor :math:`b_i`
            
  .. function:: post_process_bond_order_gradient_scaler_sqrt(raw_sum, raw_gradient, bond_params, factor_out)

    Square root scaler post process gradient
    

    Parameters:

    raw_sum: double precision  *intent(in)*    *scalar*  
        the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
    raw_gradient: double precision  *intent(in)*    *size(3)*  
        the precalculated bond order gradient sum, :math:`\nabla_\alpha \sum_j c_ij` in the above example
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` specifying the parameters
    **factor_out**: double precision  **intent(out)**    *size(3)*  
        the calculated bond order factor :math:`b_i`
            
  .. function:: post_process_bond_order_gradient_scaler_table(raw_sum, raw_gradient, bond_params, factor_out)

    Tabulated scaler post process gradient
    

    Parameters:

    raw_sum: double precision  *intent(in)*    *scalar*  
        the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
    raw_gradient: double precision  *intent(in)*    *size(3)*  
        the precalculated bond order gradient sum, :math:`\nabla_\alpha \sum_j c_ij` in the above example
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` specifying the parameters
    **factor_out**: double precision  **intent(out)**    *size(3)*  
        the calculated bond order factor :math:`b_i`
            
  .. function:: post_process_bond_order_gradient_tersoff(raw_sum, raw_gradient, bond_params, factor_out)

    Tersoff post process gradient
    

    Parameters:

    raw_sum: double precision  *intent(in)*    *scalar*  
        the precalculated bond order sum, :math:`\sum_j c_ij` in the above example
    raw_gradient: double precision  *intent(in)*    *size(3)*  
        the precalculated bond order gradient sum, :math:`\nabla_\alpha \sum_j c_ij` in the above example
    bond_params: type(bond_order_parameters)  *intent(in)*    *scalar*  
        a :data:`bond_order_parameters` specifying the parameters
    **factor_out**: double precision  **intent(out)**    *size(3)*  
        the calculated bond order factor :math:`b_i`
            
  .. function:: potential_affects_atom(interaction, atom_in, affects, position)

    Tests whether the given potential affects the specific atom.
    
    For potentials, the atoms are specified as valid targets by
    the atomic symbol, index, or tag.
    
    If position is not given, then the routine returns true if
    the atom can appear in the potential in any role.
    If position is given, then true is returned only if the atom
    is valid for that particular position.
    
    For instance, in a 3-body potential A-B-C, the potential
    May be specified so that only certain elements are valid for
    positions A and C while some other elements are valid for B.
    In a water molecule, for instance, we could have an H-O-H
    bond bending potential, but no H-H-O potentials.
    

    Parameters:

    interaction: type(potential)  *intent(in)*    *scalar*  
        the :data:`potential`
    atom_in: type(atom)  *intent(in)*    *scalar*  
        the :data:`atom`
    **affects**: logical  **intent(out)**    *scalar*  
        true if the potential affects the atom
    position: integer  *intent(in)*    *scalar*  *optional*
        specifies the particular role of the atom in the interaction
            
  .. function:: smoothening_derivative(r, hard_cut, soft_cut, factor)

    Derivative of the function for smoothening potential
    and bond order cutoffs.
    In principle any "nice" function which goes from 1 to 0
    in a finite interval could be used. Here, we choose
    
    .. math::
    
      f(r) = \frac{1}{2} ( 1 + \cos \pi \frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}} )
    
    for :math:`r \in [r_\mathrm{soft},r_\mathrm{hard}]`.
    The derivative is then
    
    .. math::
    
      f'(r) = \frac{\pi}{2 (r_\mathrm{soft}-r_\mathrm{hard})} \sin \pi \frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}}.
    
    This routine takes as arguments :math:`r`, :math:`r_\mathrm{soft}`, and :math:`r_\mathrm{hard}`, and
    returns the value of the derivative of the smoothening function.
    

    Parameters:

    r: double precision  *intent(in)*    *scalar*  
        distance :math:`r`
    hard_cut: double precision  *intent(in)*    *scalar*  
        the hard cutoff :math:`r_\mathrm{hard}`
    soft_cut: double precision  *intent(in)*    *scalar*  
        the soft cutoff :math:`r_\mathrm{soft}`
    **factor**: double precision  **intent(out)**    *scalar*  
        the calculated derivative of the smoothening factor
            
  .. function:: smoothening_factor(r, hard_cut, soft_cut, factor)

    Function for smoothening potential and bond order cutoffs.
    In principle any "nice" function which goes from 1 to 0
    in a finite interval could be used. Here, we choose
    
    .. math::
    
      f(r) = \frac{1}{2} ( 1 + \cos \pi \frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}} )
    
    for :math:`r \in [r_\mathrm{soft},r_\mathrm{hard}]`.
    
    This routine takes as arguments :math:`r`, :math:`r_\mathrm{soft}`, and :math:`r_\mathrm{hard}`, and
    returns the value of the smoothening function.
    

    Parameters:

    r: double precision  *intent(in)*    *scalar*  
        distance :math:`r`
    hard_cut: double precision  *intent(in)*    *scalar*  
        the hard cutoff :math:`r_\mathrm{hard}`
    soft_cut: double precision  *intent(in)*    *scalar*  
        the soft cutoff :math:`r_\mathrm{soft}`
    **factor**: double precision  **intent(out)**    *scalar*  
        the calculated smoothening factor
            
  .. function:: smoothening_gradient(unit_vector, r, hard_cut, soft_cut, gradient)

    Gradient of the function for smoothening potential
    and bond order cutoffs.
    In principle any "nice" function which goes from 1 to 0
    in a finite interval could be used. Here, we choose
    
    .. math::
    
      f(r) = \frac{1}{2} ( 1 + \cos \pi \frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}} )
    
    for :math:`r \in [r_\mathrm{soft},r_\mathrm{hard}]`.
    The derivative is then
    
    .. math::
    
      f'(r) = \frac{\pi}{2 (r_\mathrm{soft}-r_\mathrm{hard})} \sin \pi \frac{r-r_\mathrm{soft}}{r_\mathrm{hard}-r_\mathrm{soft}}.
    
    and the gradient with respect to :math:`r`
    
    .. math::
    
      \nabla f(r) = f'(r) \nabla r = f'(r) \hat{r}
    
    where :math:`\hat{r}` is the unit vector in the direction of :math:`\mathbf{r}`.
    
    This routine takes as arguments :math:`\hat{r}`, :math:`r`, :math:`r_\mathrm{soft}`, and :math:`r_\mathrm{hard}`, and
    returns the value of the gradient of the smoothening function.
    

    Parameters:

    unit_vector: double precision  *intent(in)*    *size(3)*  
        the vector :math:`\hat{r}`
    r: double precision  *intent(in)*    *scalar*  
        distance :math:`r`
    hard_cut: double precision  *intent(in)*    *scalar*  
        the hard cutoff :math:`r_\mathrm{hard}`
    soft_cut: double precision  *intent(in)*    *scalar*  
        the soft cutoff :math:`r_\mathrm{soft}`
    **gradient**: double precision  **intent(out)**    *size(3)*  
        the calculated derivative of the smoothening factor