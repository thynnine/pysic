
.. _pysic_interface:
        
=========================================================
pysic_interface (PyInterface.f90)
=========================================================



Pysic_interface is an interface module between the Python and Fortran
sides of pysic. When pysic is compiled, only this file is interfaced to
Python via `f2py`_ while the rest of the Fortran source files are
directly compiled to .mod Fortran modules. The main reason for this
is that F90 derived types are used extensively in the core and these
are not yet (2011) supported by `f2py`_ (although the support for
derived types is planned in third generation f2py). Because of this,
most data is passed from :mod:`pysic` to
:ref:`pysic_interface` as NumPy arrays, and conversions from objects
is needed on both sides. This is cumbersome and adds overhead, but
it is not an efficiency issue since most of the information is only
passed from Python to Fortran once and saved. Even during a molecular
dynamics simulation only the forces, coordinates and momenta
of atoms need to be communicated through the interface, which is
naturally and efficiently handled using just numeric arrays anyway.

Another limitation in current `f2py`_ is handling of string arrays.
To overcome this, string arrays are converted to integer arrays
and back using simple mapping functions in :mod:`~pysic.pysic_utility`
and :ref:`utility`.

Due to the current limitations of `f2py`_, no derived types can appear
in the module. This severly limits what the module can do, and therefore
the module has been by design made to be very light in terms of
functionality: No data is stored in the module and almost all routines
simply redirect the call to a submodule, most often :ref:`pysic_core`.
In the descriptions of the routines in this documentation,
links are provided to the submodule routines that the call is directed
to, if the routine is just a redirect of the call.

.. _f2py: http://www.scipy.org/F2py


.. only:: html


    Modules used by pysic_interface
    -------------------------------
    - :ref:`mpi`
    - :ref:`mt95`
    - :ref:`pysic_core`
    - :ref:`utility`

    List of subroutines in pysic_interface
    --------------------------------------
        
    - :func:`add_bond_order_factor`
    - :func:`add_potential`
    - :func:`allocate_bond_order_factors`
    - :func:`allocate_bond_order_storage`
    - :func:`allocate_potentials`
    - :func:`calculate_bond_order_factors`
    - :func:`calculate_bond_order_gradients`
    - :func:`calculate_bond_order_gradients_of_factor`
    - :func:`calculate_electronegativities`
    - :func:`calculate_energy`
    - :func:`calculate_forces`
    - :func:`clear_potential_multipliers`
    - :func:`create_atoms`
    - :func:`create_bond_order_factor_list`
    - :func:`create_cell`
    - :func:`create_neighbor_list`
    - :func:`create_potential_list`
    - :func:`description_of_bond_order_factor`
    - :func:`description_of_potential`
    - :func:`descriptions_of_parameters_of_bond_order_factor`
    - :func:`descriptions_of_parameters_of_potential`
    - :func:`distribute_mpi`
    - :func:`examine_atoms`
    - :func:`examine_bond_order_factors`
    - :func:`examine_cell`
    - :func:`examine_potentials`
    - :func:`finish_mpi`
    - :func:`generate_neighbor_lists`
    - :func:`get_cell_vectors`
    - :func:`get_cpu_id`
    - :func:`get_ewald_energy`
    - :func:`get_mpi_list_of_atoms`
    - :func:`get_neighbor_list_of_atom`
    - :func:`get_number_of_atoms`
    - :func:`get_number_of_cpus`
    - :func:`get_number_of_neighbors_of_atom`
    - :func:`is_bond_order_factor`
    - :func:`is_potential`
    - :func:`level_of_bond_order_factor`
    - :func:`list_valid_bond_order_factors`
    - :func:`list_valid_potentials`
    - :func:`names_of_parameters_of_bond_order_factor`
    - :func:`names_of_parameters_of_potential`
    - :func:`number_of_bond_order_factors`
    - :func:`number_of_parameters_of_bond_order_factor`
    - :func:`number_of_parameters_of_potential`
    - :func:`number_of_potentials`
    - :func:`number_of_targets_of_bond_order_factor`
    - :func:`number_of_targets_of_potential`
    - :func:`release`
    - :func:`set_ewald_parameters`
    - :func:`start_bond_order_factors`
    - :func:`start_mpi`
    - :func:`start_potentials`
    - :func:`start_rng`
    - :func:`sync_mpi`
    - :func:`update_atom_charges`
    - :func:`update_atom_coordinates`


Full documentation of subroutines in pysic_interface
----------------------------------------------------
        
        
            
  .. function:: add_bond_order_factor(n_targets, n_params, n_split, bond_name, parameters, param_split, cutoff, smooth_cut, elements, orig_elements, group_index, success)

    Creates a bond order factor in the core.
    The memory must have been allocated first using allocate_potentials.
    
    Calls :func:`core_add_bond_order_factor`
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets (interacting bodies)
    n_params: integer  *intent(in)*    *scalar*  
        number of parameters
    n_split: integer  *intent(in)*    *scalar*  
        number of subsets in the list of parameters, should equal n_targets
    bond_name: character(len=*)  *intent(in)*    *scalar*  
        bond order factor names
    parameters: double precision  *intent(in)*    *size(n_params)*  
        numeric parameters
    param_split: integer  *intent(in)*    *size(n_split)*  
        the numbers of parameters for 1-body, 2-body etc.
    cutoff: double precision  *intent(in)*    *scalar*  
        interaction hard cutoff
    smooth_cut: double precision  *intent(in)*    *scalar*  
        interaction soft cutoff
    elements: integer  *intent(in)*    *size(2, n_targets)*  
        atomic symbols specifying the elements the interaction acts on
    orig_elements: integer  *intent(in)*    *size(2, n_targets)*  
        original atomic symbols specifying the elements the interaction acts on
    group_index: integer  *intent(in)*    *scalar*  
        index denoting the potential to which the factor is connected
    **success**: logical  **intent(out)**    *scalar*  
        logical tag specifying if creation of the factor succeeded
            
  .. function:: add_potential(n_targets, n_params, pot_name, parameters, cutoff, smooth_cut, elements, tags, indices, orig_elements, orig_tags, orig_indices, pot_index, is_multiplier, success)

    Creates a potential in the core.
    The memory must have been allocated first using allocate_potentials.
    
    Calls :func:`core_add_potential`
    

    Parameters:

    n_targets: integer  *intent(in)*    *scalar*  
        number of targets (interacting bodies)
    n_params: integer  *intent(in)*    *scalar*  
        number of parameters
    pot_name: character(len=*)  *intent(in)*    *scalar*  
        potential names
    parameters: double precision  *intent(in)*    *size(n_params)*  
        numeric parameters
    cutoff: double precision  *intent(in)*    *scalar*  
        interaction hard cutoff
    smooth_cut: double precision  *intent(in)*    *scalar*  
        interaction soft cutoff
    elements: integer  *intent(in)*    *size(2, n_targets)*  
        atomic symbols specifying the elements the interaction acts on
    tags: integer  *intent(in)*    *size(n_targets)*  
        tags specifying the atoms the interaction acts on
    indices: integer  *intent(in)*    *size(n_targets)*  
        indices specifying the atoms the interaction acts on
    orig_elements: integer  *intent(in)*    *size(2, n_targets)*  
        original atomic symbols specifying the elements the interaction acts on
    orig_tags: integer  *intent(in)*    *size(n_targets)*  
        original tags specifying the atoms the interaction acts on
    orig_indices: integer  *intent(in)*    *size(n_targets)*  
        original indices specifying the atoms the interaction acts on
    pot_index: integer  *intent(in)*    *scalar*  
        index of the potential
    is_multiplier: logical  *intent(in)*    *scalar*  
        logical tag defining if the potential is a multiplier for a product potential
    **success**: logical  **intent(out)**    *scalar*  
        logical tag specifying if creation of the potential succeeded
            
  .. function:: allocate_bond_order_factors(n_bonds)

    Allocates memory for storing bond order parameters for describing the atomic interactions.
    Similar to the allocate_potentials routine.
    
    Calls :func:`core_allocate_bond_order_factors`
    

    Parameters:

    n_bonds: integer  *intent(in)*    *scalar*  
        number of bond order factors
            
  .. function:: allocate_bond_order_storage(n_atoms, n_groups, n_factors)

    Allocates memory for storing bond order factors for describing the atomic interactions.
    The difference to allocate_bond_order_factors is that this method allocates
    space for arrays used in storing actual calculated bond order factors. The other
    routine allocates space for storing the parameters used in the calculations.
    
    Calls :func:`core_allocate_bond_order_storage`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    n_groups: integer  *intent(in)*    *scalar*  
        number of bond order groups
    n_factors: integer  *intent(in)*    *scalar*  
        number of bond order parameters
            
  .. function:: allocate_potentials(n_pots)

    Allocates memory for storing potentials for describing the atomic interactions.
    It is more convenient to loop through the potentials and format them in a
    suitable way in python than in fortran. Therefore the core is first called
    through this routine in order to allocate memory for the potentials.
    Then, each potential is created individually.
    
    Calls :func:`core_allocate_potentials`
    

    Parameters:

    n_pots: integer  *intent(in)*    *scalar*  
        number of potentials
            
  .. function:: calculate_bond_order_factors(n_atoms, group_index, bond_orders)

    Returns bond order factors of the given group for all atoms.
    The group index is an identifier for the bond order parameters
    which are used for calculating one and the same factors.
    In practice, the Coordinators in pysic are indexed and this
    indexing is copied in the core. Thus the group index specifies
    the coordinator / potential.
    
    Calls :func:`core_get_bond_order_factors`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        index for the bond order factor group
    **bond_orders**: double precision  **intent(out)**    *size(n_atoms)*  
        the calculated bond order factors
            
  .. function:: calculate_bond_order_gradients(n_atoms, group_index, atom_index, gradients)

    Returns bond order factors gradients of the given group.
    The gradients of all factors are given with respect to moving the given atom.
    The group index is an identifier for the bond order parameters
    which are used for calculating one and the same factors.
    In practice, the Coordinators in pysic are indexed and this
    indexing is copied in the core. Thus the group index specifies
    the coordinator / potential.
    
    Calls :func:`core_get_bond_order_sums`
    
    and :func:`core_calculate_bond_order_gradients`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom with respect to which the factors are differentiated
    **gradients**: double precision  **intent(out)**    *size(3, n_atoms)*  
        the calculated bond order gradients
            
  .. function:: calculate_bond_order_gradients_of_factor(n_atoms, group_index, atom_index, gradients)

    Returns bond order factors gradients of the given group.
    The gradients of the given factors is given with respect to moving all atoms.
    The group index is an identifier for the bond order parameters
    which are used for calculating one and the same factors.
    In practice, the Coordinators in pysic are indexed and this
    indexing is copied in the core. Thus the group index specifies
    the coordinator / potential.
    
    Calls :func:`core_get_bond_order_sums`
    
    and :func:`core_calculate_bond_order_gradients_of_factor`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    group_index: integer  *intent(in)*    *scalar*  
        an index denoting the potential to which the factor is connected
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom whose factor is differentiated
    **gradients**: double precision  **intent(out)**    *size(3, n_atoms)*  
        the calculated bond order gradients
            
  .. function:: calculate_electronegativities(n_atoms, enegs)

    Returns electronegativities of the particles
    
    Calls :func:`core_calculate_electronegativities`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    **enegs**: double precision  **intent(out)**    *size(n_atoms)*  
        array of electronegativities on all atoms
            
  .. function:: calculate_energy(energy)

    Returns the total potential energy of the system
    
    Calls :func:`core_calculate_energy`
    

    Parameters:

    **energy**: double precision  **intent(out)**    *scalar*  
        total potential energy
            
  .. function:: calculate_forces(n_atoms, forces, stress)

    Returns forces acting on the particles and the stress tensor
    
    Calls :func:`core_calculate_forces`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    **forces**: double precision  **intent(out)**    *size(3, n_atoms)*  
        array of forces on all atoms
    **stress**: double precision  **intent(out)**    *size(6)*  
        array containing the components of the stress tensor (in order :math:`xx,yy,zz,yz,xz,xy`)
            
  .. function:: clear_potential_multipliers()

    Clears the temporary stored array of multiplier potentials

            
  .. function:: create_atoms(n_atoms, masses, charges, positions, momenta, tags, elements)

    Creates atomic particles.
    Atoms are handled as custom fortran types :data:`atom` in the core. Currently
    `f2py`_ does not support direct creation of types from Python, so instead
    all the necessary data is passed from Python as arrays and reassembled
    as types in Fortran. This is not much of an added overhead - the
    memory allocation itself already makes this a routine one does not
    wish to call repeatedly. Instead, one should call the routines
    for updating atoms whenever the actual atoms do not change
    (e.g., between MD timesteps).
    
    Calls :func:`core_generate_atoms`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    masses: double precision  *intent(in)*    *size(n_atoms)*  
        masses of atoms
    charges: double precision  *intent(in)*    *size(n_atoms)*  
        electric charges of atoms
    positions: double precision  *intent(in)*    *size(3, n_atoms)*  
        coordinates of atoms
    momenta: double precision  *intent(in)*    *size(3, n_atoms)*  
        momenta of atoms
    tags: integer  *intent(in)*    *size(n_atoms)*  
        numeric tags for the atoms
    elements: integer  *intent(in)*    *size(2, n_atoms)*  
        atomic symbols of the atoms
            
  .. function:: create_bond_order_factor_list()

    Similarly to the potential lists, also list containing all the
    bond order factors that may affect an atom are stored in a list.
    
    Calls :func:`core_assign_bond_order_factor_indices`

            
  .. function:: create_cell(vectors, inverse, periodicity)

    Creates a supercell for containing the calculation geometry
    Also the inverse cell matrix must be given,
    although it is not checked that the given inverse actually
    is the true inverse.
    
    Calls :func:`core_create_cell`
    

    Parameters:

    vectors: double precision  *intent(in)*    *size(3, 3)*  
        A 3x3 matrix containing the vectors spanning the supercell. The first index runs over xyz and the second index runs over the three vectors.
    inverse: double precision  *intent(in)*    *size(3, 3)*  
        A 3x3 matrix containing the inverse matrix of the one given in vectors, i.e. :math:`M^{-1}*M = I` for the two matrices. Since the latter represents a cell of non-zero volume, this inverse must exist. It is not tested that the given matrix actually is the inverse, the user must make sure it is.
    periodicity: logical  *intent(in)*    *size(3)*  
        A 3-element vector containing logical tags specifying if the system is periodic in the directions of the three vectors spanning the supercell.
            
  .. function:: create_neighbor_list(n_nbs, atom_index, neighbors, offsets)

    Creates neighbor lists for a single atom
    telling it which other atoms are in its
    immediate neighborhood.
    The neighbor list must be precalculated, this method only
    stores them in the core. The list must contain
    an array storing the indices of the neighboring atoms
    as well as the supercell offsets. The offsets are integer
    triplets showing how many times must the supercell vectors
    be added to the position of the neighbor to find the
    neighboring image in a periodic system.
    Note that if the system is small, one atom can in
    principle appear several times in the neighbor list.
    
    Calls :func:`core_create_neighbor_list`
    

    Parameters:

    n_nbs: integer  *intent(in)*    *scalar*  
        number of neighbors
    atom_index: integer  *intent(in)*    *scalar*  
        index of the atom for which the neighbor list is created
    neighbors: integer  *intent(in)*    *size(n_nbs)*  
        An array containing the indices of the neighboring atoms
    offsets: integer  *intent(in)*    *size(3, n_nbs)*  
        An array containing vectors specifying the offsets of the neighbors in periodic systems.
            
  .. function:: create_potential_list()

    Creates a list of indices for all atoms showing which potentials
    act on them.
    The user may define many potentials to sum up the potential energy of the
    system. However, if some potentials only act on certain atoms, they will
    be redundant for the other atoms. The potential lists are lists
    given to each atom containing the potentials which can act on the
    atom.
    
    Calls :func:`core_assign_potential_indices`

            
  .. function:: description_of_bond_order_factor(bond_name, description)

    Returns a description of the given bond order factor
    
    Calls :func:`get_description_of_bond_order_factor`
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    **description**: character(len=500)  **intent(out)**    *scalar*  
        description of the bond order actor
            
  .. function:: description_of_potential(pot_name, description)

    Returns a description of the given potential
    
    Calls :func:`get_description_of_potential`
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    **description**: character(len=500)  **intent(out)**    *scalar*  
        description of the potential
            
  .. function:: descriptions_of_parameters_of_bond_order_factor(bond_name, n_targets, param_notes)

    Lists descriptions for parameters the given bond order factor.
    Output is an array of integers. This is because `f2py`_ doesn't
    currently support string arrays. So, the characters are translated to
    integers and back in fortran and python.
    This adds a bit of overhead, but the routine is only invoked
    on user command so it doesn't matter.
    
    Calls :func:`get_descriptions_of_parameters_of_bond_order_factor`
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    **param_notes**: integer  **intent(out)**    *size(100, 12)*  
        descriptions of the parameters
            
  .. function:: descriptions_of_parameters_of_potential(pot_name, param_notes)

    Lists descriptions for parameters the given potential.
    Output is an array of integers. This is because `f2py`_ doesn't
    currently support string arrays. So, the characters are translated to
    integers and back in fortran and python.
    This adds a bit of overhead, but the routine is only invoked
    on user command so it doesn't matter.
    
    Calls :func:`get_descriptions_of_parameters_of_potential`
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    **param_notes**: integer  **intent(out)**    *size(100, 12)*  
        descriptions of the parameters
            
  .. function:: distribute_mpi(n_atoms)

    Distributes atoms among the processors.
    In the MPI scheme, atoms are distributed among
    the cpus for force and energy calculations.
    This routine initializes the arrays that
    tell each cpu which atoms it has to calculate
    interactions for. It can be called before
    the atoms are created in the core but one has to
    make sure the number of atoms specified in the last call
    matches the number of atoms in the core when a
    calculation is invoked.
    
    Calls :func:`mpi_distribute`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
            
  .. function:: examine_atoms()

    Prints some information about the atoms allocated in the core.
    This is mainly for debugging, as the python side should always
    dictate what is in the core.
    
    Calls :func:`list_atoms`

            
  .. function:: examine_bond_order_factors()

    Prints some information about the bond order factors allocated in the core.
    This is mainly for debugging, as the python side should always
    dictate what is in the core.
    
    Calls :func:`list_bonds`

            
  .. function:: examine_cell()

    Prints some information about the supercell allocated in the core.
    This is mainly for debugging, as the python side should always
    dictate what is in the core.
    
    Calls :func:`list_cell`

            
  .. function:: examine_potentials()

    Prints some information about the potential allocated in the core.
    This is mainly for debugging, as the python side should always
    dictate what is in the core.
    
    Calls :func:`list_interactions`

            
  .. function:: finish_mpi()

    Finishes MPI for parallel calculations.
    
    Calls :func:`mpi_finish`

            
  .. function:: generate_neighbor_lists(n_atoms, cutoffs)

    calculates and allocates neighbor lists

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        
    cutoffs: double precision  *intent(in)*    *size(n_atoms)*  
        
            
  .. function:: get_cell_vectors(vectors)

    Returns the vectors defining the simulation supercell.
    
    Calls :func:`core_get_cell_vectors`
    

    Parameters:

    **vectors**: double precision  **intent(out)**    *size(3, 3)*  
        A 3x3 matrix containing the vectors spanning the supercell. The first index runs over xyz and the second index runs over the three vectors.
            
  .. function:: get_cpu_id(id)

    Returns the MPI cpu id number, which is an
    integer between 0 and :math:`n_\mathrm{cpus}-1`,
    where :math:`n_\mathrm{cpus}` is the total
    number of cpus.
    

    Parameters:

    **id**: integer  **intent(out)**    *scalar*  
        cpu id number in MPI - 0 in serial mode
            
  .. function:: get_ewald_energy(real_cut, k_cut, reciprocal_cut, sigma, epsilon, energy)

    Debugging routine for Ewald

    Parameters:

    real_cut: double precision  *intent(in)*    *scalar*  
        
    k_cut: double precision  *intent(in)*    *scalar*  
        
    reciprocal_cut: integer  *intent(in)*    *size(3)*  
        
    sigma: double precision  *intent(in)*    *scalar*  
        
    epsilon: double precision  *intent(in)*    *scalar*  
        
    **energy**: double precision  **intent(out)**    *scalar*  
        
            
  .. function:: get_mpi_list_of_atoms(n_atoms, cpu_atoms)

    Returns a logical array containing true for every
    atom that is allocated to this cpu, and false
    for all other atoms.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    **cpu_atoms**: logical  **intent(out)**    *size(n_atoms)*  
        array of logical values showing which atoms are marked to be handled by this cpu
            
  .. function:: get_neighbor_list_of_atom(atom_index, n_neighbors, neighbors, offsets)

    Returns the list of neighbors for an atom

    Parameters:

    atom_index: integer  *intent(in)*    *scalar*  
        
    n_neighbors: integer  *intent(in)*    *scalar*  
        
    **neighbors**: integer  **intent(out)**    *size(n_neighbors)*  
        
    **offsets**: integer  **intent(out)**    *size(3, n_neighbors)*  
        
            
  .. function:: get_number_of_atoms(n_atoms)

    Counts the number of atoms in the current core
    
    Calls :func:`core_get_number_of_atoms`
    

    Parameters:

    **n_atoms**: integer  **intent(out)**    *scalar*  
        number of atoms
            
  .. function:: get_number_of_cpus(ncpu)

    Returns the MPI cpu count
    

    Parameters:

    **ncpu**: integer  **intent(out)**    *scalar*  
        the total number of cpus available
            
  .. function:: get_number_of_neighbors_of_atom(atom_index, n_neighbors)

    Returns the number of neighbors for an atom

    Parameters:

    atom_index: integer  *intent(in)*    *scalar*  
        
    **n_neighbors**: integer  **intent(out)**    *scalar*  
        
            
  .. function:: is_bond_order_factor(string, is_ok)

    Tells whether a given keyword defines a bond order factor or not
    
    Calls :func:`is_valid_bond_order_factor`
    

    Parameters:

    string: character(len=*)  *intent(in)*    *scalar*  
        name of a bond order factor
    **is_ok**: logical  **intent(out)**    *scalar*  
        true if string is a name of a bond order factor
            
  .. function:: is_potential(string, is_ok)

    Tells whether a given keyword defines a potential or not
    
    Calls :func:`is_valid_potential`
    

    Parameters:

    string: character(len=*)  *intent(in)*    *scalar*  
        name of a potential
    **is_ok**: logical  **intent(out)**    *scalar*  
        true if string is a name of a potential
            
  .. function:: level_of_bond_order_factor(bond_name, n_target)

    Tells the level of a bond order factor has, i.e., is it per-atom or per-pair
    
    Calls :func:`get_level_of_bond_order_factor`
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    **n_target**: integer  **intent(out)**    *scalar*  
        number of targets
            
  .. function:: list_valid_bond_order_factors(n_bonds, bond_factors)

    Lists all the keywords which define a bond order factor
    
    Calls :func:`list_bond_order_factors`
    

    Parameters:

    n_bonds: integer  *intent(in)*    *scalar*  
        number of bond order factor types
    **bond_factors**: integer  **intent(out)**    *size(11, n_bonds)*  
        names of the bond order factor types
            
  .. function:: list_valid_potentials(n_pots, potentials)

    Lists all the keywords which define a potential
    
    Calls :func:`list_potentials`
    

    Parameters:

    n_pots: integer  *intent(in)*    *scalar*  
        number of potential types
    **potentials**: integer  **intent(out)**    *size(11, n_pots)*  
        names of the potential types
            
  .. function:: names_of_parameters_of_bond_order_factor(bond_name, n_targets, param_names)

    Lists the names of parameters the given bond order factor knows.
    Output is an array of integers. This is because `f2py`_ doesn't
    currently support string arrays. So, the characters are translated to
    integers and back in fortran and python.
    This adds a bit of overhead, but the routine is only invoked
    on user command so it doesn't matter.
    
    Calls :func:`get_names_of_parameters_of_bond_order_factor`
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    **param_names**: integer  **intent(out)**    *size(10, 12)*  
        names of the parameters
            
  .. function:: names_of_parameters_of_potential(pot_name, param_names)

    Lists the names of parameters the given potential knows.
    Output is an array of integers. This is because `f2py`_ doesn't
    currently support string arrays. So, the characters are translated to
    integers and back in fortran and python.
    This adds a bit of overhead, but the routine is only invoked
    on user command so it doesn't matter.
    
    Calls :func:`get_names_of_parameters_of_potential`
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    **param_names**: integer  **intent(out)**    *size(10, 12)*  
        names of the parameters
            
  .. function:: number_of_bond_order_factors(n_bonds)

    Tells the number of differently named bond order factors the core knows
    
    Calls :func:`get_number_of_bond_order_factors`
    

    Parameters:

    **n_bonds**: integer  **intent(out)**    *scalar*  
        number of bond order factors
            
  .. function:: number_of_parameters_of_bond_order_factor(bond_name, n_targets, n_params)

    Tells how many numeric parameters a bond order factor incorporates
    
    Calls :func:`get_number_of_parameters_of_bond_order_factor`
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    n_targets: integer  *intent(in)*    *scalar*  
        number of targets
    **n_params**: integer  **intent(out)**    *scalar*  
        number of parameters
            
  .. function:: number_of_parameters_of_potential(pot_name, n_params)

    Tells how many numeric parameters a potential incorporates
    
    Calls :func:`get_number_of_parameters_of_potential`
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    **n_params**: integer  **intent(out)**    *scalar*  
        number of parameters
            
  .. function:: number_of_potentials(n_pots)

    Tells the number of differently named potentials the core knows
    
    Calls :func:`get_number_of_potentials`
    

    Parameters:

    **n_pots**: integer  **intent(out)**    *scalar*  
        number of potentials
            
  .. function:: number_of_targets_of_bond_order_factor(bond_name, n_target)

    Tells how many targets a bond order factor has, i.e., is it many-body
    
    Calls :func:`get_number_of_targets_of_bond_order_factor`
    

    Parameters:

    bond_name: character(len=*)  *intent(in)*    *scalar*  
        name of the bond order factor
    **n_target**: integer  **intent(out)**    *scalar*  
        number of targets
            
  .. function:: number_of_targets_of_potential(pot_name, n_target)

    Tells how many targets a potential has, i.e., is it a many-body potential
    
    Calls :func:`get_number_of_targets_of_potential`
    

    Parameters:

    pot_name: character(len=*)  *intent(in)*    *scalar*  
        name of the potential
    **n_target**: integer  **intent(out)**    *scalar*  
        number of targets
            
  .. function:: release()

    Deallocates all the arrays in the core
    
    Calls :func:`core_release_all_memory`

            
  .. function:: set_ewald_parameters(n_atoms, real_cut, k_radius, reciprocal_cut, sigma, epsilon, scaler)

    Sets the parameters for Ewald summation in the core.
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        
    real_cut: double precision  *intent(in)*    *scalar*  
        the real-space cutoff
    k_radius: double precision  *intent(in)*    *scalar*  
        
    reciprocal_cut: integer  *intent(in)*    *size(3)*  
        the k-space cutoffs
    sigma: double precision  *intent(in)*    *scalar*  
        the split parameter
    epsilon: double precision  *intent(in)*    *scalar*  
        electric constant
    scaler: double precision  *intent(in)*    *size(n_atoms)*  
        scaling factors for the individual charges
            
  .. function:: start_bond_order_factors()

    Initializes the bond order factors.
    A routine is called to generate descriptors for
    potentials. These descriptors are needed by the
    python interface in order to directly inquire
    the core on the types of factors available.
    
    Calls :func:`initialize_bond_order_factor_characterizers`

            
  .. function:: start_mpi()

    Initializes MPI for parallel calculations.
    
    Calls :func:`mpi_initialize`

            
  .. function:: start_potentials()

    Initializes the potentials.
    A routine is called to generate descriptors for
    potentials. These descriptors are needed by the
    python interface in order to directly inquire
    the core on the types of potentials available.
    
    Calls :func:`initialize_potential_characterizers`

            
  .. function:: start_rng(seed)

    Initialize Mersenne Twister random number generator.
    
    A seed number has to be given. In case we run in MPI
    mode, the master cpu will broadcast its seed to all other
    cpus to ensure that the random number sequences match
    in all the cpus.
    

    Parameters:

    seed: integer  *intent(in)*    *scalar*  
        a seed for the random number generator
            
  .. function:: sync_mpi()

    Syncs MPI.
    This just calls mpi_barrier, so it makes all cpus
    wait until everyone is at this particular point in
    execution.
    
    Calls :func:`mpi_sync`

            
  .. function:: update_atom_charges(n_atoms, charges)

    Updates the charges of existing atoms.
    This method does not allocate memory and so the atoms
    must already exist in the core.
    
    Calls :func:`core_update_atom_charges`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    charges: double precision  *intent(in)*    *size(n_atoms)*  
        new charges for the atoms
            
  .. function:: update_atom_coordinates(n_atoms, positions, momenta)

    Updates the positions and velocities of existing atoms.
    This method does not allocate memory and so the atoms
    must already exist in the core.
    
    Calls :func:`core_update_atom_coordinates`
    

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
    positions: double precision  *intent(in)*    *size(3, n_atoms)*  
        new coordinates for the atoms
    momenta: double precision  *intent(in)*    *size(3, n_atoms)*  
        new momenta for the atoms