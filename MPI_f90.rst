
.. _mpi:
        
=====================================
mpi (MPI.f90)
=====================================



Module for Message Parsing Interface parallellization
utilities.

This module handles the initialization of the MPI environment
and assigns the cpus their indices.
Parallellization is done by
distributing atoms on the processors and the routine for
doing this randomly is provided. Also tools for monitoring
the loads of all the cpus and redistributing them are also implemented.

.. only:: html


    Modules used by mpi
    -------------------
    - :ref:`mt95`

    List of global variables in mpi
    -------------------------------
    - :data:`all_atoms`
    - :data:`all_loads`
    - :data:`atom_buffer`
    - :data:`cpu_id`
    - :data:`is_my_atom`
    - :data:`load_length`
    - :data:`loadout`
    - :data:`loads_mask`
    - :data:`mpi_atoms_allocated`
    - :data:`mpistat`
    - :data:`mpistatus`
    - :data:`my_atoms`
    - :data:`my_load`
    - :data:`n_cpus`
    - :data:`stopwatch`
    - :data:`track_loads`

    List of subroutines in mpi
    --------------------------
        
    - :func:`balance_loads`
    - :func:`close_loadmonitor`
    - :func:`initialize_load`
    - :func:`mpi_distribute`
    - :func:`mpi_finish`
    - :func:`mpi_initialize`
    - :func:`mpi_master_bcast_int`
    - :func:`mpi_sync`
    - :func:`mpi_wall_clock`
    - :func:`open_loadmonitor`
    - :func:`record_load`
    - :func:`start_timer`
    - :func:`timer`
    - :func:`write_loadmonitor`


Full documentation of global variables in mpi
---------------------------------------------
        
        
  .. data:: all_atoms

    integer    *scalar*    
    
    the total number of atoms
    
  .. data:: all_loads

    double precision  *allocatable*  *size(:)*    
    
    list of the loads of all cpus
    
  .. data:: atom_buffer

    integer  *allocatable*  *size(:)*    
    
    list used for passing atom indices during load balancing
    
  .. data:: cpu_id

    integer    *scalar*    
    
    identification number for the cpu, from 0 to :math:`n_\mathrm{cpus}-1`
    
  .. data:: is_my_atom

    logical  *allocatable*  *size(:)*    
    
    logical array, true for the indices of the atoms that are distributed to this cpu
    
  .. data:: load_length

    integer    *scalar*    
    
    the number of times loads have been recorded
    
  .. data:: loadout

    integer    *scalar*    

    *initial value* = 2352
    
    an integer of the output channel for loads
    
  .. data:: loads_mask

    logical  *allocatable*  *size(:)*    
    
    logical array used in load rebalancing, true for cpus whose loads have not yet been balanced
    
  .. data:: mpi_atoms_allocated

    logical    *scalar*    

    *initial value* = .false.
    
    logical switch for denoting that the mpi allocatable arrays have been allocated
    
  .. data:: mpistat

    integer    *scalar*    
    
    mpi return value
    
  .. data:: mpistatus

    integer    *size(mpi_status_size)*    
    
    array for storing the mpi status return values
    
  .. data:: my_atoms

    integer    *scalar*    
    
    the number of atoms distributed to this cpu
    
  .. data:: my_load

    double precision    *scalar*    
    
    storage for the load of this particular cpu
    
  .. data:: n_cpus

    integer    *scalar*    
    
    number of cpus, :math:`n_\mathrm{cpus}`
    
  .. data:: stopwatch

    double precision    *scalar*    
    
    cpu time storage
    
  .. data:: track_loads

    logical    *scalar*  *parameter*  

    *initial value* = .false.
    
    logical switch, if true, the loads of cpus are written to a file during run
    

Full documentation of subroutines in mpi
----------------------------------------
        
        
            
  .. function:: balance_loads()

    Load balancing.
    
    The loads are gathered from all cpus and sorted. Then load (atoms)
    is passed from the most loaded cpus to the least loaded ones.

            
  .. function:: close_loadmonitor()

    Closes the output for wirting workload data

            
  .. function:: initialize_load(reallocate)

    Initializes the load monitoring arrays.
    

    Parameters:

    reallocate: logical  *intent(in)*    *scalar*  
        Logical switch for reallocating the arrays. If true, the related arrays are allocated. Otherwise only the load counters are set to zero.
            
  .. function:: mpi_distribute(n_atoms)

    distributes atoms among processors

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        number of atoms
            
  .. function:: mpi_finish()

    closes the mpi framework

            
  .. function:: mpi_initialize()

    intializes the mpi framework

            
  .. function:: mpi_master_bcast_int(sync_int)

    the master cpu broadcasts an integer value to all other cpus

    Parameters:

    **sync_int**: integer  **intent(inout)**    *scalar*  
        the broadcast integer
            
  .. function:: mpi_sync()

    syncs the cpus by calling mpi_barrier

            
  .. function:: mpi_wall_clock(clock)

    returns the global time through mpi_wtime

    Parameters:

    **clock**: double precision  **intent(out)**    *scalar*  
        the measured time
            
  .. function:: open_loadmonitor()

    Opens the output for writing workload data to a file called "mpi_load.out"

            
  .. function:: record_load(amount)

    Saves the given load.
    

    Parameters:

    amount: double precision  *intent(in)*    *scalar*  
        the load to be stored
            
  .. function:: start_timer()

    records the wall clock time to :data:`stopwatch`

            
  .. function:: timer(stamp)

    reads the elapsed wall clock time since the previous starting of the
    timer (saved in :data:`stopwatch`)
    and then restarts the timer.

    Parameters:

    **stamp**: double precision  **intent(inout)**    *scalar*  
        the elapsed real time
            
  .. function:: write_loadmonitor()

    Routine for writing force calculation workload analysis data to a file
    called "mpi_load.out"
