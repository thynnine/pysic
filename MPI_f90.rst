
.. _mpi:
        
=====================================
mpi (MPI.f90)
=====================================



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
    
    
    
  .. data:: all_loads

    double precision  *allocatable*  *size(:)*    
    
    
    
  .. data:: atom_buffer

    integer  *allocatable*  *size(:)*    
    
    
    
  .. data:: cpu_id

    integer    *scalar*    
    
    
    
  .. data:: is_my_atom

    logical  *allocatable*  *size(:)*    
    
    
    
  .. data:: load_length

    integer    *scalar*    
    
    
    
  .. data:: loadout

    integer    *scalar*    

    *initial value* = 2352
    
    
    
  .. data:: loads_mask

    logical  *allocatable*  *size(:)*    
    
    
    
  .. data:: mpi_atoms_allocated

    logical    *scalar*    

    *initial value* = .false.
    
    
    
  .. data:: mpistat

    integer    *scalar*    
    
    
    
  .. data:: mpistatus

    integer    *size(mpi_status_size)*    
    
    
    
  .. data:: my_atoms

    integer    *scalar*    
    
    
    
  .. data:: my_load

    double precision    *scalar*    
    
    
    
  .. data:: n_cpus

    integer    *scalar*    
    
    
    
  .. data:: stopwatch

    double precision    *scalar*    
    
    
    

Full documentation of subroutines in mpi
----------------------------------------
        
        
            
  .. function:: balance_loads()


            
  .. function:: close_loadmonitor()

    Closes the output for wirting workload data

            
  .. function:: initialize_load(reallocate)


    Parameters:

    reallocate: logical  *intent(in)*    *scalar*  
        
            
  .. function:: mpi_distribute(n_atoms)

    distributes atoms among processors

    Parameters:

    n_atoms: integer  *intent(in)*    *scalar*  
        
            
  .. function:: mpi_finish()

    closes the mpi framework

            
  .. function:: mpi_initialize()

    intializes the mpi framework

            
  .. function:: mpi_master_bcast_int(sync_int)

    the master cpu broadcasts an integer value to all other cpus

    Parameters:

    **sync_int**: integer  **intent(inout)**    *scalar*  
        
            
  .. function:: mpi_sync()

    syncs the cpus by calling mpi_barrier

            
  .. function:: mpi_wall_clock(clock)

    returns the global time through mpi_wtime

    Parameters:

    **clock**: double precision  **intent(out)**    *scalar*  
        
            
  .. function:: open_loadmonitor()

    Opens the output for writing workload data to a file called "mpi_load.out"

            
  .. function:: record_load(amount)


    Parameters:

    amount: double precision  *intent(in)*    *scalar*  
        
            
  .. function:: start_timer()

    records the wall clock time to stopwatch

            
  .. function:: timer(stamp)

    reads the elapsed wall clock time since the previous starting of the timer
    and then restarts the timer.

    Parameters:

    stamp: double precision  *intent()*    *scalar*  
        
            
  .. function:: write_loadmonitor()

    Routine for writing force calculation workload analysis data to a file
