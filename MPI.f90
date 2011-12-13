module mpi
  use mt95
  implicit none

#ifdef MPI
  include 'mpif.h'
#endif

  integer :: cpu_id, n_cpus, loadout = 2352
#ifdef MPI
  integer :: mpistatus(mpi_status_size), mpistat
#endif
  integer :: my_atoms, all_atoms
  integer :: load_length
  logical, allocatable :: is_my_atom(:), loads_mask(:)
  integer, allocatable :: atom_buffer(:)
  double precision :: stopwatch, my_load
  double precision, allocatable :: all_loads(:)


contains
  
  ! intializes the mpi framework
  subroutine mpi_initialize()
    implicit none
    
#ifdef MPI
    call mpi_init(mpistat)
    if(mpistat /= mpi_success)then
       write(*,*) "mpi error I1"
       stop
    end if
    call mpi_comm_size(mpi_comm_world, n_cpus, mpistat)
    if(mpistat /= mpi_success)then
       write(*,*) "mpi error I2"
       stop
    end if
    call mpi_comm_rank(mpi_comm_world, cpu_id, mpistat)
    if(mpistat /= mpi_success)then
       write(*,*) "mpi error I3"
       stop
    end if
#else
    n_cpus = 1
    cpu_id = 0
#endif
    call open_loadmonitor()

  end subroutine mpi_initialize


  ! closes the mpi framework
  subroutine mpi_finish()
    implicit none

#ifdef MPI
    call mpi_finalize(mpistat)
    if(mpistat /= mpi_success)then
       write(*,*) "mpi error F1"
       stop
    end if
#endif
    call close_loadmonitor()

  end subroutine mpi_finish


  ! syncs the cpus by calling mpi_barrier
  subroutine mpi_sync()
    implicit none

#ifdef MPI
    call mpi_barrier(mpi_comm_world, mpistat)
#endif

  end subroutine mpi_sync


  ! returns the global time through mpi_wtime
  subroutine mpi_wall_clock(clock)
    implicit none
    double precision, intent(out) :: clock

#ifdef MPI
    clock = mpi_wtime()
#else
    call cpu_time(clock)
#endif

  end subroutine mpi_wall_clock


  ! the master cpu broadcasts an integer value to all other cpus
  subroutine mpi_master_bcast_int(sync_int)
    implicit none
    integer, intent(inout) :: sync_int

#ifdef MPI
    call mpi_bcast(sync_int,1,mpi_integer,0,mpi_comm_world,mpistat)
#endif
    
  end subroutine mpi_master_bcast_int



  ! distributes atoms among processors
  subroutine mpi_distribute(n_atoms)
    implicit none
    integer, intent(in) :: n_atoms
    integer :: i, random_int
    double precision :: random_real

    allocate(is_my_atom(n_atoms))
    my_atoms = 0
    all_atoms = n_atoms
    is_my_atom = .false.
    
    do i = 1, n_atoms

       call genrand_real2(random_real)
       random_int = floor(random_real*n_cpus)
       if(random_int == cpu_id)then
          is_my_atom(i) = .true.
          my_atoms = my_atoms + 1
       end if

    end do

    call initialize_load(.true.)

  end subroutine mpi_distribute



  subroutine initialize_load(reallocate)
    implicit none
    logical, intent(in) :: reallocate
    
    my_load = 0.d0
    load_length = 0
    if(reallocate)then
       if(allocated(all_loads))then
          deallocate(all_loads)
       end if
       if(allocated(loads_mask))then
          deallocate(loads_mask)
       end if
       if(allocated(atom_buffer))then
          deallocate(atom_buffer)
       end if
       allocate(all_loads(0:n_cpus-1))
       allocate(loads_mask(0:n_cpus-1))
       allocate(atom_buffer(all_atoms))
    end if

  end subroutine initialize_load



  subroutine record_load(amount)
    implicit none
    double precision, intent(in) :: amount

    my_load = my_load + amount
    load_length = load_length + 1

  end subroutine record_load



  subroutine balance_loads()
    implicit none
    integer :: i, j, k, min_cpu, max_cpu, swap_count, random_skip
    double precision :: random

#ifdef MPI
    if(load_length < 20)then
       return
    end if

    ! gather the load estimates from all proc
    call mpi_allgather(my_load,1,mpi_double_precision,all_loads,1,mpi_double_precision,mpi_comm_world,mpistat)

    loads_mask=.true.

    ! loop over pairs of most and least loaded cpus
    do i = 1, n_cpus/2-mod(n_cpus,2)
       
       ! find the most and least loaded cpus of the ones still left
       min_cpu = minloc(all_loads,1,loads_mask)-1
       loads_mask(min_cpu) = .false.
       max_cpu = maxloc(all_loads,1,loads_mask)-1
       loads_mask(max_cpu) = .false.

       ! if load imbalance between the two processes less than 2.5%, we do not swap 
       if( (all_loads(max_cpu) - all_loads(min_cpu)) / &
            (0.5*(all_loads(max_cpu)+all_loads(min_cpu))) < 0.025) then
          exit
       end if


       ! a low load cpu takes up more tasks
       if(cpu_id == min_cpu)then

          ! get a list of atoms from a high load cpu specifying which atoms
          ! are now allocated to this cpu
          call mpi_recv(atom_buffer,all_atoms,mpi_integer,max_cpu,10,mpi_comm_world,mpistatus,mpistat)
          call mpi_get_count(mpistatus,mpi_integer,swap_count,mpistat)

          do j = 1, swap_count
             is_my_atom(atom_buffer(j)) = .true.
             my_atoms = my_atoms + 1
          end do
       end if

       ! a high load cpu gives away tasks
       if(cpu_id == max_cpu)then
          
          ! estimate the number of atoms that should be allocated to another cpu
          swap_count=max( 1, nint(0.5*my_atoms*( (all_loads(max_cpu)-all_loads(min_cpu)))/all_loads(max_cpu) ) )

          ! pick atoms to swap
          do j = 1, swap_count
             
             ! go through the list of atoms but skip some so that a random set of
             ! atoms is sent instead of a bunch
             call genrand_real1(random)
             random_skip = floor(random*my_atoms)

             do k = 1, all_atoms
                if(is_my_atom(k))then

                   ! we skipped enough atoms, so this one is chosen
                   ! for sending off - it's added to the list of atoms
                   if(random_skip == 0)then
                      atom_buffer(j) = k
                      is_my_atom(k) = .false.
                      my_atoms = my_atoms - 1
                      exit
                   end if
                   random_skip = random_skip - 1
                end if
             end do
             
          end do

          ! send the list of atoms to a low load cpu
          call mpi_send(atom_buffer,swap_count,mpi_integer,min_cpu,10,mpi_comm_world,mpistat)

       end if

    end do
    
    ! set all loads to zero
    call initialize_load(.false.)
    call write_loadmonitor()

#endif

  end subroutine balance_loads


  ! records the wall clock time to stopwatch
  subroutine start_timer()
    implicit none

    call mpi_wall_clock(stopwatch)

  end subroutine start_timer


  ! reads the elapsed wall clock time since the previous starting of the timer
  ! and then restarts the timer.
  ! *timer the elapsed real time
  subroutine timer(stamp)
    implicit none
    double precision :: stamp

    stamp = stopwatch       ! previous time
    call mpi_wall_clock(stopwatch)
    stamp = stopwatch-stamp ! time difference

  end subroutine timer






  ! Routine for writing force calculation workload analysis data to a file
  SUBROUTINE write_loadmonitor()
    IMPLICIT NONE
    INTEGER :: ii

    IF(cpu_id == 0)THEN
       DO ii = 0, n_cpus-1
          WRITE(loadout,'(I4,ES9.2,A)',advance='no') ii, all_loads(ii), ", "
       END DO
       WRITE(loadout,'(A,ES9.2,F6.1)',advance='yes') " max: ", MAXVAL(all_loads(:)), &
            n_cpus*100*MAXVAL(all_loads(:))/SUM(all_loads(:))
    END IF

  END SUBROUTINE write_loadmonitor

  ! Opens the output for writing workload data to a file called "mpi_load.out"
  SUBROUTINE open_loadmonitor()
    IMPLICIT NONE

    IF(cpu_id == 0)THEN
       OPEN(loadout,file="mpi_load.out")
    END IF

  END SUBROUTINE open_loadmonitor
  
  ! Closes the output for wirting workload data
  SUBROUTINE close_loadmonitor()
    IMPLICIT NONE

    IF(cpu_id == 0)THEN
       CLOSE(loadout)
    END IF
    
  END SUBROUTINE close_loadmonitor  




end module mpi
