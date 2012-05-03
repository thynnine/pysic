!
! Module for Message Parsing Interface parallellization
! utilities.
!
! This module handles the initialization of the MPI environment
! and assigns the cpus their indices. 
! Parallellization is done by 
! distributing atoms on the processors and the routine for
! doing this randomly is provided. Also tools for monitoring
! the loads of all the cpus and redistributing them are also implemented.
module mpi
  use mt95
  implicit none

#ifdef MPI
  include 'mpif.h'
#endif

  ! *cpu_id identification number for the cpu, from 0 to :math:`n_\mathrm{cpus}-1`
  ! *n_cpus number of cpus, :math:`n_\mathrm{cpus}`
  ! *loadout an integer of the output channel for loads
  integer :: cpu_id, n_cpus, loadout = 2352
#ifdef MPI
  ! *mpistatus array for storing the mpi status return values
  ! *mpistat mpi return value
  integer :: mpistatus(mpi_status_size), mpistat
#endif
  ! *my_atoms the number of atoms distributed to this cpu
  ! *all_atoms the total number of atoms
  integer :: my_atoms, all_atoms
  ! *load_length the number of times loads have been recorded
  integer :: load_length
  ! *is_my_atom logical array, true for the indices of the atoms that are distributed to this cpu
  ! *loads_mask logical array used in load rebalancing, true for cpus whose loads have not yet been balanced
  logical, allocatable :: is_my_atom(:), loads_mask(:)
  ! *atom_buffer list used for passing atom indices during load balancing
  integer, allocatable :: atom_buffer(:)
  ! *mpi_atoms_allocated logical switch for denoting that the mpi allocatable arrays have been allocated
  logical :: mpi_atoms_allocated = .false.
  ! *track_loads logical switch, if true, the loads of cpus are written to a file during run
  logical, parameter :: track_loads = .false.
  ! *stopwatch cpu time storage
  ! *my_load storage for the load of this particular cpu
  double precision :: stopwatch, my_load
  ! *all_loads list of the loads of all cpus
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

    if(track_loads)then
       call open_loadmonitor()
    end if

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
    
    if(track_loads)then
       call close_loadmonitor()
    end if

  end subroutine mpi_finish


  ! syncs the cpus by calling mpi_barrier
  subroutine mpi_sync()
    implicit none

#ifdef MPI
    call mpi_barrier(mpi_comm_world, mpistat)
#endif

  end subroutine mpi_sync


  ! returns the global time through mpi_wtime
  ! *clock the measured time
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
  ! *sync_int the broadcast integer
  subroutine mpi_master_bcast_int(sync_int)
    implicit none
    integer, intent(inout) :: sync_int

#ifdef MPI
    call mpi_bcast(sync_int,1,mpi_integer,0,mpi_comm_world,mpistat)
#endif
    
  end subroutine mpi_master_bcast_int



  ! distributes atoms among processors
  ! *n_atoms number of atoms
  subroutine mpi_distribute(n_atoms)
    implicit none
    integer, intent(in) :: n_atoms
    integer :: i, random_int
    double precision :: random_real

    if(.not.mpi_atoms_allocated)then
       allocate(is_my_atom(n_atoms))
       mpi_atoms_allocated = .true.
    else
       deallocate(is_my_atom)
       allocate(is_my_atom(n_atoms))
    end if
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


  ! Initializes the load monitoring arrays.
  !
  ! *reallocate Logical switch for reallocating the arrays. If true, the related arrays are allocated. Otherwise only the load counters are set to zero.
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


  ! Saves the given load.
  !
  ! *amount the load to be stored
  subroutine record_load(amount)
    implicit none
    double precision, intent(in) :: amount

    my_load = my_load + amount
    load_length = load_length + 1

  end subroutine record_load


  ! Load balancing.
  !
  ! The loads are gathered from all cpus and sorted. Then load (atoms)
  ! is passed from the most loaded cpus to the least loaded ones.
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

    if(track_loads)then
       call write_loadmonitor()
    end if

#endif

  end subroutine balance_loads


  ! records the wall clock time to :data:`stopwatch`
  subroutine start_timer()
    implicit none

    call mpi_wall_clock(stopwatch)

  end subroutine start_timer


  ! reads the elapsed wall clock time since the previous starting of the 
  ! timer (saved in :data:`stopwatch`)
  ! and then restarts the timer.
  ! *stamp the elapsed real time
  subroutine timer(stamp)
    implicit none
    double precision, intent(inout) :: stamp

    stamp = stopwatch       ! previous time
    call mpi_wall_clock(stopwatch)
    stamp = stopwatch-stamp ! time difference

  end subroutine timer






  ! Routine for writing force calculation workload analysis data to a file
  ! called "mpi_load.out"
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



#ifdef MPI
  ! stacks the "lists" from all cpus together according to the lengths given in "items"
  ! and gathers the complete list to cpu 0
  ! For example::
  ! 
  ! cpu 0          cpu 1          cpu 0
  ! abc....        12.....        abc12..
  ! de.....        3456...   ->   de3456.
  ! fghij..        78.....        fghij78
  !
  ! The stacking is done for the second array index: list(1,:,1).
  ! The stacking works so that first every cpu 2n+1 sends its data to cpu 2n,
  ! then 2*(2n+1) sends data to 2*2n, and so on, until the final cpu 2^m sends its data to cpu 0::
  !
  !  cpu
  !  0 1 2 3 4 5 6 7 8 9 10
  !  |-/ |-/ |-/ |-/ |-/ |
  !  |---/   |---/   |---/
  !  |-------/       |
  !  |---------------/
  !  x
  !
  ! *list 3d arrays containing lists to be stacked 
  ! *items the numbers of items to be stacked in each list
  ! *length the number of lists (size of list(1,1,:))
  ! *width max size of lists (size of list(1,:,1))
  ! *depth dimensionality of the stacked objects (size of list(:,1,1))
  SUBROUTINE mpi_stack(list,items,depth,length,width)
    IMPLICIT NONE
    INTEGER, POINTER :: list(:,:,:), items(:)
    INTEGER, INTENT(IN) :: length,width,depth
    INTEGER :: remainder, templist(width,length),tempitems(length), ii,level
    LOGICAL :: fine

    fine = .false.
    level = 1 ! level of communications: cpu level*(2n+1) sends to cpu level*n

    DO WHILE(.not.fine) ! "fine" is a marker for ending the receiving-sending loop for this cpu

       level = 2*level ! advance to next level
       remainder = cpu_id - (cpu_id/level)*level ! the remainder for cpu_id/level
       
       ! level gets values 2^k, so the remainders develop as follows:
       ! level  cpu/remainder
       !        0 1 2 3 4 5 6 7 8 9
       ! 2      0 1 0 1 0 1 0 1 0 1
       ! 4      0 1 2 3 0 1 2 3 0 1
       ! 8      0 1 2 3 4 5 6 7 0 1
       ! 16     0 1 2 3 4 5 6 7 8 9
       !  ...
       ! So, the cpu should receive or wait when the remainder is zero
       ! and send and finish once the remainder becomes positive

       IF(remainder == 0)THEN
          IF(n_cpus > cpu_id + level/2 )THEN ! is there a cpu that should be sending?
             ! receive the lists and the numbers of elements in the lists
             CALL MPI_RECV(templist,length*width*depth,MPI_INTEGER,cpu_id+level/2,&
                  5501+level,MPI_COMM_WORLD,mpistatus,mpistat)
             CALL MPI_RECV(tempitems,length,MPI_INTEGER,cpu_id+level/2,&
                  5502+level,MPI_COMM_WORLD,mpistatus,mpistat)

             DO ii = 1, length ! stack the lists
                IF(tempitems(ii) > 0)THEN
                   list(1:depth,items(ii)+1:items(ii)+tempitems(ii),ii) = templist(1:depth,1:tempitems(ii),ii)
                   items(ii) = items(ii) + tempitems(ii)
                END IF
             END DO             
          ELSE ! there are no more cpus that would send data to this one
             IF(cpu_id /= 0)THEN
                ! since this is not the master cpu, it must wait until it is its turn to send its data forward
                ! so, skip this level and repeat the loop (will advance the level)
             ELSE
                fine = .true. ! the master cpu has received all the data, finish now
             END IF
          END IF
       ELSE ! the first time the remainder is more than zero, send the data and finish for this cpu
          CALL MPI_SEND(list,length*width*depth,MPI_INTEGER,cpu_id-level/2,&
               5501+level,MPI_COMM_WORLD,mpistat)
          CALL MPI_SEND(items,length,MPI_INTEGER,cpu_id-level/2,&
               5502+level,MPI_COMM_WORLD,mpistat)
          fine = .true.
       END IF
    END DO

    RETURN

  END SUBROUTINE mpi_stack
#endif



end module mpi
