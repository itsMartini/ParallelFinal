program driver
  use variables_mod
  use solutions_mod
  use w_mod
  use output_mod
  use ode_solver_mod
  use mpi
  implicit none
  
  integer :: ierror, my_rank, num_cores
  integer, parameter :: MASTER = 0
  
  integer :: n_, m_
  double precision, parameter :: T_ = 4.d0
  double precision, allocatable, dimension(:) :: var_big_u_vect, x_vect
  
  double precision :: time
  
  !testing
  !  integer, parameter :: W_SIZE = 5
  !  complex(kind=8), dimension(1:W_SIZE-1) :: w_test
  !  double precision, dimension(1:W_SIZE-1) :: x_vect
  
  call mpi_init(ierror)
  
  call mpi_comm_rank(mpi_comm_world, my_rank, ierror)
  call mpi_comm_size(mpi_comm_world, num_cores, ierror)
  
  !example 1 tables
  !  call table5()
  !  call table6()
  
  !example 2 table
  !  call table7()
  
  !example 3 table
  !  call table8(my_rank, num_cores, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, MPI_SUM)
  
  !test
  !  call fdiff_ode_solve(w_test, dcmplx(0.d0, 0.d0), fn, 0.d0, 1.d0, 0.d0, 1.d0, W_SIZE, 0.d0, .FALSE., x_vect)
  !  write (*,*)
  !  write (*,*)
  !  write (*,*) realpart(w_test)
  
  !******************************************************************
  !testing parallel speedup and efficiency
  !******************************************************************
  ! 42 format(3I7, 1ES16.7)
  
  !   if (my_rank == MASTER) then
  !      open(unit = 41, file = 'config.txt', status = 'old', action = 'read')
  !      read (41,*) n_, m_
  !      close(41)
  !   end if
  
  !   call mpi_bcast(n_, 1, MPI_INTEGER,&
  !        MASTER, MPI_COMM_WORLD, ierror)
  !   call mpi_bcast(m_, 1, MPI_INTEGER,&
  !        MASTER, MPI_COMM_WORLD, ierror)
  
  !   allocate(var_big_u_vect(1:m_-1), x_vect(1:m_-1))
  
  !   time = mpi_wtime()
  
  !   var_big_u_vect = big_u_vect(T_, n_, m_, 0.5d0, 0.5d0, 1.d0, 0.5d0, x_vect,&
  !        g_3_z, 0.d0, PI, 0.d0, 0.d0, my_rank, num_cores, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, MPI_SUM)
  
  !   time = mpi_wtime()-time
  
  !   deallocate(var_big_u_vect, x_vect)
  
  !   if (my_rank == MASTER) then
  !      open(unit = 42, file = 'results/runtimes.txt', status = 'unknown', action = 'write', access = 'append')
  !      write (42,42) n_, m_, num_cores, time
  !      close(42)
  !   end if
  
  !******************************************************************
  !end testing
  !******************************************************************
  
  open(unit=50, file="results/u_values.txt", iostat=ierror, status="unknown", action="write")
  if ( ierror /= 0 ) stop "Error opening results/u_values.txt"
  
  n_ = 16
  m_ = 50
  
  allocate(var_big_u_vect(1:m_-1), x_vect(1:m_-1))
  
  do i = 1, 380, 1
     time = .2d0 + dble(i-1)*3.8d0/380.d0
     var_big_u_vect = big_u_vect(time, n_, m_, 0.5d0, 0.5d0, 1.d0, 0.5d0, x_vect, &
          g_3_z, 0.d0, PI, 0.d0, 0.d0, my_rank, num_cores, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, MPI_SUM)
     
     write(50, *) time
     write(50, *) 0.d0, x_vect, PI
     write(50, *) 0.d0, var_big_u_vect, 0.d0
  end do
  
  deallocate(var_big_u_vect, x_vect)

  call mpi_finalize(ierror)
end program driver
