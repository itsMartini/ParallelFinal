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

  integer, parameter :: N_ = 20, M_ = 80
  double precision, parameter :: T_ = 4.d0
  double precision, dimension(1:M_-1) :: var_big_u_vect, x_vect
  
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
42 format(3I6, 1ES16.7)
  open(unit = 42, file = 'results/runtimes.txt', status = 'unknown', action = 'write', access = 'append')
  
  time = mpi_wtime()
  
  var_big_u_vect = big_u_vect(T_, N_, M_, 0.5d0, 0.5d0, 1.d0, 0.5d0, x_vect,&
       g_3_z, 0.d0, PI, 0.d0, 0.d0, my_rank, num_cores, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, MPI_SUM)
  
  time = mpi_wtime()-time
  
  if (my_rank == MASTER) then
     write (42,42) N_, M_, num_cores, time
  end if
  !******************************************************************
  !end testing
  !******************************************************************

  call mpi_finalize(ierror)
end program
