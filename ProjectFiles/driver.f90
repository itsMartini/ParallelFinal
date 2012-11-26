program driver
  use variables_mod
  use solutions_mod
  use w_mod
  use output_mod
  use ode_solver_mod
  use mpi
  implicit none

  integer :: ierror, my_rank, num_cores

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
  call table8(my_rank, num_cores, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, MPI_SUM)

  !test
  !  call fdiff_ode_solve(w_test, dcmplx(0.d0, 0.d0), fn, 0.d0, 1.d0, 0.d0, 1.d0, W_SIZE, 0.d0, .FALSE., x_vect)
  !  write (*,*)
  !  write (*,*)
  !  write (*,*) realpart(w_test)

  call mpi_finalize(ierror)
end program
