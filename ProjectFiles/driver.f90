program driver
  use variables_mod
  use solutions_mod
  use w_mod
  use output_mod
  use ode_solver_mod
  implicit none
  
  double precision :: var_big_u1, little_u1, var_epsilon1
  
  !testing
  integer, parameter :: W_SIZE = 5
  complex(kind=8), dimension(1:W_SIZE-1) :: w_test
  double precision, dimension(1:W_SIZE-1) :: x_vect
  
  little_u1 = u_1(2.d0)
  var_big_u1 = big_u(2.d0, 40, 0.5d0, 0.5d0, 1.d0, 0.5d0, w_1)
  
  var_epsilon1 = abs(var_big_u1-little_u1)
  
  write (*,*) 'actual u: ', little_u1
  write (*,*) 'calculated u: ', var_big_u1
  write (*,*) 'epsilon: ', var_epsilon1
  
  !example 1 tables
  call table5()
  call table6()
  
  !example 2 table
  call table7()
  
  !example 3 table
  call table8()
  
  !test
  call fdiff_ode_solve(w_test, dcmplx(0.d0, 0.d0), fn, 0.d0, 1.d0, 0.d0, 1.d0, W_SIZE, 0.d0, .FALSE., x_vect)
  write (*,*)
  write (*,*)
  write (*,*) realpart(w_test)
end program driver
