program driver
  use variables_mod
  use solutions_mod
  use w_mod
  use output_mod
  use ode_solver_mod
  implicit none
  
  double precision :: var_big_u1, little_u1, var_epsilon1
  
  !testing
  integer, parameter :: W_SIZE = 100
  complex(kind=8), dimension(1:W_SIZE-1) :: w
  
  little_u1 = u_1(2.d0)
  var_big_u1 = big_u(2.d0, 40, 0.5d0, 0.5d0, 1.d0, 0.5d0, w_1)
  
  var_epsilon1 = abs(var_big_u1-little_u1)
  
  write (*,*) 'actual u: ', little_u1
  write (*,*) 'calculated u: ', var_big_u1
  write (*,*) 'epsilon: ', var_epsilon1
  
  call table5()
  call table6()
  
  !testing ode_solver
  call fdiff_ode_solve(w, dcmplx(0.d0,0.d0), fun, 0.d0, 1.d0, 0.d0, 4.d0, W_SIZE)
  write (*,*)
  write (*,*) 
  write (*,*) realpart(w)

end program driver
