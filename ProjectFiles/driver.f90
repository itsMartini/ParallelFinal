program driver
  use variables_mod
  use solutions_mod
  use w_mod
  implicit none
  
  double precision :: var_big_u, little_u, epsilon
  
  little_u = u_1(2.d0)
  var_big_u = big_u(2.d0, 40, 0.5d0, 0.5d0, 1.d0, 0.5d0, w_1)

  epsilon = abs(var_big_u-little_u)

  write (*,*) 'actual u: ', little_u
  write (*,*) 'calculated u: ', var_big_u
  write (*,*) 'epsilon: ', epsilon
end program driver
