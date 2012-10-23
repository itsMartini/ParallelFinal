program driver
  use variables_mod
  use solutions_mod
  use w_mod
  use output_mod
  implicit none
  
  double precision :: var_big_u1, little_u1, var_epsilon1

  little_u1 = u_1(2.d0)
  var_big_u1 = big_u(2.d0, 40, 0.5d0, 0.5d0, 1.d0, 0.5d0, w_1)

  var_epsilon1 = abs(var_big_u1-little_u1)

  write (*,*) 'actual u: ', little_u1
  write (*,*) 'calculated u: ', var_big_u1
  write (*,*) 'epsilon: ', var_epsilon1

  call table5()
end program driver
