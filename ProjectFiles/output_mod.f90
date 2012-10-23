module output_mod
  use variables_mod
  use solutions_mod
  use w_mod
  implicit none
  
  double precision :: t, little_u
  double precision, dimension(0:3) :: var_big_u, var_epsilon
  double precision, dimension(1:3) :: var_rho
  integer :: N, temp, i

contains

  subroutine table5()
1   format(1A2, 1A9, 2A11, 1A8, 1A11, 1A8, 1A11)
2   format(1F3.1, 2ES11.3, 1F8.2, 1ES11.3, 1F8.2, 1ES11.3, 1F8.2)

    open(unit=1, file='table5.txt', status='unknown', action='write')

    write (1,1) 't', 'e_10', 'e_20', 'p_20', 'e_40', 'p_40', 'e_80', 'p_80'

    do temp = 1, 12
       if (temp < 11) then
          t = dble(temp)*0.2d0
       else
          t = dble(temp-8)
       end if

       little_u = u_1(t)

       do i = 0, 3
          N = 10*2**i
          var_big_u(i) = big_u(t, N, 0.5d0, 0.5d0, 1.d0, 0.5d0, w_1)
          var_epsilon(i) = epsilon(little_u, var_big_u(i))
          if (i > 0) then
             var_rho(i) = rho(var_epsilon(i), var_epsilon(i-1))
          end if
       end do

       write (1,2) t, var_epsilon(0), var_epsilon(1), var_rho(1), var_epsilon(2),&
            var_rho(2), var_epsilon(3), var_rho(3)
    end do

    close(1)
  end subroutine table5

  subroutine table6()
3   format(1A24, 1A27)
4   format(1A2, 1A9, 2A11, 1A8, 1A11, 1A11)
5   format(1F3.1, 2ES11.3, 1F8.2, 2ES11.3, 1F8.2)
  
    open(unit=2, file='table6.txt', status='unknown', action='write')
    
    write (2,3) 'Tau = 0.25', 'Tau = 1'
    write (2,4) 't', 'e_40', 'e_80', 'p_80', 'e_40', 'e_80', 'p_80'
    
    do temp = 1, 12
       if (temp < 11) then
          t = dble(temp)*0.2d0
       else
          t = dble(temp-8)
       end if
       
       little_u = u_1(t)
       
       do i = 2, 3
          N = 10*2**i
          var_big_u(mod(i, 2)) = big_u(t, N, 0.25d0, 0.5d0, 1.d0, 0.5d0, w_1)
          var_epsilon(mod(i, 2)) = epsilon(little_u, var_big_u(mod(i, 2)))
          if (mod(i, 2) > 0) then
             var_rho(mod(i, 2)) = rho(var_epsilon(mod(i,2)),&
                  var_epsilon(mod(i, 2)-1))
          end if
          
          var_big_u(i) = big_u(t, N, 1.d0, 0.5d0, 1.d0, 0.5d0, w_1)
          var_epsilon(i) = epsilon(little_u, var_big_u(i))
          if (i > 2) then
             var_rho(i) = rho(var_epsilon(i), var_epsilon(i-1))
          end if
       end do

       write (2,5) t, var_epsilon(0), var_epsilon(1), var_rho(1),&
            var_epsilon(2), var_epsilon(3), var_rho(3)
    end do

  end subroutine table6
end module output_mod
