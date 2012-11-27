module output_mod
  use variables_mod
  use solutions_mod
  use w_mod
  use ode_solver_mod
  use mpi
  implicit none

  double precision :: t, little_u
  double precision, dimension(0:3) :: var_big_u, var_epsilon
  double precision, dimension(1:3) :: var_rho
  integer :: N, m, temp, i

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

    close(2)
  end subroutine table6

  subroutine table7()
6   format(1A2, 1A9, 2A11, 1A8, 1A11, 1A8, 1A11)
7   format(1F3.1, 2ES11.3, 1F8.2, 1ES11.3, 1F8.2, 1ES11.3, 1F8.2)

    open(unit=3, file='table7.txt', status='unknown', action='write')

    write (3,6) 't', 'e_10', 'e_20', 'p_20', 'e_40', 'p_40', 'e_80', 'p_80'

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
       
       write (3,7) t, var_epsilon(0), var_epsilon(1), var_rho(1), var_epsilon(2),&
            var_rho(2), var_epsilon(3), var_rho(3)
    end do
    
    close(3)
  end subroutine table7
  
  subroutine table8(my_rank, num_cores, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, MPI_SUM)
    integer, intent(in) :: my_rank, num_cores, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, MPI_SUM
    
    double precision, allocatable, dimension(:) :: var_big_u_vect, var_epsilon_vect
    double precision, allocatable, dimension(:) :: little_u_vect, x_vect
    integer :: j
    integer, parameter :: MASTER = 0
    
8   format(1A2, 1A9, 2A11, 1A8, 1A11, 1A8, 1A11)
9   format(1F3.1, 2ES11.3, 1F8.2, 1ES11.3, 1F8.2, 1ES11.3, 1F8.2)
    
    open(unit=4, file='table8.txt', status='unknown', action='write')
    
    write (4,8) 't', 'e_10', 'e_20', 'p_20', 'e_40', 'p_40', 'e_80', 'p_80'
    
    N = 20
    
    do temp = 1, 12
       if (temp < 11) then
          t = dble(temp)*0.2d0
       else
          t = dble(temp-8)
       end if
       
       do i = 0, 3
          m = 10*2**i
          allocate(little_u_vect(1:m-1), x_vect(1:m-1))
          allocate(var_big_u_vect(1:m-1), var_epsilon_vect(1:m-1))
          
          var_big_u_vect = big_u_vect(t, N, m, 0.5d0, 0.5d0, 1.d0, 0.5d0, x_vect,&
               g_3_z, 0.d0, PI, 0.d0, 0.d0, my_rank, num_cores, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, MPI_SUM)

          if (my_rank == MASTER) then
             little_u_vect = u_3(x_vect, t, m)

             do j = 1, m-1
                var_epsilon_vect(j) = epsilon(little_u_vect(j), var_big_u_vect(j))
             end do

             var_epsilon(i) = maxval(var_epsilon_vect)

             if (i > 0) then
                var_rho(i) = rho(var_epsilon(i), var_epsilon(i-1))
             end if
          end if

          deallocate(little_u_vect, x_vect, var_big_u_vect, var_epsilon_vect)
       end do

       if (my_rank == MASTER) then
          write (4,9) t, var_epsilon(0), var_epsilon(1), var_rho(1), var_epsilon(2),&
               var_rho(2), var_epsilon(3), var_rho(3)
       end if

    end do

    close(4)
  end subroutine table8

end module output_mod
