module solutions_mod
  use variables_mod
  use ode_solver_mod
  implicit none
  
contains
  
  !w is a function
  function summand(t, w, j, N, tau, nu, sigma, gamma)
    double precision, intent(in) :: t
    integer, intent(in) :: j, N
    complex(kind=8), external :: w
    double precision, intent(in) :: tau, nu, sigma, gamma
    complex(kind=8) :: summand

    double precision :: var_eta, var_y
    complex(kind=8) :: var_z

    var_eta = eta(j, N)
    var_y = y(var_eta, tau)
    var_z = z(phi(var_y, gamma, nu), var_y, sigma)

    summand = mu_twiddle(mu(var_eta), phi_prime(var_y, nu), sigma)*&
         exp(var_z*t)*w(var_z)
  end function summand

  !w is a function
  function big_u(t, N, tau, nu, sigma, gamma, w)
    double precision, intent(in) :: t
    integer, intent(in) :: N
    double precision, intent(in) :: tau, nu, sigma, gamma
    complex(kind=8), external :: w
    double precision :: big_u

    complex(kind=8) :: temp_u

    integer :: j

    temp_u = dcmplx(0, 0)

    temp_u = temp_u+summand(t, w, 0, N, tau, nu, sigma, gamma)/2.d0
    do j = 1, N-1
       temp_u = temp_u+summand(t, w, j, N, tau, nu, sigma, gamma)
    end do

    temp_u = temp_u/(N*tau)
    big_u = 2.d0*realpart(temp_u)
  end function big_u

  !w is a vector
  function summand_vect(t, w, j, N, m, tau, nu, sigma, gamma, var_eta, var_y, var_z)
    double precision, intent(in) :: t
    integer, intent(in) :: j, N, m
    complex(kind=8), dimension(1:m-1), intent(in) :: w
    double precision, intent(in) :: tau, nu, sigma, gamma
    double precision, intent(in) :: var_eta, var_y
    complex(kind=8), intent(in) :: var_z
    complex(kind=8), dimension(1:m-1) :: summand_vect

    summand_vect = mu_twiddle(mu(var_eta), phi_prime(var_y, nu), sigma)*&
         exp(var_z*t)*w
  end function summand_vect

  !w is a vector
  function big_u_vect(t, N, m, tau, nu, sigma, gamma, x_vect,&
       g, a, b, alpha, beta, my_rank, num_cores, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, MPI_SUM)
    double precision, intent(in) :: t
    integer, intent(in) :: N, m
    double precision, intent(in) :: tau, nu, sigma, gamma, a, b, alpha, beta
    double precision, dimension(1:m-1), intent(inout) :: x_vect
    complex(kind=8), external :: g
    double precision, dimension(1:m-1) :: big_u_vect

    !mpi variables
    integer, parameter :: MASTER = 0
    integer :: ierror
    integer, intent(in) :: my_rank, num_cores, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, MPI_SUM
    
    complex(kind=8), dimension(1:m-1) :: temp_u, temp_u_recv, w
    
    double precision :: var_eta, var_y
    complex(kind=8) :: var_z
    
    integer :: j, start, finish
    
    start = my_rank*(N-1)/num_cores + 1
    finish = (my_rank+1)*(N-1)/num_cores
    
    temp_u = dcmplx(0.d0, 0.d0)
    
    !master calculates zeroth summand
    if (my_rank == MASTER) then
       var_eta = eta(0, N)
       var_y = y(var_eta, tau)
       var_z = z(phi(var_y, gamma, nu), var_y, sigma)
       
       call fdiff_ode_solve(w, var_z, g, a, b, alpha, beta, m, t, .TRUE., x_vect)
       
       temp_u = temp_u+summand_vect(t, w, 0, N, m, tau, nu, sigma, gamma, var_eta, var_y, var_z)/2.d0
    end if
    
    !each core computes part of the total sum
    do j = start, finish
       var_eta = eta(j, N)
       var_y = y(var_eta, tau)
       var_z = z(phi(var_y, gamma, nu), var_y, sigma)
       
       call fdiff_ode_solve(w, var_z, g, a, b, alpha, beta, m, t, .TRUE., x_vect)
       
       temp_u = temp_u+summand_vect(t, w, j, N, m, tau, nu, sigma, gamma, var_eta, var_y, var_z)
    end do
    
    temp_u = temp_u/(N*tau)
    
    !master reduces the sum
    call mpi_reduce(temp_u, temp_u_recv, m-1, MPI_DOUBLE_COMPLEX,&
         MPI_SUM, MASTER, MPI_COMM_WORLD, ierror)
    
    big_u_vect = 2.d0*realpart(temp_u_recv)
  end function big_u_vect
  
  double precision function epsilon(lu, bu)
    double precision, intent(in) :: lu, bu
    epsilon = dabs(lu-bu)
  end function epsilon
  
  double precision function rho(epsilon_N, epsilon_N2)
    double precision, intent(in) :: epsilon_N, epsilon_N2
    rho = dlog(epsilon_N2/epsilon_N)/dlog(2.d0)
  end function rho
end module solutions_mod
