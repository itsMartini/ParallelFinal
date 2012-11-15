module solutions_mod
  use variables_mod
  implicit none
  
contains
  
  complex(kind=8) function summand(t, w, j, N, tau, nu, sigma, gamma)
    double precision, intent(in) :: t
    complex(kind=8), external :: w
    integer, intent(in) :: j, N
    double precision, intent(in) :: tau, nu, sigma, gamma
    double precision :: var_eta, var_y
    complex(kind=8) :: var_z
    
    var_eta = eta(j, N)
    var_y = y(var_eta, tau)
    var_z = z(phi(var_y, gamma, nu), var_y, sigma)
    
    summand = mu_twiddle(mu(var_eta), phi_prime(var_y, nu), sigma)*&
         exp(var_z*t)*w(var_z)
  end function summand
  
  double precision function big_u(t, N, tau, nu, sigma, gamma, w)
    double precision, intent(in) :: t
    integer, intent(in) :: N
    double precision, intent(in) :: tau, nu, sigma, gamma
    complex(kind=8), external :: w
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

  double precision function epsilon(lu, bu)
    double precision, intent(in) :: lu, bu
    epsilon = dabs(lu-bu)
  end function epsilon

  double precision function rho(epsilon_N, epsilon_N2)
    double precision, intent(in) :: epsilon_N, epsilon_N2
    rho = dlog(epsilon_N2/epsilon_N)/dlog(2.d0)
  end function rho
end module solutions_mod
