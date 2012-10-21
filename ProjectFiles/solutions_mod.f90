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
  end function big_u
end module solutions_mod
