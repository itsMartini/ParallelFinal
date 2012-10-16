module variables_mod
  implicit none

  double precision :: pi
  pi = acos(-1.d0)
  
contains

  double precision function eta(j, N)
    integer, intent(in) :: j, N
    eta = dble(j)/dble(N)
  end function eta

  double precision function chi(eta)
    double precision, intent(in) :: eta
    chi = log((1.d0+eta)/(1.d0-eta))
  end function chi

  double precision function mu(eta)
    double precision, intent(in) :: eta
    mu = 2.d0/(1.d0-eta**2.d0)
  end function mu

  double precision function y(eta, tau)
    double precision, intent(in) :: eta, tau
    y = 1.d0/tau*chi(eta)
  end function y

  double precision function phi(y, gamma, nu)
    double precision, intent(in) :: y, gamma, nu
    phi = gamma-sqrt(y**2.d0+nu**2.d0)
  end function phi
  
  double precision function phi_prime(y, nu)
    double precision, intent(in) :: y, nu
    phi_prime = -2.d0*y/sqrt(y**2.d0+nu**2.d0)
  end function phi_prime
  
  complex(kind=8) function z(phi, y, sigma)
    double precision, intent(in) :: phi, y, sigma
    z = (phi, sigma*y)
  end function z
  
  complex(kind=8) function mu_twiddle(mu, phi_prime, sigma)
    double precision, intent(in) :: mu, phi_prime, sigma
    mu_twiddle = (sigma*mu/(2.d0*pi), -1.d0*phi_prime*mu/(2.d0*pi))
  end function mu_twiddle
  
end module variables_mod
