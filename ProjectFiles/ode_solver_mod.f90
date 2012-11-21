module ode_solver_mod
  implicit none

contains

  subroutine generate_fdiff_matrix(m, z, DL, D, DU, ONE_OVER_H_SQUARED)
    integer, intent(in) :: m
    complex(kind=8), intent(in) :: z
    complex(kind=8), dimension(1:m-2), intent(inout) :: DL, DU
    complex(kind=8), dimension(1:m-1), intent(inout) :: D
    double precision, intent(in) :: ONE_OVER_H_SQUARED

    double precision, parameter :: ONE = 1.d0, TWO = 2.d0

    DL = -ONE*ONE_OVER_H_SQUARED
    DU = -ONE*ONE_OVER_H_SQUARED
    D = TWO*ONE_OVER_H_SQUARED+z
  end subroutine generate_fdiff_matrix

  function generate_g(g, a, H, m)
    double precision, external :: g
    double precision, intent(in) :: a, H
    integer, intent(in) :: m

    double precision, dimension(1:m-1) :: generate_g
    integer :: i

    do i = 1, m-1
       generate_g(i) = g(a+h*i)
    end do
  end function generate_g

  double precision function fun(x)
    double precision, intent(in) :: x
    fun = 0
  end function fun

  subroutine fdiff_ode_solve(w, z, g, a, b, alpha, beta, m)
    integer, intent(in) :: m
    complex(kind=8), dimension(1:m-1), intent(inout) :: w
    complex(kind=8), intent(in) :: z 
    double precision, external :: g
    double precision, intent(in) :: a, b, alpha, beta

    complex(kind=8), dimension(1:m-1) :: D
    complex(kind=8), dimension(1:m-2) :: DL, DU
    double precision :: H, ONE_OVER_H_SQUARED

    integer :: info

    H = (b-a)/dble(m)
    ONE_OVER_H_SQUARED = 1.d0/H**2.d0

    call generate_fdiff_matrix(m, z, DL, D, DU, ONE_OVER_H_SQUARED)
    w = generate_g(g, a, H, m)
    w(1) = w(1)+ONE_OVER_H_SQUARED*alpha
    w(m-1) = w(m-1)+ONE_OVER_H_SQUARED*beta

    call zgtsv(m-1, 1, DL, D, DU, w, m-1, info)
    write (*,*) 
    write (*,*) 'info:',info
  end subroutine fdiff_ode_solve

end module ode_solver_mod
