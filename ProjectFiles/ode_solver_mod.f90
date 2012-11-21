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

  function generate_g(g, a, H, m, t, t_flag)
    double precision, external :: g
    double precision, intent(in) :: a, H, t
    integer, intent(in) :: m
    logical, intent(in) :: t_flag

    double precision, dimension(1:m-1) :: generate_g
    integer :: i

    if (t_flag) then
       do i = 1, m-1
          generate_g(i) = g(a+h*i, t)
       end do
    else
       do i = 1, m-1
          generate_g(i) = g(a+h*i)
       end do
    end if
  end function generate_g

  !test function
  double precision function fun(x)
    double precision, intent(in) :: x
    fun = 0
  end function fun

  !example 2 function
  double precision function g_2(t)
    double precision, intent(in) :: t
    g_2 = (1.d0-t)**2.d0/(1.d0+t**2.d0)**2.d0
  end function g_2

  !example 3 function
  double precision function g_3(x, t)
    double precision, intent(in) :: x, t
    g_3 = exp(-t)*sin(x)+exp(-2.d0*t)*(2.d0*cos(t)-sin(t))*sin(2.d0*x)
  end function g_3

  subroutine fdiff_ode_solve(w, z, g, a, b, alpha, beta, m, t, t_flag, x_vect)
    integer, intent(in) :: m
    complex(kind=8), dimension(1:m-1), intent(inout) :: w
    double precision, dimension(1:m-1), intent(inout) :: x_vect
    complex(kind=8), intent(in) :: z 
    double precision, external :: g
    double precision, intent(in) :: a, b, alpha, beta, t
    logical, intent(in) :: t_flag
    
    complex(kind=8), dimension(1:m-1) :: D
    complex(kind=8), dimension(1:m-2) :: DL, DU
    double precision :: H, ONE_OVER_H_SQUARED
    
    integer :: info
    
    H = (b-a)/dble(m)
    ONE_OVER_H_SQUARED = 1.d0/H**2.d0
    
    do info = 1, m-1
       x_vect(info) = a+H*dble(info)
    end do

    call generate_fdiff_matrix(m, z, DL, D, DU, ONE_OVER_H_SQUARED)
    w = generate_g(g, a, H, m, t, t_flag)
    w(1) = w(1)+ONE_OVER_H_SQUARED*alpha
    w(m-1) = w(m-1)+ONE_OVER_H_SQUARED*beta
    
    call zgtsv(m-1, 1, DL, D, DU, w, m-1, info)
  end subroutine fdiff_ode_solve

end module ode_solver_mod
