module w_mod
  implicit none
  
contains
  complex(kind=8) function w_1(z)
    complex(kind=8), intent(in) :: z
    w_1 = 1.d0/((z+1)*(z+2))
  end function w_1
  
  complex(kind=8) function w_2(z)
    complex(kind=8), intent(in) :: z
    w_2 = 0
  end function w_2
  
  double precision function u_1(t)
    double precision, intent(in) :: t
    u_1 = exp(-1.d0*t)-exp(-2.d0*t)
  end function u_1
  
  function u_3(x_vect, t, m)
    double precision, intent(in) :: t
    integer :: m
    double precision, dimension(1:m-1), intent(in) :: x_vect
    double precision, dimension(1:m-1) :: u_3
    
    u_3 = (1.d0+t)*exp(-1.d0*t)*sin(x_vect)+cos(t)*exp(-2.d0*t)*sin(2.d0*x_vect)
  end function u_3
end module w_mod
