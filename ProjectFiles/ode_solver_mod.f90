module ode_solver_mod
  implicit none

contains

  function generate_fdiff_matrix(n)
    integer, intent(in) :: n

    double precision, dimension(1:n-1, 1:n-1) :: generate_fdiff_matrix
    double precision, parameter :: ONE = 1.d0, TWO = 2.d0, ZERO = 0.d0
    integer :: i, j

    do i = 1, n-1
       do j = 1, n-1
          if (i == j) then
             generate_fdiff_matrix(i, j) = TWO
          else if (abs(i-j) == 1) then
             generate_fdiff_matrix(i, j) = -ONE
          else
             generate_fdiff_matrix(i, j) = ZERO
          end if
       end do
    end do

  end function generate_fdiff_matrix

  subroutine fdiff_ode_solve(w, g, alpha, beta, n)
    integer, intent(in) :: n
    double precision, dimension(1:n-1), intent(inout) :: w
    double precision, dimension(1:n-1), intent(in) :: g
    double precision, intent(in) :: alpha, beta

    double precision, dimension(1:n-1, 1:n-1) :: fdiff_mat
    double precision :: ONE_OVER_H_SQUARED

    ONE_OVER_H_SQUARED = 1.d0/(((beta-alpha)/dble(n))**2.d0)
    fdiff_mat = ONE_OVER_H_SQUARED*generate_fdiff_matrix(n)

  end subroutine fdiff_ode_solve

end module ode_solver_mod
