! http://stackoverflow.com/questions/12523524/f2py-specifying-real-precision-in-fortran-when-interfacing-with-python

module types

implicit none
integer, parameter :: lp = selected_real_kind(3,4) ! low precision, for imputed genotypes
integer, parameter :: sp = selected_real_kind(6,37) ! single precision, eq. 32 bit precision.
integer, parameter :: dp = selected_real_kind(15,307) ! double precision, eq. 64 bit precision

real(sp) :: r_sp = 1.0
real(dp) :: r_dp = 1.0_dp

end module

!
!module input
!
!contains 
!
!subroutine input_sp(val)
!  use types
!  real(sp), intent(in) :: val
!  real(sp) :: x
!  x = val
!  write(*,*) x
!end subroutine
!
!subroutine input_dp(val)
!  use types
!  real(dp), intent(in) :: val
!  real(dp) :: x
!  x = val
!  write(*,*) dp, val, x
!end subroutine
!
!end module