MODULE random
! A module for random number generation from the following distributions:
!
!     Distribution                    Function/subroutine name
!
!     Normal (Gaussian)               random_normal
  ! https://jblevins.org/mirror/amiller/random.f90

use real_precision
IMPLICIT NONE
REAL, PRIVATE      ::  half = 0.5_dp


CONTAINS


FUNCTION random_normal() RESULT(fn_val)
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL(kind=dp) :: fn_val

!     Local variables
REAL(kind=dp) :: s = 0.449871_dp, t = -0.386595_dp, a = 0.19600_dp, &
     b = 0.25472_dp, r1 = 0.27597_dp, r2 = 0.27846_dp, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156_dp * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0_dp*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal


subroutine init_ranseed
   implicit none
   integer :: values(1:8), k
   integer, dimension(:), allocatable :: seed

   
   call date_and_time(values=values)

   call random_seed(size=k)
   allocate(seed(1:k))
   seed(:) = values(8)
   call random_seed(put=seed)
   return
end subroutine init_ranseed



end MODULE random
