
module convertUW
    
!********************************************************************************
!* Compute U from W
!*
!* ------------------------------------------------------------------------------
!*  Input:  w =    primitive variables (rho,     u,     v,     p)
!* Output:  u = conservative variables (rho, rho*u, rho*v, rho*E)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
contains

 function w2u(w) result(u)

 use edu2d_constants   , only : p2, one, half
 use prmflow

 implicit none

 real(p2), dimension(4), intent(in) :: w

!Local variables
 real(p2), dimension(4)             :: u !output

  u(1) = w(1)
  u(2) = w(1)*w(2)
  u(3) = w(1)*w(3)
  u(4) = w(4)/(gamma-one)+half*w(1)*(w(2)*w(2)+w(3)*w(3))

 end function w2u
!--------------------------------------------------------------------------------

!********************************************************************************
!* Compute U from W
!*
!* ------------------------------------------------------------------------------
!*  Input:  u = conservative variables (rho, rho*u, rho*v, rho*E)
!* Output:  w =    primitive variables (rho,     u,     v,     p)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
 function u2w(u) result(w)

 use edu2d_constants   , only : p2, one, half
 use prmflow

 implicit none

 real(p2), dimension(4), intent(in) :: u

!Local variables
 real(p2), dimension(4)             :: w !output

  w(1) = u(1)
  w(2) = u(2)/u(1)
  w(3) = u(3)/u(1)
  w(4) = (gamma-one)*( u(4) - half*w(1)*(w(2)*w(2)+w(3)*w(3)) )

 end function u2w
!--------------------------------------------------------------------------------

end module convertUW