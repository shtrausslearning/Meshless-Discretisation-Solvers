
!> Computes reference values of limited variables (density, u, v, pressure)
!! and of the control volumes. The reference values are used to normalize
!! variables within the limiter functions (Roe's upwind scheme).
!!
    subroutine LimiterRefvals

    use ModDataTypes
    use prmflow
    implicit none

    integer :: ierr
    real(rtype) :: gam1, rgas, temp, rho, cs, mach
    real(rtype) :: rvolref,limfac3

    volref=maxval(vol)
    gam1 = gamma - 1.D0
    rgas = gam1*cpgas/gamma
    limref(1) = rhoinf
    limref(2) = Sqrt(uinf*uinf+vinf*vinf)
    limref(3) = limref(2)
    limref(4) = pinf
    
    ierr=0;allocate( umin(4,allnod+igcs),stat=ierr )
    if( ierr /= 0 ) call EM("Allocation Error Umin")
    ierr=0;allocate( umax(4,allnod+igcs),stat=ierr )
    if( ierr /= 0 ) call EM("Allocation Error Umax") 
    
    umin = 0.0d0
    umax = 0.0d0
    
    limfac3 = limfac*limfac*limfac
    rvolref = 1.D0/volref**1.5D0
    eps2(1) = limfac3*limref(1)*limref(1)*rvolref
    eps2(2) = limfac3*limref(2)*limref(2)*rvolref
    eps2(3) = limfac3*limref(3)*limref(3)*rvolref
    eps2(4) = limfac3*limref(4)*limref(4)*rvolref

    return
    end subroutine LimiterRefvals
