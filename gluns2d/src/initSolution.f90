
    subroutine initialsolution
    use ModDataTypes
    use prmflow
    implicit none

    integer     :: i, j, ib, ibn, ibegn, iendn
    real(rtype) :: gam1, rgas, xmin, xmax, dx, pinl, temp, rho
    real(rtype) :: cs, mach, q, dp, dbeta, beta, p, u, v, g1cp, rhoq
! ###################################################################

    gam1 = gamma - 1.D0
    rgas = gam1*cpgas/gamma
    g1cp = gam1*cpgas

!   conservative variabls
    do i=1,allnod
    cv(1,i) = rhoinf
    cv(2,i) = rhoinf*uinf
    cv(3,i) = rhoinf*vinf
    cv(4,i) = pinf/gam1 + 0.5D0*rhoinf*qinf*qinf
    enddo
    
!   dependable variables
    do i=1,allnod
    rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
    dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
    dv(2,i) = dv(1,i)/(rgas*cv(1,i))
    dv(3,i) = Sqrt(g1cp*dv(2,i))
    dv(4,i) = gamma
    dv(5,i) = cpgas
    enddo
    
!   calculate cv,dv for wall ghost cells
    call bc_wgc
    
    end subroutine initialsolution
