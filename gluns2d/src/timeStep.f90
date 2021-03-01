

  subroutine TimeStep
! ###################################################################
  use ModDataTypes
  use prmflow
  implicit none

! local variables
    integer     :: i, j,inn
    real(rtype) :: sx, sy, ds, rrho, u, v, vc, cs, rhop, rhoT, hT, q2, &
                 theta, ra1g, a1, a4, a5, fmue, f1, f2, fac, dtv, cfac, &
                 lambdac, lambdav, tsmin, ter1, ter2, term, nx,ny
! ###################################################################
    
!   spectral radius of the inviscid flux jacobian for
!   the meshless volume scheme at node i
    do i = 1,phynod
    
     term = 0.0d0                           ! Summation of bottom term
     do j = 1,nbers(i)
     inn = conn(j,i)
     sx = aij(j,i)                     ! aij contains neighbour values for all nodes, even ghost cells
     sy = bij(j,i)                      ! bij contains neighbour values for all nodes, even ghost cells
     ds = sqrt(sx*sx+sy*sy)
     sx = sx/ds ; sy = sy/ds
     u = cv(2,i)/cv(1,i)                 ! u of MAIN NODE ONLY
     v = cv(3,i)/cv(1,i)                 ! v of MAIN NODE ONLY
     ter1 = abs(sx*u + sx*v)                      ! spectral radius part i
     ter2 = sqrt(dv(4,i)*dv(1,i)/cv(1,i))*sqrt(sx**2+sy**2) ! spectral radius part ii
     !ter2 = dv(3,i)*(sx**2+sy**2)**0.5d0  ! spectral radius part ii 
 !    term = term + (ter1+ter2)             ! semilocal, summation of all components, reset every cube     
     term = term + (ter1 + ter2)
     enddo

!    local timestep
     tstep(i) = cfl/term

    enddo 

! in case of global time-stepping - find min. time step in domain -------------
    if (ktimst == 0) then
    tsmin = 1.0d32
    do i = 1,phynod
     tsmin = min(tsmin,tstep(i))
    enddo
    do i = 1,phynod
     tstep(i) = tsmin
    enddo
    
    endif

    return
    end subroutine TimeStep
