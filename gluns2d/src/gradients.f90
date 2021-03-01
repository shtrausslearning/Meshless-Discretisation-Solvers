
    subroutine Gradients

    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none

    ! local variables
    integer     :: i, j, p, inn, m,id1,id2,jj
    real(rtype) :: dQ(4),ds,nx,ny
!   ###################################################################

    
!   Physical Nodes
    do i=1,phynod
    if( ptype(i) .eq. 1 )then
    gradx(:,i) = 0.0d0; grady(:,i) = 0.0d0
    
    do j=1,nbers(i)
    p=conn(j,i)
    
    ds = sqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
    nx = aij(j,i)
    ny = bij(j,i)
    
!   dQ
    dQ(1) = cv(1,p)         - cv(1,i)
    dQ(2) = cv(2,p)/cv(1,p) - cv(2,i)/cv(1,i)
    dQ(3) = cv(3,p)/cv(1,p) - cv(3,i)/cv(1,i)
    dQ(4) = dv(1,p)         - dv(1,i)
    
!    do jj=1,4
!    dQ(jj) = abs(dQ(jj))
!    enddo
    
!   X gradients
    do jj=1,4
    gradx(jj,i) = gradx(jj,i) + nx*dQ(jj) ! x gradients
    grady(jj,i) = grady(jj,i) + ny*dQ(jj)
    enddo

    enddo
    endif
    enddo
    
!!   Values Ghost Cells (ff,or wall ghost cells)
!!   ##########################################
!    do i=1,allnod
!    if( ntype(i) .eq. 1 )then ! if node is wall
!    
!    do j=1,nbers(i)
!    inn = conn(j,i)
!    if( ptype(inn) .eq. 2 )then ! if ghost cell
!    do m=1,4
!!    gradx(m,inn) = -gradx(m,i)
!!    grady(m,inn) = -grady(m,i)
!    gradx(m,inn) = 0.0d0
!    grady(m,inn) = 0.0d0
!    enddo
!    endif
!    enddo
!    
!    endif
!    enddo
    
!   Alternative Adoptation for Ghost Cells
    do i=1,igcs
    
    id1 = gcid(1,i) ! fluid global number
    id2 = gcid(2,i) ! ghost cell global number
    
    do m=1,4
    gradx(m,id2) = gradx(m,id1)
    grady(m,id2) = grady(m,id1)
    enddo
    
    enddo
!    !
!    open(44,file='qxqy.dat',form='formatted')
!    do i=1,allnod+igcs
!    write(44,*) gradx(1,i),grady(1,i)
!    enddo
!    close(44)
    
    
    do i=1,allnod+igcs
    if( gradx(1,i) .eq. -777.0d0 )then
    write(*,*) 'Undefined qx(1,i)'
    pause
    endif
    enddo
    
    return
    end subroutine Gradients
