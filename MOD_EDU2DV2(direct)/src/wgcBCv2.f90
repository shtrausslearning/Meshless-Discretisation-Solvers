
!   subroutine needed to set boundary values for CV of ghost nodes
!   calculate slip wall for all boundaries; farfield is overwritten later

    subroutine wgcBCv2
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    implicit none
    
    integer :: i,id1,id2
    real(p2) :: ui,vi,vnorm,ug,vg,nx,ny
    
    do i=1,igcs  ! cycle through all igcs nodes
        
    id1=gcid(1,i) ! IBW fluid node ID
    id2=gcid(2,i) ! corresp. ghost node
    
    nx = gcidnxy(1,i) ! retrieve nx
    ny = gcidnxy(2,i) ! retrieve ny

!    if( node(id2)%ptype .eq. 4 )then  ! should be wall's gc node
    
    ui = node(id1)%cv(2)/node(id1)%cv(1)  ! uip
    vi = node(id1)%cv(3)/node(id1)%cv(1)  ! vip
    vnorm = ui*nx+vi*ny
    ug = ui-2.0d0*vnorm*nx
    vg = vi-2.0d0*vnorm*ny
    
!   ghost cell node conservarive vars cv
    node(id2)%cv(1) = node(id1)%cv(1)
    node(id2)%cv(2) = ug*node(id2)%cv(1)
    node(id2)%cv(3) = vg*node(id2)%cv(1)
    node(id2)%cv(4) = node(id1)%cv(4)
    
    enddo ! i=1,igcs
    
!   overwrite farield boundary condition ; set farfield ghost nodes
    
    do i=1,nnodes+igcs
    if( node(i)%ptype .eq. 5 )then
    node(i)%cv(1) = rhoinf
    node(i)%cv(2) = rhoinf*uinf
    node(i)%cv(3) = rhoinf*vinf
    node(i)%cv(4) = pinf/(gamma-1.0d0) + 0.5d0*rhoinf*qinf*qinf
    endif
    enddo
    
    
    return
    end subroutine wgcBCv2