
!   subroutine to calculate lift and drag pressure forces

    subroutine forces
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    implicit none
    
    real(p2) :: sx,sy,pwall,nx,ny,ds,cx,cy,cp,dcy,dcx
    integer :: n1,n2,i,j
    
    cx=0.0d0
    cy=0.0d0
    
    do i=1,nbound
    do j=1,bound(i)%nbfaces
    if(trim(bound(i)%bc_type) == "slip_wall_weak") then

    n1 = bound(i)%bnode(j)
    n2 = bound(i)%bnode(j+1)
    ds = bound(i)%bfn(j)  
    sx = bound(i)%bfnx(j)*ds
    sy = bound(i)%bfny(j)*ds 
    pwall = 0.5d0*(node(n1)%dv(1)+node(n2)%dv(1))
    cp = 2.0d0*(pwall-pinf)/(rhoinf*qinf*qinf)
    dcy   = sy*cp
    dcx   = sx*cp
    cy    = cy + dcy
    cx    = cx + dcx
    
    endif
    enddo
    enddo
    
    clp = cy*cos(alpha) - cx*sin(alpha)
    cdp = cy*sin(alpha) + cx*cos(alpha)
    
    return
    end subroutine forces