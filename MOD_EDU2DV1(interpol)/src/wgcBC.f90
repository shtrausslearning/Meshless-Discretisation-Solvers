
  !   Subroutine used to generate boundary conditions for ghost cells

    subroutine wgcBC 
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    implicit none
    
    integer :: i,id1,id2,id3,ii,inn,j,iii
    real(p2) :: n_dx,n_dy,ug,vg,ui,vi,vnorm
    
    
    ii=0  ! need counter for storage
    do i=1,nbound                 ! cycle through all boundary segments
    do j=1,bound(i)%nbnodes-1     ! boundary node # exlude last as always
        
!   counter for GCID extraction ( just linear )
    ii=ii+1
    
!   (1) node relations
    id1 = gcid(1,ii)          ! 'fluid' node ID
    id2 = gcid(2,ii)          ! Ghost cell node ID    
    id3 = node(id2)%nghbr(1) ! wall node neighbour of gc
    
!   [check]
    if( id3 .ne. bound(i)%bnode(j) ) pause 'id3 /= bound(i)bnode(j)?'

!   Note: Both IP and GC were generated linearly from nbound loops, just use counter
    
!   (2) define local unit vectors : local boundary node unit vectors
    n_dx = bound(i)%bnx(j) ! boundary node nx of normal vector
    n_dy = bound(i)%bny(j) ! boundary node ny of normal vector
    
!   (3) generate id1(Image Point) CV data [averaged only here]

!   reset IP values
    node(id1)%cv = 0.0d0

    do iii=1,node(id1)%nnghbrs
    inn = node(id1)%nghbr(iii)
    node(id1)%cv(1) = node(id1)%cv(1) + node(inn)%cv(1) 
    node(id1)%cv(2) = node(id1)%cv(2) + node(inn)%cv(2) 
    node(id1)%cv(3) = node(id1)%cv(3) + node(inn)%cv(3) 
    node(id1)%cv(4) = node(id1)%cv(4) + node(inn)%cv(4) 
    enddo
    
!   set node(id1)'s / IPs VALUE
    node(id1)%cv(1) = node(id1)%cv(1)/dble(node(id1)%nnghbrs)
    node(id1)%cv(2) = node(id1)%cv(2)/dble(node(id1)%nnghbrs)
    node(id1)%cv(3) = node(id1)%cv(3)/dble(node(id1)%nnghbrs)
    node(id1)%cv(4) = node(id1)%cv(4)/dble(node(id1)%nnghbrs)
    
!   (4) generate id2(Ghost Node) CV data [ WALL ONLY]

    if( node(id2)%ptype .eq. 5 )then
!   update conservative variables of ghost cells
    ui = node(id1)%cv(2)/node(id1)%cv(1)  ! u of IP
    vi = node(id1)%cv(3)/node(id1)%cv(1)  ! v of IP
    vnorm = ui*n_dx + vi*n_dy
    
!   ghost cell values for velocity slip condition
    ug = ui-2.0d0*vnorm*n_dx
    vg = vi-2.0d0*vnorm*n_dy
    
!   ghost cell conservative variables
    node(id2)%cv(1) = 2.0d0*node(id1)%cv(1)-node(id3)%cv(1) 
    node(id2)%cv(2) = ug*node(id2)%cv(1)
    node(id2)%cv(3) = vg*node(id2)%cv(1)   
    node(id2)%cv(4) = 2.0d0*node(id1)%cv(4)-node(id3)%cv(4)
    
    endif
!   ##########################################################

    enddo
    enddo
    
!   (5) generate id2 farield boundary condition ; set farfield ghost nodes
    
    do i=1,nnodes+igcs
    if( node(i)%ptype .eq. 4 )then
    node(i)%cv(1) = rhoinf
    node(i)%cv(2) = rhoinf*uinf
    node(i)%cv(3) = rhoinf*vinf
    node(i)%cv(4) = pinf/(gamma-1.0d0) + 0.5d0*rhoinf*qinf*qinf
    endif
    enddo
    
    return
    end subroutine wgcBC