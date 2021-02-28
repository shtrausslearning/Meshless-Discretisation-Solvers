
!   subroutine to determine nx_wc,ny_wc values

    subroutine nxnywall
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    implicit none
    
    integer :: ierr,jj,j,v1
    real(p2) :: nx,ny
    
    ierr=0;allocate( nx_wc(nnodes),stat=ierr )
    if( ierr /= 0 ) pause 'allocation error nx_wc'
    ierr=0;allocate( ny_wc(nnodes),stat=ierr )
    if( ierr /= 0 ) pause 'allocation error ny_wc'
    nx_wc=0.0d0;ny_wc=0.0d0
    
    nx=-777.0d0;ny=-777.0d0
    do j=1,nbound               
    do jj=1,bound(j)%nbnodes-1
    
     v1 = bound(j)%bnode(jj  ) ! main node
     nx = 0.5d0*(bound(j)%bnx(jj)+bound(j)%bnx(jj+1))
     ny = 0.5d0*(bound(j)%bny(jj)+bound(j)%bny(jj+1))
     
     nx_wc(v1) = nx
     ny_wc(v1) = ny
     
    enddo
    enddo
    
!   [check:read] [ DONT USE ]
    nx_wc=0.0d0;ny_wc=0.0d0
    open(44,file='outnxny.dat',form='formatted')
    do j=1,nnodes
    read(44,*) jj,nx_wc(j),ny_wc(j)
    enddo
    close(44)

!   check
    open(44,file='./info/outnxny.dat',form='formatted')
    do j=1,nnodes
    write(44,*) j,nx_wc(j),ny_wc(j)
    enddo
    close(44)
    
    return
    end subroutine nxnywall

!   subroutine needed to set nx,ny components during GCID bc cycle
!   outward unit vectors are defined in boundary data only

    subroutine gcnxy
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    implicit none
    
    integer :: ierr,i,nn,id1,id2,j,k,inn,ps,ps2
    real(p2) :: nx,ny
    integer, allocatable :: wnod(:)
    
!   get nx_wc,ny_wc
    call nxnywall
    
!   allocate for all gcid, ff not used
    ierr=0;allocate( gcidnxy(2,igcs),stat=ierr )
    if( ierr /= 0 ) pause 'allocation error gcidny'
    
    do i=1,igcs
    
    id1 = gcid(1,i)
    id2 = gcid(2,i)
    
!   count how many wall neighbours
    nn=0
    do j=1,nnodes
    if( node(j)%ptype .eq. 1 )then ! if physical node is wall node
     do k=1,node(j)%nnghbrs
     inn = node(j)%nghbr(k)
     if( inn .eq. id2 )then
     nn=nn+1
     endif
     enddo
    endif
    enddo
    
    allocate(wnod(nn));wnod=0
    
!   now store
    nn=0
    do j=1,nnodes
    if( node(j)%ptype .eq. 1 )then ! if physical node is wall node
     do k=1,node(j)%nnghbrs
     inn = node(j)%nghbr(k)
     if( inn .eq. id2 )then
     nn=nn+1
     wnod(nn) = j
     endif
     enddo
    endif
    enddo
    
    if( nn .eq. 1 )then
    ps = wnod(1)
    nx = nx_wc(ps)
    ny = ny_wc(ps)
    elseif( nn .eq. 2 )then
    ps = wnod(1)
    ps2 = wnod(2)
    nx = (nx_wc(ps) + nx_wc(ps2))/2.0d0
    ny = (ny_wc(ps) + ny_wc(ps2))/2.0d0
    elseif( nn .eq. 3 )then
    ps = wnod(1)
    ps2 = wnod(2)
    nx = (nx_wc(ps) + nx_wc(ps2))/2.0d0
    ny = (ny_wc(ps) + ny_wc(ps2))/2.0d0
    elseif( nn .ne. 0 )then
    pause 'nn > 3 '
    endif
    
    gcidnxy(1,i) = nx
    gcidnxy(2,i) = ny
    deallocate(wnod)
    
    enddo
    
    return
    end subroutine gcnxy
