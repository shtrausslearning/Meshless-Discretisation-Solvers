
    Subroutine gcidarray
    use ModDataTypes
    use prmflow
    use ModInterfaces, only : EM
    implicit none
    
    integer :: ii,i,j,ic,iii,ii2,ierr

    allocate( fwnid(phynod),stat=ierr ) ! used locally only, but global variable
    if (ierr /= 0) call EM( "cannot allocate memory for fwnid()" )
  
    fwnid(:) = 0 ! reset
  
!   find fluid nodes that connect to wall nodes
    do i=1,phynod
    if( ntype(i) .eq. 3 )then  ! if node is fluid
    do j=1,nbers(i)            ! do for all neighbour nodes
    ic = conn(j,i)             ! neighbour node
    if( ntype(ic) .eq. 1 )then ! if neighbour is wall node
    fwnid(i) = 1               ! phynod marker of fluid node
    endif
    enddo
    endif
    enddo
  
!   just count how many wall dummy nodes are to be added for array allocation
!   number of ghost cells added is equivalent to the neighbour fluid nodes that 
!   are connected to the wall node
    ii=0
    do i=1,phynod
    if( fwnid(i) .eq. 1 )then
    ii=ii+1
    endif
    enddo
  
!   temporary x,y of wall ghost cell array allocation
    allocate( gcid(2,ii),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for fwnid()" )
    gcid(:,:) = 0
  
    igcs=ii     ! number of ghost cells
    iii=allnod  ! for gc numbering
    ii2 = 0
  
!   give global numbering to wall ghost cells

    do i=1,phynod
    if( fwnid(i) .eq. 1 )then ! if fluid node had wall neighbour

    iii=iii+1
    ii2=ii2+1
    gcid(1,ii2) = i   ! fluid node that has wall neighbour
    gcid(2,ii2) = iii ! wall dummy node id
  
    endif
    enddo
    
    ngcid = ii2 ! used in next subroutine
    
    return
    end subroutine gcidarray