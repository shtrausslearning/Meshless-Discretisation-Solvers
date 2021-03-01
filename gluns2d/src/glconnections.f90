
    subroutine createConnections
    use prmflow
    implicit none
    
    integer :: ierr,j,ii,ibn,i

!   main sort array allocation
    ierr=0;allocate( conn(12,phynod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for conn()" )
    conn(:,:) = 0
  
!   1. add all physical neighbours

    do j=1,phynod
    ii=0
  
    do i=1,nedint
    if( edge(1,i) .eq. j  )then
    ii=ii+1
    conn(ii,j) = edge(2,i)
    elseif( edge(2,i) .eq. j )then
    ii= ii+1
    conn(ii,j) = edge(1,i)
    endif
    enddo
    enddo
  
    ierr=0;allocate( nbers(phynod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for nbers()" )
  
    do j=1,phynod
    ii=0
    do i=1,12
    if( conn(i,j) .ne. 0 )then
    ii=ii+1
    endif
    enddo
    nbers(j) = ii  ! number of neighbours ( physical domain only )
    enddo

    !open(88,file='conn.dat',form='formatted')
    !do i=1,phynod
    !write(88,*) ( conn(j,i),j=1,12 )
    !enddo
    !close(88)
    !open(89,file='conn2.dat',form='formatted')
    !do i=1,phynod
    !write(89,*) ( conn(j,i),j=1,nbers(i) )
    !enddo
    !close(89)
    
!   Add Farfield Ghost Nodes to CONN()
  
    do j=1,phynod
    ii=nbers(j) ! current neighbour
    if( ntype(j) .eq. 2 )then ! if node is ff boundary
   
!   find dummy node of boundary node
    do ibn=1,nbnodes
    if( bnode(1,ibn) == j )then 
    conn(ii+1,j) = bnode(2,ibn)  ! bnode(2,ibn) is the farfield Ghost Node
    endif
    enddo

    endif
    enddo
  
    do j=1,phynod
    ii=nbers(j)                 ! current neighbour
    if( ntype(j) .eq. 2  )then  ! if node is ff boundary  
    nbers(j) = ii+1             ! just add extra 
    endif
    enddo
    
    return
    end subroutine createConnections