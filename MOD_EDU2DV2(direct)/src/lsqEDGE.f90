
!   Subroutine to arrange edge lsq coefficients
    subroutine lsqedge
    use edu2d_constants, only : p2
    use edu2d_my_main_data , only : node,edge,nedges
    use edu2d_my_allocation
    implicit none
    
    integer :: cnid(2) ! temporary storage  
    integer :: n1,n2,i,j, ii, inn
    real(p2) :: leaij1,leaij2,lebij1,lebij2
    
    integer :: tc1,tc2
    
!   Loop over edges and distribute the node numbers:
 
    EDGES1: do i=1,nedges

    ii=0; cnid(:) = 0
    n1 = edge(i)%n1
    n2 = edge(i)%n2
   
!   find n1
    do j=1,node(n1)%nnghbrs  ! look though all of point n1's neighbours
    inn = node(n1)%nghbr(j)
    if( inn .eq. n2 )then
    cnid(1) = j
    endif
    enddo
    
    do j=1,node(n2)%nnghbrs
    inn = node(n2)%nghbr(j)
    if( inn .eq. n1 )then
    cnid(2) = j
    endif
    enddo
    
    if( cnid(1) .eq. 0 .or. cnid(2) .eq. 0 ) pause 'cnid incorrect'
    
!   store coefficients locally
    leaij1 = node(n1)%aij(cnid(1)) ! aij edge term
    leaij2 = node(n2)%aij(cnid(2)) ! aji edge term
    lebij1 = node(n1)%bij(cnid(1)) ! bij edge term
    lebij2 = node(n2)%bij(cnid(2)) ! bji edge term
    
    call my_alloc_p2_ptr(edge(i)%eaij,2) ! local structure alloc
    call my_alloc_p2_ptr(edge(i)%ebij,2) ! local structure alloc
    edge(i)%eaij(1:2) = -777.0d0 ; edge(i)%ebij(1:2) = -777.0d0
    
!   AIJ Components -------------------------------
    edge(i)%eaij(1) = leaij1 ! add aij component 
    edge(i)%eaij(2) = leaij2 ! add aji component

!   BIJ Components ------------------------------
    edge(i)%ebij(1) = lebij1 ! add aij component 
    edge(i)%ebij(2) = lebij2 ! add aji component
    
!   check how many edges have same directivity orientation
    if( edge(i)%eaij(1) .gt. 0.0d0 .and. edge(i)%eaij(2) .gt. 0.0d0 ) tc1 = tc1 + 1
    if( edge(i)%eaij(1) .lt. 0.0d0 .and. edge(i)%eaij(2) .lt. 0.0d0 ) tc1 = tc1 + 1
    if( edge(i)%eaij(1) .gt. 0.0d0 .and. edge(i)%eaij(2) .gt. 0.0d0 ) tc2 = tc2 + 1
    if( edge(i)%eaij(1) .lt. 0.0d0 .and. edge(i)%eaij(2) .lt. 0.0d0 ) tc2 = tc2 + 1
    
    enddo EDGES1

!   if need to check same sign cases
!    write(*,*) tc1,tc2
!    pause 'tc1,tc2'
    
!   check array completeness
    do i=1,nedges
    do j=1,2
    if( edge(i)%eaij(j) .eq. -777.0d0 ) pause 'allocation error edge(i)%eaij'
    if( edge(i)%ebij(j) .eq. -777.0d0 ) pause 'allocation error edge(i)%ebij'
    enddo
    enddo
    
    return
    end subroutine lsqedge