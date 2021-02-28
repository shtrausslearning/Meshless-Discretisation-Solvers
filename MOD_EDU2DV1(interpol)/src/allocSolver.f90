
    subroutine solveralloc
    use edu2d_constants, only : p2, zero
    use prmflow
    use edu2d_grid_data, only : construct_grid_data, check_grid_data
    use edu2d_my_main_data, only : nnodes,node

    integer :: ierr,i
    
    do i=1,nnodes+2*nIPGC  ! all pys + IP/GC (1 IP/GC per boundary node) 
    allocate( node(i)%cv(nconv) )   ! conservative variable
    allocate( node(i)%cvold(nconv) )   ! previous conservative variable
    allocate( node(i)%dv(ndepv) )   ! previous conservative variable
    allocate( node(i)%du(nconv) )   ! solution change du
    allocate( node(i)%w(nconv) )   ! primitive variable
    allocate( node(i)%gradw(nconv,2) )   ! xy primitive gradients
    allocate( node(i)%res(nconv) )   ! residual 
    allocate( node(i)%rhs(nconv) )     ! RHS term 
    end do
    
!   (0) allocate global arrays once
    if( epsirs .gt. 0.0d0 )then
    ierr=0;allocate( rhsold(4,nnodes),stat=ierr )
    if( ierr /= 0 ) pause 'allocation error rhsold'
    ierr=0;allocate( rhsit(4,nnodes),stat=ierr )
    if( ierr /= 0 ) pause 'allocation error rhsit'
    ierr=0;allocate( ncontr(nnodes),stat=ierr )
    if( ierr /= 0 )pause 'allocation error ncontr'
    endif
    
    return 
    end subroutine solveralloc