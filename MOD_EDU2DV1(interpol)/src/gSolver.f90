
    subroutine gsolver
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,nnodes
    use edu2d_my_allocation
    implicit none
    
    real(p2) :: ark(5)
    integer :: nrk,irk,i, j, p
    real(p2) :: wL(4),wR(4) ! left and right states
    real(p2) :: n12(2)      ! meshless coefficient vectors
    real(p2) :: ds
    real(p2) :: fi(4),fij(4)  ! fluxes from RHS 
    real(p2) :: dwL(4),dwR(4) ! solution difference
    
!   RK Coefficients for 5 stage 
    ark(1) = 0.0695d0 ; ark(2) = 0.1602d0
    ark(3) = 0.2898d0 ; ark(4) = 0.5060d0 ; ark(5) = 1.000d0
    nrk=5
    
    do i=1,nnodes
    node(i)%cvold = node(i)%cv
    enddo
    
    call gtimestep
    
    do irk=1,nrk
    
 !  evaluate RHS ------------------------------
 
    do i=1,nnodes
    node(i)%rhs = 0.0d0     ! reset all
    do j=1,node(i)%nnghbrs  ! for all physical neighbours nnghbr -> physical nodes only

    p = node(i)%nghbr(j)         ! corresp. neighbour node ID
    ds = sqrt(node(i)%aij(j)**2+node(i)%bij(j)**2)
    n12(1)=node(i)%aij(j)/ds ; n12(2)=node(i)%bij(j)/ds
    
    dwL=0.0d0 ; dwR = 0.0d0  ! no 2nd order    
!   define states (left)
    wL(1)=node(i)%cv(1)               + dwL(1)
    wL(2)=node(i)%cv(2)/node(i)%cv(1) + dwL(2)
    wL(3)=node(i)%cv(3)/node(i)%cv(1) + dwL(3)
    wL(4)=node(i)%dv(1)               + dwL(4)
!   define states (right)
    wR(1)=node(p)%cv(1)               + dwR(1)
    wR(2)=node(p)%cv(2)/node(p)%cv(1) + dwR(2)
    wR(3)=node(p)%cv(3)/node(p)%cv(1) + dwR(3)
    wR(4)=node(p)%dv(1)               + dwR(4)
    
!   flux evaluation (fij/fL)
    if(trim(inviscid_flux)=='roe')then
    call groe(wL,wR,n12,fij,fi)   ! updated
    elseif(trim(inviscid_flux)=='rhll')then
    call grhll(wL,wR,n12,fij,fi)
    endif
    
   !residual
    node(i)%rhs(1) = node(i)%rhs(1) + 2.0d0*( fij(1) - fi(1) )*ds
    node(i)%rhs(2) = node(i)%rhs(2) + 2.0d0*( fij(2) - fi(2) )*ds
    node(i)%rhs(3) = node(i)%rhs(3) + 2.0d0*( fij(3) - fi(3) )*ds
    node(i)%rhs(4) = node(i)%rhs(4) + 2.0d0*( fij(4) - fi(4) )*ds
    
    enddo
    enddo
 
 !  --------------------------------------------
 
 !  calculate timestep
!    call gtimestep ! makes abit slower
!    
    do i=1,nnodes
    node(i)%rhs = ark(irk)*node(i)%dt*node(i)%rhs
    enddo
    
    if( epsirs .gt. 0.0d0 ) call cirs
    
!   update solution
    do i=1,nnodes
    node(i)%cv = node(i)%cvold - node(i)%rhs
    enddo
    
!   update solution of BC (conservative variables)
    call wgcBC    
    call dependentvarsall

!   [check]
!    open(44,file='./info/checkME.dat',form='formatted')    
!    do i=1,nnodes+igcs
!    write(44,*) ( node(i)%cv(j), j=1,4 )
!    enddo
!    pause

    
    enddo
    
    return
    end subroutine gsolver