
    subroutine lsnPLS
!   ###################################################################
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,nnodes
    use edu2d_my_allocation , only : my_alloc_p2_ptr
    implicit none
    
    integer, parameter :: ndim=2 ! Dimension of the Space
    integer, parameter :: mtot=1 ! Number of Test points
    integer, parameter :: htot=mtot*ndim ! Number of derivatives predicted at test points
    integer, parameter :: stot=3 ! Number of Terms in the regression (stot=1 gives constnt mean function)
    
    integer :: i, j, inn,ierr,jj
    
!   not varying definitions
    integer :: ntot
    real(p2) :: dr,A(stot,stot)

!   local varying arrays
    real(p2),allocatable :: Kmat(:,:)    
    real(p2),allocatable :: LX(:,:)      ! ndim,ntot
    real(p2),allocatable :: HT(:,:)
    real(p2),allocatable :: H(:,:)       ! Polynomial dx,dy    
    real(p2),allocatable :: C(:,:)
    real(p2),allocatable :: N(:,:)
   
    do i=1,nnodes
       
     ntot = node(i)%nnghbrs            
     call my_alloc_p2_ptr(node(i)%aij,ntot) ! dynamically allocate array aij
     call my_alloc_p2_ptr(node(i)%bij,ntot) ! dynamically allocate array bij
     
!    collocation array
     ierr=0;allocate( lx(2,ntot+1),stat=ierr )
     if (ierr /= 0) pause 'lspls alloc error'
     lx=0.0d0
     
     lx(1,1) = node(i)%x
     lx(2,1) = node(i)%y
     
     jj=1 ! central node
     do j=1,ntot   
      jj=jj+1
      inn = node(i)%nghbr(j)   ! neighbour node id
      LX(1,jj) = node(inn)%x
      LX(2,jj) = node(inn)%y 
     enddo  
     
     ierr=0;allocate( C(stot,ntot+1),stat=ierr )
     if (ierr /= 0) pause 'lspls alloc error'
     ierr=0;allocate( N(stot,ntot+1),stat=ierr )
     if (ierr /= 0) pause 'lspls alloc error'
     ierr=0;allocate( H(stot,ntot+1),stat=ierr )
     if (ierr /= 0) pause 'lspls alloc error'
     ierr=0;allocate( HT(ntot+1,stot),stat=ierr )
     if (ierr /= 0) pause 'lspls alloc error'
     ierr=0;allocate( Kmat(ntot+1,ntot+1),stat=ierr ) ! must include centre
     if (ierr /= 0) pause 'lspls alloc error'
     C=0.0d0;N=0.0d0;H=0.0d0;HT=0.0d0;Kmat=0.0d0
     
!    Polynomial Matrix & its transpose matrix HT

!    Centre Node Information
     H(1,1) = 1.0d0
     H(2,1) = node(i)%x
     H(3,1) = node(i)%y
     
     jj=1
     do j=1,ntot 
      jj=jj+1
      H(1,jj) = 1.0d0
      H(2,jj) = lx(1,jj)  ! x components
      H(3,jj) = lx(2,jj)  ! y components
     enddo
     HT=transpose(H)  ! reset & define

!    WEIGHTING MATRIX ( Diagonal )

     Kmat(1,1) = 1.0d0

     jj =1
     do j=1,ntot
      jj = jj + 1
      Kmat(jj,jj)=1.0d0  ! weighted diagonal components of neighbours
     enddo
     
!    Matrix Multiplications to get LS coefficients
     call matrixmult(stot,ntot+1,H,ntot+1,Kmat,N)  ! Matrix N = H*W [stot,ntot]
     call matrixmult(stot,ntot+1,N,3,HT,A)                  ! N*H^t ( Least Square Matrix ) [stot,stot]
     call invertsymmetric(stot,A)                            ! Invert Least Square Matrix
     
!    Least Squares Coefficients
     N=0.0d0 ; call symmatrixmult(stot,A,ntot+1,H,N)  ! A^-1 * H           ! New N (stot,ntot)
     call matrixmult(stot,ntot+1,N,ntot+1,Kmat,C)       ! A^-1 * H^t * K^-1

!    store meshless coefficients
     node(i)%aij=0.0d0;node(i)%bij=0.0d0;jj=1
     do j=1,ntot
     jj=jj+1
     node(i)%aij(j) = C(2,jj)
     node(i)%bij(j) = C(3,jj)
     enddo
    
!    deallocate local arrays
     deallocate(C);deallocate(N);deallocate(H);deallocate(Lx)
     deallocate(HT);deallocate(Kmat)

    enddo
    
      
    return
!   ###################################################################
    end subroutine lsnPLS