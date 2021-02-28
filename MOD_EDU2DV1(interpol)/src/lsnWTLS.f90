
!   [node based WTLS]
!   subroutine to generate WTLS meshless coefficients with weighting

    subroutine lsnWTLS
!   ###################################################################
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,nnodes
    use edu2d_my_allocation , only : my_alloc_p2_ptr
    implicit none

    integer, parameter :: ndim=2 ! Dimension of the Space
    integer, parameter :: mtot=1 ! Number of Test points
    integer, parameter :: htot=mtot*ndim ! Number of derivatives predicted at test points
    integer, parameter :: stot=2 ! Number of Terms in the regression (stot=1 gives constnt mean function)

    integer :: i, j, inn,ierr, ntot
    real(p2) :: x0,y0,dr
    real(p2) :: A(stot,stot)

!   local varying arrays
    real(p2),allocatable :: Kmat(:,:)    
    real(p2),allocatable :: LX(:,:)      ! ndim,ntot
    real(p2),allocatable :: HT(:,:)
    real(p2),allocatable :: H(:,:)       ! Polynomial dx,dy    
    real(p2),allocatable :: C(:,:)
    real(p2),allocatable :: N(:,:)
    real(p2),allocatable :: w(:)
    
!   calculate meshless coefficients for tls
    do i=1,nnodes
       
!    define centre node data
     x0 = node(i)%x
     y0 = node(i)%y      
     ntot = node(i)%nnghbrs            
     call my_alloc_p2_ptr(node(i)%aij,ntot) ! dynamically allocate array aij
     call my_alloc_p2_ptr(node(i)%bij,ntot) ! dynamically allocate array bij

!    store local neighbour node placements
     ierr=0;allocate( lx(2,ntot),stat=ierr )
     if (ierr /= 0) pause 'lstls alloc error'
     lx=0.0d0
     
     do j=1,ntot           
      inn = node(i)%nghbr(j)   ! neighbour node id
      lx(1,j) = node(inn)%x      ! x coordinate of neighbour
      lx(2,j) = node(inn)%y      ! y coordinate of neighbour
     enddo  
     
!    allocate local arrays
     ierr=0 ; allocate( w(ntot),stat=ierr )
     if (ierr /= 0) pause 'lstls alloc error'
     ierr=0 ; allocate( C(stot,ntot),stat=ierr )
     if (ierr /= 0) pause 'lstls alloc error'
     ierr=0 ; allocate( N(stot,ntot),stat=ierr )
     if (ierr /= 0) pause 'lstls alloc error'
     ierr=0 ; allocate( H(stot,ntot),stat=ierr )
     if (ierr /= 0) pause 'lstls alloc error'
     ierr=0 ; allocate( HT(ntot,stot),stat=ierr )
     if (ierr /= 0) pause 'lstls alloc error'
     ierr=0 ; allocate( Kmat(ntot,ntot),stat=ierr ) 
     if (ierr /= 0) pause 'lstls alloc error'
     C=0.0d0;N=0.0d0;H=0.0d0;H=0.0d0;HT=0.0d0;Kmat=0.0d0;w=0.0d0
     
!    Polynomial Matrix & its transpose matrix HT
     do j=1,ntot 
      H(1,j) = LX(1,j) - x0
      H(2,j) = LX(2,j) - y0
      dr = dsqrt( H(1,j)**2 + H(2,j)**2 )
      w(j) = 1.0d0/(dr**dble(gradient_weight_p))
     enddo
     HT=transpose(H)  ! reset & define
     
!    WEIGHTING MATRIX ( Diagonal )

     do j=1,ntot
      Kmat(j,j)=w(j)   ! Weighted Matrix
     enddo
     
     call matrixmult(stot,ntot,H,ntot,Kmat,N)      ! H*K^-1 [stot,ntot]
     call matrixmult(stot,ntot,N,stot,HT,A)        ! H*K^-1*H^t ( Least Square Matrix ) [stot,stot]
     call invertsymmetric(stot,A)                  ! Cholesky Inversion
     N=0.0d0 ; call symmatrixmult(stot,A,ntot,H,N)  ! A^-1 * H           ! New N (stot,ntot)
!     N=0.0d0 ; call matrixmult(stot,stot,A,ntot,H,N) ! A^-1 * H 
     call matrixmult(stot,ntot,N,ntot,Kmat,C)       ! A^-1 * H^t * K^-1
     
!    store meshless coefficients
     node(i)%aij=0.0d0;node(i)%bij=0.0d0
     do j=1,ntot
     node(i)%aij(j) = C(1,j)
     node(i)%bij(j) = C(2,j)
     enddo
    
!    deallocate local arrays
     deallocate(C);deallocate(N);deallocate(H);deallocate(lx)
     deallocate(HT);deallocate(Kmat);deallocate(w)

    enddo

    return
!   ###################################################################
    end subroutine lsnWTLS