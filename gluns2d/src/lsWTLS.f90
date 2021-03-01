
    subroutine lsWTLS
!   ###################################################################
    use prmflow
    implicit none

    integer :: i, j, inn,ierr,jj
    
!   not varying definitions
    integer :: ntot
    real(rtype) :: x0,y0,dr
    real(rtype) :: A(2,2)   ! Least Square Matrix

!   local varying arrays
    real(rtype),allocatable :: Kmat(:,:)    
    real(rtype),allocatable :: LX(:,:)      ! ndim,ntot
    real(rtype),allocatable :: HT(:,:)
    real(rtype),allocatable :: H(:,:)       ! Polynomial dx,dy    
    real(rtype),allocatable :: C(:,:)
    real(rtype),allocatable :: N(:,:)
    real(rtype),allocatable :: w(:)
   
    ierr=0;allocate( aij(12,phynod),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for bij()" )
    ierr=0;allocate( bij(12,phynod),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for bij()" )
    aij(:,:) = 0.0d0 ; bij(:,:) = 0.0d0 
    
    do i=1,phynod
       
     x0 = x(i)
     y0 = y(i)              ! centre point
     ntot = nbers(i)    ! number of neighbours of current node

     ierr=0;allocate( LX(2,ntot),stat=ierr );if(ierr /= 0) call EM( "cannot allocate memory for LX()" )
     LX=0.0d0
     
     do j=1,ntot           
      inn = conn(j,i)
      LX(1,j) = x(inn)
      LX(2,j) = y(inn)   
     enddo  
     
     ierr=0 ; allocate( w(ntot),stat=ierr );if( ierr /= 0 ) call EM( "cannot allocate memory for w()" )
     w=0.0d0
     
!    allocate local arrays
     ierr=0;allocate( C(2,ntot),stat=ierr );if(ierr /= 0) call EM( "cannot allocate memory for C()" )
     ierr=0;allocate( N(2,ntot),stat=ierr );if(ierr /= 0) call EM( "cannot allocate memory for N()" )
     ierr=0;allocate( H(2,ntot),stat=ierr );if(ierr /= 0) call EM( "cannot allocate memory for H()" ) 
     ierr=0;allocate( HT(ntot,2),stat=ierr );if(ierr /= 0) call EM( "cannot allocate memory for HT()" )
     ierr=0;allocate( Kmat(ntot,ntot),stat=ierr );if( ierr /= 0) call EM( "cannot allocate memory for Kmat()" )
     C=0.0d0;N=0.0d0;H=0.0d0;HT=0.0d0;Kmat=0.0d0
     
!    Polynomial Matrix & its transpose matrix HT
     do j=1,ntot 
      H(1,j) = LX(1,j) - x0
      H(2,j) = LX(2,j) - y0
      dr = dsqrt( H(1,j)**2 + H(2,j)**2 )
   !   dr = H(1,j)**2 + H(2,j)**2
      w(j) = 1.0d0/(dr**3.0d0)
     enddo
     HT=transpose(H)  ! reset & define

     do j=1,ntot
!      Kmat(j,j)=w(j)   ! Weighted Matrix
      Kmat(j,j)=1.0d0
     enddo
     
!    Matrix Multiplications to get LS coefficients
     call matrixmult(2,ntot,H,ntot,Kmat,N)  ! Matrix N = H*W [stot,ntot]
     call matrixmult(2,ntot,N,2,HT,A)          ! N*H^t ( Least Square Matrix ) [stot,stot]
     call invertsymmetric(2,A)                     ! Invert Least Square Matrix
     call symmatrixmult(2,A,ntot,N,C)         ! A^-1 * (HW)
     
     aij(:,i) = 0.0d0 ; bij(:,i) = 0.0d0
     do j=1,ntot
      aij(j,i) = C(1,j)
      bij(j,i) = C(2,j)
     enddo
    
!    deallocate local arrays
     deallocate(C);deallocate(N);deallocate(H);deallocate(Lx)
     deallocate(HT);deallocate(Kmat);deallocate(w)

    enddo
      
    return
!   ###################################################################
    end subroutine lsWTLS