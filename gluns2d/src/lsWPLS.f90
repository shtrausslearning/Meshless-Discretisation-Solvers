
    subroutine lsWPLS
!   ###################################################################
    use prmflow
    implicit none

    integer :: i, j, inn,ierr,jj
    
!   not varying definitions
    integer :: ntot
    real(rtype) :: x0,y0,dr
    real(rtype) :: A(3,3)   ! Least Square Matrix

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
       
     ntot = nbers(i)   ! number of neighbour
     
     ierr=0;allocate( LX(2,ntot+1),stat=ierr );if(ierr /= 0) call EM( "cannot allocate memory for LX()" )
     LX=0.0d0
     
     LX(1,1) = x(i) ; LX(2,1) = y(i) ; jj=1 ! central node
     do j=1,ntot   
      jj=jj+1 ; inn = conn(j,i)
      LX(1,jj) = x(inn); LX(2,jj) = y(inn)  ! satellite nodes
     enddo  
     
     ierr=0 ; allocate( w(ntot+1),stat=ierr );if( ierr /= 0 ) call EM( "cannot allocate memory for w()" )
     w=0.0d0
     
!    allocate local arrays
     ierr=0;allocate( C(3,ntot+1),stat=ierr );if(ierr /= 0) call EM( "cannot allocate memory for C()" )
     ierr=0;allocate( N(3,ntot+1),stat=ierr );if(ierr /= 0) call EM( "cannot allocate memory for N()" )
     ierr=0;allocate( H(3,ntot+1),stat=ierr );if(ierr /= 0) call EM( "cannot allocate memory for H()" ) 
     ierr=0;allocate( HT(ntot+1,3),stat=ierr );if(ierr /= 0) call EM( "cannot allocate memory for HT()" )
     ierr=0;allocate( Kmat(ntot+1,ntot+1),stat=ierr );if( ierr /= 0) call EM( "cannot allocate memory for Kmat()" )
     C=0.0d0;N=0.0d0;H=0.0d0;HT=0.0d0;Kmat=0.0d0
     
!    Polynomial Matrix & its transpose matrix HT
     H(1,1) = 1.0d0 ; H(2,1) = x(i) ; H(3,1) = y(i) ; jj=1 ; w(1) = 1.0d0
     do j=1,ntot 
      jj=jj+1
      H(1,jj) = 1.0d0 ; H(2,jj) = LX(1,jj)  ; H(3,jj) = LX(2,jj) 
      dr = dsqrt( H(1,jj)**2 + H(2,jj)**2 ) ; w(jj) = 1.0d0/(dr)
     enddo
     HT=transpose(H)  ! reset & define

     jj =1 ; Kmat(1,1) = w(1) ! first diagonal component
     do j=1,ntot
      jj = jj + 1
      Kmat(jj,jj)=w(jj)  ! weighted diagonal components
!      Kmat(jj,jj)=1.0d0
     enddo
     
!    Matrix Multiplications to get LS coefficients
     call matrixmult(3,ntot+1,H,ntot+1,Kmat,N)  ! Matrix N = H*W [stot,ntot]
     call matrixmult(3,ntot+1,N,3,HT,A)                  ! N*H^t ( Least Square Matrix ) [stot,stot]
     call invertsymmetric(3,A)                            ! Invert Least Square Matrix
     call symmatrixmult(3,A,ntot+1,N,C)                 ! A^-1 * (HW)
     
     aij(:,i) = 0.0d0 ; bij(:,i) = 0.0d0 ; jj=1
     do j=1,ntot
      jj = jj + 1 
      aij(j,i) = C(2,jj)
      bij(j,i) = C(3,jj)
     enddo
    
!    deallocate local arrays
     deallocate(C);deallocate(N);deallocate(H);deallocate(Lx)
     deallocate(HT);deallocate(Kmat);deallocate(w)

    enddo
      
    return
!   ###################################################################
    end subroutine lsWPLS