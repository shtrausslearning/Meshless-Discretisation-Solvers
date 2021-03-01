    subroutine lsTLS
!   ###################################################################
    use prm_flow
    implicit none

    integer :: i, j, inn,ierr,jj
    
!   not varying definitions
    integer :: ntot
    double precision :: x0,y0,dr
    double precision :: A(2,2)   ! Least Square Matrix

!   local varying arrays
    double precision,allocatable :: Kmat(:,:)    
    double precision,allocatable :: LX2(:,:)      ! ndim,ntot
    double precision,allocatable :: HT(:,:)
    double precision,allocatable :: H(:,:)       ! Polynomial dx,dy    
    double precision,allocatable :: C(:,:)
    double precision,allocatable :: N(:,:)
    double precision,allocatable :: w(:)
   
    ax(:,:) = 0.0d0 ; ay(:,:) = 0.0d0 
    
    do i=1,npts
    if( ptype(i) .eq. 1 )then
       
     x0 = coord(1,i)
     y0 = coord(2,i)              ! centre point
     ntot = nnbr(i)    ! number of neighbours of current node

     ierr=0;allocate( LX2(2,ntot),stat=ierr );if(ierr /= 0) pause 'alloc error'
     LX2=0.0d0
     
     do j=1,ntot           
      inn = conn(j,i)
      LX2(1,j) = coord(1,inn)
      LX2(2,j) = coord(2,inn)   
     enddo  
     
     ierr=0 ; allocate( w(ntot),stat=ierr )
     w=0.0d0
     
!    allocate local arrays
     ierr=0;allocate( C(2,ntot),stat=ierr )
     ierr=0;allocate( N(2,ntot),stat=ierr )
     ierr=0;allocate( H(2,ntot),stat=ierr )
     ierr=0;allocate( HT(ntot,2),stat=ierr )
     ierr=0;allocate( Kmat(ntot,ntot),stat=ierr )
     C=0.0d0;N=0.0d0;H=0.0d0;HT=0.0d0;Kmat=0.0d0
     
!    Polynomial Matrix & its transpose matrix HT
     do j=1,ntot 
      H(1,j) = LX2(1,j) - x0
      H(2,j) = LX2(2,j) - y0
      dr = dsqrt( H(1,j)**2 + H(2,j)**2 )
!      dr = H(1,j)**2 + H(2,j)**2
      w(j) = 1.0d0/(dr**1.0d0)
     enddo
     HT=transpose(H)  ! reset & define

     do j=1,ntot
!      Kmat(j,j)=w(j)   ! Weighted Matrix
      Kmat(j,j) = 1.0d0
     enddo
     
!    Matrix Multiplications to get LS coefficients
     call matrixmult(2,ntot,H,ntot,Kmat,N)  ! Matrix N = H*W [stot,ntot]
     call matrixmult(2,ntot,N,2,HT,A)          ! N*H^t ( Least Square Matrix ) [stot,stot]
     
!     if( i .eq. 101 )then
!     do j=1,2
!     write(*,*) ( A(j,jj),jj=1,2 )
!     enddo
!     pause
!     endif
     
     call invertsymmetric(2,A)                     ! Invert Least Square Matrix
     call symmatrixmult(2,A,ntot,N,C)         ! A^-1 * (HW)
     
     ax(:,i) = 0.0d0 ; ay(:,i) = 0.0d0
     do j=1,ntot
      ax(j,i) = C(1,j)
      ay(j,i) = C(2,j)
     enddo
     
!     if( i .eq. 101 )then
!     do j=1,ntot
!     write(*,*) ax(j,i)
!     write(*,*) ''
!     write(*,*) ay(j,i)
!     enddo
!     pause
!     endif
    
!    deallocate local arrays
     deallocate(C);deallocate(N);deallocate(H);deallocate(Lx2)
     deallocate(HT);deallocate(Kmat);deallocate(w)

    endif
    enddo
      
    return
!   ###################################################################
    end subroutine lsTLS