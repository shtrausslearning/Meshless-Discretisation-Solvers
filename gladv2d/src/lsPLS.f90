
    subroutine lsPLS
!   ###################################################################
    use prm_flow
    implicit none

    integer :: i, j, inn,ierr,jj
    
!   not varying definitions
    integer :: ntot
    double precision :: x0,y0,dr
    double precision :: A(3,3)   ! Least Square Matrix

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
       
     ntot = nnbr(i)   ! number of neighbour
     
     ierr=0;allocate( LX2(2,ntot+1),stat=ierr );if(ierr /= 0) pause 'allocation error LX2'
     LX2=0.0d0
     
     LX2(1,1) = coord(1,i) ; LX2(2,1) = coord(2,i) ; jj=1 ! central node
     do j=1,ntot   
      jj=jj+1 ; inn = conn(j,i)
      LX2(1,jj) = coord(1,inn); LX2(2,jj) = coord(2,inn)  ! satellite nodes
     enddo  
     
     ierr=0 ; allocate( w(ntot+1),stat=ierr );if( ierr /= 0 ) pause 'allocation error LX2'
     w=0.0d0
     
!    allocate local arrays
     ierr=0;allocate( C(3,ntot+1),stat=ierr );if(ierr /= 0) pause 'allocation error LX2'
     ierr=0;allocate( N(3,ntot+1),stat=ierr );if(ierr /= 0) pause 'allocation error LX2'
     ierr=0;allocate( H(3,ntot+1),stat=ierr );if(ierr /= 0) pause 'allocation error LX2' 
     ierr=0;allocate( HT(ntot+1,3),stat=ierr );if(ierr /= 0) pause 'allocation error LX2'
     ierr=0;allocate( Kmat(ntot+1,ntot+1),stat=ierr );if( ierr /= 0) pause 'allocation error LX2'
     C=0.0d0;N=0.0d0;H=0.0d0;HT=0.0d0;Kmat=0.0d0
     
!    Polynomial Matrix & its transpose matrix HT
     H(1,1) = 1.0d0 ; H(2,1) = coord(1,i) ; H(3,1) = coord(2,i) ; jj=1 ; w(1) = 1.0d0
     do j=1,ntot 
      jj=jj+1
      H(1,jj) = 1.0d0 ; H(2,jj) = LX2(1,jj)  ; H(3,jj) = LX2(2,jj) 
      dr = dsqrt( H(1,jj)**2 + H(2,jj)**2 ) ; w(jj) = 1.0d0/(dr**3.0d0)
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
     
     ax(:,i) = 0.0d0 ; ay(:,i) = 0.0d0 ; jj=1
     do j=1,ntot
      jj = jj + 1 
      ax(j,i) = C(2,jj)
      ay(j,i) = C(3,jj)
     enddo
    
!    deallocate local arrays
     deallocate(C);deallocate(N);deallocate(H);deallocate(Lx2)
     deallocate(HT);deallocate(Kmat);deallocate(w)

    endif
    enddo
      
    return
!   ###################################################################
    end subroutine lsPLS