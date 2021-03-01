    subroutine lsCWTLS
!   ###################################################################
    use prm_flow
    use covars
    implicit none

    integer :: i, j, nn,ierr,ii,jj
    integer, parameter :: ndim=2 ! Dimension of the Space
    integer, parameter :: mtot=1 ! Number of Test points
    integer, parameter :: htot=mtot*ndim ! Number of derivatives predicted at test points
    integer, parameter :: stot=2 ! Number of Terms in the regression (stot=1 gives constnt mean function)
    
!   not varying definitions
    integer :: ntot
    double precision :: x0,y0
    double precision :: theta(ndim),sigma,sigmaN
    double precision :: A(stot,stot)
    
!   local varying arrays
    double precision,allocatable :: KijT(:,:)    ! training covariance
    double precision,allocatable :: LX2(:,:)      ! ndim,ntot
    double precision,allocatable :: HT(:,:)      ! Polynomial 
    double precision,allocatable :: H(:,:)       ! Polynomial dx,dy    
    double precision,allocatable :: C(:,:)
    double precision,allocatable :: N(:,:)
  
    ax(:,:) = 0.0d0 ; ay(:,:) = 0.0d0 
    
    open(5,file='kriging.dat',form='formatted')
    read(5,*) ( theta(i),i=1,ndim )
    read(5,*) sigma
    read(5,*) sigmaN
    read(5,*) covarflag
    close(5)

    do i=1,npts
    if( ptype(i) .eq. 1 )then
        
     x0 = coord(1,i) ; y0 = coord(2,i)              ! centre point
     ntot = nnbr(i)                    ! number of neighbours of current node

     ierr=0;allocate( LX2(ndim,ntot),stat=ierr )
     if( ierr /= 0) pause 'allocation error'
     
     do j=1,ntot           
      nn    = conn(j,i) ! neighbour id   
      LX2(1,j) = coord(1,nn)   ! neighbour global X
      LX2(2,j) = coord(2,nn)   ! neighbour global Y
     enddo  
        
!    allocate local arrays
     ierr=0 ; allocate( C(stot,ntot),stat=ierr )
     if( ierr /= 0) pause 'allocation error'
     ierr=0 ; allocate( N(stot,ntot),stat=ierr )
     if( ierr /= 0) pause 'allocation error'
     ierr=0 ; allocate( H(stot,ntot),stat=ierr )
     if( ierr /= 0) pause 'allocation error'
     ierr=0 ; allocate( HT(ntot,stot),stat=ierr )
     if( ierr /= 0) pause 'allocation error'
     ierr=0 ; allocate( KijT(ntot,ntot),stat=ierr ) 
     if( ierr /= 0) pause 'allocation error'
     C=0.0d0;N=0.0d0;H=0.0d0;H=0.0d0;HT=0.0d0;KijT=0.0d0
    
!    Polynomial Matrix & its transpose matrix HT
     do j=1,nnbr(i) 
      H(1,j) = LX2(1,j) - x0
      H(2,j) = LX2(2,j) - y0
     enddo
     HT=transpose(H)  ! reset & define
     
!    Covariance Matrix Kij
     call CovarT(ndim,ntot,LX2,theta,KijT)  !  Calculate Covariance Matrix 
     KijT=sigma**2*KijT                            !  Add in Noise and magnitude  
     do j=1,ntot
      KijT(j,j)=KijT(j,j)+sigmaN**2              ! diagonal modification only (noise)
     enddo
     
!    Least Square Matrix (A)
!     call invertsymmetric(ntot,KijT)               ! K^-1   [ntot,ntot]
     
     do j=1,ntot
     KijT(j,j) = KijT(j,j)**2.0d0
     enddo
     
     call matrixmult(stot,ntot,H,ntot,KijT,N)   ! H*K^-1 [stot,ntot]
     call matrixmult(stot,ntot,N,stot,HT,A)     ! H*K^-1*H^t ( Least Square Matrix ) [stot,stot]
     if( i .eq. 101 )then
     do j=1,stot
     write(*,*) ( A(j,jj),jj=1,stot )
     enddo
     pause
     endif
     call invertsymmetric(stot,A)                   ! Cholesky Inversion

!    Least Squares Coefficients
     N=0.0d0 ; call symmatrixmult(stot,A,ntot,H,N)       ! A^-1 * H           ! New N (stot,ntot)
     call matrixmult(stot,ntot,N,ntot,KijT,C)                    ! A^-1 * H^t * K^-1
     
     ax(:,i) = 0.0d0 ; ay(:,i) = 0.0d0
     do j=1,ntot
      ax(j,i) = C(1,j)
      ay(j,i) = C(2,j)
     enddo
     
!     if( i .eq. 100 )then
!     do j=1,ntot
!     write(*,*) ax(j,i)
!     write(*,*) ''
!     write(*,*) ay(j,i)
!     enddo
!     pause
!     endif
        
!    deallocate local arrays
     deallocate(C);deallocate(N);deallocate(H);deallocate(Lx2)
     deallocate(HT);deallocate(KijT)

    endif
    enddo
      
    return
!   ###################################################################
    end subroutine lsCWTLS