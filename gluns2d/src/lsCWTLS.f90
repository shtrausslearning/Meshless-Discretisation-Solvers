    subroutine lsCWTLS
!   ###################################################################
    use prmflow
    use covars
    implicit none

    integer :: i, j, nn,ierr,ii
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
    double precision,allocatable :: LX(:,:)      ! ndim,ntot
    double precision,allocatable :: HT(:,:)      ! Polynomial 
    double precision,allocatable :: H(:,:)       ! Polynomial dx,dy    
    double precision,allocatable :: C(:,:)
    double precision,allocatable :: N(:,:)
  
    ierr=0;allocate( aij(12,phynod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for bij()" )
    ierr=0;allocate( bij(12,phynod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for bij()" )
    aij(:,:) = 0.0d0 ; bij(:,:) = 0.0d0 
    
    open(5,file='kriging.dat',form='formatted')
    read(5,*) ( theta(i),i=1,ndim )
    read(5,*) sigma
    read(5,*) sigmaN
    read(5,*) covarflag
    close(5)

    do i=1,phynod
       
     x0 = x(i) ; y0 = y(i)              ! centre point
     ntot = nbers(i)                    ! number of neighbours of current node

     ierr=0;allocate( LX(ndim,ntot),stat=ierr )
     if (ierr /= 0) call EM( "cannot allocate memory for LX()" )
     
     do j=1,ntot           
      nn    = conn(j,i) ! neighbour id   
      LX(1,j) = x(nn)   ! neighbour global X
      LX(2,j) = y(nn)   ! neighbour global Y
     enddo  
        
!    allocate local arrays
     ierr=0 ; allocate( C(stot,ntot),stat=ierr )
     if (ierr /= 0) call EM( "cannot allocate memory for C()" )
     ierr=0 ; allocate( N(stot,ntot),stat=ierr )
     if (ierr /= 0) call EM( "cannot allocate memory for N()" )
     ierr=0 ; allocate( H(stot,ntot),stat=ierr )
     if (ierr /= 0) call EM( "cannot allocate memory for H()" ) 
     ierr=0 ; allocate( HT(ntot,stot),stat=ierr )
     if (ierr /= 0) call EM( "cannot allocate memory for HT()" )
     ierr=0 ; allocate( KijT(ntot,ntot),stat=ierr ) 
     if( ierr /= 0) call EM( "cannot allocate memory for Kmat()" )
     C=0.0d0;N=0.0d0;H=0.0d0;H=0.0d0;HT=0.0d0;KijT=0.0d0
    
!    Polynomial Matrix & its transpose matrix HT
     do j=1,nbers(i) 
      H(1,j) = LX(1,j) - x0
      H(2,j) = LX(2,j) - y0
     enddo
     HT=transpose(H)  ! reset & define
     
!    Covariance Matrix Kij
     call CovarT(ndim,ntot,LX,theta,KijT)  !  Calculate Covariance Matrix 
     KijT=sigma**2*KijT                            !  Add in Noise and magnitude  
     do j=1,ntot
      KijT(j,j)=KijT(j,j)+sigmaN**2              ! diagonal modification only (noise)
     enddo
     
!    Least Square Matrix (A)
     call invertsymmetric(ntot,KijT)               ! K^-1   [ntot,ntot]
     call matrixmult(stot,ntot,H,ntot,KijT,N)   ! H*K^-1 [stot,ntot]
     call matrixmult(stot,ntot,N,stot,HT,A)     ! H*K^-1*H^t ( Least Square Matrix ) [stot,stot]
     call invertsymmetric(stot,A)                   ! Cholesky Inversion

!    Least Squares Coefficients
     N=0.0d0 ; call symmatrixmult(stot,A,ntot,H,N)       ! A^-1 * H           ! New N (stot,ntot)
     call matrixmult(stot,ntot,N,ntot,KijT,C)                    ! A^-1 * H^t * K^-1
     
     aij(:,i) = 0.0d0 ; bij(:,i) = 0.0d0
     do j=1,ntot
      aij(j,i) = C(1,j);bij(j,i) = C(2,j)
     enddo
     
!    deallocate local arrays
     deallocate(C);deallocate(N);deallocate(H);deallocate(Lx)
     deallocate(HT);deallocate(KijT)

    enddo
      
    return
!   ###################################################################
    end subroutine lsCWTLS