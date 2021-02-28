
!   [edge based WTLS]
!   subroutine to generate WTLS meshless coefficients with weighting

    subroutine eOWTLS
!   ###################################################################
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,nnodes,edge,nedges
    use edu2d_my_allocation , only : my_alloc_p2_ptr
    implicit none
    
    integer :: OMPiD(2) ! omp node information
    integer :: cnID(2),lcn,i,j
    integer :: ie,nbrM,ierr,oP,idL,ntot,jj,inn
    real(p2) :: c1,c2,x0,y0,dr
    real(p2) :: mAij,mBij,nAij,nBij
    real(p2) :: WFpwrDP,xl1,xl2,yl1,yl2
    real(p2) :: mNx,mNy,nNx,nNy
    real(p2) :: x1,y1,nds,mds
    real(p2) :: A(2,2)   ! Least Square Matrix

!   local varying arrays
    real(p2),allocatable :: Kmat(:,:)    
    real(p2),allocatable :: LX(:,:)      ! ndim,ntot
    real(p2),allocatable :: HT(:,:)
    real(p2),allocatable :: H(:,:)       ! Polynomial dx,dy    
    real(p2),allocatable :: C(:,:)
    real(p2),allocatable :: N(:,:)
    real(p2),allocatable :: w(:)
    

!   (1) first allocate all necessary aij/bij components
    do i=1,nnodes
     ntot = node(i)%nnghbrs            
     call my_alloc_p2_ptr(node(i)%aij,ntot) ! dynamically allocate array aij
     call my_alloc_p2_ptr(node(i)%bij,ntot) ! dynamically allocate array bij
    enddo
    
!   (2) Edge based coefficient evaluation for physical domain only

    OMPiD=0 ; cnID=0
    do ie=1,nedges ! do for all physical edges
    
!    physical nodes
     ompID(1) = edge(ie)%n1   ! 1st node of edge
     ompID(2) = edge(ie)%n2   ! 2nd node of edge
     
!    find and store 

     do jj=1,node(ompID(1))%nnghbrs  ! look though all of point n1's neighbours
     inn = node(ompID(1))%nghbr(jj)
     if( inn .eq. edge(ie)%n2 )then
     cnid(1) = jj
     endif
     enddo
    
     do jj=1,node(ompID(2))%nnghbrs
     inn = node(ompID(2))%nghbr(jj)
     if( inn .eq. edge(ie)%n1 )then
     cnid(2) = jj
     endif
     enddo

     do oP = 1,2  ! if a coefficient hasn't been calculated yet
     
!     Local main stuff
      idL = ompID(oP)    ! local node in question 
      ntot = node(idL)%nnghbrs    
      lcn = cnID(oP)
      x0 = node(idL)%x  ! x coordinate of node in question
      y0 = node(idL)%y  ! y coordinate of node in question
      
!     allocate local arrays
      ierr=0;allocate( C(2,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( N(2,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( H(2,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( HT(ntot,2),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( Kmat(ntot,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( w(ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      C=0.0d0;N=0.0d0;H=0.0d0;HT=0.0d0;Kmat=0.0d0;w(:)=0.0d0
     
!     define the polynomial matrix/collocation matrix
      do jj=1,ntot
      inn = node(idL)%nghbr(jj)
      H(1,jj) = node(inn)%x - x0
      H(2,jj) = node(inn)%y - y0
      dr = dsqrt( H(1,jj)**2 + H(2,jj)**2 ) ! local dr
      w(jj) = 1.0d0/(dr**dble(gradient_weight_p))  ! local diagonal weight
      enddo
!     transpose the collocation matrix
      HT=transpose(H)  ! reset & define
      
!     define the weighting matrix
      do jj=1,ntot
       Kmat(jj,jj) = 1.0d0
      enddo
      
!     Matrix Multiplications to get LS coefficients
      call matrixmult(2,ntot,H,ntot,Kmat,N)     ! Matrix N = H*W [stot,ntot]
      call matrixmult(2,ntot,N,2,HT,A)          ! N*H^t ( Least Square Matrix ) [stot,stot]
      call invertsymmetric(2,A)                 ! Invert Least Square Matrix (cholesky)
      call symmatrixmult(2,A,ntot,N,C)          ! A^-1 * (HW) ! gridless coefficients (aij,bij)
      
!     store edge coefficients
      node(idL)%aij(lcn) = C(1,lcn)
      node(idL)%bij(lcn) = C(2,lcn)
      
!     deallocate local arrays
      deallocate(C);deallocate(N);deallocate(H)
      deallocate(HT);deallocate(Kmat);deallocate(w)
    
     enddo
     
    enddo  ! do ie=1,nedint

!   (3) generate additional components for ghost nodes 

    call lsTLSffGC
    call lsTLSwgc 

    return
!   ###################################################################
    end subroutine eOWTLS
    
!   Additional Components needed for Ghost Nodes

    subroutine lsTLSwgc
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
    if( node(i)%ptype .eq. 1 )then    ! if wall node
       
!    define centre node data
     x0 = node(i)%x
     y0 = node(i)%y      
     ntot = node(i)%nnghbrs            

!    store local neighbour node placements
     ierr=0;allocate( lx(2,ntot),stat=ierr )
     if (ierr /= 0) pause 'lstls alloc error'
     lx=0.0d0
     
     do j=1,ntot           
      inn = node(i)%nghbr(j)   ! neighbour node id
      lx(1,j) = node(inn)%x      ! x coordinate of neighbour
      lx(2,j) = node(inn)%y      ! y coordinate of neighbour
     enddo  
     
!     allocate local arrays
      ierr=0;allocate( C(2,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( N(2,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( H(2,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( HT(ntot,2),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( Kmat(ntot,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( w(ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      C=0.0d0;N=0.0d0;H=0.0d0;HT=0.0d0;Kmat=0.0d0;w(:)=0.0d0
     
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
     call matrixmult(stot,ntot,N,ntot,Kmat,C)       ! A^-1 * H^t * K^-1
     
     do j=1,ntot
      inn = node(i)%nghbr(j)
      if( node(inn)%ptype .eq. 4 )then ! if neighbour is wall ghost node
      node(i)%aij(j) = C(1,j)
      node(i)%bij(j) = C(2,j)
      endif
     enddo
    
!    deallocate local arrays
     deallocate(C);deallocate(N);deallocate(H);deallocate(lx)
     deallocate(HT);deallocate(Kmat);deallocate(w)

    endif
    enddo

    return
!   ###################################################################
    end subroutine lsTLSwgc
    
    subroutine lsTLSffGC
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
    if( node(i)%ptype .eq. 2 )then    ! if fardield wall node
       
!    define centre node data
     x0 = node(i)%x
     y0 = node(i)%y      
     ntot = node(i)%nnghbrs            

!    store local neighbour node placements
     ierr=0;allocate( lx(2,ntot),stat=ierr )
     if (ierr /= 0) pause 'lstls alloc error'
     lx=0.0d0
     
     do j=1,ntot           
      inn = node(i)%nghbr(j)   ! neighbour node id
      lx(1,j) = node(inn)%x      ! x coordinate of neighbour
      lx(2,j) = node(inn)%y      ! y coordinate of neighbour
     enddo  
     
!     allocate local arrays
      ierr=0;allocate( C(2,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( N(2,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( H(2,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( HT(ntot,2),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( Kmat(ntot,ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      ierr=0;allocate( w(ntot),stat=ierr )
      if (ierr /= 0) pause 'lsetls alloc error'
      C=0.0d0;N=0.0d0;H=0.0d0;HT=0.0d0;Kmat=0.0d0;w(:)=0.0d0
     
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
     call matrixmult(stot,ntot,N,ntot,Kmat,C)       ! A^-1 * H^t * K^-1
     
     do j=1,ntot
      inn = node(i)%nghbr(j)
      if( node(inn)%ptype .eq. 5 )then ! if neighbour is farfield ghost node
      node(i)%aij(j) = C(1,j)
      node(i)%bij(j) = C(2,j)
      endif
     enddo
    
!    deallocate local arrays
     deallocate(C);deallocate(N);deallocate(H);deallocate(lx)
     deallocate(HT);deallocate(Kmat);deallocate(w)

    endif
    enddo

    return
!   ###################################################################
    end subroutine lsTLSffGC