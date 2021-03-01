
    Module prm_flow
    implicit none 
    
    integer :: veltype    
    integer :: npts         ! Number of global points
    integer :: makeanim     ! write output I/O
    integer :: order        ! Order
    integer :: nx           ! number of points in x direction
    integer :: ny           ! number of points in y direction
    integer :: np0          ! number of inner calculating cubes / not with bc cubes
    integer :: perturb      ! add grid perturbation I/O
    integer :: animinterval 
    integer :: ictype       ! initial condition type [ SetInitCond.f90 ]
    
    double precision :: pertx  ! perturbation x direction
    double precision :: perty  ! perturbation y direction
    double precision :: dt     ! dt timestep
    double precision :: Ttime
    
    double precision :: ymin
    double precision :: ymax
    double precision :: xmin
    double precision :: xmax
    double precision :: lx    ! lx=xmax-xmin
    double precision :: ly
    
    double precision :: kkk  ! MUSCL coefficient
    double precision :: varmax0
    double precision :: varmin0
    
    double precision :: vdrho
    
!   Some Parameters
    integer :: nxmax
    integer :: nxmin
    integer :: nymin
    integer :: nymax
    integer :: npmax
    integer :: ntmax
    integer :: nnbrmax 
    
!   allocatables
    double precision,allocatable :: coord(:,:)
    double precision,allocatable :: var(:)
    double precision,allocatable :: var_old(:)
    double precision,allocatable :: ax(:,:)
    double precision,allocatable :: ay(:,:)
    double precision,allocatable :: res(:)
    double precision,allocatable :: varx(:)
    double precision,allocatable :: vary(:)
    double precision,allocatable :: varini(:)
    
    integer,allocatable :: nnum(:,:)   ! not important
    integer,allocatable :: nnbr(:)     ! number of neighbours eg. 4(face) or 8(diagonal)
    integer,allocatable :: ptype(:)    ! point type ( interior: 1 | exterior bc : 2  
    integer,allocatable :: conn(:,:)   ! neighbour cube index conn(j,ic) where i is main node, j is neighbour
    

    end module
    
    subroutine alloc
    use prm_flow
    implicit none
    
    allocate( var(npmax)        )
    allocate( var_old(npmax)    )
    allocate( ax(nnbrmax,npmax) )
    allocate( ay(nnbrmax,npmax) )
    allocate( res(npmax)        )
    allocate( varx(npmax)       )
    allocate( vary(npmax)       ) 
    
    allocate( coord(2,npmax)    )
    allocate( nnum(0:nxmax,0:nymax) ) 
    allocate(  nnbr(npmax)          ) 
    allocate( ptype(npmax)          ) 
    allocate( conn(nnbrmax,npmax)   ) 

    
    end subroutine alloc