

    module prmflow
    use edu2d_constants
    implicit none

    !Ratio of specific heats = 1.4 fpr air
    real(p2) :: gamma
    real(p2) :: M_inf
    real(p2) :: alphIn

    !Inout data files
    character(80) :: datafile_grid_in  !Grid file
    character(80) :: datafile_bcmap_in !Boundary condition specification file
    character(80) :: datafile_tec      !Tecplot file for viewing the result
    character(80) :: gridInF            !Grid Input File (.grid)
    character(80) :: bcInF

    !Scheme parameters
    character(80) :: inviscid_flux !Numerial flux for the inviscid terms (Euler)
    character(80) :: limiter_type  ! 2nd order limiter

    character(80) ::    gradient_weight ! "none" or "inverse_distance"
    real(p2)  :: gradient_weight_p  !  1.0  or any other real value

    !Parameter for explicit scheme
    real(p2) :: CFLexp              ! Input CFL number for explicit RK2 scheme
    real(p2) :: CFL                 ! Actual CFL number (set to CFLexp) - for future implicit 
    real(p2) :: tolerance           ! L1 Tolerance for steady convergence

    integer :: ngc                  ! number of gcs generated, selected w/ [nnodes+1,nnodes+ngc] 
    integer :: nIP                  ! number of IPs generated, selected w/ [nnodes+1,ngc,nnodes+ngc+nIP]
    !   specific boundary cases are selected with node(i)%ptype values

    integer :: tbcn                    ! total number of boundary nodes (preread grid) -> needed to know how many node()% are needed
    integer :: nIPGC                   ! total number of IP/GC ( 1 GC/IP per boundary node ) : tbcn - nboundaries
    integer :: igcs                    ! gcid array maximum value : provided nIP,nGC is identical
    integer, allocatable :: gcid(:,:)  ! gcid array (1: fluid,2:GC node) 
    
    integer :: outstep                 ! step iteration at which to output something (eg. flowfield)
    integer :: iter                    ! global iteration count
    integer :: maxiter                 ! self defined max number of iterations
    
    real(p2) :: Rdrho                  ! convergence L2R
    
    real(p2) :: cpgas
    real(p2) :: pinf
    real(p2) :: qinf
    real(p2) :: rhoinf
    real(p2) :: uinf
    real(p2) :: vinf
    real(p2) :: alpha
    real(p2) :: tinf
    real(p2) :: machinf
    real(p2) :: xref,yref,cref
    integer :: nconv,ndepv
    
    real(p2) :: clp,cdp ! forces
    real(p2) :: rad
    
!   cirs related
    real(p2), allocatable :: rhsold(:,:)
    real(p2), allocatable :: rhsit(:,:)
    integer, allocatable :: ncontr(:)
    integer :: nitirs
    real(p2) :: epsirs
    
    character(80) :: casename
    

    

    end module prmflow