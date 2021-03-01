    module prmflow
    use ModDataTypes
    
    real(rtype) :: hinf
    character(chrlen) :: fnGrid    !< grid and topology data
    character(chrlen) :: solin     !< restart solution - input
    character(chrlen) :: solout    !< restart solution - output
    character(chrlen) :: fnConv    ! convergence file name
    character(chrlen) :: fnSurf    ! surface data file name
    
    integer :: lrest !< use of previous solution for restart ("Y"=yes, "N"=no)
    integer :: maxiter    !< max. number of iterations
    integer :: outstep    !< number of iterations between solution dumps
    integer :: iter       !< actual iteration number

    real(rtype) :: convtol !< convergence criterion (2-norm of density change for
                           !! which the iteration process is stopped)
    
  integer :: allnod !< number of grid nodes (including dummy nodes)\n
  integer :: phynod, &  !< number of physical grid nodes (allnod - dummy nodes)
             nt, &   !< number of triangles
             nsegs, &   !< number of boundary segments (composed of boundary faces)
             nbfaces, & !< number of all boundary faces
             nbnodes    !< total number of boundary nodes

  character(chrlen), allocatable :: bname(:) !< names of boundary segments (used for plots)

  integer, allocatable :: btype(:) !< types of boundary conditions:\n
    !! 100-199 = inflow
    !! 200-299 = outflow
    !! 300-399 = viscous wall
    !! 400-499 = inviscid wall
    !! 500-599 = symmetry line
    !! 600-699 = far-field
    !! 700-799 = periodic boundary
  integer, allocatable :: bface(:,:) !< indexes of two nodes defining a face (NOT used for periodic boundaries!)
  integer, allocatable :: bnode(:,:) !< data related to boundary nodes\n
  integer, allocatable :: ibound(:,:) !< pointer from boundary segment to boundary faces and nodes\n

  real(rtype) :: xref, & !< x-coordinate of the reference point
                 yref, & !< y-coordinate of the reference point
                 cref    !< reference length or airfoil chord

  integer, allocatable :: elem(:,:) !< node indexes of triangle elements

  real(rtype), allocatable :: x(:), &  !< x-coordinates of grid points
                              y(:)     !< y-coordinates of grid points
  real(rtype), allocatable :: sij(:,:) !< x,y-components of the face vector (n*dS)\n
  real(rtype), allocatable :: sbf(:,:) !< normal vector of boundary face (outward pointing,
  real(rtype), allocatable :: sproj(:,:) !< projections of control volumes on the x- and y-axis

  
  integer, allocatable :: wbfnode(:)   ! wall node only array
  integer :: nwbfnode                  ! total number of wall nodes
  integer, allocatable :: ntype(:)        ! node type
  integer, allocatable :: ptype(:)
  
  integer,allocatable :: conn(:,:)     ! neighbours of node
  integer,allocatable :: nbers(:)      ! number of neighbours
  
  integer, allocatable :: ffbfnode(:)   ! farfield node only array
  integer :: nffbfnode
  integer :: igcs                     ! number of wall ghost cells
  integer, allocatable :: gcid(:,:)   ! [1] fluid node [2] ghost cell node
  integer :: ngcid                    ! gcid(1/2,:) global counter
  
  real(rtype), allocatable :: aij(:,:)   
  real(rtype), allocatable :: bij(:,:)
  
  real(rtype), allocatable :: nx_wc(:) 
  real(rtype), allocatable :: ny_wc(:)
  integer, allocatable :: fwnid(:) ! local only
  
  !integer, allocatable :: gcid_ibm(:,:) ! ibm instead of katz gcid
  !integer :: igcs_ibm                   ! imb number of ghost cells
  
  integer      :: ktimst    !< switch between local (="1") and global (="0") time-stepping
  integer      :: lvort     !< far-field vortex correction ("1"=yes, "0"=no)

  integer :: nedges, & !< total number of edges (including edges between boundary and dummy nodes)
             nedint, & !< number of edges excluding those to dummy nodes
             iorder, & !< order of Roe's upwind scheme (1 or 2)
             iflux,  &
             nitirs, & !< number of Jacobi iterations (implicit residual smoothing)
             nrk,    & !< number of stages (Runge-Kutta scheme); max. = 5
             ldiss(5)  !< dissipation evaluation per stage (0=no, 1=yes)

  integer, allocatable :: edge(:,:) !< edge list (node i, node j)\n


  real(rtype) :: cfl,      & !< CFL-number
                 epsirs,   & !< coefficient of implicit residual smoothing
                 epsentr,  & !< entropy correction coefficient (Roe's upwind scheme)
                 ark(5),   & !< stage coefficients
                 betrk(5)    !< dissipation-blending coefficients

  real(rtype) :: volref     ! reference volume
  real(rtype) :: limref(4)  ! reference values in limiter
  real(rtype),allocatable :: vol(:) ! FVM based Volume
  real(rtype),allocatable :: umin(:,:) !min./max. values around a node
  real(rtype),allocatable :: umax(:,:)
  real(rtype) :: limfac
  real(rtype) :: eps2(4)

  real(rtype), allocatable :: cvold(:,:), & !< conservative variables from previous time step
                              diss(:,:),  & !< artificial dissipation
                              lim(:,:),   & ! limiter
                              rhs(:,:),   & !< residual (right-hand side)
                              tstep(:)      !< time steps (without the CFL-number)

  real(rtype), allocatable :: gradx(:,:) 
  real(rtype), allocatable :: grady(:,:)
  real(rtype), allocatable :: tgcid(:,:)

  real(rtype) :: pi, & !< 3.14...
                 rad   !< 180./pi
                 
  character(1) :: iflow, & !< equations solved ("E"=Euler, "N"=Navier-Stokes)
                  soltype    !< type of flow ("E"=external, "I"=internal)

! reference values

  real(rtype) :: gamma,  & !< ratio of specific heat coefficients
                 cpgas,  & !< specific heat coefficient at constant pressure
                 prlam,  & !< laminar Prandtl number
                 renum,  & !< Reynolds number
                 refvel, & !< reference velocity (internal flow only; for external flow computed from the far-field boundary)
                 refrho, & !< reference density (internal flow only; for external flow computed from the far-field boundary)
                 refvisc   !< reference dynamic viscosity coefficient (computed from renum, refvel, cref and refrho)

! boundary conditions - external flow

  real(rtype) :: machinf, & !< Mach-number at infinity
                 alpha,   & !< angle of attack
                 pinf,    & !< static pressure at infinity
                 tinf,    & !< static temperature at infinity
                 rhoinf,  & !< density at infinity
                 uinf,    & !< u-component of velocity vector at infinity
                 vinf,    & !< v-component of velocity vector at infinity
                 qinf       !< total velocity (= SQRT(uinf**2+vinf**2))

! flow variables
  integer :: nconv, & !< number of conservative variables (cv)
             ndepv    !< number of dependent variables (dv)

  real(rtype), allocatable :: cv(:,:) !< conservative variables\n
  real(rtype), allocatable :: dv(:,:) !< dependent variables\n
  
  integer, parameter :: mxquant =13, & !< total number of plot variables
                        mxqfield=11    !< no. of plot variables in the field (cf and Cp only at the boundaries)

  character(chrlen) :: title           !< title of the simulation case

  real(rtype) :: drho    !< change of the density residual (convergence criterion)
  real(rtype) :: drho1   !< initial change of the density residual (used for normalization)
  real(rtype) :: cl      !< lift coefficient (pressure forces only; external flow)
  real(rtype) :: cd      !< drag coefficient (pressure forces only; external flow)
  real(rtype) :: cm      !< pitching moment coefficient wrp. to the reference point
                 
  integer :: vtkout    
  integer :: cpcoff
  real(rtype) :: Hfs
  
  end module 