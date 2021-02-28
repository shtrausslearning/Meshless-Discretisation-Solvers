
!* 2. module grid_data_type
!*
!* This module defines custom grid data types for unstructured grids
!*
!* NOTE: These data types are designed to make it easier to understand the code.
!*       They may not be the best in terms of efficiency.
!*
!* NOTE: Custom grid data types (derived types) are very useful.
!*       For example, if I declare a variable, "a", by the statemant:
!*           type(node_type), dimension(100) :: a
!*       The variable, a, is a 1D array each component of which contains all data
!*       defined as below. These data can be accessed by %, e.g.,
!*           a(1)%x, a(1)%y, a(1)%nghbr(1:nnghbrs), etc.
!*       In C-programming, this type of data is called "structure", I think.
!*
!*
!* NOTE: Not all data types below are used in EDU2D-Euler-Implicit.

 module edu2d_grid_data_type

  use edu2d_constants, only : p2

  implicit none

  private

  public ::  node_type
  public ::   elm_type
  public ::  edge_type
  public :: bgrid_type
  public ::  face_type
  public ::   jac_type

!----------------------------------------------------------
! Data type for nodal quantities (used for node-centered schemes)
! Note: Each node has the following data.
!----------------------------------------------------------
  type node_type
!  to be read from a grid file
   real(p2)                          :: x, y      !nodal coordinates
!  to be constructed in the code
   integer                           :: nnghbrs   !number of neighbors
   integer,   dimension(:), pointer  :: nghbr     !list of neighbors
   integer                           :: nelms     !number of elements
   integer,   dimension(:), pointer  :: elm       !list of elements
   real(p2)                          :: vol       !dual-cell volume
   integer                           :: bmark     !Boundary mark
   integer                           :: nbmarks   !# of boundary marks
!  to be computed in the code
   !Below are arrays always allocated.
   real(p2), dimension(:)  , pointer :: cv         !conservative variables n
   real(p2), dimension(:)  , pointer :: cvold      !conservative variables n-1
   real(p2), dimension(:)  , pointer :: dv         !dependent variables n
   real(p2), dimension(:)  , pointer :: uexact    !conservative variables
   real(p2), dimension(:,:), pointer :: gradu     !gradient of u
   real(p2), dimension(:)  , pointer :: res       !residual (rhs)
   real(p2), dimension(:)  , pointer :: diss       !residual (rhs)
   real(p2)                          :: ar        ! Control volume aspect ratio
   real(p2), dimension(:)  , pointer :: aij       !    Linear LSQ coefficient for wx
   real(p2), dimension(:)  , pointer :: bij       !     Linear LSQ coefficient for wy
   real(p2), dimension(:)  , pointer :: hoaij     ! Quadratic LSQ coefficient for wx
   real(p2), dimension(:)  , pointer :: hobij     ! Quadratic LSQ coefficient for wy
   real(p2), dimension(:  ), pointer :: dx, dy    ! Extra data used by Quadratic LSQ
   real(p2), dimension(:,:), pointer :: dw        ! Extra data used by Quadratic LSQ
   integer                           :: ptype     ! node ptype

   !Below are optional: Pointers need to be allocated in the main program if necessary.
   real(p2), dimension(:)  , pointer :: du        !change in conservative variables
   real(p2), dimension(:)  , pointer :: w         !primitive variables(optional)
   real(p2), dimension(:,:), pointer :: gradw     !gradient of w
   real(p2)                          :: phi       !limiter function (0 <= phi <= 1)
   real(p2)                          :: dt        !local time step
   real(p2)                          :: wsn       !Half the max wave speed at face
   real(p2), dimension(:), pointer   :: r_temp    ! For GCR implementation
   real(p2), dimension(:), pointer   :: u_temp    ! For GCR implementation
   real(p2), dimension(:), pointer   :: w_temp    ! For GCR implementation
   real(p2), dimension(:), pointer   :: rhs       ! rhs term

   end type node_type

!----------------------------------------------------------
! Data type for element/cell quantities (used for cell-centered schemes)
! Note: Each element has the following data.
!----------------------------------------------------------
  type elm_type
!  to be read from a grid file
   integer                           :: nvtx     !number of vertices
   integer,   dimension(:), pointer  :: vtx      !list of vertices
!  to be constructed in the code
   integer                           :: nnghbrs  !number of neighbors
   integer,   dimension(:), pointer  :: nghbr    !list of neighbors
   real(p2)                          :: x, y     !cell center coordinates
   real(p2)                          :: vol      !cell volume

   integer,  dimension(:)  , pointer :: edge     !list of edges
   real(p2), dimension(:)  , pointer :: u        !conservative variables
   real(p2), dimension(:)  , pointer :: uexact   !conservative variables
!NotUsed   real(p2), dimension(:)  , pointer :: du       !change in conservative variables
   real(p2), dimension(:,:), pointer :: gradu    !gradient of u
   real(p2), dimension(:)  , pointer :: res      !residual (rhs)
   real(p2)                          :: dt       !local time step
   real(p2)                          :: wsn      !
   integer                           :: bmark    !Boundary mark
   integer                           :: nvnghbrs !number of vertex neighbors
   integer,  dimension(:), pointer   :: vnghbr   !list of vertex neighbors
   real(p2)                          :: ar       !Element volume aspect ratio
   real(p2), dimension(:) , pointer  :: aij      !Linear LSQ coefficient for ux
   real(p2), dimension(:) , pointer  :: bij      !Linear LSQ coefficient for uy

  end type elm_type

!----------------------------------------------------------
! Data type for edge quantities (used for node-centered scheemes)
! Note: Each edge has the following data.
!----------------------------------------------------------
  type edge_type
!  to be constructed in the code
   integer                          :: n1, n2 !associated nodes
   integer                          :: e1, e2 !associated elements
   real(p2),           dimension(2) :: dav    !unit directed-area vector
   real(p2)                         :: da     !magnitude of the directed-area vector
   real(p2),           dimension(2) :: ev     !unit edge vector
   real(p2)                         :: e      !magnitude of the edge vector
   integer                          :: kth_nghbr_of_1 !neighbor index
   integer                          :: kth_nghbr_of_2 !neighbor index
   real(p2), dimension(:), pointer  :: eaij    ! edge based aij grouping
   real(p2), dimension(:), pointer  :: ebij    ! edge based bij grouping
   
   
  end type edge_type

!----------------------------------------------------------
! Data type for boundary quantities (for both node/cell-centered schemes)
! Note: Each boundary segment has the following data.
!----------------------------------------------------------
  type bgrid_type
!  to be read from a boundary grid file
   character(80)                    :: bc_type !type of boundary condition
   integer                          :: nbnodes !# of boundary nodes
   integer,   dimension(:), pointer :: bnode   !list of boundary nodes
!  to be constructed in the code
   integer                          :: nbfaces !# of boundary faces
   real(p2),  dimension(:), pointer :: bfnx    !x-component of the face outward normal
   real(p2),  dimension(:), pointer :: bfny    !y-component of the face outward normal
   real(p2),  dimension(:), pointer :: bfn     !magnitude of the face normal vector
   real(p2),  dimension(:), pointer :: bnx     !x-component of the outward normal
   real(p2),  dimension(:), pointer :: bny     !y-component of the outward normal
   real(p2),  dimension(:), pointer :: bn      !magnitude of the normal vector
   integer ,  dimension(:), pointer :: belm    !list of elm adjacent to boundary face
   integer ,  dimension(:), pointer :: kth_nghbr_of_1
   integer ,  dimension(:), pointer :: kth_nghbr_of_2
  end type bgrid_type

!----------------------------------------------------------
! Data type for face quantities (used for cell-centered schemes)
!
! A face is defined by a line segment connecting two nodes.
! The directed area is defined as a normal vector to the face,
! pointing in the direction from e1 to e2.
!
!      n2
!       o------------o
!     .  \         .
!    .    \   e2  .
!   .  e1  \    .
!  .        \ .         Directed area is positive: n1 -> n2
! o----------o         e1: left element
!             n1       e2: right element (e2 > e1 or e2 = 0)
!
! Note: Each face has the following data.
!----------------------------------------------------------
  type face_type
! to be constructed in the code (NB: boundary faces are excluded.)
   integer                         :: n1, n2 !associated nodes
   integer                         :: e1, e2 !associated elements
   real(p2),          dimension(2) :: dav    !unit directed-area vector
   real(p2)                        :: da     !magnitude of the directed-area vector
  end type face_type


  type jac_type
!  to be constructed in the code
   real(p2),  dimension(:,:)  , pointer :: diag ! diagonal block of Jacobian matrix
   real(p2),  dimension(:,:,:), pointer :: off  ! off-diagonal block of Jacobian matrix
  end type jac_type


 end module edu2d_grid_data_type
!********************************************************************************