
!* 3. module my_main_data 
!*
!* This module defines the main data that will be used in the code.
!*
!* The main data include parameters and data arrays. They can be accessed by
!* other routines via the use statement: 'use my_main_data'.

 module edu2d_my_main_data

  use edu2d_constants     , only : p2, one
  use edu2d_grid_data_type, only : node_type, elm_type, edge_type, bgrid_type, face_type, jac_type

  implicit none

  private

  public :: nnodes, node
  public :: ntria, nquad, nelms, elm
  public :: nedges, edge
  public :: nbound, bound
  public :: nfaces, face

  public :: nq
  public :: gradient_type
  public :: inviscid_jac 
  public :: rho_inf, u_inf, v_inf, p_inf 
  public :: jac
  public :: CFL1, CFL2, CFL_ramp_steps
  public :: sweeps

!  Parameters

   !Number of equtaions/variables in the target equtaion.
   integer       :: nq 

   !LSQ gradient related parameteres:
   character(80) ::     gradient_type  ! "linear"; or for node-centered schemes can use "quadratic2"

   !Reference quantities
   real(p2) :: rho_inf, u_inf, v_inf, p_inf

   !Implicit solver parameters
    character(80) :: inviscid_jac        ! Inviscid flux for Jacobian
    real(p2)      :: CFL1, CFL2          ! Initial and terminal CFL number for ramping
    integer       :: CFL_ramp_steps      ! Number of iterations to reach CFL2 from CFL1
    integer       :: sweeps              ! Number of GS relaxation

!  Node data
   integer                                 :: nnodes !total number of nodes
   type(node_type), dimension(:), pointer  :: node   !array of nodes

!  Element data (element=cell)
   integer                                 :: ntria  !total number of triangler elements
   integer                                 :: nquad  !total number of quadrilateral elements
   integer                                 :: nelms  !total number of elements
   type(elm_type),  dimension(:), pointer  :: elm    !array of elements

!  Edge data
   integer                                 :: nedges !total number of edges
   type(edge_type), dimension(:), pointer  :: edge   !array of edges

!  Boundary data
   integer                                 :: nbound !total number of boundary types
   type(bgrid_type), dimension(:), pointer :: bound  !array of boundary segments

!  Face data (cell-centered scheme only)
   integer                                 :: nfaces !total number of cell-faces
   type(face_type), dimension(:), pointer  :: face   !array of cell-faces

!  Matrix data
   type(jac_type), dimension(:)  , pointer  :: jac   !array of edges

 end module edu2d_my_main_data
!********************************************************************************
