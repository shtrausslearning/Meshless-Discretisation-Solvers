
!********************************************************************************
!* Read the grid and the exact solution.
!* ------------------------------------------------------------------------------
!*  Input: datafile_grid_in  = filename of the grid file
!*         datafile_bcmap_in = filename of the bc file
!*
!* Output: nnodes, ncells, node(:), elm(:), bound(:) = data used in the solver
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
!* 1. "datafile_grid_in" is assumed to have been written in the following format:
!*
!*   -----------------------------------------------------------------------
!*    write(*,*) nnodes, ntria, nquad !Numbers of nodes, triangles and quads
!*
!*   do i = 1, nnodes
!*    write(*,*) x(i), y(i) !(x,y) coordinates of each node
!*   end do
!*
!*   do i = 1, ntria        !Nodes of triangles ordered counterclockwise
!*    write(*,*) node_1(i), node_2(i), node_3(i)
!*   end do
!*
!*   do i = 1, nquad        !Nodes of quadrilaterals ordered counterclockwise
!*    write(*,*) node_1(i), node_2(i), node_3(i), node_4(i)
!*   end do
!* 
!*    write(*,*) nbound     !Number of boundary segments
!*
!*   do i = 1, nbound
!*    write(*,*) nbnodes(i) !Number of nodes on each segment
!*   end do
!*
!*   do i = 1, nbound
!*    do j = 1, nbnodes(i)
!*     write(*,*) bnode(j)  !Node number of each node j in segment i
!*    end do
!*   end do
!*   -----------------------------------------------------------------------
!*
!*   NOTE: Add the first node to the end if the segment is closed
!*         (e.g., airfoil) The number of nodes will be the actual number + 1
!*         in that case.
!*
!*   NOTE: Boundary nodes must be ordered such that the domain is on the left.
!*
!********************************************************************************
!*
!* 2. "datafile_bcmap_in" is assumed have been written in the following format:
!*
!*   -----------------------------------------------------------------------
!*    write(*,*) "Boundary Segment              Boundary Condition"
!*   do i = 1, nbound
!*    write(*,*) i, bc_name
!*   end do
!*   -----------------------------------------------------------------------
!*
!*   NOTE: bc_name is the name of the boundary condition, e.g.,
!*
!*         1. "freestream"
!*             Roe flux with freestream condition on the right state.
!*
!*         2. "slip_wall"
!*             Solid wall condition. Mass flux through the boundary is set zero.
!*
!*         3. "outflow_supersonic"
!*             Just compute the boundary flux by the physical Euler flux
!*             (equivalent to the interior-extrapolation condition.)
!*
!*         4. "outflow_back_pressure"
!*             Fix the back pressure. This should work for subsonic flows in a
!*             large enough domain.
!*
!*         Something like the above needs to be implemented in a solver.
!*
!********************************************************************************
!* Data to be read and stored:
!*
!* 1. Some numbers
!*    nnodes        = Number of nodes
!*    ntria         = Number of triangular elements
!*    nquad         = Number of quadrilateral elements
!*    nelms         = Total number of elements (=ntria+nquad)
!*
!* 2. Element data:
!*    elm(1:nelms)%nvtx   =  Number of vertices of each element
!*    elm(1:nelms)%vtx(:) = Pointer to vertices of each element
!*
!* 3. Node data: nodes are stored in a 1D array
!*    node(1:nnodes)%x     = x-coordinate of the nodes
!*    node(1:nnodes)%y     = y-coordinate of the nodes
!*
!* 4. Boundary Data:
!*    nbound                   = Number of boundary segments
!*    bound(1:nbound)%nbnodes  = Number of nodes in each segment
!*    bound(1:nbound)%bnode(:) = List of node numbers for each segment
!*    bound(1:nbound)%bc_type  = Boundary condition name for each segment
!*    bound(1:nbound)%bc_type  = Boundary condition name for each segment
!*
!********************************************************************************

 subroutine read_grid

 use edu2d_my_main_data, only : nnodes, node, ntria, nquad, nelms, elm, nbound, bound
 use prmflow

 implicit none

!Local variables
 integer  :: i, j, os, dummy_int

!--------------------------------------------------------------------------------
! 1. Read grid file>: datafile_grid_in

  write(101,*) "Reading the grid file....", gridInF
  
!  Open the input file.
   open(unit=1, file='./grid/'//gridInF, status="unknown", iostat=os)

! READ: Get the size of the grid.
  read(1,*) nnodes, ntria, nquad
  nelms = ntria + nquad   ! total number of elements

!  Allocate node and element arrays.
   allocate(node(nnodes+500))
   allocate(elm(  nelms))

! READ: Read the nodal coordinates
  do i = 1, nnodes
   read(1,*) node(i)%x, node(i)%y
  end do

! Read element-connectivity information

! Triangles: assumed that the vertices are ordered counterclockwise
!
!         v3
!         /\
!        /  \
!       /    \
!      /      \
!     /        \
!    /__________\
!   v1           v2

! READ: read connectivity info for triangles
  if ( ntria > 0 ) then
   do i = 1, ntria
    elm(i)%nvtx = 3
    allocate(elm(i)%vtx(3))
    read(1,*) elm(i)%vtx(1), elm(i)%vtx(2), elm(i)%vtx(3)
   end do
  endif

! Quads: assumed that the vertices are ordered counterclockwise
!
!        v4________v3
!         /        |
!        /         |
!       /          |
!      /           |
!     /            |
!    /_____________|
!   v1             v2

! READ: read connectivity info for quadrilaterals
  if ( nquad > 0 ) then
   do i = 1, nquad
    elm(ntria+i)%nvtx = 4
    allocate( elm(ntria+i)%vtx(4))
    read(1,*) elm(ntria+i)%vtx(1), elm(ntria+i)%vtx(2), &
              elm(ntria+i)%vtx(3), elm(ntria+i)%vtx(4)
   end do
  endif

!  Write out the grid data.

   write(101,*)
   write(101,*) " Total numbers:"
   write(101,*) "      nodes = ", nnodes
   write(101,*) "  triangles = ", ntria
   write(101,*) "      quads = ", nquad
   write(101,*)

! Read the boundary grid data

! READ: Number of boundary condition types
  read(1,*) nbound
  allocate(bound(nbound))

! READ: Number of Boundary nodes (including the starting one at the end if
! it is closed such as an airfoil.)
  do i = 1, nbound
   read(1,*) bound(i)%nbnodes
   allocate(bound(i)%bnode(bound(i)%nbnodes))
  end do

! READ: Read boundary nodes
  do i = 1, nbound
   do j = 1, bound(i)%nbnodes
   read(1,*) bound(i)%bnode(j)
   end do
  end do

!  Print the boundary grid data.
   write(101,*) " Boundary nodes:"
   write(101,*) "    segments = ", nbound
    do i = 1, nbound
     write(101,'(a9,i3,2(a11,i5))') " boundary", i, "  bnodes = ", bound(i)%nbnodes, &
                                                  "  bfaces = ", bound(i)%nbnodes-1
    end do
   write(101,*) ''

  close(1)
  
! End of Read grid file>: datafile_grid_in
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
! 2. Read the boundary condition data file

   write(101,*)
   write(101,*) "Reading the boundary condition file....", bcInF
   write(101,*)

! Open the input file.
  open(unit=2, file='./grid/'//bcInF, status="unknown", iostat=os)

    read(2,*) 

! READ: Read the boundary condition type
  do i = 1, nbound
    read(2,*) dummy_int, bound(i)%bc_type
   end do

!  Print the data
    write(101,*) " Boundary conditions:"
   do i = 1, nbound
    write(101,'(a9,i3,a12,a35)') " boundary", i, "  bc_type = ", trim(bound(i)%bc_type)
   end do

    i = dummy_int !Never mind. Just to avoid a compilation warning.

    write(101,*)

  close(2)

! End of Read the boundary condition data file
!--------------------------------------------------------------------------------

 end subroutine read_grid