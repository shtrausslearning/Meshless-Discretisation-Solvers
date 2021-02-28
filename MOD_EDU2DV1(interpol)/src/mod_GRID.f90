
 module edu2d_grid_data

 private
 public :: construct_grid_data
 public :: check_grid_data

 contains

!********************************************************************************
!* Construct the grid data:
!*
!* The following data, needed for NCFV method, will be constructed based on the
!* data read from the grid file.
!*
!* 1. Element data:
!*    elm(:)%nnghbrs  = Number of element neighbors of each element
!*    elm(:)%nghbr(:) = List of element neighbors of each element
!*    elm(:)%x        = x-coordinate of the centroid
!*    elm(:)%y        = y-coordinate of the centroid
!*    elm(:)%vol      = Volume of the element
!*
!*
!* 2. Node data:
!*    node(:)%nnghbrs = Number of node neighbors of each node
!*    node(:)%nghbr(:)= List of node neighbors of each node
!*    node(:)%nelms   = Number of adjacent elements of each node
!*    node(:)%elm     = List of adjacent elements of each node
!*    node(:)%vol     = Volume of the dual volume around each node
!*
!* 3. Edge data:
!*    edge(:)%n1, n2  = End nodes of each edge (edge points n1 -> n2)
!*    edge(:)%e1, e2  = Left and right elements of each edge
!*    edge(:)%dav     = Unit directed area vector of each edge
!*    edge(:)%da      = Magnitude of the directed area vector for each edge
!*    edge(:)%ev      = Unit edge vector of each edge (vector n1 -> n2)
!*    edge(:)%e       = Magnitude of the edge vector for each edge
!*
!*
!* 4. Boudnary data
!*    bound(:)%bnx    = Outward normal at boundary nodes (x-component of unit vector)
!*    bound(:)%bny    = Outward normal at boundary nodes (y-component of unit vector)
!*    bound(:)%bn     = Magnitude of (bnx,bny)
!*    NOTE: In this code, the above normal vector at boundary nodes is computed by
!*          a quadratic fit. It is sufficiently accuarte for 3rd-order schemes.
!*          See http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2015v281pp518-555_preprint.pdf
!*          for details on the quadratic approximation for computing more accurate normals.
!*    bound(:)%bfnx   = Outward normal at boundary nodes (x-component of unit vector)
!*    bound(:)%bfny   = Outward normal at boundary nodes (y-component of unit vector)
!*    bound(:)%bfn    = Magnitude of (bfnx,bfny)
!*    bound(:)%belm   = Element to which the boundary face belongs
!*
!********************************************************************************
 subroutine construct_grid_data

 use edu2d_my_main_data , only : nnodes, node, nelms, elm, nedges, edge, nbound, bound, face, nfaces
 use edu2d_constants    , only : p2, zero, half, third
 use edu2d_my_allocation, only : my_alloc_int_ptr, my_alloc_p2_ptr, my_alloc_p2_matrix_ptr

 implicit none

!Local variables
 integer  ::  i, j, k, ii, in, im, jelm, v1, v2, v3, v4
 real(p2) :: x1, x2, x3, x4, y1, y2, y3, y4, xm, ym, xc, yc
 real(p2) :: xj, yj, xm1, ym1, xm2, ym2, dsL,dsR,dx,dy
 logical  :: found
 integer  :: vL, vR, n1, n2, e1, e2
 integer  :: vt1, vt2, ielm

 integer  :: ave_nghbr, min_nghbr, max_nghbr, imin, imax

 integer :: iedge

 real(p2)                          :: ds

! Some initialization
 v2 = 0
 vL = 0
 im = 0
 jelm = 0

  write(101,*) "Constructing grid data...."

! Initializations
  do i = 1, nnodes
   node(i)%nelms = 0
  end do
   nedges = 0

!--------------------------------------------------------------------------------
! Loop over elements and construct the fololowing data.
!
! 1. Surrounding elements: node(:)%nelms, node(:)%elm(:)
!
!    Example: Node i is surrounded by the eleemnts, 23, 101, 13, 41.
!             node(i)%nelms = 4
!             node(i)%elm(1) = 23
!             node(i)%elm(2) = 13
!             node(i)%elm(3) = 41
!             node(i)%elm(4) = 101
!
!        o-------o-------------o
!       /        |   .         |
!      /    23   |      41     |
!     o----------o-------------o
!      \        i \            |
!       \   101    \     13    |
!        \          \          | 
!         o----------o---------o
!
! 2. Element quantities  : elm(:)%x,elm(:)%y,elm(:)%vol
!
!  o-----------o            
!   \          |            o
!    \    (x,y)|           / \
!     \   .    |          /   \
!      \       |         /  .  \    (x,y): centroid coordinates
!       \      |        / (x,y) \     vol: volume of element
!        o-----o       o---------o

  elements : do i = 1, nelms

   v1 = elm(i)%vtx(1)
   v2 = elm(i)%vtx(2)
   v3 = elm(i)%vtx(3)

   x1 = node(v1)%x
   x2 = node(v2)%x
   x3 = node(v3)%x

   y1 = node(v1)%y
   y2 = node(v2)%y
   y3 = node(v3)%y

! Distribute the element index to nodes.

   node(v1)%nelms = node(v1)%nelms + 1
   call my_alloc_int_ptr(node(v1)%elm, node(v1)%nelms)
   node(v1)%elm(node(v1)%nelms) = i

   node(v2)%nelms = node(v2)%nelms + 1
   call my_alloc_int_ptr(node(v2)%elm, node(v2)%nelms)
   node(v2)%elm(node(v2)%nelms) = i

   node(v3)%nelms = node(v3)%nelms + 1
   call my_alloc_int_ptr(node(v3)%elm, node(v3)%nelms)
   node(v3)%elm(node(v3)%nelms) = i

! Compute the cell center and cell volume.
   tri_or_quad : if (elm(i)%nvtx==3) then

!   Triangle centroid and volume
    elm(i)%x   = third*(x1+x2+x3)
    elm(i)%y   = third*(y1+y2+y3)
    elm(i)%vol = tri_area(x1,x2,x3,y1,y2,y3)

   elseif (elm(i)%nvtx==4) then

!   OK, this is a quad. Get the 4th vertex.
    v4 = elm(i)%vtx(4)
    x4 = node(v4)%x
    y4 = node(v4)%y
!   Centroid: median dual
!   (Note: There is an alternative. See Appendix B in Nishikawa AIAA2010-5093.)
    xm1 = half*(x1+x2)
    ym1 = half*(y1+y2)
    xm2 = half*(x3+x4)
    ym2 = half*(y3+y4)
    elm(i)%x   = half*(xm1+xm2)
    elm(i)%y   = half*(ym1+ym2)
!   Volume is computed as a sum of two triangles: 1-2-3 and 1-3-4.
    elm(i)%vol = tri_area(x1,x2,x3,y1,y2,y3) + tri_area(x1,x3,x4,y1,y3,y4)

     xc = elm(i)%x
     yc = elm(i)%y
    if (tri_area(x1,x2,xc,y1,y2,yc)<zero) then
     write(101,*) " Centroid outside the quad element 12c: i=",i
     write(101,'(a10,2es10.2)') "  (x1,y1)=",x1,y1
     write(101,'(a10,2es10.2)') "  (x2,y2)=",x2,y2
     write(101,'(a10,2es10.2)') "  (x3,y3)=",x3,y3
     write(101,'(a10,2es10.2)') "  (x4,y4)=",x4,y4
     write(101,'(a10,2es10.2)') "  (xc,yc)=",xc,yc
     stop
    endif

    if (tri_area(x2,x3,xc,y2,y3,yc)<zero) then
     write(101,*) " Centroid outside the quad element 23c: i=",i
     write(101,'(a10,2es10.2)') "  (x1,y1)=",x1,y1
     write(101,'(a10,2es10.2)') "  (x2,y2)=",x2,y2
     write(101,'(a10,2es10.2)') "  (x3,y3)=",x3,y3
     write(101,'(a10,2es10.2)') "  (x4,y4)=",x4,y4
     write(101,'(a10,2es10.2)') "  (xc,yc)=",xc,yc
     stop
    endif

    if (tri_area(x3,x4,xc,y3,y4,yc)<zero) then
     write(101,*) " Centroid outside the quad element 34c: i=",i
     write(101,'(a10,2es10.2)') "  (x1,y1)=",x1,y1
     write(101,'(a10,2es10.2)') "  (x2,y2)=",x2,y2
     write(101,'(a10,2es10.2)') "  (x3,y3)=",x3,y3
     write(101,'(a10,2es10.2)') "  (x4,y4)=",x4,y4
     write(101,'(a10,2es10.2)') "  (xc,yc)=",xc,yc
     stop
    endif

    if (tri_area(x4,x1,xc,y4,y1,yc)<zero) then
     write(101,*) " Centroid outside the quad element 41c: i=",i
     write(101,'(a10,2es10.2)') "  (x1,y1)=",x1,y1
     write(101,'(a10,2es10.2)') "  (x2,y2)=",x2,y2
     write(101,'(a10,2es10.2)') "  (x3,y3)=",x3,y3
     write(101,'(a10,2es10.2)') "  (x4,y4)=",x4,y4
     write(101,'(a10,2es10.2)') "  (xc,yc)=",xc,yc
     stop
    endif

!  Distribution of element number to the 4th node of the quadrilateral
   node(v4)%nelms = node(v4)%nelms + 1
   call my_alloc_int_ptr(node(v4)%elm, node(v4)%nelms)
   node(v4)%elm(node(v4)%nelms) = i

   endif tri_or_quad

  end do elements

! Median dual volume

  do i = 1, nnodes
   node(i)%vol = zero
  end do

  elementsv : do i = 1, nelms

   v1 = elm(i)%vtx(1)
   v2 = elm(i)%vtx(2)
   v3 = elm(i)%vtx(3)

   tri_or_quadv : if (elm(i)%nvtx==3) then
!   Dual volume is exactly 1/3 of the volume of the triangle.
    node(v1)%vol = node(v1)%vol + third*elm(i)%vol
    node(v2)%vol = node(v2)%vol + third*elm(i)%vol
    node(v3)%vol = node(v3)%vol + third*elm(i)%vol

   elseif (elm(i)%nvtx==4) then
    v4 = elm(i)%vtx(4)

    x1 = node(v1)%x
    x2 = node(v2)%x
    x3 = node(v3)%x
    x4 = node(v4)%x
    xc = elm(i)%x

    y1 = node(v1)%y
    y2 = node(v2)%y
    y3 = node(v3)%y
    y4 = node(v4)%y
    yc = elm(i)%y

! - Vertex 1
     xj = node(v1)%x
     yj = node(v1)%y
    xm1 = half*(xj+x2)
    ym1 = half*(yj+y2)
    xm2 = half*(xj+x4)
    ym2 = half*(yj+y4)

!   Median volume is computed as a sum of two triangles.
    node(v1)%vol = node(v1)%vol + & 
                   tri_area(xj,xm1,xc,yj,ym1,yc) + tri_area(xj,xc,xm2,yj,yc,ym2)

! - Vertex 2
     xj = node(v2)%x
     yj = node(v2)%y
    xm1 = half*(xj+x3)
    ym1 = half*(yj+y3)
    xm2 = half*(xj+x1)
    ym2 = half*(yj+y1)

!   Median volume is computed as a sum of two triangles.
    node(v2)%vol = node(v2)%vol + &
                   tri_area(xj,xm1,xc,yj,ym1,yc) + tri_area(xj,xc,xm2,yj,yc,ym2)

! - Vertex 3
     xj = node(v3)%x
     yj = node(v3)%y
    xm1 = half*(xj+x4)
    ym1 = half*(yj+y4)
    xm2 = half*(xj+x2)
    ym2 = half*(yj+y2)

!   Median volume is computed as a sum of two triangles.
    node(v3)%vol = node(v3)%vol + &
                   tri_area(xj,xm1,xc,yj,ym1,yc) + tri_area(xj,xc,xm2,yj,yc,ym2)

! - Vertex 4
     xj = node(v4)%x
     yj = node(v4)%y
    xm1 = half*(xj+x1)
    ym1 = half*(yj+y1)
    xm2 = half*(xj+x3)
    ym2 = half*(yj+y3)

!   Median volume is computed as a sum of two triangles.
    node(v4)%vol = node(v4)%vol + &
                   tri_area(xj,xm1,xc,yj,ym1,yc) + tri_area(xj,xc,xm2,yj,yc,ym2)
 
   endif tri_or_quadv

  end do elementsv

!--------------------------------------------------------------------------------
! Loop over elements 2
!
!  Allocate elm(:)%nghbr(:) : elm(:)%nnghrs, elm(:)%nghr(:)
!  Construct element nghbr data: elm(:)%nghbr(:)
!  Order of neighbor elements [e1,e2,e3,..] are closely related to
!  the order of vertices [v1,v2,v3,..] (see below).
!
!          o------o
!          |      |                
!        v4|  e1  |v3                     v3
!    o-----o------o------o      o---------o------------o
!    |     |      |      |       .      .   .        .
!    | e2  |      |  e4  |        . e2 .     . e1  .
!    o-----o------o------o         .  .       .  .
!       v1 |     .v2              v1 o---------o v2   
!          | e3 .                     .   e3  .
!          |   .                        .    .
!          |  .                           . .
!          | .                             o
!          o
!

! Allocate the neighbor array

  do i = 1, nelms

!  3 neighbors for triangle
   if (elm(i)%nvtx==3) then

    elm(i)%nnghbrs = 3
    allocate(elm(i)%nghbr(3))

!  4 neighbors for quadrilateral
   elseif (elm(i)%nvtx==4) then

    elm(i)%nnghbrs = 4
    allocate(elm(i)%nghbr(4))

   endif

  end do

! Begin constructing the element-neighbor data

  elements2 : do i = 1, nelms

   elm_vertex : do k = 1, elm(i)%nvtx

!   Get the face of the element i:
!
!             vL      vR
!              o------o
!             /       |
!            /        |
!           o---------o
!
    if (k  < elm(i)%nvtx) vL = elm(i)%vtx(k+1)
    if (k == elm(i)%nvtx) vL = elm(i)%vtx(1)     
    vR = elm(i)%vtx(k)

!   Loop over the surrounding elements of the node vR,
!   and find the element neighbor from them.
    found = .false.
    elms_around_vR : do j = 1, node(vR)%nelms
    jelm = node(vR)%elm(j)

     edge_matching : do ii = 1, elm(jelm)%nvtx
                   v1 = elm(jelm)%vtx(ii)
      if (ii  > 1) v2 = elm(jelm)%vtx(ii-1)
      if (ii == 1) v2 = elm(jelm)%vtx(elm(jelm)%nvtx)

      if (v1==vR .and. v2==vL) then
       found = .true.
       im = ii+1
       if (im > elm(jelm)%nvtx) im = im - elm(jelm)%nvtx
       exit edge_matching
      endif
     end do edge_matching

     if (found) exit elms_around_vR

    end do elms_around_vR

     in = k + 2
     if (in > elm(i)%nvtx) in = in - elm(i)%nvtx

    if (found) then
     elm(   i)%nghbr(in) = jelm
     elm(jelm)%nghbr(im) = i
    else
     elm(   i)%nghbr(in) = 0
    endif

   end do elm_vertex

  end do elements2

!--------------------------------------------------------------------------------
! Edge-data for node-centered (edge-based) scheme.
!
! Loop over elements 3
! Construct edge data: edge(:)%n1, n2, e1, e2.
! Edge points from node n1 to node n2.
!
!      n2
!       o------------o
!     .  \         .
!    .    \   e2  .
!   .  e1  \    .
!  .        \ .         Directed area is positive: n1 -> n2
! o----------o         e1: left element
!             n1       e2: right element (e2 > e1 or e2 = 0)

! First count the number of edges.
!
! NOTE: Count edges only if the neighbor element number is
!       greater than the current element (i) to avoid double
!       count. Zero element number indicates that it is outside
!       the domain (boundary face).

  elements0 : do i = 1, nelms

   v1 = elm(i)%vtx(1)
   v2 = elm(i)%vtx(2)
   v3 = elm(i)%vtx(3)

   tri_quad0 : if (elm(i)%nvtx==3) then

    if ( elm(i)%nghbr(3) > i  .or. elm(i)%nghbr(3)==0 ) then
     nedges = nedges + 1
    endif

    if ( elm(i)%nghbr(1) > i .or. elm(i)%nghbr(1)==0 ) then
     nedges = nedges + 1
    endif

    if ( elm(i)%nghbr(2) > i .or. elm(i)%nghbr(2)==0 ) then
     nedges = nedges + 1
    endif

   elseif (elm(i)%nvtx==4) then

    v4 = elm(i)%vtx(4)

    if ( elm(i)%nghbr(3) > i .or. elm(i)%nghbr(3) ==0 ) then
     nedges = nedges + 1
    endif

    if ( elm(i)%nghbr(4) > i .or. elm(i)%nghbr(4) ==0 ) then
     nedges = nedges + 1
    endif

    if ( elm(i)%nghbr(1) > i .or. elm(i)%nghbr(1) ==0 ) then
     nedges = nedges + 1
    endif

    if ( elm(i)%nghbr(2) > i .or. elm(i)%nghbr(2) ==0 ) then
     nedges = nedges + 1
    endif

   endif tri_quad0

  end do elements0

! Allocate the edge array.
  allocate(edge(nedges))
  nedges = 0
  edge(:)%e1 = 0
  edge(:)%e2 = 0

! Construct the edge data:
!  two end nodes (n1, n2), and left and right elements (e1, e2)

  elements3 : do i = 1, nelms

   v1 = elm(i)%vtx(1)
   v2 = elm(i)%vtx(2)
   v3 = elm(i)%vtx(3)

! Triangular element
   tri_quad2 : if (elm(i)%nvtx==3) then

    if ( elm(i)%nghbr(3) > i  .or. elm(i)%nghbr(3)==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v1
     edge(nedges)%n2 = v2
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(3)
    endif

    if ( elm(i)%nghbr(1) > i .or. elm(i)%nghbr(1)==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v2
     edge(nedges)%n2 = v3
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(1)
    endif

    if ( elm(i)%nghbr(2) > i .or. elm(i)%nghbr(2)==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v3
     edge(nedges)%n2 = v1
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(2)
    endif

!  Quadrilateral element
   elseif (elm(i)%nvtx==4) then

    v4 = elm(i)%vtx(4)

    if ( elm(i)%nghbr(3) > i .or. elm(i)%nghbr(3) ==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v1
     edge(nedges)%n2 = v2
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(3)
    endif

    if ( elm(i)%nghbr(4) > i .or. elm(i)%nghbr(4) ==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v2
     edge(nedges)%n2 = v3
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(4)
    endif

    if ( elm(i)%nghbr(1) > i .or. elm(i)%nghbr(1) ==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v3
     edge(nedges)%n2 = v4
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(1)
    endif

    if ( elm(i)%nghbr(2) > i .or. elm(i)%nghbr(2) ==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v4
     edge(nedges)%n2 = v1
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(2)
    endif

   endif tri_quad2

  end do elements3

! Loop over edges
! Construct edge vector and directed area vector.
!
! Edge vector is a simple vector pointing froom n1 to n2.
! For each edge, add the directed area vector (dav) from
! the left and right elements.
!
!              n2
!   o-----------o-----------o
!   |     dav   |  dav      |
!   |       ^   |   ^       |
!   |       |   |   |       |
!   |   c - - - m - - -c    |
!   |           |           |
!   |           |           |    m: edge midpoint
!   |           |           |    c: element centroid
!   o-----------o-----------o
!                n1
!
  edges : do i = 1, nedges

   n1 = edge(i)%n1
   n2 = edge(i)%n2
   e1 = edge(i)%e1
   e2 = edge(i)%e2
   xm = half*( node(n1)%x + node(n2)%x )
   ym = half*( node(n1)%y + node(n2)%y )

   edge(i)%dav = zero

! Contribution from the left element
  if (e1 > 0) then
   xc = elm(e1)%x
   yc = elm(e1)%y
   edge(i)%dav(1) = -(ym-yc)
   edge(i)%dav(2) =   xm-xc
  endif

! Contribution from the right element
  if (e2 > 0) then
   xc = elm(e2)%x
   yc = elm(e2)%y
   edge(i)%dav(1) = edge(i)%dav(1) -(yc-ym)
   edge(i)%dav(2) = edge(i)%dav(2) + xc-xm
  endif

  if (e1 < 0 .and. e2 < 0) then
   write(101,*) "!!!!! e1 and e2 are both negative... No way..."
  endif

! Magnitude and unit vector
   edge(i)%da  = sqrt( edge(i)%dav(1)**2 + edge(i)%dav(2)**2 )
   edge(i)%dav = edge(i)%dav / edge(i)%da

! Edge vector

  edge(i)%ev(1) = node(n2)%x - node(n1)%x
  edge(i)%ev(2) = node(n2)%y - node(n1)%y
  edge(i)%e     = sqrt( edge(i)%ev(1)**2 + edge(i)%ev(2)**2 )
  edge(i)%ev    = edge(i)%ev / edge(i)%e

  end do edges

!--------------------------------------------------------------------------------
! Construct node neighbor data:
!  pointers to the neighbor nodes(o)
!
!        o     o
!         \   / 
!          \ /
!     o-----*-----o
!          /|
!         / |
!        /  o        *: node in interest
!       o            o: neighbors (edge-connected nghbrs)
!

  do i = 1, nnodes
   node(i)%nnghbrs = 0
  end do

! Loop over edges and distribute the node numbers:

  edges4 : do i = 1, nedges

   n1 = edge(i)%n1
   n2 = edge(i)%n2

! (1) Add node1 to the neighbor list of n2
   node(n1)%nnghbrs = node(n1)%nnghbrs + 1
   call my_alloc_int_ptr(node(n1)%nghbr, node(n1)%nnghbrs)
   node(n1)%nghbr(node(n1)%nnghbrs) = n2

! (2) Add node2 to the neighbor list of n1
   node(n2)%nnghbrs = node(n2)%nnghbrs + 1
   call my_alloc_int_ptr(node(n2)%nghbr, node(n2)%nnghbrs)
   node(n2)%nghbr(node(n2)%nnghbrs) = n1

  end do edges4

!--------------------------------------------------------------------------------
! Boundary normal at nodes constructed by accumulating the contribution
! from each boundary face normal. This vector will be used to enforce
! the tangency condition, for example.
!
!
!        Interior domain      /
!                            o
!                  .        /
!                  .       /
! --o-------o-------------o
!           j   |  .  |   j+1
!               v  .  v
!
!        Left half added to the node j, and
!       right half added to the node j+1.
!

! Allocate and initialize the normal vector arrays
  do i = 1, nbound

   allocate(bound(i)%bnx(bound(i)%nbnodes))
   allocate(bound(i)%bny(bound(i)%nbnodes))
   allocate(bound(i)%bn( bound(i)%nbnodes))

   do j = 1, bound(i)%nbnodes
    bound(i)%bnx(j) = zero
    bound(i)%bny(j) = zero
    bound(i)%bn( j) = zero
   end do

  end do

! Normal vector at boundary nodes
! Note: Below it describes normals of linear approximation.
!       We will overwrite it by a quadratic approximation.
!
! Linear approximation:
!
! Step 1. Compute the outward normals
  do i = 1, nbound
   do j = 1, bound(i)%nbnodes-1

    x1 = node(bound(i)%bnode(j  ))%x
    y1 = node(bound(i)%bnode(j  ))%y

    x2 = node(bound(i)%bnode(j+1))%x
    y2 = node(bound(i)%bnode(j+1))%y

!   Normal vector pointing into the domain at this point.
    bound(i)%bnx(j) = bound(i)%bnx(j) + half*( -(y2-y1) )
    bound(i)%bny(j) = bound(i)%bny(j) + half*(   x2-x1  )

    bound(i)%bnx(j+1) = bound(i)%bnx(j+1) + half*( -(y2-y1) )
    bound(i)%bny(j+1) = bound(i)%bny(j+1) + half*(   x2-x1  )

   end do
  end do

! Step 2. Compute the magnitude and turn (bnx,bny) into a unit vector
  do i = 1, nbound
   do j = 1, bound(i)%nbnodes

    bound(i)%bn(j)  = sqrt( bound(i)%bnx(j)**2 + bound(i)%bny(j)**2 )
!   Minus sign to make it pont out towards the outside of the domain.
    bound(i)%bnx(j) =  - bound(i)%bnx(j) / bound(i)%bn(j)
    bound(i)%bny(j) =  - bound(i)%bny(j) / bound(i)%bn(j)

   end do
  end do

! Now, ignore the linear approximation, and let us construct
! more accurate surfae normal vectors and replace the linear ones.
! So, we will overwrite the unit normal vectors: bnx, bny.
! Note: We keep the magnitude of the normal vector.
!
! Quadratic approximation:
! See http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2015v281pp518-555_preprint.pdf
! for details on the quadratic approximation for computing more accurate normals.
! 
!  boundary_type0 : do i = 1, nbound
!   boundary_nodes0 : do j = 1, bound(i)%nbnodes
!
!     if (j==1) then
!      v1 = bound(i)%bnode(j  )
!      v2 = bound(i)%bnode(j+1)
!      v3 = bound(i)%bnode(j+2)
!     elseif (j==bound(i)%nbnodes) then
!      v1 = bound(i)%bnode(j-2)
!      v2 = bound(i)%bnode(j-1)
!      v3 = bound(i)%bnode(j  )
!     else
!      v1 = bound(i)%bnode(j-1)
!      v2 = bound(i)%bnode(j)
!      v3 = bound(i)%bnode(j+1)
!     endif
!
!     x1 = node(v1)%x
!     x2 = node(v2)%x
!     x3 = node(v3)%x
!
!     y1 = node(v1)%y
!     y2 = node(v2)%y
!     y3 = node(v3)%y
!
!!----------------------------------------------------------------------
!!   Fit a quadratic over 3 nodes
!
!!    Skip the last one if the boundary segment is a closed boundary 
!!    in which case the last node is the same as the first one.
!     if (j==bound(i)%nbnodes .and. bound(i)%bnode(j)==bound(i)%bnode(1) ) then
!      bound(i)%bn(j)  = bound(i)%bn(1)
!      bound(i)%bnx(j) = bound(i)%bnx(1)
!      bound(i)%bny(j) = bound(i)%bny(1)
!      cycle
!     endif
!
!     dsL = sqrt( (x2-x1)**2 + (y2-y1)**2 )
!     dsR = sqrt( (x3-x2)**2 + (y3-y2)**2 )
!      dx = dsR*x1/(dsL*(-dsL-dsR))-x2/dsR+x2/dsL+dsL*x3/((dsR+dsL)*dsR)
!      dy = dsR*y1/(dsL*(-dsL-dsR))-y2/dsR+y2/dsL+dsL*y3/((dsR+dsL)*dsR)
!
!     ds  = sqrt( dx**2 + dy**2 )
!     bound(i)%bnx(j) = -( -dy / ds )
!     bound(i)%bny(j) = -(  dx / ds )
!
!   end do boundary_nodes0
!  end do boundary_type0

!--------------------------------------------------------------------------------
! Construct neighbor index over edges
!
!  Example:
!
!        o     o
!         \   / 
!          \j/       k-th neighbor
!     o-----*----------o
!          /|  edge i
!         / |
!        /  o        Note: k-th neighbor is given by "node(j)%nghbr(k)"
!       o
!
!  Consider the edge i
!
!   node j        k-th neighbor
!       *----------o
!      n1  edge i  n2
!
!   We store "k" in the edge data structure as
!
!    edge(i)%kth_nghbr_of_1: n2 is the "edge(i)%kth_nghbr_of_1"-th neighbor of n1
!    edge(i)%kth_nghbr_of_2: n1 is the "edge(i)%kth_nghbr_of_3"-th neighbor of n2
!
!   That is,  we have
!
!    n2 = node(n1)%nghbr(edge(i)%kth_nghbr_of_1)
!    n1 = node(n2)%nghbr(edge(i)%kth_nghbr_of_2)
!
!   We make use of this data structure to access off-diagonal entries in Jacobian matrix.
!

! Loop over edges

  edges5 : do i = 1, nedges

   n1 = edge(i)%n1
   n2 = edge(i)%n2

   do k = 1, node(n2)%nnghbrs

    if ( n1 == node(n2)%nghbr(k) ) then
     edge(i)%kth_nghbr_of_2 = k
    endif

   end do

   do k = 1, node(n1)%nnghbrs

    if ( n2 == node(n1)%nghbr(k) ) then
     edge(i)%kth_nghbr_of_1 = k
    endif

   end do

  end do edges5

! Boundary mark: It should be an array actually because some nodes are associated with
!                more than one boundaries.
  do i = 1, nnodes
   node(i)%bmark   = 0
   node(i)%nbmarks = 0
  end do

  do i = 1, nbound
   do j = 1, bound(i)%nbnodes
    node( bound(i)%bnode(j) )%bmark   = i
    node( bound(i)%bnode(j) )%nbmarks = node( bound(i)%bnode(j) )%nbmarks + 1
   end do
  end do

!--------------------------------------------------------------------------------
! Boundary face data
!
!      |     Domain      |
!      |                 |
!      o--o--o--o--o--o--o  <- Boundary segment
!   j= 1  2  3  4  5  6  7
!
!   In the above case, nbnodes = 7, nbfaces = 6
!

  do i = 1, nbound
   bound(i)%nbfaces = bound(i)%nbnodes-1
   allocate(bound(i)%bfnx(    bound(i)%nbfaces   ))
   allocate(bound(i)%bfny(    bound(i)%nbfaces   ))
   allocate(bound(i)%bfn(     bound(i)%nbfaces   ))
   allocate(bound(i)%belm(    bound(i)%nbfaces   ))
   allocate(bound(i)%kth_nghbr_of_1(    bound(i)%nbfaces   ))
   allocate(bound(i)%kth_nghbr_of_2(    bound(i)%nbfaces   ))
  end do

! Boundary face vector: outward normal
  do i = 1, nbound
   do j = 1, bound(i)%nbfaces

    x1 = node(bound(i)%bnode(j  ))%x
    y1 = node(bound(i)%bnode(j  ))%y
    x2 = node(bound(i)%bnode(j+1))%x
    y2 = node(bound(i)%bnode(j+1))%y

    bound(i)%bfn(j)  =  sqrt( (x1-x2)**2 + (y1-y2)**2 )
    bound(i)%bfnx(j) = -(y1-y2) / bound(i)%bfn(j)
    bound(i)%bfny(j) =  (x1-x2) / bound(i)%bfn(j)

   end do
  end do

! Boundary normal vector at nodes: outward normal
  do i = 1, nbound
   do j = 1, bound(i)%nbfaces

    x1 = node(bound(i)%bnode(j  ))%x
    y1 = node(bound(i)%bnode(j  ))%y
    x2 = node(bound(i)%bnode(j+1))%x
    y2 = node(bound(i)%bnode(j+1))%y

    bound(i)%bfn(j)  =  sqrt( (x1-x2)**2 + (y1-y2)**2 )
    bound(i)%bfnx(j) = -(y1-y2) / bound(i)%bfn(j)
    bound(i)%bfny(j) =  (x1-x2) / bound(i)%bfn(j)

   end do
  end do

! Neighbor index over boundary edges (faces)

  do i = 1, nbound
   do j = 1, bound(i)%nbfaces

    n1 = bound(i)%bnode(j  )  !Left node
    n2 = bound(i)%bnode(j+1)  !Right node

    do k = 1, node(n2)%nnghbrs
     if ( n1 == node(n2)%nghbr(k) ) then
      bound(i)%kth_nghbr_of_2(j) = k
     endif
    end do

    do k = 1, node(n1)%nnghbrs
     if ( n2 == node(n1)%nghbr(k) ) then
      bound(i)%kth_nghbr_of_1(j) = k
     endif
    end do

   end do
  end do

! Find element adjacent to the face: belm
!
!  NOTE: This is useful to figure out what element
!        each boundary face belongs to. Boundary flux needs
!        special weighting depending on the element.
!
!      |_________|_________|________|
!      |         |         |        | 
!      |         |         |        | 
!      |_________|_________|________|
!      |         |         |        |     <- Grid (e.g., quads)
!      |         | elmb(j) |        |
!   ---o---------o---------o--------o---  <- Boundary segment
!                 j-th face
!
! elmb(j) is the element number of the element having the j-th boundary face.
!

  do i = 1, nbound
   do j = 1, bound(i)%nbfaces

!   bface is defined by the nodes v1 and v2.
    v1 = bound(i)%bnode(j  )
    v2 = bound(i)%bnode(j+1)

    found = .false.
!   Find the element having the bface from the elements
!   around the node v1.
    do k = 1, node(v1)%nelms
     ielm = node(v1)%elm(k)
     do ii = 1, elm(ielm)%nvtx
      in = ii
      im = ii+1
      if (im > elm(ielm)%nvtx) im = im - elm(ielm)%nvtx !return to 1
      vt1 = elm(ielm)%vtx(in)
      vt2 = elm(ielm)%vtx(im)
       if (vt1 == v1 .and. vt2 == v2) then
        found = .true.
        exit
       endif
     end do
     if (found) exit
    end do

    if (found) then
     bound(i)%belm(j) = ielm
    else
     write(101,*) " Boundary-adjacent element not found. Error..."
     stop
    endif

   end do
  end do

!--------------------------------------------------------------------------------
! Construct least-squares matrix for node-centered schemes.
!
!        o     o
!         \   / 
!          \ /
!     o-----*-----o
!          /|
!         / |
!        /  o        *: node in interest
!       o            o: neighbors (edge-connected nghbrs)
!

! Check the number of neighbor nodes (must have at least 2 neighbors)
  write(101,*) " --- Node neighbor data:"

  ave_nghbr = node(1)%nnghbrs
  min_nghbr = node(1)%nnghbrs
  max_nghbr = node(1)%nnghbrs
       imin = 1
       imax = 1
   if (node(1)%nnghbrs==2) then
    write(101,*) "--- 2 neighbors for the node = ", 1
   endif

  do i = 2, nnodes
   ave_nghbr = ave_nghbr + node(i)%nnghbrs
   if (node(i)%nnghbrs < min_nghbr) imin = i
   if (node(i)%nnghbrs > max_nghbr) imax = i
   min_nghbr = min(min_nghbr, node(i)%nnghbrs)
   max_nghbr = max(max_nghbr, node(i)%nnghbrs)
   if (node(i)%nnghbrs==2) then
    write(101,*) "--- 2 neighbors for the node = ", i
   endif
  end do

  write(101,*) "      ave_nghbr = ", ave_nghbr/nnodes
  write(101,*) "      min_nghbr = ", min_nghbr, " at node ", imin
  write(101,*) "      max_nghbr = ", max_nghbr, " at node ", imax
  write(101,*)

!--------------------------------------------------------------------------------
! Cell centered scheme data
!--------------------------------------------------------------------------------

  write(101,*) "Generating CC scheme data......"

  do i = 1, nelms   
   allocate(elm(i)%edge( elm(i)%nnghbrs ) )
  end do

  edges3 : do i = 1, nedges

   e1 = edge(i)%e1
   e2 = edge(i)%e2

! Left element
  if (e1 > 0) then
   do k = 1, elm(e1)%nnghbrs
    if (elm(e1)%nghbr(k)==e2) elm(e1)%edge(k) = i
   end do
  endif

! Right element
  if (e2 > 0) then
   do k = 1, elm(e2)%nnghbrs
    if (elm(e2)%nghbr(k)==e1) elm(e2)%edge(k) = i
   end do
  endif

  end do edges3

! Face-data for cell-centered (edge-based) scheme.
!
! Loop over elements 4
! Construct face data:
! face is an edge across elements pointing
! element e1 to element e2 (e2 > e1):
!
!       e2
!        \    
!         \ face: e1 -> e2 
!          \
!  n1 o--------------o n2 <-- face
!            \
!             \          n1, n2: end nodes of the face
!              \         e1: element 1
!              e1        e2: element 2  (e2 > e1)
!
! Note: Face data is dual to the edge data.
!       It can be trivially constructed from the edge data, but
!       here the face data is constructed by using the element
!       neighbor data just for an educational purpose.

  nfaces = 0
  elements4 : do i = 1, nelms
   do k = 1, elm(i)%nnghbrs
   jelm = elm(i)%nghbr(k)
    if (jelm > i) then
     nfaces = nfaces + 1
    endif
   end do
  end do elements4

  allocate(face(nfaces))

  nfaces = 0

  elements5 : do i = 1, nelms
   do k = 1, elm(i)%nnghbrs
   jelm = elm(i)%nghbr(k)

    if (jelm > i) then

     nfaces = nfaces + 1

     face(nfaces)%e1 = i
     face(nfaces)%e2 = jelm

     iedge = elm(i)%edge(k)
     v1 = edge(iedge)%n1
     v2 = edge(iedge)%n2

     if (edge(iedge)%e1 == jelm) then
      face(nfaces)%n1 = v1
      face(nfaces)%n2 = v2
     else
      face(nfaces)%n1 = v2
      face(nfaces)%n2 = v1
     endif
   
    elseif (jelm == 0) then
!    Skip boundary faces.
    endif

   end do
  end do elements5

! Loop over faces
! Construct directed area vector.

  faces : do i = 1, nfaces

   n1 = face(i)%n1
   n2 = face(i)%n2
   e1 = face(i)%e1
   e2 = face(i)%e2

! Face vector
  face(i)%dav(1) = -( node(n2)%y - node(n1)%y )
  face(i)%dav(2) =    node(n2)%x - node(n1)%x
  face(i)%da     = sqrt( face(i)%dav(1)**2 + face(i)%dav(2)**2 )
  face(i)%dav    = face(i)%dav / face(i)%da

  end do faces

! Construct vertex-neighbor data for cell-centered scheme.
!
! For each element, i, collect all elements sharing the nodes
! of the element, i, including face-neighors.
!
!      ___________
!     |     |     |
!     |  o  |  o  |
!     |_____|_____|
!    /\    / \    \
!   / o\ o/ i \  o \
!  /____\/_____\____\
!  \    /      /\    \
!   \o /  o   / o\ o  \
!    \/______/____\____\
!
!          i: Element of interest
!          o: Vertex neighbors (k = 1,2,...,9)

  write(101,*) " --- Vertex-neighbor data:"

  do i = 1, nelms
   elm(i)%nvnghbrs = 1
   call my_alloc_int_ptr(elm(i)%vnghbr, 1)
  end do

  ave_nghbr = 0
  min_nghbr = 10000
  max_nghbr =-10000
       imin = 1
       imax = 1

! Initialization
  elements6 : do i = 1, nelms
   elm(i)%nvnghbrs = 0
  end do elements6

! Collect vertex-neighbors
  elements7 : do i = 1, nelms

! (1)Add face-neighbors
   do k = 1, elm(i)%nnghbrs
    if ( elm(i)%nghbr(k) > 0 ) then
     elm(i)%nvnghbrs = elm(i)%nvnghbrs + 1
     call my_alloc_int_ptr(elm(i)%vnghbr, elm(i)%nvnghbrs)
     elm(i)%vnghbr(elm(i)%nvnghbrs) = elm(i)%nghbr(k)
    endif
   end do

! (2)Add vertex-neighbors
   do k = 1, elm(i)%nvtx
    v1 = elm(i)%vtx(k)

    velms : do j = 1, node(v1)%nelms
     e1 = node(v1)%elm(j)
     if (e1 == i) cycle velms

!    Check if the element is already added.
       found = .false.
     do ii = 1, elm(i)%nvnghbrs
      if ( e1 == elm(i)%vnghbr(ii) ) then
       found = .true.
       exit
      endif
     end do

!    Add the element, e1, if not added yet.
     if (.not.found) then
      elm(i)%nvnghbrs = elm(i)%nvnghbrs + 1
      call my_alloc_int_ptr(elm(i)%vnghbr, elm(i)%nvnghbrs)
      elm(i)%vnghbr(elm(i)%nvnghbrs) = e1
     endif
    end do velms

   end do

   ave_nghbr = ave_nghbr + elm(i)%nvnghbrs
   if (elm(i)%nvnghbrs < min_nghbr) imin = i
   if (elm(i)%nvnghbrs > max_nghbr) imax = i
   min_nghbr = min(min_nghbr, elm(i)%nvnghbrs)
   max_nghbr = max(max_nghbr, elm(i)%nvnghbrs)
   if (elm(i)%nvnghbrs < 3) then
    write(101,*) "--- Not enough neighbors: elm = ", i, &
               "elm(i)%nvnghbrs=",elm(i)%nvnghbrs
   endif

  end do elements7

  write(101,*) "      ave_nghbr = ", ave_nghbr/nelms
  write(101,*) "      min_nghbr = ", min_nghbr, " elm = ", imin
  write(101,*) "      max_nghbr = ", max_nghbr, " elm = ", imax
  write(101,*)


  do i = 1, nelms
   elm(i)%bmark = 0
  end do

  bc_loop : do i = 1, nbound
   if (trim(bound(i)%bc_type) == "dirichlet") then
    do j = 1, bound(i)%nbfaces
     elm( bound(i)%belm(j) )%bmark = 1
    end do
   endif
  end do bc_loop

!--------------------------------------------------------------------------------

 return

 end subroutine construct_grid_data

!********************************************************************************

!********************************************************************************
!* Compute the area of the triangle defined by the nodes, 1, 2, 3.
!*
!*              3 (x3,y3)
!*              o 
!*             / \ 
!*            /   \
!* (x1,y1) 1 o-----o 2 (x2,y2)
!*
!* Nodes must be ordered counterclockwise (otherwise it gives negative area)
!*
!********************************************************************************
 function tri_area(x1,x2,x3,y1,y2,y3) result(area)
 use edu2d_constants, only : p2, half
 implicit none
 real(p2), intent(in) :: x1,x2,x3,y1,y2,y3
 real(p2) :: area

  area = half*( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) )

 end function tri_area


!********************************************************************************
!* Check the grid data.
!*
!* 1. Directed area must sum up to zero around every node.
!* 2. Directed area must sum up to zero over the entire grid.
!* 3. Global sum of the boundary normal vectors must vanish.
!* 4. Global sum of the boundary face normal vectors must vanish.
!* 5. Check element volumes which must be positive.
!* 6. Check dual volumes which must be positive.
!* 7. Global sum of the dual volumes must be equal to the sum of element volumes.
!*
!* Add more tests you can think of.
!*
!********************************************************************************
 subroutine check_grid_data

 use edu2d_my_main_data  , only : nnodes, node,  nelms,   elm, nedges, edge, &
                            nbound, bound
 use edu2d_constants     , only : p2, zero, half

 implicit none
!Local variables
 integer  :: i, j, n1, n2, ierr, k
 real(p2), dimension(nnodes,2) :: sum_dav_i
 real(p2), dimension(2) :: sum_dav, sum_bn
 real(p2), dimension(2) :: sum_bfn
 real(p2)               :: sum_volc, sum_vol
 real(p2)               :: mag_dav, mag_bn
 real(p2)               :: vol_min, vol_max, vol_ave

  write(101,*) "Checking grid data...."

  mag_dav = zero
  mag_bn  = zero

!--------------------------------------------------------------------------------
! Directed area sum check
!--------------------------------------------------------------------------------

! Compute the sum of the directed area for each node.

   sum_dav_i = zero
  do i = 1, nedges
   n1 = edge(i)%n1
   n2 = edge(i)%n2
   sum_dav_i(n1,:) = sum_dav_i(n1,:) + edge(i)%dav(:)*edge(i)%da
   sum_dav_i(n2,:) = sum_dav_i(n2,:) - edge(i)%dav(:)*edge(i)%da
   mag_dav = mag_dav + edge(i)%da
  end do
   mag_dav = mag_dav/real(nedges,p2)

! Add contribution from boundary edges.
  do i = 1, nbound
   do j = 1, bound(i)%nbfaces

     n1 = bound(i)%bnode(j)
     n2 = bound(i)%bnode(j+1)

     sum_dav_i(n1,1) = sum_dav_i(n1,1) + half*bound(i)%bfnx(j)*bound(i)%bfn(j)
     sum_dav_i(n1,2) = sum_dav_i(n1,2) + half*bound(i)%bfny(j)*bound(i)%bfn(j)

     sum_dav_i(n2,1) = sum_dav_i(n2,1) + half*bound(i)%bfnx(j)*bound(i)%bfn(j)
     sum_dav_i(n2,2) = sum_dav_i(n2,2) + half*bound(i)%bfny(j)*bound(i)%bfn(j)

   end do
  end do

! Compute also the sum of the boundary normal vector (at nodes).

  sum_bn = 0
  do i = 1, nbound
   do j = 1, bound(i)%nbnodes
     k = bound(i)%bnode(j)
     if (j > 1 .and. k==bound(i)%bnode(1)) cycle !Skip if the last node is equal to the first node).
    sum_bn(1)      = sum_bn(1)      + bound(i)%bnx(j)*bound(i)%bn(j)
    sum_bn(2)      = sum_bn(2)      + bound(i)%bny(j)*bound(i)%bn(j)
    mag_bn = mag_bn + abs(bound(i)%bn(j))
   end do
    mag_bn = mag_bn/real(bound(i)%nbnodes,p2)
  end do

! Global sum of boundary normal vectors must vanish.

!  if (sum_bn(1) > 1.0e-12_p2*mag_bn .and. sum_bn(2) > 1.0e-12_p2*mag_bn) then
!   write(*,*) "--- Global sum of the boundary normal vector:"
!   write(*,'(a19,es10.3)') "    sum of bn_x = ", sum_bn(1)
!   write(*,'(a19,es10.3)') "    sum of bn_y = ", sum_bn(2)
!   write(*,*) "Error: boundary normal vectors do not sum to zero..."
!   stop
!  endif

! Sum of the directed area vectors must vanish at every node.

  do i = 1, nnodes
   if (abs(sum_dav_i(i,1))>1.0e-12_p2*mag_dav .or. abs(sum_dav_i(i,2))>1.0e-12_p2*mag_dav) then
   write(101,'(a11,i5,a7,2es10.3,a9,2es10.3)') &
    " --- node=", i, " (x,y)=", node(i)%x, node(i)%y, " sum_dav=",sum_dav_i(i,:)
   endif
  end do

   write(101,*) "--- Max sum of directed area vector around a node:"
   write(101,*) "  max(sum_dav_i_x) = ", maxval(sum_dav_i(:,1))
   write(101,*) "  max(sum_dav_i_y) = ", maxval(sum_dav_i(:,2))

  if (maxval(abs(sum_dav_i(:,1)))>1.0e-12_p2*mag_dav .or. &
      maxval(abs(sum_dav_i(:,2)))>1.0e-12_p2*mag_dav) then
   write(101,*) "--- Max sum of directed area vector around a node:"
   write(101,*) "  max(sum_dav_i_x) = ", maxval(sum_dav_i(:,1))
   write(101,*) "  max(sum_dav_i_y) = ", maxval(sum_dav_i(:,2))
   write(101,*) "Error: directed area vectors do not sum to zero..."
   stop
  endif

! Of course, the global sum of the directed area vector sum must vanish.
   sum_dav = zero
  do i = 1, nnodes
   sum_dav = sum_dav + sum_dav_i(i,:)
  end do

   write(101,*) "--- Global sum of the directed area vector:"
   write(101,'(a19,es10.3)') "    sum of dav_x = ", sum_dav(1)
   write(101,'(a19,es10.3)') "    sum of dav_y = ", sum_dav(2)

  if (sum_dav(1) > 1.0e-12_p2*mag_dav .and. sum_dav(2) > 1.0e-12_p2*mag_dav) then
   write(101,*) "Error: directed area vectors do not sum globally to zero..."
   write(101,*) "--- Global sum of the directed area vector:"
   write(101,'(a19,es10.3)') "    sum of dav_x = ", sum_dav(1)
   write(101,'(a19,es10.3)') "    sum of dav_y = ", sum_dav(2)
   stop
  endif

!--------------------------------------------------------------------------------
! Global sum check for boundary face vector
!--------------------------------------------------------------------------------
  sum_bfn = 0
  do i = 1, nbound
   do j = 1, bound(i)%nbfaces
     sum_bfn(1) =  sum_bfn(1) + bound(i)%bfnx(j)*bound(i)%bfn(j)
     sum_bfn(2) =  sum_bfn(2) + bound(i)%bfny(j)*bound(i)%bfn(j)
   end do
  end do

   write(101,*) "--- Global sum of the boundary face vector:"
   write(101,'(a19,es10.3)') "    sum of bfn_x = ", sum_bfn(1)
   write(101,'(a19,es10.3)') "    sum of bfn_y = ", sum_bfn(2)

  if (sum_bfn(1) > 1.0e-12_p2*mag_bn .and. sum_bfn(2) > 1.0e-12_p2*mag_bn) then
   write(101,*) "Error: boundary face normals do not sum globally to zero..."
   write(101,*) "--- Global sum of the boundary face normal vector:"
   write(101,'(a19,es10.3)') "    sum of bfn_x = ", sum_bfn(1)
   write(101,'(a19,es10.3)') "    sum of bfn_y = ", sum_bfn(2)
   stop
  endif

!--------------------------------------------------------------------------------
! Volume check
!--------------------------------------------------------------------------------
! (1)Check the element volume: make sure there are no zero or negative volumes

   vol_min =  1.0e+15
   vol_max = -1.0
   vol_ave =  zero

       ierr = 0
   sum_volc = zero
  do i = 1, nelms

      vol_min = min(vol_min,elm(i)%vol)
      vol_max = max(vol_max,elm(i)%vol)
      vol_ave = vol_ave + elm(i)%vol

   sum_volc = sum_volc + elm(i)%vol

   if (elm(i)%vol < zero) then
     write(101,*) "Negative volc=",elm(i)%vol, " elm=",i, " stop..."
     ierr = ierr + 1
   endif

   if (abs(elm(i)%vol) < 1.0e-14_p2) then
     write(101,*) "Vanishing volc=",elm(i)%vol, " elm=",i, " stop..."
     ierr = ierr + 1
   endif

  end do

   vol_ave = vol_ave / real(nelms)

   write(101,*)
   write(101,'(a30,es25.15)') "    minimum element volume = ", vol_min
   write(101,'(a30,es25.15)') "    maximum element volume = ", vol_max
   write(101,'(a30,es25.15)') "    average element volume = ", vol_ave
   write(101,*)

!--------------------------------------------------------------------------------
! (2)Check the dual volume (volume around a node)

   vol_min =  1.0e+15
   vol_max = -1.0
   vol_ave =  zero

      ierr = 0
   sum_vol = zero
  do i = 1, nnodes

      vol_min = min(vol_min,node(i)%vol)
      vol_max = max(vol_max,node(i)%vol)
      vol_ave = vol_ave + node(i)%vol

   sum_vol = sum_vol + node(i)%vol

   if (node(i)%vol < zero) then
     write(101,*) "Negative vol=",node(i)%vol, " node=",i, " stop..."
     ierr = ierr + 1
   endif

   if (abs(node(i)%vol) < 1.0e-14_p2) then
     write(101,*) "Vanishing vol=",node(i)%vol, " node=",i, " stop..."
     ierr = ierr + 1
   endif

  end do

   vol_ave = vol_ave / real(nnodes)

   write(101,*)
   write(101,'(a30,es25.15)') "    minimum dual volume = ", vol_min
   write(101,'(a30,es25.15)') "    maximum dual volume = ", vol_max
   write(101,'(a30,es25.15)') "    average dual volume = ", vol_ave
   write(101,*)


  if (ierr > 0) stop

  if (abs(sum_vol-sum_volc) > 1.0e-08_p2*sum_vol) then
   write(101,*) "--- Global sum of volume: must be the same"
   write(101,'(a19,es10.3)') "    sum of volc = ", sum_volc
   write(101,'(a19,es10.3)') "    sum of vol  = ", sum_vol
   write(101,'(a22,es10.3)') " sum_vol-sum_volc  = ", sum_vol-sum_volc
   write(101,*) "Error: sum of dual volumes and cell volumes do not match..."
   stop
  endif

  call check_skewness_nc
  call compute_ar

  write(101,*)
  write(101,*) "Grid data look good!"

 end subroutine check_grid_data

!*******************************************************************************
!* Skewness computation for edges.
!*******************************************************************************
 subroutine check_skewness_nc
 use edu2d_my_main_data  , only : nedges, edge
 use edu2d_constants     , only : p2, zero

 implicit none
 integer :: i
 real(p2) :: e_dot_n, e_dot_n_min, e_dot_n_max, alpha

     e_dot_n = zero
 e_dot_n_min = 100000.0_p2
 e_dot_n_max =-100000.0_p2

  edges : do i = 1, nedges

   alpha = edge(i)%ev(1)*edge(i)%dav(1) + edge(i)%ev(2)*edge(i)%dav(2)
   e_dot_n     = e_dot_n + abs(alpha)
   e_dot_n_min = min(e_dot_n_min, abs(alpha))
   e_dot_n_max = max(e_dot_n_max, abs(alpha))

  end do edges

  e_dot_n = e_dot_n / real(nedges,p2)

 write(101,*)
 write(101,*) " ------ Skewness check (NC control volume) ----------"
 write(101,*) "   L1(e_dot_n) = ", e_dot_n
 write(101,*) "  Min(e_dot_n) = ", e_dot_n_min
 write(101,*) "  Max(e_dot_n) = ", e_dot_n_max
 write(101,*) " ----------------------------------------------------"

 end subroutine check_skewness_nc


!*******************************************************************************
!* Control volume aspect ratio
!*******************************************************************************
 subroutine compute_ar
 use edu2d_my_main_data  , only : node, nnodes, elm, nelms
 use edu2d_constants     , only : p2, zero, one, half, two

 implicit none
 integer :: i, n1, n2
 real(p2) :: ar, ar_min, ar_max, nnodes_eff

 integer  :: k
 real(p2) :: side_max, side(4), side_mid, side_min, height

! Initialization

  node1 : do i = 1, nnodes
   node(i)%ar     = zero
  end do node1

! Compute element aspect-ratio: longest_side^2 / vol

  elm1: do i = 1, nelms

   side_max = -one
   
   do k = 1, elm(i)%nvtx

     n1 = elm(i)%vtx(k)
    if (k == elm(i)%nvtx) then
     n2 = elm(i)%vtx(1)
    else
     n2 = elm(i)%vtx(k+1)
    endif

     side(k) = sqrt( (node(n2)%x-node(n1)%x)**2 + (node(n2)%y-node(n1)%y)**2 )
    side_max =  max(side_max, side(k))

   end do

   if (elm(i)%nvtx == 3) then

 ! AR for triangle:  Ratio of a half of a square with side_max to volume
    elm(i)%ar = (half*side_max**2) / elm(i)%vol

    if     (side(1) >= side(2) .and. side(1) >= side(3)) then

       side_max = side(1)
      if (side(2) >= side(3)) then
       side_mid = side(2); side_min = side(3)
      else
       side_mid = side(3); side_min = side(2)
      endif

    elseif (side(2) >= side(1) .and. side(2) >= side(3)) then

       side_max = side(2)
      if (side(1) >= side(3)) then
       side_mid = side(1); side_min = side(3)
      else
       side_mid = side(3); side_min = side(1)
      endif

    else

       side_max = side(3)
      if (side(1) >= side(2)) then
       side_mid = side(1); side_min = side(2)
      else
       side_mid = side(2); side_min = side(1)
      endif

    endif

       height = two*elm(i)%vol/side_mid
    elm(i)%ar = side_mid/height

   else

  ! AR for quad: Ratio of a square with side_max to volume
    elm(i)%ar = side_max**2 / elm(i)%vol

   endif

  end do elm1

! Compute the aspect ratio:
  node2 : do i = 1, nnodes

    node(i)%ar = zero
   do k = 1, node(i)%nelms
    node(i)%ar = node(i)%ar + elm(node(i)%elm(k))%ar
   end do

    node(i)%ar = node(i)%ar / real(node(i)%nelms, p2)

  end do node2

! Compute the min/max and L1 of AR

  nnodes_eff= zero
         ar = zero
     ar_min = 100000.0_p2
     ar_max =-100000.0_p2

  node3: do i = 1, nnodes
   if (node(i)%bmark /= 0) cycle node3
   ar     = ar + abs(node(i)%ar)
   ar_min = min(ar_min, abs(node(i)%ar))
   ar_max = max(ar_max, abs(node(i)%ar))
   nnodes_eff = nnodes_eff + one
  end do node3

  ar = ar / nnodes_eff

 write(101,*)
 write(101,*) " ------ Aspect ratio check (NC control volume) ----------"
 write(101,*) " Interior nodes only"
 write(101,*) "   L1(AR) = ", ar
 write(101,*) "  Min(AR) = ", ar_min
 write(101,*) "  Max(AR) = ", ar_max

  nnodes_eff= zero
         ar = zero
     ar_min = 100000.0_p2
     ar_max =-100000.0_p2

  node4: do i = 1, nnodes
   if (node(i)%bmark == 0) cycle node4
   ar     = ar + abs(node(i)%ar)
   ar_min = min(ar_min, abs(node(i)%ar))
   ar_max = max(ar_max, abs(node(i)%ar))
   nnodes_eff = nnodes_eff + one
  end do node4

  ar = ar / nnodes_eff

 write(101,*)
 write(101,*) " Boundary nodes only"
 write(101,*) "   L1(AR) = ", ar
 write(101,*) "  Min(AR) = ", ar_min
 write(101,*) "  Max(AR) = ", ar_max
 write(101,*) " --------------------------------------------------------"

 end subroutine compute_ar



 end module edu2d_grid_data