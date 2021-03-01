!     Generate points and connectivity
      subroutine GenPtsConn
      use prm_flow
      implicit none

      integer          :: c, i, j, exists, p0, p1, p2, iseed, ids
      double precision :: x, y, drand, rdx, rdy
      double precision :: dx,dy

!     seed for randum number generator
      iseed = 29374

      do i=1,npmax
       ptype(i)   = -1
       nnbr(i)    = -1
       coord(1,i) = 2.0d0*dmax1(xmin, xmax)
       coord(2,i) = 2.0d0*dmax1(ymin, ymax)
      enddo

      lx = xmax - xmin
      ly = ymax - ymin
      dx = lx/(nx-1)
      dy = ly/(ny-1)
      write(*,*) 'lx, ly =',lx,ly
      write(*,*) 'dx, dy =',dx,dy

      c  = 0
!     Generate Point Type, X,Y coordinates
      do i=1,nx                      ! first to last x
      do j=1,ny                      ! first to last y
      
       c          = c + 1            ! global node count
       x          = xmin + (i-1)*dx  ! spacing of x
       y          = ymin + (j-1)*dy  ! spacing of y
       
       if(perturb.eq.1)then ! add perturbation to all internal nodes except outer
        if( i.ne.1 .and. i.ne.nx .and. j.ne.1 .and. j.ne.ny )then
        rdx = drand(iseed)
        rdy = drand(iseed)
        x   = x + rdx*pertx*dx
        y   = y + rdy*perty*dy
        endif
       endif
       
       coord(1,c) = x  ! global coordinates
       coord(2,c) = y
       ptype(c)   = 1  ! internal node 
       
       if(i.eq.nx .or. j.eq.ny)then
       ptype(c)= 2   ! external node
       endif
       
       nnum(i,j) = c ! global number of i,j
       
      enddo
      enddo

!     These are points in the computational domain
      np0 = c

!     Generate Ghost Node Identity
!     ############################
      
!     bottom layer
      do i=1,nx
      c         = c + 1
      nnum(i,0) = c
      ptype(c)  = 2
      enddo

!     left layer
      do j=1,ny
      c         = c + 1
      nnum(0,j) = c
      ptype(c)  = 2
      enddo

!     Left-bottom corner point
      c         = c + 1
      nnum(0,0) = c
      ptype(c)  = 2

      npts = c   ! total number of points, including ghost layers
      print*,'Number of points =',npts

      
!     Define the neighbours of all bounded nodes
!     ###########################################

      do i=1,nx-1
      do j=1,ny-1
      
      c  = 0          ! neighbour counter
      p0 = nnum(i,j)  ! global number of i,j
      
!     Face Neighbours

!     -x neighbour node
      if( exists(i-1,j) .eq. 1 )then
      c          = c + 1
      conn(c,p0) = nnum(i-1,j)
      endif

!     +x neighbour node
      if( exists(i+1,j).eq.1 )then
      c          = c + 1
      conn(c,p0) = nnum(i+1,j)
      endif

!     -y neighbour node
      if(exists(i,j-1).eq.1)then
      c          = c + 1
      conn(c,p0) = nnum(i,j-1)
      endif

!     +y neighbour node
      if(exists(i,j+1).eq.1)then
      c          = c + 1
      conn(c,p0) = nnum(i,j+1)
      endif
      
!     Diagonal Neighbours 
      
!!     -x,-y neighbour node
!      if( exists(i-1,j-1) .eq. 1 )then
!      c = c + 1
!      conn(c,p0) = nnum(i-1,j-1)
!      endif
!      
!!     -x,+y neighbour node
!      if( exists(i-1,j+1) .eq. 1 )then
!      c = c + 1
!      conn(c,p0) = nnum(i-1,j+1)
!      endif
!      
!!     +x,-y neighbour node
!      if( exists(i+1,j-1) .eq. 1 )then
!      c = c + 1
!      conn(c,p0) = nnum(i+1,j-1)
!      endif
!      
!!     +x,+y neighbour node
!      if( exists(i+1,j+1) .eq. 1 )then
!      c = c + 1
!      conn(c,p0) = nnum(i+1,j+1)
!      endif
      
      nnbr(p0) = c ! total neighbours [ 4 face 4 diagonal ]
      !if( nnbr(p0) < 8 )then
      ! pause 'nnbr(p0) < 8'
      !endif 

      enddo
      enddo

!     Bottom ghost layer [ except bot right ] ///////////////////////////////
      do i=1,nx-1
       p0         = nnum(i,0)         ! bottom node itself
       p1         = nnum(i,ny-1)      ! second from top node 
       nnbr(p0)   = 1                 ! bot node has 1 connection
       conn(1,p0) = p1                ! the connection point is p1
       coord(1,p0)= coord(1,p1)       ! x coordinate of bot ghost point
       coord(2,p0)= coord(2,p1) - ly  ! y coordinate of bot ghost point
      enddo
!     only the bot right node case
      p0         = nnum(nx,0)         ! bot right corner node only
      p1         = nnum(1,ny-1)       ! top left of ptype(i)=1
      nnbr(p0)   = 1                  ! also only 1 connection
      conn(1,p0) = p1                 ! its connection is p1
      p2         = nnum(nx,ny-1)      ! use to define poistion of p0 only
      coord(1,p0)= coord(1,p2)        ! x coordinate of bot right ghost point
      coord(2,p0)= coord(2,p2) - ly   ! y coordinate of bot right ghost point
      
!     Left Ghost Layer //////////////////////////////////////////////////////
      do j=1,ny-1
       p0         = nnum(0,j)         ! left node itself
       p1         = nnum(nx-1,j)      ! last node of ptype(i)=1
       nnbr(p0)   = 1                 ! it has one connection
       conn(1,p0) = p1                ! which is p1
       coord(1,p0)= coord(1,p1) - lx  ! x coordinate of left ghost point
       coord(2,p0)= coord(2,p1)       ! y coordinaet of left ghost point
      enddo
!     only the top left node case
      p0         = nnum(0,ny)         ! top left node itself
      p1         = nnum(nx-1,1)       ! bot right ptype(i) = 1 node
      nnbr(p0)   = 1                  ! only 1 connection  
      conn(1,p0) = p1                 ! which is p1 
      p2         = nnum(nx-1,ny)      ! p2 only used to define poision of p0
      coord(1,p0)= coord(1,p2) - lx
      coord(2,p0)= coord(2,p2)

!     Left-botton corner point  ////////////////////////////////////////////
      p0         = nnum(0,0)          ! bot left corner ghost node itself
      p1         = nnum(nx-1,ny-1)    ! top right node of ptype(i)=1  
      nnbr(p0)   = 1                  ! only one connection
      conn(1,p0) = p1                 ! which is p1
      coord(1,p0)= coord(1,p1) - lx   ! x coordinate of bot left ghost node
      coord(2,p0)= coord(2,p1) - ly   ! y coordinate of bot right ghost node

!     Top periodic layer [ from original point creation ] /////////////////
      do i=1,nx-1                     
       p0 = nnum(i,ny)                ! corresponding top node
       nnbr(p0) = 1                   ! one neighbour only
       conn(1,p0) = nnum(i,1)         ! which is the bot ghost 
      enddo

!     Right periodic layer [ from original point creation ] //////////////
      do j=1,ny-1
         p0 = nnum(nx,j)              ! corresponding right node
         nnbr(p0) = 1                 ! which has one neighbour
         conn(1,p0) = nnum(1,j)       ! most left ptype(i)=1 node
      enddo

!     Top-Right corner point [ from original point creation ] ////////////
      p0 = nnum(nx,ny)          ! top right node itself ( ptype(i) = 2 )
      nnbr(p0) = 1              ! which has one node
      conn(1,p0) = nnum(1,1)    ! most bot left ptype(i) = 1 node

!     check
      do i=1,npts
      if(ptype(i).eq.-1)then
       write(*,*) '!!! ptype not set for node =',i
      endif
      if(nnbr(i).eq.-1)then
       write(*,*) '!!! nnbr not set for node =',i
      endif
      enddo
      
!      ids=101
!      write(*,*) ptype(ids)
!      open(333,file='test.dat',form='formatted')
!!      write(333,*) coord(1,230), coord(2,230)
!      do j=1,nnbr(ids)
!      write(333,*) coord(1,conn(j,ids)), coord(2,conn(j,ids))
!      enddo
!      close(333)
!      pause

      return
      end


!   point (i,j) exists in the grid
    integer function exists(i, j)
    use prm_flow

    implicit none
    integer :: i, j

    exists = 1

    if(i.lt.0)  exists = 0
    if(i.gt.nx) exists = 0
    if(j.lt.0)  exists = 0
    if(j.gt.ny) exists = 0

    return
    end


    real*8 function drand(iseed)
    use prm_flow

!                  mod(iseed*7141 + 54773, 259200)
!         ran = -----------------------------------
!                            259200

    integer :: iseed
    integer :: ia, ic, im
    parameter(ia = 7141, ic = 54773, im = 259200)
!
    iseed    = abs(mod(iseed*ia+ic, im))
!
    drand    = dble(iseed)/dble(im)
!
    return
    end