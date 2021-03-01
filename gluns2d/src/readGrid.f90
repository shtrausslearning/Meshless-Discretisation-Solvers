
    
!! Input of grid coordinates, elements and boundaries.
!> Reads in grid data and boundary segments.
  subroutine ReadGrid
  use prmflow
  use ModInterfaces
  implicit none

! local variables
  integer :: ierr, i, ib, ibn, ibf, ibegf, iendf, ibegn, iendn
  integer :: tid,tid2,j,k,m,ii,wtid,fftid,ic
  
! temporary arrays
  integer, allocatable :: wbfconn(:) ! wall bc from bface
  integer, allocatable :: ffbfconn(:)  ! ff bc from bface
  integer, allocatable :: res(:)       ! output for bsort
  
  real(rtype) :: dxf,dyf,dsf,unx,uny,x0,y0,x1,y1
  integer :: p1,p2

! ###################################################################

  open(20, file='./grids/'//fnGrid, status="old", action="read", iostat=ierr)
  if (ierr /= 0) call EM( "cannot open grid file" )

  read(20,"(1X)")
  read(20,"(1X)")
  read(20,"(1X)")
  read(20,*) phynod,nt,nsegs ! numbers of physical nodes, triangles and boundary segments

! boundary type, no. of boundary faces & nodes, boundary name
  allocate( btype(nsegs),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for btype()" )
  allocate( bname(nsegs),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for bname()" )
  allocate( ibound(2,nsegs),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for ibound()" )

! read boundary elements
  read(20,"(1X)")
  do ib=1,nsegs
  read(20,  *  ) btype(ib),ibound(1,ib),ibound(2,ib)
  read(20,"(A)") bname(ib)
  enddo

  nbfaces = ibound(1,nsegs) ! final boundary left number 
  nbnodes = ibound(2,nsegs) ! final boundary right number

! definition of boundary faces / periodic nodes
  allocate( bnode(3,nbnodes),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for bnode()" )
  allocate( bface(2,nbfaces),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for bface()" )

  do ibn=1,nbnodes
    bnode(1,ibn) = -777
    bnode(2,ibn) = -777      ! set in DummyNodes
    bnode(3,ibn) = -777      ! set in EdgesFinalize
  enddo
  do ibf=1,nbfaces
    bface(1,ibf) = -777
    bface(2,ibf) = -777
  enddo

  tid = 0  ! counter for bnode
  tid2 = 0 ! counter for bface
  
  read(20,"(1X)")
  ibegf = 1
  ibegn = 1
  do ib=1,nsegs
  
    iendf = ibound(1,ib)
    iendn = ibound(2,ib)
  
!   periodic nodes always written first
    if (btype(ib)>=700 .and. btype(ib)<800) then   ! periodic nodes
     do ibn=ibegn,iendn
     tid = tid + 1                                 ! just counter
     read(20,*) bnode(1,ibn),bnode(2,ibn)
     enddo
    else                                           ! boundary faces
    do ibf=ibegf,iendf
     tid2 = tid2 + 1
     read(20,*) bface(1,ibf),bface(2,ibf)
     enddo
    endif
    
    ibegf = iendf + 1 ! update to next line
    ibegn = iendn + 1 ! update to next line
    
  enddo
  
! confirm, that bface array is completely defined
  do ibf=1,nbfaces
  if (bface(1,ibf)<0 .or. bface(2,ibf)<0) then
  call EM( "array bface() not completely defined" )
  endif
  enddo
    
  call DummyNodes ! generate dummy nodes
  
! grid nodes 1-allnod only
  ierr=0;allocate( x(allnod),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for x()" )
  ierr=0;allocate( y(allnod),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for y()" )

  read(20,"(1X)")
  do i=1,phynod
  read(20,*) x(i),y(i)
  enddo
  
! Nx,Ny at Wall Cells (all)
  
  ierr=0;allocate( nx_wc(phynod),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for local nx_wc()" )
  ierr=0;allocate( ny_wc(phynod),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for local ny_wc()" )
  nx_wc(:) = -777.0d0 ; ny_wc(:) = -777.0d0 ! reset
  
  do ibf=1,nbfaces
  
    p1 = bface(1,ibf)
    p2 = bface(2,ibf)

    dxf = x(p2) - x(p1)
    dyf = y(p2) - y(p1)
    dsf = dsqrt(dxf**2+dyf**2)
    unx = dyf/dsf  ! original
    uny = -dxf/dsf 
   ! unx = -dyf/dsf ! alternative
    !uny = dxf/dsf
    
    nx_wc(p2) = unx 
    ny_wc(p2) = uny
    
  enddo  
  

! ######################### wall boundary #########################
  wtid  = 0
  ibegf = 1
  do ib=1,nsegs
   iendf = ibound(1,ib)
   if (btype(ib)>=300 .and. btype(ib)<500) then   ! wall
   do ibf=ibegf,iendf
   wtid = wtid + 1
   enddo
   endif
   ibegf = iendf + 1 ! update to next line
  enddo
  
  allocate( wbfconn(2*wtid),stat=ierr )
  if (ierr /= 0) call EM( "cannot allocate memory for wbfconn()" )
  
  ii=0
  ibegf = 1
  do ib=1,nsegs
   iendf = ibound(1,ib)
   if (btype(ib)>=300 .and. btype(ib)<500) then   ! wall
   do ibf=ibegf,iendf
   ii=ii+1
   wbfconn(ii) = bface(1,ibf)
   enddo
   endif
   ibegf = iendf + 1 ! update to next line
  enddo
  
  ibegf = 1
  do ib=1,nsegs
   iendf = ibound(1,ib)
   if (btype(ib)>=300 .and. btype(ib)<500) then   ! wall
   do ibf=ibegf,iendf
   ii=ii+1
   wbfconn(ii) = bface(2,ibf)
   enddo
   endif
   ibegf = iendf + 1 ! update to next line
  enddo

! remove duplicates of wbfsort array
  allocate( res(2*wtid),stat=ierr )
  if (ierr /= 0) call EM( "cannot allocate memory for res()" )
  
  k = 1
  res(1) = wbfconn(1)
  outer: do i=2,2*wtid
     do j=1,k
        if (res(j) == wbfconn(i)) then
           ! Found a match so start looking again
           cycle outer
        end if
     end do
     ! No match found so add it to the output
     k = k + 1
     res(k) = wbfconn(i)
     
  end do outer
  
! number of wall boundary nodes
  nwbfnode = k
    
  allocate( wbfnode(nwbfnode),stat=ierr )
  if (ierr /= 0) call EM( "cannot allocate memory for wbfnode()" )
  
  do i=1,k
  wbfnode(i) = res(i)
  enddo
  
  deallocate( res )
  deallocate( wbfconn )
  
! ########################  farfield #############################

  fftid = 0
  ibegf = 1
  do ib=1,nsegs
   iendf = ibound(1,ib)
   if (btype(ib) .eq. 600 ) then   ! farfield
   do ibf=ibegf,iendf
   fftid = fftid + 1
   enddo
   endif
   ibegf = iendf + 1 ! update to next line
  enddo
  
  allocate( ffbfconn(2*fftid),stat=ierr )
  if (ierr /= 0) call EM( "cannot allocate memory for wbfconn()" )
  
  ii=0
  ibegf = 1
  do ib=1,nsegs
   iendf = ibound(1,ib)
   if ( btype(ib) .eq. 600 ) then   ! ff
   do ibf=ibegf,iendf
   ii=ii+1
   ffbfconn(ii) = bface(1,ibf)
   enddo
   endif
   ibegf = iendf + 1 ! update to next line
  enddo
  
  ibegf = 1
  do ib=1,nsegs
   iendf = ibound(1,ib)
   if ( btype(ib) .eq. 600 ) then   ! ff
   do ibf=ibegf,iendf
   ii=ii+1
   ffbfconn(ii) = bface(2,ibf)
   enddo
   endif
   ibegf = iendf + 1 ! update to next line
  enddo
  
! remove duplicates of wbfsort array
  allocate( res(2*fftid),stat=ierr )
  if (ierr /= 0) call EM( "cannot allocate memory for res()" )
  
  k = 1
  res(1) = ffbfconn(1)
  outer2: do i=2,2*fftid
     do j=1,k
        if (res(j) == ffbfconn(i)) then
           ! Found a match so start looking again
           cycle outer2
        end if
     end do
     ! No match found so add it to the output
     k = k + 1
     res(k) = ffbfconn(i)
     
  end do outer2
  
! number of farfield boundary nodes
  nffbfnode = k
  
  allocate( ffbfnode(nffbfnode),stat=ierr )
  if (ierr /= 0) call EM( "cannot allocate memory for ffbfnode()" )
  
  do i=1,k
  ffbfnode(i) = res(i)
  enddo
  
! confirm
  open(933,file='ffbfnode.dat',form='formatted')
  do i=1,k
  write(933,*) ffbfconn(i)
  enddo
  close(933)
  
  deallocate( res )
  deallocate( ffbfconn )

! ##################################################################
  
! triangles
  ierr=0
  allocate( elem(3,nt),stat=ierr )
  if (ierr /= 0) call EM( "cannot allocate memory for elem()" )

  read(20,"(1X)")
  do i=1,nt
  read(20,*) elem(1,i),elem(2,i),elem(3,i)
  enddo
  close(20)
  
! node type allocation
  ierr = 0
  allocate( ntype(allnod),stat=ierr )
  if (ierr /= 0) call EM( "cannot allocate memory for ntype()" )
  
  do i=1,allnod
  ntype(i) = -777
  enddo
  
! initialisation : all fluid nodes
  do i=1,phynod
  ntype(i) = 3                  ! fluid node
  enddo

! physical farfield nodes 
  do i=1,nffbfnode
  ic = ffbfnode(i)
  do j=1,phynod
  if( j .eq. ic )then
  ntype(j) = 2 
  endif
  enddo
  enddo
  
! physical wallnodes
  do i=1,nwbfnode
  ic = wbfnode(i)
  do j=1,phynod
  if( j .eq. ic )then
  ntype(j) = 1 
  endif
  enddo
  enddo
  
! enter p

  return
  end subroutine ReadGrid  

  
  !open(332,file='form.dat',form='formatted')
  !open(333,file='form2.dat',form='formatted')
  !do i=1,phynod
  !if( ntype(i) .eq. 1 )then
  !x0=x(i)
  !y0=y(i)
  !x1= x0 + nx_wc(i)*0.01d0
  !y1= y0 + ny_wc(i)*0.01d0
  !write(332,*) x0,y0
  !write(333,*) x1,y1
  !endif
  !enddo
  !close(332)
  !close(333)
  !
  !pause '333'
  !
  !open(321,file='bface1.dat',form='formatted')
  !do ibf=1,nbfaces
  !!write(321,*) x( bface(1,ibf) ),y( bface(1,ibf) )
  !write(321,*) bface(1,ibf)
  !enddo
  !close(321)
  !
  !open(322,file='bface2.dat',form='formatted')
  !do ibf=1,nbfaces
  !!write(322,*) x( bface(2,ibf) ),y( bface(2,ibf) )
  !write(322,*) bface(2,ibf)
  !enddo
  !close(322)
  !
  !pause '322'