
    Subroutine make_wgc
    use ModDataTypes
    use prmflow
    use ModInterfaces, only : EM
    implicit none
    
    integer :: nidd,ii,i,zz,ierr,itemp,ic,j,wnbr,iii,ix,ibb,nibb,ix2,ibb2,jj,p1,p2,inn
    real(rtype) :: dxx,dyy,mid,cid,x0,y0,x_gc,y_gc
    
    integer, allocatable :: wnp(:)     ! loop wall node allocation
    real(rtype), allocatable :: xta(:) ! temporary x,y arrays
    real(rtype), allocatable :: yta(:) !
    integer, allocatable :: tntype(:)  ! temporary ntype array
!   #########################################################################
    
!   TEMPORARY X,Y ARRAYS WITH ADDED SLOT FOR RESET WALL GHOST CELLS X,Y
    
!   Now redefine, xy array to include additional gc components of wall
    ierr=0;allocate( xta(allnod+igcs),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for xtemp()" )
    ierr=0;allocate( yta(allnod+igcs),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for ytemp()" )
    xta(:) = -777.0d0 ; yta(:) = -777d0
  
    do i=1,allnod
    xta(i) = x(i);yta(i) = y(i)
    enddo
  
    deallocate( x );deallocate( y )
    
    !write(*,*) phynod,allnod
    
    !open(54,file='confirmmee.dat',form='formatted')
    !do i=1,ngcid
    !write(54,*) gcid(1,i),gcid(2,i)
    !enddo
    !close(54)
    !
    !pause '54'
    
!   DEFINE WALL GHOST CELL X,Y COORDINATES FOR GHOST CELLS [ MAIN LOOP ] 
    
!   Main loop ( uses xta,yta )
    nidd = 0
    ii=0  ! reset
    do i=1,phynod 
    if( fwnid(i) .eq. 1 )then   ! if fluid node was identified as wall neighbour (through edge connection)
  
      itemp=-777;nidd = nidd + 1
    ! find ghost cell id
      do j=1,ngcid
      if( gcid(1,j) .eq. i )then ! run through all 
      itemp = gcid(2,j)          ! fluid related ghost cell id ( itemp = gcid(2,i) )       
      endif
      enddo
      if( itemp .eq. -777 ) then 
      pause 'itemp = -777' 
      endif
  
      wnbr=0 ! look through conn and count wall neighbours [ for wnp local allocation ]
      do j=1,nbers(i)              ! do for all neighbours of current fluid node fwnid(i) = 1
      inn = conn(j,i)              ! corresponding neighbour node
      if( ntype(inn) .eq. 1 )then  ! if neighbour is a wall node; count it
      wnbr=wnbr+1
      endif
      enddo
      if( wnbr .eq. 0 )then
      pause 'fwnid found no wall neighbours'
      endif
  
    ! allocate/deallocate within loop
      ierr=0;allocate( wnp(wnbr),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for wnp()" )
  
      ii=0; ! Store the wall neighbour nodes into wnp() array
      do j=1,nbers(i)             ! do for all neighbours
      inn=conn(j,i)                ! neighbour
      if( ntype(inn) .eq. 1 )then  ! if neighbour is wall
      ii=ii+1
      wnp(ii) = inn                ! corresponding wall node array for fluid node
      endif
      enddo
  
    ! local wnp() array contains all wall nodes that are connected to the fluid node [ from conn ]
      
    ! Add ghost cell neighbour to all wnp() wall nodes
      do j=1,wnbr                  ! for all wall neighbour nodes
      
!      Main manipulator
       ibb  = wnp(j)                ! global node number of wall [1,phynod]
       nibb = nbers(ibb)            ! current number of neighbours of this wall node
      
       do ii=1,igcs                  ! search for main node
        if( gcid(1,ii) .eq. i )then  ! find the fwnid node in gcid(1,:)
        ibb2 = gcid(2,ii)            ! identify the ghost cell id [ ibb2 ]
        endif
       enddo
       
!      update neighbour info for main manipulator, by updating w/ ghost cell id
       nbers(ibb) = nibb+1           ! add additional connection
       conn(nibb+1,ibb) = ibb2       ! add ghost cell id
       
      enddo
  
    ! select corresponding scenario for fluid node
      
      if( wnbr .eq. 1 )then          ! if only 1 wall node neighbour    
    ! #########################################################################
  
      do jj=1,nbfaces
      if( bface(1,jj) .eq. wnp(1) )then  
      p1 = wnp(1)                        ! corresponding wall neighbour
      p2 = bface(2,jj)                   ! bface edge partner
      endif
      enddo
  
      dxx = xta(p2) - xta(p1)
      dyy = yta(p2) - yta(p1)
      if( dxx == 0.0d0 )then
      dxx = 1d-20
      endif
  
      mid = dyy/dxx
      cid = yta(p1) - mid*xta(p1)
  
      x0 = xta(i) ! actual main fluid node
      y0 = yta(i) 
  
    ! Reflection Point X,Y Coordinates 
      x_gc = ((1.0d0-mid**2)*x0 + 2.0d0*mid*(y0-cid))/(1.0d0+mid**2)
      y_gc = ( 2.0d0*mid*x0 - (1.0d0-mid**2)*y0 + 2.0d0*cid )/(1.0d0+mid**2)
    ! storage of gc value
      xta(itemp) = x_gc
      yta(itemp) = y_gc
  
      elseif( wnbr .eq. 2 )then ! if 2 wall node neighbours
    ! #########################################################################
  
      p1 = wnp(1)
      p2 = wnp(2)
  
      dxx = xta(p2) - xta(p1)
      dyy = yta(p2) - yta(p1)
      if( dxx == 0.0d0 )then
      dxx = 1d-20
      endif
  
      mid = dyy/dxx
      cid = yta(p1) - mid*xta(p1)
  
      x0 = xta(i) ! actual main fluid node
      y0 = yta(i) 
  
    ! Reflection Point X,Y Coordinates 
      x_gc = ((1.0d0-mid**2)*x0 + 2.0d0*mid*(y0-cid))/(1.0d0+mid**2)
      y_gc = ( 2.0d0*mid*x0 - (1.0d0-mid**2)*y0 + 2.0d0*cid )/(1.0d0+mid**2)
    ! storage of gc value
      xta(itemp) = x_gc
      yta(itemp) = y_gc
  
      elseif( wnbr .eq. 3 )then ! if 3 wall node neighbours exists
    ! #########################################################################
  
      p1 = wnp(2) ! second in order
      p2 = wnp(3) ! third in order
  
      dxx = xta(p2) - xta(p1)
      dyy = yta(p2) - yta(p1)
      if( dxx == 0.0d0 )then
      dxx = 1d-20
      endif
  
      mid = dyy/dxx
      cid = yta(p1) - mid*xta(p1)
  
      x0 = xta(i) ! actual main fluid node
      y0 = yta(i) 
  
    ! Reflection Point X,Y Coordinates 
      x_gc = ((1.0d0-mid**2)*x0 + 2.0d0*mid*(y0-cid))/(1.0d0+mid**2)
      y_gc = ( 2.0d0*mid*x0 - (1.0d0-mid**2)*y0 + 2.0d0*cid )/(1.0d0+mid**2)
    ! storage of gc value
      xta(itemp) = x_gc
      yta(itemp) = y_gc
  
      elseif( wnbr .eq. 4 )then ! if 3 wall node neighbours exists
    ! #########################################################################

      p1 = wnp(2)
      p2 = wnp(3) 
  
      dxx = xta(p2) - xta(p1)
      dyy = yta(p2) - yta(p1)
      if( dxx == 0.0d0 )then
      dxx = 1d-20
      endif
  
      mid = dyy/dxx
      cid = yta(p1) - mid*xta(p1)
  
      x0 = xta(i) ! actual main fluid node
      y0 = yta(i) 
  
    ! Reflection Point X,Y Coordinates 
      x_gc = ((1.0d0-mid**2)*x0 + 2.0d0*mid*(y0-cid))/(1.0d0+mid**2)
      y_gc = ( 2.0d0*mid*x0 - (1.0d0-mid**2)*y0 + 2.0d0*cid )/(1.0d0+mid**2)
    ! storage of gc value
      xta(itemp) = x_gc
      yta(itemp) = y_gc
  
      endif
    ! #########################################################################
  
      deallocate( wnp ) ! local allocate/reallocate
  
    endif
    enddo ! End main loop
      
!   now redefine x,y main array
    ierr=0;allocate( x(allnod+igcs),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for x() again" )
    ierr=0;allocate( y(allnod+igcs),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for y() again" )
  
!   add all physical,ff/w gcs to main xy array
    do i=1,allnod+igcs
    x(i) = xta(i);y(i) = yta(i)
    enddo
    
    do j=i,allnod+igcs
    if( x(i) .eq. -777.0d0 .or. y(i) .eq. -777.0d0 )then
    pause 'allocation error in make_wgc'
    endif
    enddo

    deallocate( xta );deallocate( yta ) ! deallocate temporary arrays
     
!   Update ntype ( add ghost cells )
    ierr=0;allocate( tntype(allnod+igcs),stat=ierr );if (ierr /= 0) call EM( "cannot allocate memory for tntype() again" )
  
    do i=1,allnod
    tntype(i) = ntype(i)
    enddo
    do i=allnod+1,allnod+igcs
    tntype(i) = 4 ! wall ghost cells
    enddo
  
!   redefine ntype (point type)
    deallocate( ntype )
    allocate( ntype(allnod+igcs),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for ntype again() again" )
  
    do i=1,allnod+igcs
    ntype(i) = tntype(i)
    enddo
  
    deallocate( tntype );deallocate( fwnid  ) ! deallocate temporary array
    
    return
    End subroutine make_wgc