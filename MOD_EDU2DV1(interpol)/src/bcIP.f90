

!   Subroutine to create ghost cells for all boundaries
!   IP is added to node(), given x,y and neighbour data for interpolation

    subroutine bcIP
    use edu2d_constants, only : p2
    use prmflow, only : nIP,ngc
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    use edu2d_my_allocation
    implicit none
    
    logical, parameter :: confirmID = .true.  ! show confirmation 
    integer :: i,j,ii,iii,nid,inn,v,vv
    real(p2) :: n_dx,n_dy,dr,drav,x0,x1,y0,y1
    character(256) :: linebuf
    integer, allocatable :: IPnbr(:) ! local nodes used to determine IP's neighbours
    integer :: iiii,ierr
    
    iii=0
    ii=nnodes+ngc                 ! starting point after which IPs are generated in nnodes()%...
    do i=1,nbound                 ! cycle through all boundary segments
    do j=1,bound(i)%nbnodes-1     ! boundary node #
            
     iii=iii+1                               ! general counter to determine final node value after IP generation
     ii=ii+1                                 ! counter over nnodes
     n_dx = bound(i)%bnx(j)                  ! boundary node nx of normal vector
     n_dy = bound(i)%bny(j)                  ! boundary node ny of normal vector
     nid = node(bound(i)%bnode(j))%nnghbrs-1 ! current neighbour count for boundary node (exclude GC, should be last, need check)
     
!    (1) first determine IP's X,Y location -----------------------------------------------------

!    (1a) using average radius of boundary node's physical neighbours
     drav=0.0d0
     do v=1,nid
     inn = node(bound(i)%bnode(j))%nghbr(v) ! neighbour node of boundary node
     x1 = node(node(bound(i)%bnode(j))%nghbr(v))%x ; y1 = node(node(bound(i)%bnode(j))%nghbr(v))%y
     x0 = node(bound(i)%bnode(j))%x                ; y0 = node(bound(i)%bnode(j))%y
     dr = sqrt( (x1-x0)**2 + (y1-y0)**2 ) ! local edge dr
     drav = drav + dr
     enddo
     drav = drav/dble(nid) ! local edge radius
    
!    (1b) using the boundary node's unit directivity find x,y 
     node(ii)%x = node(bound(i)%bnode(j))%x - n_dx*drav
     node(ii)%y = node(bound(i)%bnode(j))%y - n_dy*drav
     
!    (2) now lets find its neighbour connectivity ( simplest case first using wall node's fluid neighbours )

!    1) Add IBW fluid nodes to neighbour of IP -------------------------------------------------

     do v=1,nid                             ! search through all neighbours of wall node
     inn = node(bound(i)%bnode(j))%nghbr(v) ! neighbour node of boundary node
     if( node(inn)%ptype .eq. 3 )then       ! if the boundary node's neighbour is IBW node           
      node(ii)%nnghbrs = node(ii)%nnghbrs + 1                    ! add number of neighours of IP 
      call my_alloc_int_ptr(node(ii)%nghbr, node(ii)%nnghbrs)    ! expand array dynamically
      node(ii)%nghbr(node(ii)%nnghbrs) = inn                     ! add corresponding node
     endif
     enddo
     
!    check
     if( node(ii)%nnghbrs .eq. 0 ) pause 'no neighbours found'

!!    2) Add IBW's neighbour nodes -------------------------------------------------------------
!
!!    (2a1) lets first count how many nodes there are for subsequent array allocation 
!     iiii=0
!     do v=1,nid                             ! search through all neighbours of wall node
!     inn = node(bound(i)%bnode(j))%nghbr(v) ! neighbour node of boundary node
!     if( node(inn)%ptype .eq. 3 )then       ! if the boundary node's neighbour is IBW node           
!     iiii=iiii + node(inn)%nnghbrs          ! number of neighbours which IBW has
!     endif
!     enddo
!     
!!    (2a2) allocate array to store all of IP's neighbours
!     ierr=0;allocate( IPnbr(iiii),stat=ierr )
!     if( ierr /= 0 ) pause 'allocation error IPnbr'
!     IPnbr=0 ! reset
!     
!!    Array IPnbr will contain overlaps, but not translated into %nghbr
!     
!!    and add to the IPnbr array of course
!     iiii=0
!     do v=1,nid                             ! search through all neighbours of wall node
!     inn = node(bound(i)%bnode(j))%nghbr(v) ! neighbour node of boundary node
!     if( node(inn)%ptype .eq. 3 )then       ! if the boundary node's neighbour is IBW node           
!           
!      do vv=1,node(inn)%nnghbrs             ! go through all neighbours of IBW's node
!      iiii=iiii+1       
!      IPnbr(iiii) = node(inn)%nghbr(vv)     ! add the ID to the IPnbr array
!      enddo
!      
!     endif
!     enddo
!    
!!    check the array is correctly allocated/full
!     do v=1,iiii
!     if( IPnbr(v) .eq. 0 ) pause 'allocation error'
!     enddo
!
!!    (2a3) Add gc node to the neighbor list of boundary node, node(bound(i)%bnode(j))
!     do v=1,iiii                                                 ! add all nodes of IPnbr array to IP's neighbours
!      if( ANY( node(ii)%nghbr .eq. IPnbr(v) )) cycle             ! add only to nghbr if not already included
!      node(ii)%nnghbrs = node(ii)%nnghbrs + 1                    ! add number of neighours of IP 
!      call my_alloc_int_ptr(node(ii)%nghbr, node(ii)%nnghbrs)    ! expand array dynamically
!      node(ii)%nghbr(node(ii)%nnghbrs) = IPnbr(v)                ! add corresponding node
!     enddo
!     
!!    (2a4) deallocate array for subsequent reuse 
!     deallocate( IPnbr ) ! deallocate locally
!     
!!   -------------------------------------------------------------------------------------------

!    added two new ptypes (might be useful) : 
!    ptype( node(:)... ) = 6 ( farfield IP node )
!    ptype( node(:)... ) = 7 ( wall IP node )

!    (3) define local ghost cell ptype from boundary type
     if (trim(bound(i)%bc_type) == "freestream") then
     node(ii)%ptype = 6
     elseif(trim(bound(i)%bc_type) == "slip_wall_weak")then
     node(ii)%ptype = 7  
     endif
     
    enddo
    enddo
    
!   define the subsequent number of IPs generated
    nIP = iii
    
!   check correct allocation in python  #################
    if( confirmID .eqv. .true. )then
    
    open(44,file='./info/xyIPall.dat',form='formatted')
    do i=nnodes+ngc+1,nnodes+ngc+nIP
    write(linebuf,*) node(i)%x,',',node(i)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
!    write(44,*) node(i)%x,node(i)%y
    enddo
    close(44)
    
    open(44,file='./info/xyIPFF.dat',form='formatted')
    do i=1,nnodes+ngc+nIP
    if( node(i)%ptype .eq. 6 )then
    write(linebuf,*) node(i)%x,',',node(i)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
!    write(44,*) node(i)%x,node(i)%y
    endif
    enddo
    close(44)
    
    open(44,file='./info/xyIPW.dat',form='formatted')
    do i=1,nnodes+ngc+nIP
    if( node(i)%ptype .eq. 7 )then
    write(linebuf,*) node(i)%x,',',node(i)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
!    write(44,*) node(i)%x,node(i)%y
    endif
    enddo
    close(44)
    
    endif
!   ########################################################

!   check 2 : IPs neighbours check
    open(44,file='./info/IPsNEIGH.dat',form='formatted')
    i=nnodes+ngc+90
    do v=1,node(i)%nnghbrs
    inn = node(i)%nghbr(v) ! neighbour node of boundary node
    write(linebuf,*) node(inn)%x,',',node(inn)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
    enddo
    close(44)
    
!    open(44,file='temp.dat',form='formatted')
!    do i=nnodes+ngc+1,nnodes+ngc+nIP
!    write(44,*) node(i)%nghbr(:)
!    enddo
!    close(44)
    
    return
    end subroutine bcIP
    
    
    
    
     
    