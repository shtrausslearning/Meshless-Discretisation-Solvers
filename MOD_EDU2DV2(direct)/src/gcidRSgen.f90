
!   GCID array based GC generation ; 
!   fluid nodes with wall neighbours are identified
!                               

    subroutine gcidRS
    use edu2d_constants
    use prmflow
    use edu2d_my_main_data
    use edu2d_grid_data_type
    use edu2d_my_allocation
    implicit none
    
    integer :: i,inn,j,iii,ierr, id1,id2, ii2, iiff,iiw, wnbr, iwn, wn
    real(p2) :: dxx,dyy,x0,y0,mid,cid
    character(256) :: linebuf
    integer :: fwnid(nnodes)          ! mask to identify IBW with wall neighbours
    integer, allocatable :: wnbrID(:) ! array to store ID1's wall neighbours
    
    fwnid=0 ! reset
!   (1) identify IBW nodes that have boundary neighbours
    do i=1,nnodes
    if( node(i)%ptype .eq. 3 )then ! if IBW fluid node
    do j=1,node(i)%nnghbrs         ! for all physical neighbours nnghbr -> physical nodes only
    inn = node(i)%nghbr(j)         ! corresp. neighbour node ID
    
!   identify wall neighbours
    if( node(inn)%ptype .eq. 1 )then ! if IDW has wall neighbour
    fwnid(i) = 1 
    endif
    if( node(inn)%ptype .eq. 2 )then ! if farfield node
    fwnid(i) = 2 
    endif
    
    enddo
    endif
    enddo
    
!   (2) now count how many IDW were seleted ( includes all boundaries)
    iiW=0;iiFF=0
    do i=1,nnodes
    if( fwnid(i) .eq. 1 )then ! if IDW node has wall neighbour
    iiW=iiW+1
    elseif( fwnid(i) .eq. 2 )then ! if IFW node has ff neighbour
    iiFF=iiFF+1
    endif
    enddo
    
!   (3) allocate gcid global array (both wall and ff)
    ierr=0;allocate( gcid(2,iiw+iiff),stat=ierr )
    if( ierr /= 0 ) pause 'allocation error gcid'
    gcid=0 ! reset
    
    igcs = iiw+iiff ! set the global gcid array maximum
    
!   (4) fill in the gcid array
    iii=nnodes  
    ii2=0 ! gcid counter only
    
!   write wall nodes first (more significant in BC)
    do i=1,nnodes               ! go through all phy nodes
    if( fwnid(i) .eq. 1 )then   ! if IDW node has wall neighbour
    iii=iii+1;ii2=ii2+1
    gcid(1,ii2) = i
    gcid(2,ii2) = iii
    endif
    enddo
    
!   now add the farfield nodes
    do i=1,nnodes
    if( fwnid(i) .eq. 2 )then ! if IDW has a ff neighbour boundary
    iii=iii+1 ! continue counting
    ii2=ii2+1 ! continue counting   
    gcid(1,ii2) = i
    gcid(2,ii2) = iii
    endif
    enddo
    
!   check
    if( iii .ne. nnodes+igcs ) pause 'iii/=nnodes+igcs'

!   (5) generate ptype info for node(i)%ptype of gc nodes
    do i=1,igcs
    id1=gcid(1,i);id2=gcid(2,i)
    if( fwnid(id1) .eq. 1 ) then ! if wall GC node
    node(id2)%ptype = 4
    elseif( fwnid(id1) .eq. 2 )then ! if ff gc node
    node(id2)%ptype = 5
    endif
    enddo
    
!   (6) generate node(i)%x,y data for all ghost nodes
    do i=1,igcs
   
     id1=gcid(1,i)  ! fluid IBW node   
     id2=gcid(2,i)  ! ghost node pair

!    count how many wall nodes are connected to ID1 node
     wnbr=0
     do j=1,node(id1)%nnghbrs         ! for all physical neighbours nnghbr -> physical nodes only
     inn = node(id1)%nghbr(j)         ! corresp. neighbour node ID
     if( node(inn)%ptype .eq. 1 .or. node(inn)%ptype .eq. 2  )then 
     wnbr=wnbr+1 ! count how many wall node neighbours ID1 has
     endif
     enddo
     if( wnbr == 0 ) pause 'gcid1 has no wall neighbours?'
     allocate(wnbrID(wnbr));wnbrID=0 ! allocate for wallnode storage,reset
    
     iwn=0
     do j=1,node(id1)%nnghbrs         ! for all physical neighbours nnghbr -> physical nodes only
     inn = node(id1)%nghbr(j)         ! corresp. neighbour node ID
     if( node(inn)%ptype .eq. 1 .or. node(inn)%ptype .eq. 2  )then 
     iwn=iwn+1 ! storage counter
     wnbrID(iwn) = inn
     endif
     enddo
     
!    if fluid IDW node ID1, has more than 1 wall node neighbour
     if( iwn .gt. 1 )then 
    
!    if more than 2 neighbours
     dxx = node( wnbrID(1) )%x - node( wnbrID(2) )%x
     dyy = node( wnbrID(1) )%y - node( wnbrID(2) )%y
     if( dxx .eq. 0.0d0 ) dxx = 1.0d-10 ! in case
     if( dyy .eq. 0.0d0 ) dyy = 1.0d-10 ! in case
     x0 = node(id1)%x ; y0 = node(id1)%y ! x,y coordinate of IDW fluid node
    
     mid=dyy/dxx
     cid= node( wnbrID(1) )%y - mid * node( wnbrID(1) )%x 
     node(id2)%x = ((1.0d0-mid**2.0d0)*x0 + 2.0d0*mid*(y0-cid))/(1.0d0+mid**2.0d0)
     node(id2)%y = ( 2.0d0*mid*x0 - (1.0d0-mid**2.0d0)*y0 + 2.0d0*cid )/(1.0d0+mid**2.0d0)

!    if fluid IDW node ID1 has only 1 neighbour
     elseif( iwn == 1 )then ! if one node only ( lets simply extend edge )
     
     dxx = node(wnbrID(1))%x - node(id1)%x
     dyy = node(wnbrID(1))%y - node(id1)%y
     node(id2)%x = node(id1)%x + 2.0d0*dxx
     node(id2)%y = node(id1)%y + 2.0d0*dyy
    
     endif
     
!    (7) lets also add fluid IDW nodes to neighbours of wall nodes

     do j=1,iwn ! iwn is the max neighbour count; same as wnbr
         
      wn = wnbrID(j) ! wall node ID number
         
      node(wn)%nnghbrs = node(wn)%nnghbrs + 1                    ! add extra connection 
      call my_alloc_int_ptr(node(wn)%nghbr, node(wn)%nnghbrs)    ! expand array dynamically
      node(wn)%nghbr(node(wn)%nnghbrs) = id2                     ! add corresponding node
      
    enddo
     
!   ######################################################################################
    
!   deallocate local
    deallocate(wnbrID)

    enddo
    
!   check (output gcid1)
    open(44,file='./info/gcIDNEW1.dat',form='formatted')
    do i=1,igcs
    id1 = gcid(1,i)
    id2 = gcid(2,i)
    write(linebuf,*) node(id1)%x,',',node(id1)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
    enddo
    close(44)
    
!   check (output gcid2)
    open(44,file='./info/gcIDNEW2.dat',form='formatted')
    do i=1,igcs
    id1 = gcid(1,i)
    id2 = gcid(2,i)
    write(linebuf,*) node(id2)%x,',',node(id2)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
    enddo
    close(44)
    
!!   [check ] check wall neighbour
!    i=167
!    do j=1,node(i)%nnghbrs
!    inn = node(i)%nghbr(j)
!    write(*,*) node(inn)%ptype,inn
!    enddo
!    
!   check (output gcid2)
    open(44,file='./info/wallneighb.dat',form='formatted')
    i=64
    do j=1,node(i)%nnghbrs
    inn = node(i)%nghbr(j)
    write(linebuf,*) node(inn)%x,',',node(inn)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
    enddo
    close(44)
    
!    pause 'please confirm ghost nodes'

    return
    end subroutine gcidRS