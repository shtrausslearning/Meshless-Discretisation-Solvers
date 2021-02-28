
!   Subroutine to create ghost cells for all boundaries
!   ghost cells are generated for each boundary node, as opposed to old approach
!   cell values will be averaged based on fluid neighbour nodes of wall

    subroutine bcgc
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    use edu2d_my_allocation
    implicit none
    
    logical, parameter :: confirmID = .true.  ! show confirmation 
    integer :: i,j,ii,iii,nid,inn,v
    real(p2) :: n_dx,n_dy,dr,drav,x0,x1,y0,y1
    character(256) :: linebuf

    iii=0
    ii=nnodes                   ! starting node value ( generate IP after GC only )
    do i=1,nbound               ! cycle through all boundary segments
    do j=1,bound(i)%nbnodes-1   ! cycle though all nods in the particular segment (-1 needed to not include 1st node again)
        
     iii=iii+1 ! counter to determien how many ghost cells are generated  
     ii=ii+1   ! counter over nnodes
     n_dx = bound(i)%bnx(j) ! boundary node nx of normal vector
     n_dy = bound(i)%bny(j) ! boundary node ny of normal vector
     nid  = node(bound(i)%bnode(j))%nnghbrs ! current neighbour count for boundary node
     
     drav=0.0d0
     do v=1,nid
     inn = node(bound(i)%bnode(j))%nghbr(v) ! neighbour node of boundary node
     x1 = node(node(bound(i)%bnode(j))%nghbr(v))%x ; y1 = node(node(bound(i)%bnode(j))%nghbr(v))%y
     x0 = node(bound(i)%bnode(j))%x                ; y0 = node(bound(i)%bnode(j))%y
     dr = sqrt( (x1-x0)**2 + (y1-y0)**2 ) ! local edge dr
     drav = drav + dr
     enddo
     drav = drav/dble(nid) ! local edge radius
     
!    [1] average radius is used to determine the length of the ghost cell connection
!    [2] direction of the ghost cell placement is determined from the node unit vector
     
!    (1) Add gc node to the neighbor list of boundary node, node(bound(i)%bnode(j))
     node(bound(i)%bnode(j))%nnghbrs = node(bound(i)%bnode(j))%nnghbrs + 1
     call my_alloc_int_ptr(node(bound(i)%bnode(j))%nghbr, node(bound(i)%bnode(j))%nnghbrs)
     node(bound(i)%bnode(j))%nghbr(node(bound(i)%bnode(j))%nnghbrs) = ii
   
!    (2) define its node x,y location
     node(ii)%x = node(bound(i)%bnode(j))%x + n_dx*drav
     node(ii)%y = node(bound(i)%bnode(j))%y + n_dy*drav
     
!    (3) define local ghost cell ptype from boundary type
     if (trim(bound(i)%bc_type) == "freestream") then
     node(ii)%ptype = 4
     elseif(trim(bound(i)%bc_type) == "slip_wall_weak")then
     node(ii)%ptype = 5  
     endif
     
!   (4) Add corresponding boundary node to GC's neighbours ( for boundary treatment )
!       note : Only required for wall ghost cells but generated for all

     node(ii)%nnghbrs = node(ii)%nnghbrs + 1
     call my_alloc_int_ptr(node(ii)%nghbr,node(ii)%nnghbrs)
     node(ii)%nghbr(node(ii)%nnghbrs) = bound(i)%bnode(j)   ! add the boundary node to GC's neighbours
     
    enddo
    enddo
    
 !  set the number of ghost IP for future use
    ngc = iii
    
!   [check] only
!    do i=1,nnodes+ngc
!    iii=0
!    do j=1,node(i)%nnghbrs
!    if( node(node(i)%nghbr(j))%ptype .eq. 4 .and. node(i)%ptype .eq. 2 )then
!!    if( node(node(i)%nghbr(j))%ptype .eq. 5 .and. node(i)%ptype .eq. 1 )then
!    iii=iii+1
!    endif
!    enddo
!    if( iii .gt. 1 ) write(*,*) 'found'
!    enddo
    
!    i=1
!    write(*,*) node(i)%ptype
!    write(*,*) node(i)%nnghbrs
!    write(*,*) node(i)%nghbr(:)
!    write(*,*) node(node(i)%nghbr(:))%ptype
    
!   [check] correct allocation in python  #################
    if( confirmID .eqv. .true. )then
    
    open(44,file='./info/xyGCall.dat',form='formatted')
    do i=nnodes+1,nnodes+ngc
    write(linebuf,*) node(i)%x,',',node(i)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
!    write(44,*) node(i)%x,node(i)%y
    enddo
    close(44)

    open(44,file='./info/xyGCFF.dat',form='formatted')
    do i=1,nnodes+ngc
    if( node(i)%ptype .eq. 4 )then
    write(linebuf,*) node(i)%x,',',node(i)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
!    write(44,*) node(i)%x,node(i)%y
    endif
    enddo
    close(44)
    
    open(44,file='./info/xyGCW.dat',form='formatted')
    do i=1,nnodes+ngc
    if( node(i)%ptype .eq. 5 )then
    write(linebuf,*) node(i)%x,',',node(i)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
!    write(44,*) node(i)%x,node(i)%y
    endif
    enddo
    close(44)
    
    endif
!   ########################################################

!!   [check] if wall gc's neighbour is boundary node and has only 1
!    do i=nnodes+1,nnodes+ngc
!    if( node(i)%ptype .eq. 4 )then     
!     if( node(i)%nnghbrs .ne. 2 ) pause 'gc error'
!    do j=1,node(i)%nnghbrs
!!    write(*,*) node(node(i)%nghbr(j))%ptype
!    if( node(node(i)%nghbr(j))%ptype .ne. 1 )then
!    pause 'gc error2'
!    write(*,*) node(node(i)%nghbr(j))%ptype
!    endif
!    enddo
!    
!    endif
!    enddo
    

    return
    end subroutine bcgc
    
    
     
    