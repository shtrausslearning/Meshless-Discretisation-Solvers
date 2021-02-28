
!   subroutine to allocate physical ptype to node(i)%ptype
    subroutine ptypePHYS
    use edu2d_constants, only : p2
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    use edu2d_my_allocation
    implicit none
    
    logical, parameter :: confirmID = .true.  ! show confirmation 
    integer :: i,j
    character(256) :: linebuf
 
 !  Ptype Summary
 !  ---------------
!   [1] physical wall [2] outer boundary [3] fluid in between [4] wall gc [5] farfield gc
 
!   (1) set physical ptype for all nodes in domain to fluid ptype=3

    NODE1: do i=1,nnodes
    node(i)%ptype=3
    enddo NODE1
    
!   (2) Identify physical boundary nodes ( only ff/wall for now )
!       Wall takes over ff if corner node
    BC_LOOP: do i = 1, nbound  ! cycle though wall boundaries

     if (trim(bound(i)%bc_type) == "freestream") then
      do j = 1, bound(i)%nbfaces
       node(bound(i)%bnode(j  ))%ptype = 2  !Left node
       node(bound(i)%bnode(j+1))%ptype = 2  !Right node
      enddo
     elseif(trim(bound(i)%bc_type) == "slip_wall_weak")then
      do j = 1, bound(i)%nbfaces
       node(bound(i)%bnode(j  ))%ptype = 1  !Left node
       node(bound(i)%bnode(j+1))%ptype = 1  !Right node
      enddo
     endif
     
     enddo BC_LOOP
     
!   check correct allocation in python  #################
    if( confirmID .eqv. .true. )then

    open(44,file='./info/xyPall.dat',form='formatted')
    do i=1,nnodes
    write(linebuf,*) node(i)%x,',',node(i)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
!    write(44,*) node(i)%x,node(i)%y
    enddo
    close(44)
    
    open(44,file='./info/xyPwall.dat',form='formatted')
    do i=1,nnodes
    if( node(i)%ptype .eq. 1 )then
    write(linebuf,*) node(i)%x,',',node(i)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
!    write(44,*) node(i)%x,node(i)%y
    endif
    enddo
    close(44)
    
    open(44,file='./info/xyPff.dat',form='formatted')
    do i=1,nnodes
    if( node(i)%ptype .eq. 2 )then
    write(linebuf,*) node(i)%x,',',node(i)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
!    write(44,*) node(i)%x,node(i)%y
    endif
    enddo
    close(44)
    
    open(44,file='./info/xyPIBW.dat',form='formatted')
    do i=1,nnodes
    if( node(i)%ptype .eq. 3 )then
    write(linebuf,*) node(i)%x,',',node(i)%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
!    write(44,*) node(i)%x,node(i)%y
    endif
    enddo
    close(44)
    
    endif
!   ########################################################
    
    return
    end subroutine ptypePHYS