
    subroutine plotsurfaces
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    implicit none
    
    character(80) :: fname
    character(80) :: dname
    integer :: i
    real(p2) :: hl,lu,lv,lmach,lc,lent,lcp
    
!   Entropy Output
    open(45,file='./output/surfentropy.dat',form='formatted')
    
    do i=1,nnodes
    if( node(i)%ptype .eq. 1 )then ! if wall node
    lent = (node(i)%dv(1)/pinf)/( (node(i)%cv(1)/rhoinf)**gamma ) - 1.0d0
    write(dname,*) node(i)%x,',',lent
    call del_spaces(dname)
    write(45,*) trim(dname)

    endif
    enddo
    close(45)
    
!   Cp output 
    open(45,file='./output/surfCp.dat',form='formatted')

    do i=1,nnodes
    if( node(i)%ptype .eq. 1 )then ! if wall node
    lcp = 2.0d0*(pinf-node(i)%dv(1))/(rhoinf*qinf*qinf)
    write(dname,*) node(i)%x,',',lcp
    call del_spaces(dname)
    write(45,*) trim(dname)
    endif
    enddo
    close(45)

!   cv(1) output 
    open(45,file='./output/surfcv1.dat',form='formatted')

    do i=1,nnodes
    if( node(i)%ptype .eq. 1 )then ! if wall node
    write(dname,*) node(i)%x,',',node(i)%cv(1)
    call del_spaces(dname)
    write(45,*) trim(dname)
    endif
    enddo
    close(45)
    
!   cv(2) output 
    open(45,file='./output/surfcv2.dat',form='formatted')

    do i=1,nnodes
    if( node(i)%ptype .eq. 1 )then ! if wall node
    write(dname,*) node(i)%x,',',node(i)%cv(2)
    call del_spaces(dname)
    write(45,*) trim(dname)
    endif
    enddo
    close(45)
    
!   cv(3) output 
    open(45,file='./output/surfcv3.dat',form='formatted')

    do i=1,nnodes
    if( node(i)%ptype .eq. 1 )then ! if wall node
    write(dname,*) node(i)%x,',',node(i)%cv(3)
    call del_spaces(dname)
    write(45,*) trim(dname)
    endif
    enddo
    close(45)
    
!   cv(4) output 
    open(45,file='./output/surfcv4.dat',form='formatted')

    do i=1,nnodes
    if( node(i)%ptype .eq. 1 )then ! if wall node
    write(dname,*) node(i)%x,',',node(i)%cv(4)
    call del_spaces(dname)
    write(45,*) trim(dname)
    endif
    enddo
    close(45)
    
!   enthalpy output 
    open(45,file='./output/surfenthalpy.dat',form='formatted')

    do i=1,nnodes
    if( node(i)%ptype .eq. 1 )then ! if wall node
    lu = node(i)%cv(2)/node(i)%cv(1)
    lv = node(i)%cv(3)/node(i)%cv(1)
    hl = ((gamma-1.0d0)/gamma)*(node(i)%dv(1)/node(i)%cv(1)) + 0.5d0*( lu**2 + lv**2 )
    write(dname,*) node(i)%x,',',hl
    call del_spaces(dname)
    write(45,*) trim(dname)
    endif
    enddo
    close(45)
    
!   enthalpy output 
    open(45,file='./output/surfmach.dat',form='formatted')

    do i=1,nnodes
    if( node(i)%ptype .eq. 1 )then ! if wall node
    lu = node(i)%cv(2)/node(i)%cv(1)
    lv = node(i)%cv(3)/node(i)%cv(1)
    lc = node(i)%dv(3)
    lmach = sqrt( lu*lu + lv*lv )/lc
    write(dname,*) node(i)%x,',',lmach
    call del_spaces(dname)
    write(45,*) trim(dname)
    endif
    enddo
    close(45)
    
    return
    end subroutine plotsurfaces