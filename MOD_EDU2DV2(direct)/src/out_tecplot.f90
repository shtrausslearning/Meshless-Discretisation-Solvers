&
    
!   WRITE TECPLOT FILE
    subroutine write_tecplot_file
    use prmflow
    use edu2d_my_main_data, only : nnodes, node, elm, nelms
    implicit none
    integer :: i, k, os, ierr
    real(p2) :: lu,lv,lc,lmach,lpress,ggm1
    real(p2), allocatable :: lvars(:,:) ! local variable storage array

    open(unit=1, file='./output/'//datafile_tec, status="unknown", iostat=os)
    write(1,*) 'title = "grid"'
    write(1,'(a100)') 'variables = "x","y","rho","u","v","mach","Cp","Total Enthalpy","Entropy"'
    write(1,*) 'zone n=',nnodes,' e =', nelms,' et=quadrilateral, f=fepoint'

    ierr=0;allocate( lvars(7,nnodes),stat=ierr )
    if( ierr /= 0 ) pause 'allocation error lvars'
    
    do i=1,nnodes
    lu = node(i)%cv(2)/node(i)%cv(1) 
    lv = node(i)%cv(3)/node(i)%cv(1)
    lvars(1,i) = node(i)%cv(1)/rhoinf
    lvars(2,i) = (node(i)%cv(2)/node(i)%cv(1))/uinf
    lvars(3,i) = (node(i)%cv(3)/node(i)%cv(1))/vinf
    enddo

!   additional output variables
!   1. mach number
    do i=1,nnodes
    lu = node(i)%cv(2)/node(i)%cv(1) 
    lv = node(i)%cv(3)/node(i)%cv(1)
    lc = node(i)%dv(3)
    lvars(4,i) = sqrt(lu*lu+lv*lv)/lc
    enddo
    
!   2. cp (nondim)
    do i=1,nnodes
    lpress = node(i)%dv(1)
    lvars(5,i) = 2.0d0*(pinf-lpress)/(rhoinf*qinf*qinf)
    enddo
    
!   3. Total Enthalpy, H
    do i=1,nnodes
    ggm1 = gamma/(gamma-1.0d0)
    lu = node(i)%cv(2)/node(i)%cv(1) 
    lv = node(i)%cv(3)/node(i)%cv(1)
    lvars(6,i) = ggm1*node(i)%dv(1)/node(i)%cv(1) + 0.5d0*( lu*lu + lv*lv )
    enddo
    
!   4. entropy ( -1 form ) Kuya said to use rho?
    do i=1,nnodes
    lvars(7,i) = (node(i)%dv(1)/pinf)/( (node(i)%cv(1)/rhoinf)**gamma ) -1.0d0
    enddo
    
!   Store all data now
    do i = 1, nnodes
    write(1,*) node(i)%x,node(i)%y,(lvars(k,i),k=1,7)
    end do

!   write geometric stuff
    do i = 1, nelms
     if(elm(i)%nvtx == 3)then
      write(1,*) elm(i)%vtx(1), elm(i)%vtx(2), elm(i)%vtx(3), elm(i)%vtx(3)
     elseif (elm(i)%nvtx == 4) then
      write(1,*) elm(i)%vtx(1), elm(i)%vtx(2), elm(i)%vtx(3), elm(i)%vtx(4)
     else
      write(101,*) " Error in elm%vtx data... Stop..: elm(i)%nvtx=",elm(i)%nvtx
      stop
     endif
    enddo

    close(1)
    deallocate(lvars)
    
    end subroutine write_tecplot_file
