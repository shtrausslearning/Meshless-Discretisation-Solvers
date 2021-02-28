
!   Subroutine needed to know how many ghost cells and IPs are to be allocated
!   The reason this is required is that array node() has dependencies so it cannot
!   be reallocated without deleting all of its dependenices ( as far as Im aware )

    subroutine preRgrid

    use edu2d_my_main_data
    use prmflow

    implicit none

    integer  :: i, os, ii,iii
  
    !  Open the input file.
    open(unit=1, file='./grid/'//gridInF, status="unknown", iostat=os)

    ! READ: Get the size of the grid.
    read(1,*) nnodes, ntria, nquad
    nelms = ntria + nquad   ! total number of elements

    do i = 1, nnodes
    read(1,*) 
    end do

    if ( ntria > 0 ) then
    do i = 1, ntria
    read(1,*) 
    end do
    endif

    if ( nquad > 0 ) then
    do i = 1, nquad
    read(1,*) 
    end do
    endif

    read(1,*) nbound

    iii=0
    do i = 1, nbound
    read(1,*) ii
    iii=iii+ii
    end do

    close(1)
  
    tbcn = iii
    nIPGC = tbcn - nbound

    return
    end subroutine preRgrid