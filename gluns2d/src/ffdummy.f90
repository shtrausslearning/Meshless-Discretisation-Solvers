
    subroutine DummyNodes
    use prmflow
    use ModInterfaces, only : EM
    implicit none

!   local variables
    character(chrlen) :: msg
    logical :: flag
    integer :: ierr, ibegf, iendf, ibegn, iendn, itype
    integer :: i, ib, ibf, idn
    integer, allocatable :: marker(:)  ! node marker
! ###################################################################

! temporary array
    allocate( marker(phynod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate temporary marker" )

!   loop over boundary segments

    ibegf = 1
    ibegn = 1
    idn   = 0  ! counter of dummy nodes

    do ib=1,nsegs
    iendf = ibound(1,ib)
    iendn = ibound(2,ib)
    itype = btype(ib)
    flag  = .false.  ! true for inlet/outlet/far-field
    
    if (itype>=100 .and. itype<200) flag = .true. ! inlet nodes
    if (itype>=200 .and. itype<300) flag = .true. ! outlet nodes
    if (itype>=600 .and. itype<700) flag = .true. ! farfield nodes
    
    if (itype<700)then                            ! NOT a periodic boundary

      marker(:) = -777 ! reset

! --- loop over faces of boundary "ib" and mark nodes
      do ibf=ibegf,iendf
      marker(bface(1,ibf)) = 1 ! left node of bc column
      marker(bface(2,ibf)) = 1 ! right node of bc column
      enddo

! --- store node indexes in "bnode";
!     count dummy nodes (idn) for inlet/outlet/far-field boundary
      
      do i=1,phynod                      ! cycle through all nodes
      if (marker(i) == 1) then           ! if bc node of current itype
      
!      check only
       if (ibegn > nbnodes) then        ! check dimension
       call EM( "max. no. of boundary nodes exceeded" )
       endif
       
       if(flag)then                   ! *** inlet/outlet/far-field
        idn            = idn + 1
        bnode(1,ibegn) = i             ! actual bc node
        bnode(2,ibegn) = phynod + idn  ! new dummy node
        bnode(3,ibegn) = -1            ! set in "EdgesFinalize"
       else                             ! *** other boundary type
        bnode(1,ibegn) = i             ! index of boundary node
       endif
       
        ibegn = ibegn + 1                ! count boundary nodes
        
      endif
      enddo

! --- check no. of boundary nodes
      if (ibegn-1 .ne. iendn) then
      write(msg,1000) itype,iendn,ibegn-1
      call EM( msg )
      endif

    endif  ! periodic?

! - reset pointers to faces and to nodes
    ibegf = iendf + 1
    ibegn = iendn + 1

    enddo  ! ib

    deallocate( marker )

!   set total number of nodes (add dummy nodes)
    allnod = phynod + idn

1000 format("no. of nodes for boundary ",I3," is wrong. It should be ",I5, &
            " but it is ",I5)

    end subroutine DummyNodes
