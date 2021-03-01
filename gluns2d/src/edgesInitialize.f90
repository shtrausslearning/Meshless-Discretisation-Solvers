
  subroutine EdgesInitialize( niedge,iedge )

    use prmflow
    use ModInterfaces, only : EM
    implicit none

    ! parameters
    integer, intent(out) :: niedge(:), iedge(:,:)

    ! local variables
    integer :: d, cedge, cedge2, mxedges
    integer :: i, j, ic, ie, n

!   ###################################################################

    mxedges = ubound(iedge,2)  ! max. possible number of edges

!   reset all pointers

    do i=1,phynod
    niedge(i) = -777
    enddo
    do ie=1,mxedges
    iedge(1,ie) = -777
    iedge(2,ie) = -777
    iedge(3,ie) = -777
    enddo

! loop over nodes of all triangles
  nedint = 0

  do n=1,3

! - loop over triangles

    do ic=1,nt

      i = elem(n,ic)
      if (n < 3) then
        j = elem(n+1,ic)
      else
        j = elem(1,ic)
      endif
      if (i > j) then  ! lower index first
        d = i
        i = j
        j = d
      endif

      if (niedge(i) < 0) then

! ----- define new edge
        nedint = nedint + 1
        if (nedint > mxedges) then
          call EM( "max. number of edges reached" )
        endif
        niedge(i)       = nedint
        iedge(1,nedint) = j

      else

! ----- insert node "j" into list of adjacent nodes
        cedge = niedge(i)
10      continue
          if (iedge(1,cedge) == j) goto 20
          cedge2 = iedge(2,cedge)
          if (cedge2 < 0) then
            nedint = nedint + 1
            if (nedint > mxedges) then
              call EM( "max. number of edges reached" )
            endif
            iedge(2,cedge ) = nedint
            iedge(1,nedint) = j
            goto 20
          endif
          cedge = cedge2
        goto 10
20      continue
      endif

    enddo ! loop over triangles

  enddo   ! loop over nodes of triangles

! set total no. of edges (add edges to dummy nodes)
  nedges = nedint + (allnod-phynod)

  end subroutine EdgesInitialize

  
    subroutine EdgesFinalize( niedge,iedge )

    use prmflow
    use ModInterfaces, only : EM
    implicit none

!   parameters
    integer :: niedge(:), iedge(:,:)

!   local variables
    integer :: ierr, i, ibn, ie, cedge
    integer :: ii,j,req,ic, inn
!   ###################################################################

    allocate( edge(2,nedges),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for edge()" )

    ie = 0  ! edge counter

    do i=1,phynod

! - loop over all grid nodes
    cedge = niedge(i)
    if (cedge > 0) then
10    continue

! ----- loop over all edges adjacent to node "i"
        ie             = ie + 1
        edge(1,ie)     = i
        edge(2,ie)     = iedge(1,cedge)
        iedge(3,cedge) = ie                 ! we need it in InitMetrics
        cedge          = iedge(2,cedge)     ! next adjacent edge
        if (cedge < 0) goto 20
        
      goto 10
20    continue
    endif

    enddo

    if (ie /= nedint) then
    call EM( "did not get the correct number of interior edges" )
    endif
   
!   add edges to dummy nodes;
!   store 'dummy' edges in "bnode(3,*)"
    do ibn=1,nbnodes
    if (bnode(3,ibn) == -1) then      ! dummy node here (see DummyNodes)
    ie           = ie + 1
    edge(1,ie)   = bnode(1,ibn)     ! boundary node first
    edge(2,ie)   = bnode(2,ibn)     ! dummy node second
    bnode(3,ibn) = ie
    endif
    enddo

!   check
    if (ie /= nedges) then
    call EM( "did not get the correct number of dummy edges" )
    endif
    
    write(*,*) ''
    write(*,*) ' Grid Information for ',trim(fnGrid)
    write(*,*) '######################################################'
    write(*,1000) phynod,allnod-phynod,nt,nedint,nedges,nbfaces,nbnodes
    write(*,*) '######################################################'

1000 format(" No. of interior nodes: ",I8,/, &
            " No. of dummy nodes   : ",I8,/, &
            " No. of grid cells    : ",I8,/, &
            " No. of interior edges: ",I8,/, &
            " Total number of edges: ",I8,/, &
            " No. of boundary faces: ",I8,/, &
            " No. of boundary nodes: ",I8)
    
    return
    end subroutine EdgesFinalize