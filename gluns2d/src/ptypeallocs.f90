
    Subroutine ptypealloc
    use ModDataTypes
    use prmflow
    
!   point type
    allocate( ptype(allnod+igcs),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate ntype " )
  
    do i=1,phynod
    ptype(i) = 1            ! all physical nodes
    enddo                   ! [1] wall boundary [2] outer boundary [3] fluid nodes
    do i=phynod+1,allnod+igcs
    ptype(i) = 2            ! all ghost cell nodes
    enddo                   ! [5] farfield dummy [4] wall dummy
    
    return
    end subroutine ptypealloc