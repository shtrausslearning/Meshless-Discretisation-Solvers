

    subroutine alloc
    use ModDataTypes
    use prmflow
    use ModInterfaces, only : EM
    implicit none

!   local variables
    integer :: ierr,totalnod

!   ###################################################################
  
    totalnod = allnod+igcs
  
!   base flow variables
    ierr=0;allocate( cv(nconv,totalnod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for cv()" )
    ierr=0;allocate( dv(ndepv,totalnod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for dv()" )

!   general numerical variables
    ierr=0;allocate( cvold(nconv,totalnod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for cvold()" )
    ierr=0;allocate( rhs(nconv,totalnod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for rhs()" )
    ierr=0;allocate( tstep(totalnod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for tstep()" )
    tstep=0;rhs=0.0;cvold=0.0;cv=0.0;dv=0.0
    
    ierr=0;allocate( gradx(4,totalnod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for gradx()" )
    ierr=0;allocate( grady(4,totalnod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for grady()" )
    gradx=0;grady=0
    
    ierr=0;allocate( lim(4,totalnod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for lim()" )
    lim=0
  
    end subroutine alloc
