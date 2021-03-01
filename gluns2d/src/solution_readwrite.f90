
  subroutine ReadSolution
! Binary Format
  use prmflow
  use ModInterfaces, only : EM
  implicit none

! local variables
  integer :: ierr, i, n, chk1, chk2

! ###################################################################

  open(unit=60, file=solin, status="old", action="read", &
       form="unformatted", iostat=ierr)
  if (ierr /= 0) call EM( "cannot open solution file" )

! Check
  read(60) chk1,chk2
  if (chk1 /= allnod)  &
  call EM( "no. of nodes differs from the grid file" )
  if (chk2 /= nconv)  &
  call EM( "different number of conservative variables" )
! read
  read(60) drho1,iter ! start residual, iteration
  read(60) ((cv(n,i), i=1,allnod), n=1,nconv)
  close(60)

  end subroutine ReadSolution

  subroutine WriteSolution

  use prmflow
  use ModInterfaces, only : EM
  implicit none

! local variables
  integer :: ierr, i, n

! ###################################################################

  open(unit=70, file=solout, status="unknown", action="write", &
       form="unformatted", iostat=ierr)
  if (ierr /= 0) call EM( "cannot open solution file" )

! dimensions
  write(70) allnod,nconv
  write(70) drho1,iter ! initial residual, iteration # and solution
  write(70) ((cv(n,i), i=1,allnod), n=1,nconv)
  close(70)

end subroutine WriteSolution
