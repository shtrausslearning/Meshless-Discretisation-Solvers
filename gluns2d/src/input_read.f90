
    subroutine ReadParams
    use prmflow
    use ModInterfaces, only : EM, ReadChar
    implicit none

!!   parameters
!    character(*), intent(in) :: fname

!   local variables
    character(1) :: ch
    integer :: ierr, i

! ###################################################################

    !open(10, file=fname, status="old", action="read", iostat=ierr)
    !if (ierr /= 0) call EM( "cannot open input file" )
    
    open(10,file='input.dat',form='formatted')
    read(10,"(A)") title
    read(10,"(A)") fnGrid
    read(10,"(A)") solin
    read(10,"(A)") solout
    read(10,*) lrest
    read(10,*) machinf
    read(10,*) alpha
    read(10,*) xref
    read(10,*) yref
    read(10,*) cref
    read(10,*) cpcoff
    read(10,*) maxiter
    read(10,*) outstep
    read(10,*) vtkout
    read(10,*) convtol
    read(10,*) cfl
    read(10,*) ktimst
    read(10,*) iorder
    read(10,*) limfac
    read(10,*) iflux
    read(10,*) epsirs
    read(10,*) nitirs
    read(10,*) lvort
    read(10,*) nrk
    read(10,*) (ark  (i), i=1,nrk)
    close(10)
  
    soltype = "E" ! external flow
    gamma = 1.4d0
    cpgas = 1004.5d0 ! specific heat coeff. at constant pressure [J/kgK]
    renum = 5000.0d0
    iflow = "E"
    prlam = 0.72d0
    pinf = 1.01E5   ! static pressure at infinity [Pa]
    tinf = 288.0d0 ! static temperature at infinity [K]
    epsentr = 0.000001d0 ! entropy correction coefficient
!   smaller values move shock backward
!   too high reduces Cp significantly

    return
    end subroutine ReadParams


    function ReadChar( iunit )
    use ModDataTypes
    implicit none
    
    integer, intent(in) :: iunit
    character(1) :: ReadChar
    character(chrlen) :: str
    integer :: i
! ###################################################################

    read(iunit,"(A)") str

    do i=1,Len_trim(str)
    ReadChar = str(i:i)
    if (ReadChar /= " ") exit  ! first non-empty char should be the option ...
    enddo

    end function ReadChar
