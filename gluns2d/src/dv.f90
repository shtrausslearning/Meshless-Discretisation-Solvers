
!!
    subroutine dvall
    use ModDataTypes
    use prmflow
    implicit none

!   local variables
    integer     :: i
    real(rtype) :: gam1, rgas, g1cp, rhoq, s1, s2, s12, rat, cppr
!   ###################################################################

    gam1 = gamma - 1.D0
    rgas = gam1*cpgas/gamma
    g1cp = gam1*cpgas

!   Euler equations
    do i=1,phynod
    rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
    dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
    dv(2,i) = dv(1,i)/(rgas*cv(1,i))
    dv(3,i) = Sqrt(g1cp*dv(2,i))
    dv(4,i) = gamma
    dv(5,i) = cpgas
    enddo

    end subroutine dvall

    
subroutine DependentVarsOne( i )

  use ModDataTypes
  use prmflow
  implicit none

! parameters
  integer, intent(in) :: i

! local variables
  real(rtype) :: gam1, rgas, g1cp, rhoq, s1, s2, s12, rat

! ###################################################################

  gam1 = gamma - 1.D0
  rgas = gam1*cpgas/gamma
  g1cp = gam1*cpgas

! Euler equations

    rhoq    = cv(2,i)*cv(2,i) + cv(3,i)*cv(3,i)
    dv(1,i) = gam1*(cv(4,i)-0.5D0*rhoq/cv(1,i))
    dv(2,i) = dv(1,i)/(rgas*cv(1,i))
    dv(3,i) = Sqrt(g1cp*dv(2,i))
    dv(4,i) = gamma
    dv(5,i) = cpgas

end subroutine DependentVarsOne
