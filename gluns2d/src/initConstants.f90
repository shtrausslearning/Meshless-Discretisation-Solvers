
    subroutine InitConstants
    use ModDataTypes
    use prmflow
    implicit none

!   local variables
    real(rtype) :: gam1, rgas

    pi  = 4.D0*Atan(1.D0)
    rad = 180.D0/pi
   
    gam1 = gamma - 1.D0
    rgas = gam1*cpgas/gamma

    alpha   = alpha/rad
    rhoinf  = pinf/(rgas*tinf)
    qinf    = machinf * Sqrt(gam1*cpgas*tinf)
    uinf    = qinf*Cos(alpha)
    vinf    = qinf*Sin(alpha)
    refrho  = rhoinf
    refvel  = qinf

    return
    end subroutine InitConstants
