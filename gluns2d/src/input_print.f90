
    subroutine PrintParams
    use prmflow
    implicit none

!   local variables
    integer :: i
  
    write(*,*) ' Solution Information : '
    write(*,*) '##################################################################'
    
    ! physics - general
    write(*,1040) soltype,"E=external flow, I=internal flow"
    !write(*,1040) iflow,"E=Euler, N=Navier-Stokes"
    !write(*,1045) gamma,"ratio of specific heats"
    !write(*,1045) cpgas,"specific heat coefficient (p=const.)"
    write(*,1045) renum,"Reynolds number"
    !write(*,1045) refvel,"reference velocity"
    !write(*,1045) refrho,"reference density"
    !write(*,1045) refvisc,"laminar viscosity"
    !write(*,1045) prlam,"laminar Prandtl number"

    write(*,1045) machinf,"Mach-number at infinity"
    write(*,1045) alpha*rad,"angle of attack [deg]"
    write(*,1045) pinf,"static pressure at infinity [Pa]"
    write(*,1045) tinf,"static temperature at infinity [K]"
    write(*,1050) maxiter,"max. number of iterations"
    write(*,1050) outstep,"number of iterations between solution dumps"
    write(*,1045) convtol,"2-norm of density change to stop the iteration"
    write(*,1055) lrest,"use previous solution for restart (Y=1, N=0)"
    write(*,1045) cfl,"CFL-number"
    write(*,1055) ktimst,"1=local, 0=global time-stepping"
    write(*,1055) iorder,"1st-order (1) / 2nd-order (2) Roe scheme"
    write(*,1045) epsentr,"entropy correction coefficient (Roe scheme only)"
    write(*,1045) epsirs, "IRS Coefficient (0:Not Used)"
    write(*,1050) nitirs, "IRS Iterations if used "
    write(*,1055) lvort,"correction of far-field due to single vortex (external flow)"
    write(*,1055) nrk,"number of Runge-Kutta stages (steady flow only)"
    write(*,1060) (ark  (i), i=1,nrk)
    write(*,*) '##################################################################'

1000  format(/,A,/)
1005  format("#",/,"# Physics - general",/,"# ",17("-"))
1010  format("#",/,"# Physics - external flow",/,"# ",23("-"))
1020  format("#",/,"# Geometrical reference values",/,"# ",28("-"))
1025  format("#",/,"# Iteration control",/,"# ",17("-"))
1030  format("#",/,"# Numerical parameters",/,"# ",20("-"))
1035  format("#",/,"# Quantities to plot",/,"# ",18("-"))
1040  format(2X,A1,12X,"# ",A)
1041  format(2X,A1,2X,"# ",A)
1045  format(1X,1PE11.4,3X,"# ",A)
1050  format(2X,I6,7X,"# ",A)
1055  format(2X,I1,12X,"# ",A)
1060  format(5(2X,F6.4))
1065  format(5(2X,I6))

    end subroutine PrintParams
