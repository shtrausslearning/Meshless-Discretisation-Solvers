
!> Monitors the convergence, prints it out and stores it in a file. For external
!! flow, it also prints out the lift, the drag and the moment coefficients. For
!! internal flow, it prints out the mass flow and the mass flow ratio.
!!
    subroutine Convergence
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none

!   local variables
    integer     :: i, idr
    real(rtype) :: dr, drmax

! ###################################################################
! compute the residual

    drho  = 0.D0
    drmax = 0.D0
    do i=1,phynod
  
    dr   = cv(1,i) - cvold(1,i)
    drho = drho + dr*dr
    if (Abs(dr) >= drmax) then
      drmax = Abs(dr)
      idr   = i
    endif
    
    enddo

    if (iter == 1) then
    drho1 = Sqrt(drho) + 1.D-32
    drho  = 1.D0
    else
    drho  = Sqrt(drho)/drho1
    endif

    call Forces

    if(mod(iter,10)==1)then
    write(50,1000) iter,Log10(drho),drmax,idr,cl,cd
    write(51,1006) iter,log10(drho)
    write(52,1006) iter,cl
    write(53,1006) iter,cd
    write(*,1000) iter,Log10(drho),drmax,idr,cl,cd
    endif
    
    if( isnan(Log10(drho)) )then
    write(*,*) 'drho erro'
    pause
    stop
    endif

1000  format(I7,1P,2X,E12.5,2X,E12.5,0P,I8,1P,2(2X,E12.5))
1005  format(I7,1P,2X,E12.4,2X,E10.4,0P)
1006  format(I0,1P,2X,F0.4)

    end subroutine Convergence