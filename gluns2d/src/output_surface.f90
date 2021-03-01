
  subroutine PlotSurfaces
  use ModDataTypes
  use prmflow
  use ModInterfaces, only : EM
  implicit none

! local variables
  character(chrlen) :: fname
  integer     :: errFlag, itype, ibegf, iendf, ibegn, iendn, ibf1, ibf2, &
                 nquant, nsurfs,ierr
  integer     :: i, ib, ibf, ibn, m
  real(rtype) :: rrho, u, v, e, press, temp, c, ptot, ttot, mach, machis, &
                 ptloss, pratio, ptotinf, gam1, ggm1
  real(rtype) :: cf, cp, visc, sx, sy, ds, sxn, syn, grdnx, grdny, grdnn, &
                 dvdnx, dvdny, dvdna, sgn

! *****************************************************************************

  write(fnSurf,"(A,A,A)") 'surfp_',Trim(title),".dat"
  open(44, file=fnSurf, status="unknown", action="write", iostat=ierr)
  if (ierr /= 0) call EM( "cannot open plot file Surface.dat" )

  do i=1,phynod
  if( ntype(i) .eq. 1 )then ! if wall node
  if( i .le. cpcoff )then

    rrho      = 1.D0/cv(1,i)
    u          = cv(2,i)*rrho
    v          = cv(3,i)*rrho
    e          = cv(4,i)*rrho
    press   = dv(1,i)
    c           = dv(3,i)
    cp         = 2.0d0*(pinf-press)/(rhoinf*qinf*qinf)
    write(44,1006) x(i),cp

  endif  
  endif
  enddo 

  close(44)

1020  format(20E16.8,20E16.8)
1006  format(F0.5,1X,F0.5)

end subroutine PlotSurfaces

