
subroutine BoundaryConditions2( work )

  use ModDataTypes
  use prmflow
  use ModInterfaces, only : EM
  implicit none

! parameters
  real(rtype) :: work(:)
  integer :: ib, ibegn, iendn, itype, wdim

  ibegn = 1

  do ib=1,nsegs

    itype = btype(ib)
    iendn = ibound(2,ib)

    if (itype>=600 .and. itype<700) then

      wdim = Ubound(work,1)
      if ((4*nbnodes) > wdim) then
        call EM( "insufficient work space in BoundaryConditions" )
      endif
      call BcondFarfield2( ibegn,iendn, &
                          work( 1           :  nbnodes), &
                          work((1+  nbnodes):2*nbnodes), &
                          work((1+2*nbnodes):3*nbnodes), &
                          work((1+3*nbnodes):4*nbnodes) )
    endif

    ibegn = iendn + 1

  enddo ! ib
  
  call bc_wgc

end subroutine BoundaryConditions2
