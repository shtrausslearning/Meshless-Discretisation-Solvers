
  subroutine output_vtk
! ##########################################################################################
! WRITING VTK V2.0 FILE FORMAT 
! ##########################################################################################
  use modDataTypes
  use prmflow
  use modInterfaces
  implicit none

  character(chrlen) :: fname
  integer           :: ierr, nquant, i, m, CELLTYPE
  real(rtype)       :: rrho, u, v, e, press, temp, c, mach, ttot, ptot
  real(rtype)       :: ptloss, pratio, ptotinf, gam1, ggm1, visc, machis
  real(rtype)       :: varout(mxqfield+2)
! ######################################################################################

  ptotinf = 0.D0
  if (soltype == "E") then
    gam1    = gamma - 1.D0
    ggm1    = gamma/gam1
    ptotinf = pinf*(1.D0+0.5D0*gam1*machinf*machinf)**ggm1
  endif

  write(fname,"(A,I5.5,A)") Trim(title),iter,".vtk"
  
  open(30, file="./vtk/"//fname, status="unknown", action="write", iostat=ierr)
  if (ierr /= 0) call EM( "cannot open plot file (field)" )

  write(30,1006)
  write(30,1007)
  write(30,1008) 
  write(30,*)""
  write(30,1009)
  write(30,1004) phynod
  
  DO i=1,phynod
    varout(1) = x(i)
    varout(2) = y(i)
    varout(3) = 0.0d0
    WRITE(30,1011) (varout(m), m=1,3)
  enddo
  
  CELLTYPE = 4*nt
      
  write(30,*) ""
  write(30,1013) nt,CELLTYPE
      
  do i=1,nt
    write(30,1040) 3,elem(1,i)-1,elem(2,i)-1,elem(3,i)-1
  enddo
  
  write(30,*) ""
  write(30,1014) nt
  
  do i=1,nt
    write(30,1015) 5
  enddo 

   write(30,*) ''
   write(30,1012) phynod
      
!  LOCAL MACH NUMBER OUTPUT
!  ###########################################################################
    write(30,*) 'SCALARS LOCAL-MACH float'  
    write(30,*) 'LOOKUP_TABLE default'
      
    DO i=1,phynod
       
     rrho  = 1.D0/cv(1,i)
     u     = cv(2,i)*rrho
     v     = cv(3,i)*rrho
     c     = dv(3,i)
      
     mach  = SQRT(u*u+v*v)/c
       
     WRITE(30,1020) mach
    ENDDO  

    close(30)

1000  format(A,/,"1",/,"Flow Field",/,"1 ",I2,/,"x [m]",/,"y [m]")
1010  format("0 0"/,I6,I6," 0",/,"Unstructured")
1020  FORMAT(1P,20E14.6)
1040  FORMAT(1I2,3I7)

1006  FORMAT("# vtk DataFile Version 2.0" )
1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")
1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I7,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
1013  FORMAT('CELLS ',I7,I7)
1014  FORMAT('CELL_TYPES ',I7)
1015  FORMAT(I1)
1012  FORMAT('POINT_DATA ',I7)

! ######################################################################################
  end subroutine output_vtk
  
  subroutine output_vtk2(input)
! ##########################################################################################
! WRITING VTK V2.0 FILE FORMAT 
! ##########################################################################################
  use modDataTypes
  use prmflow
  use modInterfaces
  implicit none

  character(chrlen) :: fname
  integer           :: ierr, nquant, i, m, CELLTYPE
  real(rtype)       :: rrho, u, v, e, press, temp, c, mach, ttot, ptot
  real(rtype)       :: ptloss, pratio, ptotinf, gam1, ggm1, visc, machis
  real(rtype)       :: varout(mxqfield+2)
  integer, intent(in) :: input
! ######################################################################################


  write(fname,"(A,I5.5,A)") Trim(title),iter,".vtk"
  
  open(30, file="./vtk/"//fname, status="unknown", action="write", iostat=ierr)
  if (ierr /= 0) call EM( "cannot open plot file (field)" )

  write(30,1006)
  write(30,1007)
  write(30,1008) 
  write(30,*)""
  write(30,1009)
  write(30,1004) input
  
  DO i=1,allnod+igcs
  if( ntype(i) .eq. 4 )then ! if wall ghost cell
    varout(1) = x(i)
    varout(2) = y(i)
    varout(3) = 0.0d0
    WRITE(30,1011) (varout(m), m=1,3)
  endif  
  enddo
  
  close(30)

1000  format(A,/,"1",/,"Flow Field",/,"1 ",I2,/,"x [m]",/,"y [m]")
1010  format("0 0"/,I6,I6," 0",/,"Unstructured")
1020  FORMAT(1P,20E14.6)
1040  FORMAT(1I2,3I7)

1006  FORMAT("# vtk DataFile Version 2.0" )
1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")

1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I6,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
1013  FORMAT('CELLS ',I5,I7)
1014  FORMAT('CELL_TYPES ',I5)
1015  FORMAT(I1)
1012  FORMAT('POINT_DATA ',I5)

! ######################################################################################
  end subroutine output_vtk2
  
  subroutine output_vtk3
! ##########################################################################################
! WRITING VTK V2.0 FILE FORMAT 
! ##########################################################################################
  use modDataTypes
  use prmflow
  use modInterfaces
  implicit none

  character(chrlen) :: fname
  integer           :: ierr, nquant, i, m, CELLTYPE,id1,id2
  real(rtype)       :: rrho, u, v, e, press, temp, c, mach, ttot, ptot
  real(rtype)       :: ptloss, pratio, ptotinf, gam1, ggm1, visc, machis
  real(rtype)       :: varout(mxqfield+2)
! ######################################################################################

    open(30, file="ghostcells.vtk",form='formatted')

    write(30,1006)
    write(30,1007)
    write(30,1008) 
    write(30,*)""
    write(30,1009)
    !write(30,1004) igcs
    write(30,1004) 1
  
    DO i=1,igcs
  
    id1 = gcid(1,i) ! fluid node
    id2 = gcid(2,i) ! ghost cell node
    if( id2 .eq. 13913 )then

    varout(1) = x(id2)
    varout(2) = y(id2)
    varout(3) = 0.0d0
    WRITE(30,1011) (varout(m), m=1,3)
  
    endif
    enddo
  
    close(30)

1000  format(A,/,"1",/,"Flow Field",/,"1 ",I2,/,"x [m]",/,"y [m]")
1010  format("0 0"/,I6,I6," 0",/,"Unstructured")
1020  FORMAT(1P,20E14.6)
1040  FORMAT(1I2,3I7)

1006  FORMAT("# vtk DataFile Version 2.0" )
1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")

1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I6,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
1013  FORMAT('CELLS ',I5,I7)
1014  FORMAT('CELL_TYPES ',I5)
1015  FORMAT(I1)
1012  FORMAT('POINT_DATA ',I5)

! ######################################################################################
  end subroutine output_vtk3
  
  subroutine output_vtk4
  use modDataTypes
  use prmflow
  use modInterfaces
  implicit none

  character(chrlen) :: fname
  integer           :: ierr, nquant, i, m, CELLTYPE,id1,id2
  real(rtype)       :: rrho, u, v, e, press, temp, c, mach, ttot, ptot
  real(rtype)       :: ptloss, pratio, ptotinf, gam1, ggm1, visc, machis
  real(rtype)       :: varout(mxqfield+2)
! ######################################################################################

    open(30, file="fluidnode.vtk",form='formatted')

    write(30,1006)
    write(30,1007)
    write(30,1008) 
    write(30,*)""
    write(30,1009)
    write(30,1004) igcs
  
    DO i=1,igcs
  
    id1 = gcid(1,i) ! fluid node
    id2 = gcid(2,i) ! ghost cell node

    varout(1) = x(id1)
    varout(2) = y(id1)
    varout(3) = 0.0d0
    WRITE(30,1011) (varout(m), m=1,3)
  
    enddo
  
    close(30)

1000  format(A,/,"1",/,"Flow Field",/,"1 ",I2,/,"x [m]",/,"y [m]")
1010  format("0 0"/,I6,I6," 0",/,"Unstructured")
1020  FORMAT(1P,20E14.6)
1040  FORMAT(1I2,3I7)

1006  FORMAT("# vtk DataFile Version 2.0" )
1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")

1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I6,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
1013  FORMAT('CELLS ',I5,I7)
1014  FORMAT('CELL_TYPES ',I5)
1015  FORMAT(I1)
1012  FORMAT('POINT_DATA ',I5)

! ######################################################################################
  end subroutine output_vtk4
  
  
    
    subroutine output_vtk5
    use modDataTypes
    use prmflow
    use modInterfaces
    implicit none

    character(chrlen) :: fname
    integer           :: ierr, nquant, i, m, CELLTYPE,id1,id2,id,inn,j
    real(rtype)       :: rrho, u, v, e, press, temp, c, mach, ttot, ptot
    real(rtype)       :: ptloss, pratio, ptotinf, gam1, ggm1, visc, machis
    real(rtype)       :: varout(mxqfield+2)
  
    id=302

    open(38, file="requiredID.vtk",form='formatted')

    write(38,1006)
    write(38,1007)
    write(38,1008) 
    write(38,*)""
    write(38,1009)
    write(38,1004) nbers(id)
  
    do j=1,nbers(id)
    inn = conn(j,id)
    write(38,1011) x(inn),y(inn),0.0d0
    enddo
  
    close(38)

1000  format(A,/,"1",/,"Flow Field",/,"1 ",I2,/,"x [m]",/,"y [m]")
1010  format("0 0"/,I6,I6," 0",/,"Unstructured")
1020  FORMAT(1P,20E14.6)
1040  FORMAT(1I2,3I7)

1006  FORMAT("# vtk DataFile Version 2.0" )
1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")

1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I6,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
1013  FORMAT('CELLS ',I5,I7)
1014  FORMAT('CELL_TYPES ',I5)
1015  FORMAT(I1)
1012  FORMAT('POINT_DATA ',I5)

  end subroutine output_vtk5
  
    subroutine previewgrid
    use modDataTypes
    use prmflow
    use modInterfaces
    implicit none

    character(chrlen) :: fname
    integer           :: ierr, nquant, i, m, CELLTYPE
    real(rtype)       :: rrho, u, v, e, press, temp, c, mach, ttot, ptot
    real(rtype)       :: ptloss, pratio, ptotinf, gam1, ggm1, visc, machis
    real(rtype)       :: varout(mxqfield+2)


    write(fname,"(A,I5.5,A)") Trim(title),"_preview.vtk"
  
    open(30, file=fname, status="unknown", action="write", iostat=ierr)
    if (ierr /= 0) call EM( "cannot open plot file (field)" )

    write(30,1006)
    write(30,1007)
    write(30,1008) 
    write(30,*)""
    write(30,1009)
    write(30,1004) phynod
  
    DO i=1,phynod
    varout(1) = x(i)
    varout(2) = y(i)
    varout(3) = 0.0d0
    WRITE(30,1011) (varout(m), m=1,3)
    enddo
  
    CELLTYPE = 4*nt
      
    write(30,*) ""
    write(30,1013) nt,CELLTYPE
      
    do i=1,nt
    write(30,1040) 3,elem(1,i)-1,elem(2,i)-1,elem(3,i)-1
    enddo
  
    write(30,*) ""
    write(30,1014) nt
  
    do i=1,nt
    write(30,1015) 5
    enddo 

    close(30)

1000  format(A,/,"1",/,"Flow Field",/,"1 ",I2,/,"x [m]",/,"y [m]")
1010  format("0 0"/,I6,I6," 0",/,"Unstructured")
1020  FORMAT(1P,20E14.6)
1040  FORMAT(1I2,3I7)

1006  FORMAT("# vtk DataFile Version 2.0" )
1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")

1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I7,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
1013  FORMAT('CELLS ',I7,I7)
1014  FORMAT('CELL_TYPES ',I7)
1015  FORMAT(I1)
1012  FORMAT('POINT_DATA ',I7)

    end subroutine previewgrid
    
  subroutine output_vtkv3
! ##########################################################################################
! WRITING VTK V3.0 FILE FORMAT 
! ##########################################################################################
  use modDataTypes
  use prmflow
  use modInterfaces
  implicit none

  character(chrlen) :: fname3
  integer           :: ierr, nquant, i, m, CELLTYPE
  double precision       :: rrho, u, v, e, press, temp, c, mach, ttot, ptot
  double precision       :: ptloss, pratio, ptotinf, gam1, ggm1, visc, machis
  double precision       :: varout(mxqfield+2)
  double precision       :: lpress,cp
! ######################################################################################

  write(fname3,"(A,A,I5.5,A)") Trim(title),"_",iter,".vtk"
  open(30, file="./output/"//fname3, status="unknown", action="write", iostat=ierr)
  
  write(30,1006)
  write(30,111) machinf,alpha
  write(30,1008) 
  write(30,1009)
  write(30,1004) phynod
  
  DO i=1,phynod
    varout(1) = x(i)
    varout(2) = y(i)
    varout(3) = 0.0d0
    WRITE(30,1011) (varout(m), m=1,3)
  enddo
  
  CELLTYPE = 4*nt
      
  write(30,*) ""
  write(30,1013) nt,CELLTYPE
      
  do i=1,nt
    write(30,'(i4,3i10)') 3,elem(1,i)-1,elem(2,i)-1,elem(3,i)-1
  enddo
  
  write(30,*) ""
  write(30,1014) nt
  
  do i=1,nt
    write(30,1015) 5
  enddo 

   write(30,*) ''
   write(30,1012) phynod
      
!  WRITE OUTPUTS
    write(30,*) 'SCALARS M float'  
    write(30,*) 'LOOKUP_TABLE default'
      
    DO i=1,phynod
     rrho  = 1.D0/cv(1,i)
     u     = cv(2,i)*rrho
     v     = cv(3,i)*rrho
     c     = dv(3,i)
     mach  = dsqrt(u*u+v*v)/c
     write(30,1020) mach
    ENDDO 
!
    write(30,*) 'SCALARS Rho float'  
    write(30,*) 'LOOKUP_TABLE default'
    do i=1,phynod
    write(30,1020) cv(1,i)
    enddo
    write(30,*) 'SCALARS Rho-U float'  
    write(30,*) 'LOOKUP_TABLE default'
    do i=1,phynod
    write(30,1020) cv(2,i)
    enddo
    write(30,*) 'SCALARS Rho-V float'  
    write(30,*) 'LOOKUP_TABLE default'
    do i=1,phynod
    write(30,1020) cv(3,i)
    enddo
    write(30,*) 'SCALARS Rho-E float'  
    write(30,*) 'LOOKUP_TABLE default'
    do i=1,phynod
    write(30,1020) cv(4,i)
    enddo
    write(30,*) 'SCALARS Cp float'  
    write(30,*) 'LOOKUP_TABLE default'
    do i=1,phynod
     lpress = dv(1,i)
     cp         = 2.0d0*(pinf-lpress)/(rhoinf*qinf*qinf)
     write(30,1020) cp
    enddo
!    write(30,25) 'Velocity'
!    do i=1,phynod
!    write(30,'(3e18.8)') cv(2,i)/cv(1,i),cv(3,i)/cv(1,i),0.0
!    enddo

    close(30)

1000  format(A,/,"1",/,"Flow Field",/,"1 ",I2,/,"x [m]",/,"y [m]")
1010  format("0 0"/,I6,I6," 0",/,"Unstructured")
1020  FORMAT(1P,20E14.6)
!1040  FORMAT(1I2,3I7)

1006  FORMAT("# vtk DataFile Version 3.0" )
!1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")

1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I7,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
1013  FORMAT('CELLS ',I7,I7)
1014  FORMAT('CELL_TYPES ',I7)
1015  FORMAT(I1)
1012  FORMAT('POINT_DATA ',I7)
25    format('VECTORS ', a10, '   float')
111   format('Mach =', f6.3, 2x, ' AOA = ', f6.3)
112   format('Mach =', f6.3, 2x, ' AOA = ', f6.3, ' Reynolds = ', e10.4)

! ######################################################################################
  end subroutine output_vtkv3