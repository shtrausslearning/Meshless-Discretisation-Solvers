
      subroutine output_vtk(iter)
      use prm_flow
      implicit none

    
      integer :: mt,ic,id,is,errFlag
      real :: u,v,cc,mach,rrho
      character(6) :: c
      
      integer, intent(in) :: iter  
    
      mt = 1
      
      write(c, '(I6)' ) iter
      c = adjustr(c)
      do is = 1,6
       if( ichar( c(is:is) ) == 32 ) c(is:is) = '0'
      enddo
     
      open(mt, file="./output/var_"//c//".vtk",form='formatted')
      write(mt,1006)
      write(mt,1007)
      write(mt,1008)
      write(mt,*)
      write(mt,1009)
      write(mt,1004) np0
      do ic = 1,np0
      if( var(ic) < 1.0d-10 )then
      write(mt,1011) coord(1,ic),coord(2,ic),0.0d0
      else
!      write(mt,1011) coord(1,ic),coord(2,ic),0.0d0
      write(mt,1011) coord(1,ic),coord(2,ic),var(ic)
      endif
      enddo
      
      write(mt,*) ''
      write(mt,1012) np0
      write(mt,*) 'SCALARS VarU float'  
      write(mt,*) 'LOOKUP_TABLE default'
      
      DO ic=1,np0
       if( var(ic) < 1d-10 )then
        write(mt,1020) 0.0d0
       else
        WRITE(mt,1020) var(ic)
       endif
      ENDDO  

      close(mt)

1006  FORMAT("# vtk DataFile Version 2.0" )
1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")
1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I6,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
1020  FORMAT(1P,20E14.6)
1012  FORMAT('POINT_DATA ',I5)

      end subroutine output_vtk
      
      subroutine output_vtkdiff(iter)
      use prm_flow
      implicit none

    
      integer :: mt,ic,id,is,errFlag
      real :: u,v,cc,mach,rrho
      character(6) :: c
      
      integer, intent(in) :: iter  
    
      mt = 1
      
      write(c, '(I6)' ) iter
      c = adjustr(c)
      do is = 1,6
       if( ichar( c(is:is) ) == 32 ) c(is:is) = '0'
      enddo
     
      open(mt, file="./output/error_"//c//".vtk",form='formatted')
      write(mt,1006)
      write(mt,1007)
      write(mt,1008)
      write(mt,*)
      write(mt,1009)
      write(mt,1004) np0
      do ic = 1,np0
      if( var(ic) < 1.0d-10 )then
      write(mt,1011) coord(1,ic),coord(2,ic),0.0d0
      else
      write(mt,1011) coord(1,ic),coord(2,ic),0.0d0
      endif
      enddo
      
      write(mt,*) ''
      write(mt,1012) np0
      write(mt,*) 'SCALARS difference float'  
      write(mt,*) 'LOOKUP_TABLE default'
      
      DO ic=1,np0
       if( var(ic) < 1d-10 )then
        write(mt,1020) 0.0d0
       else
        WRITE(mt,1020) var(ic)-varini(ic)
       endif
      ENDDO  

      close(mt)

1006  FORMAT("# vtk DataFile Version 2.0" )
1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")
1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I6,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
1020  FORMAT(1P,20E14.6)
1012  FORMAT('POINT_DATA ',I5)

      end subroutine output_vtkdiff
      
!      subroutine output_vtk2(iter)
!      use prm_flow
!      implicit none
!
!    
!      integer :: mt,ic,id,is,errFlag
!      real :: u,v,cc,mach,rrho
!      character(6) :: c
!      
!      integer, intent(in) :: iter  
!    
!      mt = 1
!      
!      write(c, '(I6)' ) iter
!      c = adjustr(c)
!      do is = 1,6
!       if( ichar( c(is:is) ) == 32 ) c(is:is) = '0'
!      enddo
!     
!      open(mt, file="./output/vtkoutputB"//c//".vtk",form='formatted')
!      write(mt,1006)
!      write(mt,1007)
!      write(mt,1008)
!      write(mt,*)
!      write(mt,1009)
!      write(mt,1004) np0
!      do ic = 1,np0
!      if( var(ic) < 1.0d-10 )then
!      write(mt,1011) coord(1,ic),coord(2,ic),0.0d0
!      else
!      write(mt,1011) coord(1,ic),coord(2,ic),0.0d0
!      endif
!      enddo
!      
!      write(mt,*) ''
!      write(mt,1012) np0
!      write(mt,*) 'SCALARS LOCAL-MACH float'  
!      write(mt,*) 'LOOKUP_TABLE default'
!      
!      DO ic=1,np0
!       if( var(ic) < 1d-10 )then
!        write(mt,1020) 0.0d0
!       else
!        WRITE(mt,1020) var(ic)
!       endif
!      ENDDO  
!
!      close(mt)
!
!1006  FORMAT("# vtk DataFile Version 2.0" )
!1007  FORMAT("VTK FORMAT")
!1008  FORMAT("ASCII")
!1009  FORMAT("DATASET UNSTRUCTURED_GRID")
!1004  FORMAT('POINTS',I6,' float ')
!1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
!1020  FORMAT(1P,20E14.6)
!1012  FORMAT('POINT_DATA ',I5)
!
!      end subroutine output_vtk2

      
!      subroutine output_vtkbcm(iter)
!      use prm_flow
!      implicit none
!
!    
!      integer :: mt,ic,id,is,errFlag
!      real :: u,v,cc,mach,rrho
!      character(6) :: c
!      
!      integer, intent(in) :: iter  
!    
!      mt = 1
!      
!      write(c, '(I6)' ) iter
!      c = adjustr(c)
!      do is = 1,6
!       if( ichar( c(is:is) ) == 32 ) c(is:is) = '0'
!      enddo
!     
!      open(mt, file="./output/vtkoutput"//c//".vtk",form='formatted')
!      write(mt,1006)
!      write(mt,1007)
!      write(mt,1008)
!      write(mt,*)
!      write(mt,1009)
!      write(mt,1004) np0
!      do ic = 1,np0
!      if( var(ic) < 1.0d-10 )then
!      write(mt,1011) coord(1,ic),coord(2,ic),0.0d0
!      else
!      write(mt,1011) coord(1,ic),coord(2,ic),var(ic)
!      endif
!      enddo
!      
!      write(mt,*) ''
!      write(mt,1012) np0
!      write(mt,*) 'SCALARS LOCAL-MACH float'  
!      write(mt,*) 'LOOKUP_TABLE default'
!      
!      DO ic=1,np0
!       if( var(ic) < 1d-10 )then
!        write(mt,1020) 0.0d0
!       else
!        WRITE(mt,1020) var(ic)
!       endif
!      ENDDO  
!
!      close(mt)
!
!1006  FORMAT("# vtk DataFile Version 2.0" )
!1007  FORMAT("VTK FORMAT")
!1008  FORMAT("ASCII")
!1009  FORMAT("DATASET UNSTRUCTURED_GRID")
!1004  FORMAT('POINTS',I6,' float ')
!1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
!1020  FORMAT(1P,20E14.6)
!1012  FORMAT('POINT_DATA ',I5)
!
!      end subroutine output_vtkbcm

