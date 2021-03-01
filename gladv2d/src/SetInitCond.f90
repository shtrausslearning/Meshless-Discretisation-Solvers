
!     Set initial condition
      subroutine SetInitCond
      use prm_flow
      implicit none

      integer          :: i, p
      double precision :: x, y, InitCond

      varmin0 = 1.0d20
      varmax0 =-1.0d20

      do i=1,np0
        x       = coord(1,i)
        y       = coord(2,i)
        var(i)  = InitCond(x, y)
        varmin0 = dmin1(varmin0, var(i))
        varmax0 = dmax1(varmax0, var(i))
        varini(i) = var(i)
      enddo
      
!     periodic condition
      do i=np0+1,npts
       p = conn(1,i)   ! neighbour point
       var(i) = var(p)
      enddo

       write(*,*) 'Initial min,max =',varmin0,varmax0
       
       do i=1,np0
        varini(i) = var(i)
       enddo 

      return
      end


!     Initial condition
      double precision function InitCond(x,y)
      use prm_flow
      implicit none

      double precision :: x, y
      double precision :: r, r2

      if(ictype.eq.1)then
!        use xmin=0, xmax=1, ymin=0, ymax=1
         r2 = (x-0.5d0)**2 + (y-0.5d0)**2
         r  = dsqrt(r2)
         if(r.le.0.5d0)then
!            InitCond = (0.25d0**2 - r2)/0.25d0**2
            InitCond = 1.0d0 - 4.0d0**2*r2
            if(InitCond.lt.0.0d0)InitCond = 0.0d0
         else
            InitCond = 0.0d0
         endif
      elseif(ictype.eq.2) then
       InitCond = sin(x)*cos(y)
          
      else      
         print*,'ictype is not known, ictype =',ictype
         stop
      endif

      return
      end
