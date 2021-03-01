
      subroutine AdvecVel(x,y,vx,vy)
      use prm_flow
      implicit none
      double precision :: x, y, vx, vy

      if(veltype.eq.1)then
         vx = 1.0d0
         vy = 0.0d0
      elseif(veltype.eq.2)then
         vx = 1.0d0
         vy = 1.0d0
      elseif(veltype.eq.3)then
         vx = y
         vy = -x
      else
         write(*,*) 'Advec velocity type veltype is unknown'
         write(*,*) 'veltype =',veltype
         stop
      endif

      return
      end
