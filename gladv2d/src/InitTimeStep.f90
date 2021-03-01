!c     find initial global time step allowed by stability condition
      subroutine InitTimeStep
      use prm_flow
      implicit none

      integer          :: i, j
      double precision :: x, y, vx, vy, v, vsum, dtlocal, dtglobal

      dtglobal = 1.0d20

      do i=1,np0
         vsum = 0.0d0
         x    = coord(1,i)
         y    = coord(2,i)
         call AdvecVel(x,y,vx,vy)
         do j=1,nnbr(i)
            v    = ax(j,i)*vx + ay(j,i)
            vsum = vsum + dabs(v)
         enddo
         vsum    = vsum + 1.0d-16
         dtlocal = 1.0d0/vsum
         dtglobal= dmin1(dtglobal, dtlocal)
      enddo
    
      write(*,*) 'Chosen time step =',dt
      write(*,*) 'CFL time step    =',dtglobal

      if(dt.gt.dtglobal+1.0d-15)then
         pause 'Time step is greater than stability limit'
      endif

      return
      end
