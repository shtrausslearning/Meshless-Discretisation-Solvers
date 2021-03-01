      subroutine Update(rks)
      use prm_flow
      implicit none

      integer          :: i, rks
      double precision :: aa(10), bb(10)

!c    3-stage TVD-RK 
      !aa(1) = 0.0d0
      !aa(2) = 3.0d0/4.0d0
      !aa(3) = 1.0d0/3.0d0
      !bb(1) = 1.0d0
      !bb(2) = 1.0d0/4.0d0
      !bb(3) = 2.0d0/3.0d0
!     Standard 3-Stage RK
      aa(1) = 0.1481d0
      aa(2) = 0.4d0
      aa(3) = 1.0d0

      do i=1,npts
      if(ptype(i).eq.1)then
      !var(i) = aa(rks)*var_old(i) + bb(rks)*(var(i) - dt*res(i))
      var(i) = var_old(i) - aa(rks)*dt*res(i)
      endif
      enddo

!c     apply periodic condition
      do i=1,npts
      if(ptype(i).eq.2)then
       var(i) = var(conn(1,i))
      endif
      enddo

      return
      end
