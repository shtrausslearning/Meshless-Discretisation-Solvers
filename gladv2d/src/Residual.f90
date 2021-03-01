
    subroutine Residual
    use prm_flow
    implicit none

    integer          :: i, j, p
    double precision :: sx, sy, du, dm, dp, sl, sr, vl, vr, fxi, fyi 
    double precision :: vxi, vyi, xij, yij, fij, vxij, vyij, vn, valimit

    do i=1,npts
    res(i) = 0.0d0
    if(ptype(i).eq.1)then
    call AdvecVel(coord(1,i), coord(2,i), vxi, vyi)
    fxi   = vxi*var(i)
    fyi   = vyi*var(i)
    do j=1,nnbr(i)
            
        p  = conn(j,i)                              ! satellite
        sx = coord(1,p) - coord(1,i)                ! dx
        sy = coord(2,p) - coord(2,i)                ! dy 
        du = var(p) - var(i)                        ! check
        dm = 2.0d0*(sx*varx(i) + sy*vary(i)) - du   ! du or dQ
        dp = 2.0d0*(sx*varx(p) + sy*vary(p)) - du
        sl = valimit(dm, du)
        sr = valimit(dp, du)

        vl = var(i) + 0.25d0*sl*( (1.0d0-kkk*sl)*dm + (1.0d0+kkk*sl)*du)
        vr = var(p) - 0.25d0*sr*( (1.0d0-kkk*sr)*dp + (1.0d0+kkk*sr)*du)

        xij = 0.5d0*(coord(1,i) + coord(1,p)) ! midpoint location x
        yij = 0.5d0*(coord(2,i) + coord(2,p)) ! midpoint location y 
        call AdvecVel(xij, yij, vxij, vyij)   ! midpoint U      
        vn = ax(j,i)*vxij + ay(j,i)*vyij      ! interface velocity ijcoord

        if(vn .gt. 0.0d0)then
         fij = vn*vl
        else
         fij = vn*vr
        endif

        res(i) = res(i) + 2.0d0*(fij - (ax(j,i)*fxi + ay(j,i)*fyi))

        
    enddo
    endif
    enddo

    return
    End Subroutine Residual

      double precision function valimit(a, b)
      use prm_flow
      implicit none
      double precision :: a, b

      double precision num, den, epsi
      parameter(epsi=1.0d-15)

      if(order.eq.1)then
         valimit = 0.0d0
      else if(order.eq.2)then
         valimit = 1.0d0
      else if(order.eq.3)then
         num     = 2.0d0*a*b + epsi
         den     = a**2 + b**2 + epsi
         valimit = dmax1(0.0d0, num/den)
      else
         print*,'order is not properly set'
         print*,'I dont know order =',order
         stop
      endif

      return
      end
