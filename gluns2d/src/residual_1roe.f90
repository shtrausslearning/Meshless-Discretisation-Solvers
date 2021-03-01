

    Subroutine residual_roe
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    
    integer :: i,j,p,m,ic,mm,jjj,jj
    real(rtype) :: ds,nx,ny,rrho,rl,ul,vl,pl,hl,gaml,rr,ur,vr,pr,hr,gamr
    real(rtype) :: rav,gam1,dd,dd1,uav,hav,vav,q2a,c2a,cav,uv,du,encorf
    real(rtype) :: h1,h2,h4,delta,eabs1,eabs2,eabs4,fd(4),h3,h5
    real(rtype) :: rhoul,rhovl,rhl,ql,rhour,rhovr,rhr,qr,pav,fi(4),fc(4),fij(4)
    real(rtype) :: fx(4),fy(4),ggm1
    real(rtype) :: dup(4),dm(4),dp(4),sl(4),sr(4),kkk
    
    real(rtype) :: voll,volr,rx,ry,limval,eps2nr,eps2nl,d2r,d2l,d1minr,d1minl,d1maxr
    real(rtype) :: d1maxl

    do i=1,phynod 
    rhs(:,i) = 0.0d0
    if( ptype(i) .eq. 1 )then
    do j=1,nbers(i)
    
     p = conn(j,i)  
     ds = dsqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
     nx = aij(j,i)/ds
     ny = bij(j,i)/ds

!    LEFT STATE DEFINITION
     rl   = cv(1,i)
     ul   = cv(2,i)/cv(1,i)
     vl   = cv(3,i)/cv(1,i)
     pl   = dv(1,i)
     ql = (ul*nx + vl*ny)   ! Convariant Velocity
     gam1 = dv(4,i) - 1.0d0
     ggm1 = dv(4,i)/gam1
     hl  = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)

!    RIGHT STATE DEFINITION
     rr   = cv(1,p)
     ur   = cv(2,p)/cv(1,p)
     vr   = cv(3,p)/cv(1,p)
     pr   = dv(1,p)
     qr =  (ur*nx+ vr*ny)
     gam1 = dv(4,p) - 1.0d0
     ggm1 = dv(4,p)/gam1
     hr  = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
     
!    roe averaging
     rav  = dsqrt(rl*rr)
     gam1 = 0.5d0*(dv(4,i)+dv(4,p)) - 1.0d0
     dd   = rav/rl
     dd1  = 1.0d0/(1.0d0+dd)
     uav  = (ul+dd*ur)*dd1
     vav  = (vl+dd*vr)*dd1
     hav  = (hl+dd*hr)*dd1
     q2a  = 0.5d0*(uav*uav+vav*vav)
     c2a  = gam1*(hav-q2a)
     
     cav  = dsqrt(c2a)
     uv   = uav*nx + vav*ny
     du   = (ur-ul)*nx + (vr-vl)*ny
     
!    eigenvalues
     h1 = abs(uv - cav)
     h2 = abs(uv)
     h4 = abs(uv + cav)
     delta =  epsentr*h4
     
     eabs1 = EntropyCorr1( h1,delta )
     eabs2 = EntropyCorr1( h2,delta )
     eabs4 = EntropyCorr1( h4,delta )
    
     h1 = rav*cav*du
     h2 = eabs1*(pr-pl - h1)/(2.D0*c2a)
     h3 = eabs2*(rr-rl - (pr-pl)/c2a)
     h4 = eabs2*rav
     h5 = eabs4*(pr-pl + h1)/(2.D0*c2a)
     
     fd(1) = h2 + h3 + h5
     fd(2) = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-du*nx) + h5*(uav+cav*nx)
     fd(3) = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-du*ny) + h5*(vav+cav*ny)
     fd(4) = h2*(hav-cav*uv) + h3*q2a + h4*(uav*(ur-ul)+vav*(vr-vl)-uv*du) + h5*(hav+cav*uv)
     
!   individual flux components of cloud components
    pav = 0.5d0*(pl+pr)
    fc(1) = 0.5d0*( ql*rl     + qr*rr        ) 
    fc(2) = 0.5d0*( ql*rl*ul  + qr*rr*ur  ) + nx*pav 
    fc(3) = 0.5d0*( ql*rl*vl  + qr*rr*vr  ) + ny*pav 
    fc(4) = 0.5d0*( ql*rl*hl  + qr*rr*hr   )
    
   !individual flux components of cloud components
    fi(1) = ql*rl 
    fi(2) = ql*rl*ul + nx*pl 
    fi(3) = ql*rl*vl + ny*pl
    fi(4) = ql*rl*hl  
    
    fij(1) = fc(1) - 0.5d0*fd(1)
    fij(2) = fc(2) - 0.5d0*fd(2)
    fij(3) = fc(3) - 0.5d0*fd(3)
    fij(4) = fc(4) - 0.5d0*fd(4)

   !residual
    rhs(1,i) = rhs(1,i) + 2.0d0*( fij(1) - fi(1) )
    rhs(2,i) = rhs(2,i) + 2.0d0*( fij(2) - fi(2) )
    rhs(3,i) = rhs(3,i) + 2.0d0*( fij(3) - fi(3) )
    rhs(4,i) = rhs(4,i) + 2.0d0*( fij(4) - fi(4) )
    
    enddo
    endif
    enddo

    
contains

  real(rtype) function EntropyCorr1( z,d )
    implicit none
    real(rtype) :: z, d

    if (z > d) then
      EntropyCorr1 = z
    else
      EntropyCorr1 = 0.5D0*(z*z+d*d)/d
    endif
  end function EntropyCorr1
  
    real(rtype) function valimit(a,b,epsi)
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    real(rtype) :: a, b
    double precision num, den, epsi

    if(iorder.eq.1)then
    valimit = 0.0d0
    elseif(iorder.eq.2)then
    valimit = 1.0d0
    elseif(iorder.eq.3)then ! Van Albada 
    num     = 2.0d0*a*b + epsi
    den     = a**2 + b**2 + epsi
    valimit = dmax1(0.0d0, num/den)
    endif

    return
    end function valimit
    
    real(rtype) function valimit2(a, b, ds)
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    real(rtype) :: a, b
    double precision num, den, epsi,ds
    parameter(epsi=1.0d-15)

    valimit2  = 1 - abs( ( a-b )/( abs(a) + abs(b) + epsi*ds ) )**3.0d0

    return
    end function valimit2
    
    real(rtype) function Venkat( d2,d1min,d1max,eps2 )
    implicit none
    real(rtype) :: d2, d1min, d1max, eps2
    real(rtype) :: num, den
    
    eps2 = 1d20

    Venkat = 1.D0
    if (d2 > 1.D-12) then
      num    = (d1max*d1max+eps2)*d2 + 2.D0*d2*d2*d1max
      den    = d2*(d1max*d1max+2.D0*d2*d2+d1max*d2+eps2)
      Venkat = num/den
    else if (d2 < -1.D-12) then
      num    = (d1min*d1min+eps2)*d2 + 2.D0*d2*d2*d1min
      den    = d2*(d1min*d1min+2.D0*d2*d2+d1min*d2+eps2)
      Venkat = num/den
    endif
    end function Venkat
    
    end subroutine residual_roe
    


    Subroutine residual_roer
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    
    integer :: i,j,p,m,ic,mm,jjj,jj
    real(rtype) :: ds,nx,ny,rrho,rl,ul,vl,pl,hl,gaml,rr,ur,vr,pr,hr,gamr
    real(rtype) :: rav,gam1,dd,dd1,uav,hav,vav,q2a,c2a,cav,uv,du,encorf
    real(rtype) :: h1,h2,h4,delta,eabs1,eabs2,eabs4,fd(4),h3,h5
    real(rtype) :: rhoul,rhovl,rhl,ql,rhour,rhovr,rhr,qr,pav,fi(4),fc(4),fij(4)
    real(rtype) :: fx(4),fy(4),ggm1
    real(rtype) :: dup(4),dm(4),dp(4),sl(4),sr(4),kkk
    
    real(rtype) :: voll,volr,rx,ry,limval,eps2nr,eps2nl,d2r,d2l,d1minr,d1minl,d1maxr
    real(rtype) :: d1maxl

    do i=1,phynod 
    rhs(:,i) = 0.0d0
    if( ptype(i) .eq. 1 )then
    do j=1,nbers(i)
    
     p = conn(j,i)  
     ds = dsqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
     nx = aij(j,i)/ds
     ny = bij(j,i)/ds

    rx   = (x(p)-x(i)); ry   = (y(p)-y(i))

    voll = vol(i)**1.5D0    ;volr = vol(p)**1.5D0

    eps2nl   = eps2(1)*voll
    eps2nr   = eps2(1)*volr    
    d1minl   = umin(1,i) - cv(1,i)
    d1maxl   = umax(1,i) - cv(1,i)
    d1minr   = umin(1,p) - cv(1,p)
    d1maxr   = umax(1,p) - cv(1,p)
    d2l      =  gradx(1,i)*rx + grady(1,i)*ry
!    d2r      =  -gradx(1,p)*rx - grady(1,p)*ry
    d2r      =   gradx(1,p)*rx + grady(1,p)*ry
    limval   = Venkat( d2l,d1minl,d1maxl,eps2nl )
    sl(1) = Min(limval,lim(1,i))
    limval   = Venkat( d2r,d1minr,d1maxr,eps2nr )
    sr(1) = Min(limval,lim(1,p))

    ul       = cv(2,i)/cv(1,i)
    ur       = cv(2,p)/cv(1,p)
    eps2nl   = eps2(2)*voll
    eps2nr   = eps2(2)*volr
    d1minl   = umin(2,i) - ul
    d1maxl   = umax(2,i) - ul
    d1minr   = umin(2,p) - ur
    d1maxr   = umax(2,p) - ur
    d2l      =  gradx(2,i)*rx + grady(2,i)*ry
!    d2r      =  - gradx(2,p)*rx - grady(2,p)*ry
    d2r      =  gradx(2,p)*rx + grady(2,p)*ry
    limval   = Venkat( d2l,d1minl,d1maxl,eps2nl )
    sl(2)    = Min(limval,lim(2,i))
    limval   = Venkat( d2r,d1minr,d1maxr,eps2nr )
    sr(2)    = Min(limval,lim(2,p))

    vl       = cv(3,i)/cv(1,i)
    vr       = cv(3,p)/cv(1,p)
    eps2nl   = eps2(3)*voll
    eps2nr   = eps2(3)*volr
    d1minl   = umin(3,i) - vl
    d1maxl   = umax(3,i) - vl
    d1minr   = umin(3,p) - vr
    d1maxr   = umax(3,p) - vr
    d2l      =  gradx(3,i)*rx + grady(3,i)*ry
!    d2r      =  - gradx(3,p)*rx - grady(3,p)*ry
    d2r      =  gradx(3,p)*rx + grady(3,p)*ry
    limval   = Venkat( d2l,d1minl,d1maxl,eps2nl )
    sl(3)    = Min(limval,lim(3,i))
    limval   = Venkat( d2r,d1minr,d1maxr,eps2nr )
    sr(3)    = Min(limval,lim(3,p))

! - pressure

    eps2nl   = eps2(4)*voll
    eps2nr   = eps2(4)*volr
    d1minl   = umin(4,i) - dv(1,i)
    d1maxl   = umax(4,i) - dv(1,i)
    d1minr   = umin(4,p) - dv(1,p)
    d1maxr   = umax(4,p) - dv(1,p)
    d2l      =  gradx(4,i)*rx + grady(4,i)*ry
!    d2r      =  - gradx(4,p)*rx - grady(4,p)*ry
    d2r      =  gradx(4,p)*rx + grady(4,p)*ry
    limval   = Venkat( d2l,d1minl,d1maxl,eps2nl )
    sl(4)    = Min(limval,lim(4,i))
    limval   = Venkat( d2r,d1minr,d1maxr,eps2nr )
    sr(4)    = Min(limval,lim(4,p))

!!    LEFT STATE DEFINITION
!     rl   = cv(1,i)
!     ul   = cv(2,i)/cv(1,i)
!     vl   = cv(3,i)/cv(1,i)
!     pl   = dv(1,i)
!     ql = (ul*nx + vl*ny)   ! Convariant Velocity
!     gam1 = dv(4,i) - 1.0d0
!     ggm1 = dv(4,i)/gam1
!     hl  = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)
!
!!    RIGHT STATE DEFINITION
!     rr   = cv(1,p)
!     ur   = cv(2,p)/cv(1,p)
!     vr   = cv(3,p)/cv(1,p)
!     pr   = dv(1,p)
!     qr =  (ur*nx+ vr*ny)
!     gam1 = dv(4,p) - 1.0d0
!     ggm1 = dv(4,p)/gam1
!     hr  = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
     
!    LEFT STATE DEFINITION
     rl   = cv(1,i)         + sl(1)*(gradx(1,i)*rx+grady(1,i)*ry)
     ul   = cv(2,i)/cv(1,i) + sl(2)*(gradx(2,i)*rx+grady(2,i)*ry)
     vl   = cv(3,i)/cv(1,i) + sl(3)*(gradx(3,i)*rx+grady(3,i)*ry)
     pl   = dv(1,i)         + sl(4)*(gradx(4,i)*rx+grady(4,i)*ry)
     ql = (ul*nx + vl*ny)   ! Convariant Velocity
     gam1 = dv(4,i) - 1.0d0
     ggm1 = dv(4,i)/gam1
     hl  = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)

!    RIGHT STATE DEFINITION
     rr   = cv(1,p)         - sr(1)*(gradx(1,p)*rx+grady(1,p)*ry)
     ur   = cv(2,p)/cv(1,p) - sr(2)*(gradx(2,p)*rx+grady(2,p)*ry)
     vr   = cv(3,p)/cv(1,p) - sr(3)*(gradx(3,p)*rx+grady(3,p)*ry)
     pr   = dv(1,p)         - sr(4)*(gradx(4,p)*rx+grady(4,p)*ry)
     qr =  (ur*nx+ vr*ny)
     gam1 = dv(4,p) - 1.0d0
     ggm1 = dv(4,p)/gam1
     hr  = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
     
!    roe averaging
     rav  = dsqrt(rl*rr)
     gam1 = 0.5d0*(dv(4,i)+dv(4,p)) - 1.0d0
     dd   = rav/rl
     dd1  = 1.0d0/(1.0d0+dd)
     uav  = (ul+dd*ur)*dd1
     vav  = (vl+dd*vr)*dd1
     hav  = (hl+dd*hr)*dd1
     q2a  = 0.5d0*(uav*uav+vav*vav)
     c2a  = gam1*(hav-q2a)
     
     cav  = dsqrt(c2a)
     uv   = uav*nx + vav*ny
     du   = (ur-ul)*nx + (vr-vl)*ny
     
!    eigenvalues
     h1 = abs(uv - cav)
     h2 = abs(uv)
     h4 = abs(uv + cav)
     delta =  epsentr*h4
     
     eabs1 = EntropyCorr1( h1,delta )
     eabs2 = EntropyCorr1( h2,delta )
     eabs4 = EntropyCorr1( h4,delta )
    
     h1 = rav*cav*du
     h2 = eabs1*(pr-pl - h1)/(2.D0*c2a)
     h3 = eabs2*(rr-rl - (pr-pl)/c2a)
     h4 = eabs2*rav
     h5 = eabs4*(pr-pl + h1)/(2.D0*c2a)
     
     fd(1) = h2 + h3 + h5
     fd(2) = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-du*nx) + h5*(uav+cav*nx)
     fd(3) = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-du*ny) + h5*(vav+cav*ny)
     fd(4) = h2*(hav-cav*uv) + h3*q2a + h4*(uav*(ur-ul)+vav*(vr-vl)-uv*du) + h5*(hav+cav*uv)
             
!   individual flux components of cloud components
    pav = 0.5d0*(pl+pr)
    fc(1) = 0.5d0*( ql*rl     + qr*rr        ) 
    fc(2) = 0.5d0*( ql*rl*ul  + qr*rr*ur  ) + nx*pav 
    fc(3) = 0.5d0*( ql*rl*vl  + qr*rr*vr  ) + ny*pav 
    fc(4) = 0.5d0*( ql*rl*hl  + qr*rr*hr   )
    
   !individual flux components of cloud components
    fi(1) = ql*rl 
    fi(2) = ql*rl*ul + nx*pl 
    fi(3) = ql*rl*vl + ny*pl
    fi(4) = ql*rl*hl  
    
    fij(1) = fc(1) - 0.5d0*fd(1)
    fij(2) = fc(2) - 0.5d0*fd(2)
    fij(3) = fc(3) - 0.5d0*fd(3)
    fij(4) = fc(4) - 0.5d0*fd(4)

   !residual
    rhs(1,i) = rhs(1,i) + 2.0d0*( fij(1) - fi(1) )
    rhs(2,i) = rhs(2,i) + 2.0d0*( fij(2) - fi(2) )
    rhs(3,i) = rhs(3,i) + 2.0d0*( fij(3) - fi(3) )
    rhs(4,i) = rhs(4,i) + 2.0d0*( fij(4) - fi(4) )
    
    enddo
    endif
    enddo

    
contains

  real(rtype) function EntropyCorr1( z,d )
    implicit none
    real(rtype) :: z, d

    if (z > d) then
      EntropyCorr1 = z
    else
      EntropyCorr1 = 0.5D0*(z*z+d*d)/d
    endif
  end function EntropyCorr1
  
    real(rtype) function valimit(a,b,epsi)
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    real(rtype) :: a, b
    double precision num, den, epsi

    if(iorder.eq.1)then
    valimit = 0.0d0
    elseif(iorder.eq.2)then
    valimit = 1.0d0
    elseif(iorder.eq.3)then ! Van Albada 
    num     = 2.0d0*a*b + epsi
    den     = a**2 + b**2 + epsi
    valimit = dmax1(0.0d0, num/den)
    endif

    return
    end function valimit
    
    real(rtype) function valimit2(a, b, ds)
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    real(rtype) :: a, b
    double precision num, den, epsi,ds
    parameter(epsi=1.0d-15)

    valimit2  = 1 - abs( ( a-b )/( abs(a) + abs(b) + epsi*ds ) )**3.0d0

    return
    end function valimit2
    
    real(rtype) function Venkat( d2,d1min,d1max,eps2 )
    implicit none
    real(rtype) :: d2, d1min, d1max, eps2
    real(rtype) :: num, den

    Venkat = 1.D0
    if (d2 > 1.D-12) then
      num    = (d1max*d1max+eps2)*d2 + 2.D0*d2*d2*d1max
      den    = d2*(d1max*d1max+2.D0*d2*d2+d1max*d2+eps2)
      Venkat = num/den
    else if (d2 < -1.D-12) then
      num    = (d1min*d1min+eps2)*d2 + 2.D0*d2*d2*d1min
      den    = d2*(d1min*d1min+2.D0*d2*d2+d1min*d2+eps2)
      Venkat = num/den
    endif
    end function Venkat
    
    end subroutine residual_roer
    
    Subroutine residual_roeMU
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    
    integer :: i,j,m,ic,mm,jjj,ii,p,jj
    real(rtype) :: ds,nx,ny,rrho,rl,ul,vl,pl,hl,gaml,rr,ur,vr,pr,hr,gamr
    real(rtype) :: rav,gam1,dd,dd1,uav,hav,vav,q2a,c2a,cav,uv,dup,encorf
    real(rtype) :: h1,h2,h4,delta,eabs1,eabs2,eabs4,fd(4),h3,h5
    real(rtype) :: rhoul,rhovl,ql,rhour,rhovr,qr,pav,fi(4),fc(4),fij(4)
    real(rtype) :: fx(4),fy(4),ggm1
    real(rtype) :: sx,sy, dp(4), dm(4), du(4), sl(4), sr(4), eps, kkk
    real(rtype) :: voll,volr,eps2nl(4),eps2nr(4)
    
    do i=1,phynod 
    rhs(:,i) = 0.0d0
    if( ptype(i) .eq. 1 )then
    do j=1,nbers(i)
    
     p = conn(j,i)  
     sx = x(p) - x(i)
     sy = y(p) - y(i)
     ds = sqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
     nx = aij(j,i)/ds
     ny = bij(j,i)/ds
     
     du(1) = cv(1,p)         - cv(1,i)
     du(2) = cv(2,p)/cv(1,p) - cv(2,i)/cv(1,i)
     du(3) = cv(3,p)/cv(1,p) - cv(3,i)/cv(1,i)
     du(4) = dv(1,p)         - dv(1,i)
     
     voll = vol(i)**1.5D0 ;volr = vol(p)**1.5D0
     do jj=1,4
     eps2nl(jj)   = eps2(jj)*voll  ! activation limit eps from fvm
     eps2nr(jj)   = eps2(jj)*volr   
     enddo
     
     do jj=1,4
     dm(jj) = (sx*gradx(jj,i) + sy*grady(jj,i)) - du(jj)
     dp(jj) = (sx*gradx(jj,p) + sy*grady(jj,p)) - du(jj)
     sl(jj) = valimit(dm(jj),du(jj),eps2nl(jj))
     sr(jj) = valimit(dp(jj),du(jj),eps2nr(jj))
     enddo

!     kkk = 1.0d0/3.0d0 
     kkk = 1.0d0
!     kkk = 0.0d0

!    1st iorder -> sl/sr = 0
!    2nd iorder -> sl/sr = 1
!    3rd iorder -> limiter on sl/sr
     
!    LEFT STATE DEFINITION
     rl   = cv(1,i)         + 0.25d0*sl(1)*( (1.0d0-kkk*sl(1))*dm(1) + (1.0d0+kkk*sl(1))*du(1) )
     ul   = cv(2,i)/cv(1,i) + 0.25d0*sl(2)*( (1.0d0-kkk*sl(2))*dm(2) + (1.0d0+kkk*sl(2))*du(2) )
     vl   = cv(3,i)/cv(1,i) + 0.25d0*sl(3)*( (1.0d0-kkk*sl(3))*dm(3) + (1.0d0+kkk*sl(3))*du(3) )
     pl   = dv(1,i)         + 0.25d0*sl(4)*( (1.0d0-kkk*sl(4))*dm(4) + (1.0d0+kkk*sl(4))*du(4) )
     ql = ul*aij(j,i) + vl*bij(j,i)   ! Convariant Velocity
     gam1 = dv(4,i) - 1.0d0
     ggm1 = dv(4,i)/gam1
     hl  = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)
     
!    RIGHT STATE DEFINITION
     rr   = cv(1,p)           - 0.25d0*sr(1)*( (1.0d0-kkk*sr(1))*dp(1) + (1.0d0+kkk*sr(1))*du(1))
     ur   = cv(2,p)/cv(1,p)   - 0.25d0*sr(2)*( (1.0d0-kkk*sr(2))*dp(2) + (1.0d0+kkk*sr(2))*du(2))
     vr   = cv(3,p)/cv(1,p)   - 0.25d0*sr(3)*( (1.0d0-kkk*sr(3))*dp(3) + (1.0d0+kkk*sr(3))*du(3))
     pr   = dv(1,p)           - 0.25d0*sr(4)*( (1.0d0-kkk*sr(4))*dp(4) + (1.0d0+kkk*sr(4))*du(4))
     qr =  ur*aij(j,i)+ vr*bij(j,i)
     gam1 = dv(4,p) - 1.0d0
     ggm1 = dv(4,p)/gam1
     hr  = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
     
!    roe averaging
     rav  = sqrt(rl*rr)
     gam1 = 0.5d0*(dv(4,i)+dv(4,p)) - 1.0d0
     dd   = rav/rl
     dd1  = 1.0d0/(1.0d0+dd)
     uav  = (ul+dd*ur)*dd1
     vav  = (vl+dd*vr)*dd1
     hav  = (hl+dd*hr)*dd1
     q2a  = 0.5d0*(uav*uav+vav*vav)
     c2a  = gam1*(hav-q2a)
     
     cav  = sqrt(c2a)
     uv   = uav*nx + vav*ny
     dup   = (ur-ul)*nx + (vr-vl)*ny
     
!    eigenvalues
     h1 = abs(uv - cav)
     h2 = abs(uv)
     h4 = abs(uv + cav)
     delta =  epsentr*h4

!    entropy correction
     eabs1 = EntropyCorr2( h1,delta )
     eabs2 = EntropyCorr2( h2,delta )
     eabs4 = EntropyCorr2( h4,delta )
     
     h1 = rav*cav*dup
     h2 = eabs1*(pr-pl - h1)/(2.D0*c2a)
     h3 = eabs2*(rr-rl - (pr-pl)/c2a)
     h4 = eabs2*rav
     h5 = eabs4*(pr-pl + h1)/(2.D0*c2a)
     
     fd(1) = h2 + h3 + h5
     fd(2) = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-dup*nx) + h5*(uav+cav*nx)
     fd(3) = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-dup*ny) + h5*(vav+cav*ny)
     fd(4) = h2*(hav-cav*uv) + h3*q2a + h4*(uav*(ur-ul)+vav*(vr-vl)-uv*dup) + h5*(hav+cav*uv)

!   individual flux components of cloud components
    pav = 0.5d0*(pl + pr)
    fc(1) = 0.5d0*( ql*rl     + qr*rr     ) 
    fc(2) = 0.5d0*( ql*rl*ul  + qr*rr*ur  ) + aij(j,i)*pav 
    fc(3) = 0.5d0*( ql*rl*vl  + qr*rr*vr  ) + bij(j,i)*pav 
    fc(4) = 0.5d0*( ql*rl*hl  + qr*rr*hr  )
    
   !individual flux components of cloud components
    fi(1) = ql*rl   
    fi(2) = ql*rl*ul + aij(j,i)*pl 
    fi(3) = ql*rl*vl + bij(j,i)*pl 
    fi(4) = ql*rl*hl
    
    fij(1) = fc(1)/ds - 0.5d0*fd(1)
    fij(2) = fc(2)/ds - 0.5d0*fd(2)
    fij(3) = fc(3)/ds - 0.5d0*fd(3)
    fij(4) = fc(4)/ds - 0.5d0*fd(4)

   !residual
    rhs(1,i) = rhs(1,i) + 2.0d0*( fij(1) - fi(1)/ds )
    rhs(2,i) = rhs(2,i) + 2.0d0*( fij(2) - fi(2)/ds )
    rhs(3,i) = rhs(3,i) + 2.0d0*( fij(3) - fi(3)/ds )
    rhs(4,i) = rhs(4,i) + 2.0d0*( fij(4) - fi(4)/ds )
    
    enddo
    endif
    enddo
    
contains

    real(rtype) function valimit(a, b, epsi)
    implicit none
    real(rtype) :: a, b, epsi
    double precision num, den

    if(iorder.eq.1)then
    valimit = 0.0d0
    elseif(iorder.eq.2)then
    valimit = 1.0d0
    elseif(iorder.eq.3)then
    num     = 2.0d0*a*b + epsi
    den     = a**2 + b**2 + epsi
    valimit = dmax1(0.0d0, num/den)
    endif

    return
    end function valimit
    
    real(rtype) function valimit2(a, b, epsi)
    implicit none
    real(rtype) :: a, b, epsi
    double precision num, den,ds

    valimit2  = 1 - abs( ( a-b + epsi )/( abs(a) + abs(b) + epsi ) )**3.0d0

    return
    end function valimit2
    

  real(rtype) function EntropyCorr2( z,d )
    implicit none
    real(rtype) :: z, d

    if (z > d) then
      EntropyCorr2 = z
    else
      EntropyCorr2 = 0.5D0*(z*z+d*d)/d
    endif
  end function EntropyCorr2
    
    end subroutine residual_roeMU