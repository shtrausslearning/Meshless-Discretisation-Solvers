!
!    
!    Subroutine residual_roe2
!    use ModDataTypes
!    use prmflow
!    use ModInterfaces
!    implicit none
!    
!    integer :: i,j,p,m,ic,mm,jjj,ii
!    real(rtype) :: ds,nx,ny,rrho,rl,ul,vl,pl,hl,gaml,rr,ur,vr,pr,hr,gamr
!    real(rtype) :: rav,gam1,dd,dd1,uav,hav,vav,q2a,c2a,cav,uv,du,encorf
!    real(rtype) :: h1,h2,h4,delta,eabs1,eabs2,eabs4,fd(4),h3,h5
!    real(rtype) :: rhoul,rhovl,rhl,ql,rhour,rhovr,rhr,qr,pav,fi(4),fc(4),fij(4)
!    real(rtype) :: fx(4),fy(4),ggm1
!    real(rtype) :: sx,sy, dQ(4), dWjrj(4), dWirj(4), sl(4), sr(4), eps, kkk
!    
!
!    do i=1,phynod 
!    rhs(:,i) = 0.0d0
!    if( ptype(i) .eq. 1 )then
!    do j=1,nbers(i)
!    
!     p = conn(j,i)  
!     sx = -0.5d0*(x(i) - x(p))
!     sy = 0.5d0*(y(p) - y(i))
!     
!     ds = sqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
!     nx = aij(j,i)/ds ; ny = bij(j,i)/ds
!     
!     eps = 1.0d-12
!    
!     dQ(1) = cv(1,p) - cv(1,i)
!     dQ(2) = cv(2,p)/cv(1,p) - cv(2,i)/cv(1,i)
!     dQ(3) = cv(3,p)/cv(1,p) - cv(3,i)/cv(1,i)
!     dQ(4) = dv(1,p) - dv(1,i)
!
!     dWirj(1) = sx*gradx(1,i) + sy*grady(1,i)
!     dWirj(2) = sx*gradx(2,i) + sy*grady(2,i)
!     dWirj(3) = sx*gradx(3,i) + sy*grady(3,i)
!     dWirj(4) = sx*gradx(4,i) + sy*grady(4,i)
!     
!     dWjrj(1) = sx*gradx(1,p) + sy*grady(1,p)
!     dWjrj(2) = sx*gradx(2,p) + sy*grady(2,p)
!     dWjrj(3) = sx*gradx(3,p) + sy*grady(3,p)
!     dWjrj(4) = sx*gradx(4,p) + sy*grady(4,p)
!     
!     sl(1) = ( dWirj(1)*dQ(1) + dabs( dWirj(1)*dQ(1) ) + eps )/( dWirj(1)**2 + dQ(1)**2 + eps )
!     sl(2) = ( dWirj(2)*dQ(2) + dabs( dWirj(2)*dQ(2) ) + eps )/( dWirj(2)**2 + dQ(2)**2 + eps )
!     sl(3) = ( dWirj(3)*dQ(3) + dabs( dWirj(3)*dQ(3) ) + eps )/( dWirj(3)**2 + dQ(3)**2 + eps )
!     sl(4) = ( dWirj(4)*dQ(4) + dabs( dWirj(4)*dQ(4) ) + eps )/( dWirj(4)**2 + dQ(4)**2 + eps )
!     
!     sr(1) = ( dWjrj(1)*dQ(1) + dabs( dWjrj(1)*dQ(1) ) + eps )/( dWjrj(1)**2 + dQ(1)**2 + eps )
!     sr(2) = ( dWjrj(2)*dQ(2) + dabs( dWjrj(2)*dQ(2) ) + eps )/( dWjrj(2)**2 + dQ(2)**2 + eps )
!     sr(3) = ( dWjrj(3)*dQ(3) + dabs( dWjrj(3)*dQ(3) ) + eps )/( dWjrj(3)**2 + dQ(3)**2 + eps )
!     sr(4) = ( dWjrj(4)*dQ(4) + dabs( dWjrj(4)*dQ(4) ) + eps )/( dWjrj(4)**2 + dQ(4)**2 + eps )
!     
!!    LEFT STATE DEFINITION
!     rl   = cv(1,i)                 + 0.5d0*sl(1)*dWirj(1)
!     ul   = (cv(2,i)/cv(1,i)        + 0.5d0*sl(2)*dWirj(2))
!     vl   = (cv(3,i)/cv(1,i)        + 0.5d0*sl(3)*dWirj(3) )
!     pl   = dv(1,i)                 + 0.5d0*sl(4)*dWirj(4)
!     ql = ul*aij(j,i) + vl*bij(j,i)   ! Convariant Velocity
!
!     gam1 = dv(4,i) - 1.0d0
!     ggm1 = dv(4,i)/gam1
!     hl  = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)
!     rhl = hl*rl
!     
!!    RIGHT STATE DEFINITION
!     rr   = cv(1,p)  - 0.5d0*sr(1)*dWjrj(1)
!     ur   = (cv(2,p)/cv(1,p) - 0.5d0*sr(2)*dWjrj(2))
!     vr   = (cv(3,p)/cv(1,p) - 0.5d0*sr(3)*dWjrj(3))
!     pr   = dv(1,p)  - 0.5d0*sr(4)*dWjrj(4)
!     qr =  ur*aij(j,i) + vr*bij(j,i)
!
!     gam1 = dv(4,p) - 1.0d0
!     ggm1 = dv(4,p)/gam1
!     hr  = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
!     rhr = rr*hr
!     
!!    roe averaging
!     rav  = sqrt(rl*rr)
!     gam1 = 0.5d0*(dv(4,i)+dv(4,p)) - 1.0d0
!     dd   = rav/rl
!     dd1  = 1.0d0/(1.0d0+dd)
!     uav  = (ul+dd*ur)*dd1
!     vav  = (vl+dd*vr)*dd1
!     hav  = (hl+dd*hr)*dd1
!     q2a  = 0.5d0*(uav*uav+vav*vav)
!     c2a  = gam1*(hav-q2a)
!     
!     cav  = sqrt(c2a)
!     uv   = uav*nx + vav*ny
!     du   = (ur-ul)*nx + (vr-vl)*ny
!     
!!    eigenvalues
!     h1 = abs(uv - cav)
!     h2 = abs(uv)
!     h4 = abs(uv + cav)
!     delta =  epsentr*h4
!
!!    entropy correction
!     eabs1 = EntropyCorr2( h1,delta )
!     eabs2 = EntropyCorr2( h2,delta )
!     eabs4 = EntropyCorr2( h4,delta )
!     
!     h1 = rav*cav*du
!     h2 = eabs1*(pr-pl - h1)/(2.D0*c2a)
!     h3 = eabs2*(rr-rl - (pr-pl)/c2a)
!     h4 = eabs2*rav
!     h5 = eabs4*(pr-pl + h1)/(2.D0*c2a)
!     
!     fd(1) = h2 + h3 + h5
!     fd(2) = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-du*nx) + h5*(uav+cav*nx)
!     fd(3) = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-du*ny) + h5*(vav+cav*ny)
!     fd(4) = h2*(hav-cav*uv) + h3*q2a + h4*(uav*(ur-ul)+vav*(vr-vl)-uv*du) + &
!             h5*(hav+cav*uv)
!
!
!!   individual flux components of cloud components
!    pav = 0.5d0*(pl + pr)
!    fc(1) = 0.5d0*( ql*rl     + qr*rr     ) 
!    fc(2) = 0.5d0*( ql*rl*ul  + qr*rr*ur  ) + aij(j,i)*pav 
!    fc(3) = 0.5d0*( ql*rl*vl  + qr*rr*vr  ) + bij(j,i)*pav 
!    fc(4) = 0.5d0*( ql*rhl    + qr*rhr    )
!    
!   !individual flux components of cloud components
!    fi(1) = ql*rl   
!    fi(2) = ql*rl*ul + aij(j,i)*pl 
!    fi(3) = ql*rl*vl + bij(j,i)*pl 
!    fi(4) = ql*rhl  
!    
!    fij(1) = fc(1) - 0.5d0*fd(1)*ds
!    fij(2) = fc(2) - 0.5d0*fd(2)*ds
!    fij(3) = fc(3) - 0.5d0*fd(3)*ds
!    fij(4) = fc(4) - 0.5d0*fd(4)*ds
!
!   !residual
!    rhs(1,i) = rhs(1,i) + 2.0d0*( fij(1) - fi(1) )
!    rhs(2,i) = rhs(2,i) + 2.0d0*( fij(2) - fi(2) )
!    rhs(3,i) = rhs(3,i) + 2.0d0*( fij(3) - fi(3) )
!    rhs(4,i) = rhs(4,i) + 2.0d0*( fij(4) - fi(4) )
!    
!    enddo
!    endif
!    enddo
!    
!contains
!
!    real(rtype) function valimit(a, b)
!    use ModDataTypes
!    use prmflow
!    use ModInterfaces
!    implicit none
!    real(rtype) :: a, b
!    double precision num, den, epsi
!    parameter(epsi=1.0d-15)
!
!    if(iorder.eq.1)then
!    valimit = 0.0d0
!    elseif(iorder.eq.2)then
!    valimit = 1.0d0
!    elseif(iorder.eq.3)then
!    num     = 2.0d0*a*b + epsi
!    den     = a**2 + b**2 + epsi
!    valimit = dmax1(0.0d0, num/den)
!    endif
!
!    return
!    end function valimit
!
!  real(rtype) function EntropyCorr2( z,d )
!    implicit none
!    real(rtype) :: z, d
!
!    if (z > d) then
!      EntropyCorr2 = z
!    else
!      EntropyCorr2 = 0.5D0*(z*z+d*d)/d
!    endif
!  end function EntropyCorr2
!    
!    end subroutine residual_roe2
!
!    Subroutine residual_roe2b
!    use ModDataTypes
!    use prmflow
!    use ModInterfaces
!    implicit none
!    
!    integer :: i,j,m,ic,mm,jjj,ii,p
!    real(rtype) :: ds,nx,ny,rrho,rl,ul,vl,pl,hl,gaml,rr,ur,vr,pr,hr,gamr
!    real(rtype) :: rav,gam1,dd,dd1,uav,hav,vav,q2a,c2a,cav,uv,dup,encorf
!    real(rtype) :: h1,h2,h4,delta,eabs1,eabs2,eabs4,fd(4),h3,h5
!    real(rtype) :: rhoul,rhovl,rhl,ql,rhour,rhovr,rhr,qr,pav,fi(4),fc(4),fij(4)
!    real(rtype) :: fx(4),fy(4),ggm1
!    real(rtype) :: sx,sy, dp(4), dm(4), du(4), sl(4), sr(4), eps, kkk
!    
!    do i=1,phynod 
!    rhs(:,i) = 0.0d0
!    if( ptype(i) .eq. 1 )then
!    do j=1,nbers(i)
!    
!     p = conn(j,i)  
!     sx = x(p) - x(i)
!     sy = y(p) - y(i)
!     ds = sqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
!     nx = aij(j,i)/ds
!     ny = bij(j,i)/ds
!     
!     du(1) = cv(1,p)         - cv(1,i)
!     du(2) = cv(2,p)/cv(1,p) - cv(2,i)/cv(1,i)
!     du(3) = cv(3,p)/cv(1,p) - cv(3,i)/cv(1,i)
!     du(4) = dv(1,p)         - dv(1,i)
!     dm(1) = 2.0d0*(sx*gradx(1,i) + sy*grady(1,i)) - du(1)
!     dm(2) = 2.0d0*(sx*gradx(2,i) + sy*grady(2,i)) - du(2)
!     dm(3) = 2.0d0*(sx*gradx(3,i) + sy*grady(3,i)) - du(3)
!     dm(4) = 2.0d0*(sx*gradx(4,i) + sy*grady(4,i)) - du(4)
!     dp(1) = 2.0d0*(sx*gradx(1,p) + sy*grady(1,p)) - du(1)
!     dp(2) = 2.0d0*(sx*gradx(2,p) + sy*grady(2,p)) - du(2)
!     dp(3) = 2.0d0*(sx*gradx(3,p) + sy*grady(3,p)) - du(3)
!     dp(4) = 2.0d0*(sx*gradx(4,p) + sy*grady(4,p)) - du(4)
!     
!     sl(1) = valimit(dm(1),du(1))
!     sl(2) = valimit(dm(2),du(2))
!     sl(3) = valimit(dm(3),du(3))
!     sl(4) = valimit(dm(4),du(4))
!     
!     sr(1) = valimit(dp(1),du(1))
!     sr(2) = valimit(dp(2),du(2))
!     sr(3) = valimit(dp(3),du(3))
!     sr(4) = valimit(dp(4),du(4))
!
!
!     kkk = 1.0d0/3.0d0
!!!     kkk = 1.0d0
!!     kkk = 0.0d0
!     
!!    LEFT STATE DEFINITION
!     rl   = cv(1,i)         + 0.25d0*sl(1)*( (1.0d0-kkk*sl(1))*dm(1) + (1.0d0+kkk*sl(1))*du(1) )
!     ul   = cv(2,i)/cv(1,i) + 0.25d0*sl(2)*( (1.0d0-kkk*sl(2))*dm(2) + (1.0d0+kkk*sl(2))*du(2) )
!     vl   = cv(3,i)/cv(1,i) + 0.25d0*sl(3)*( (1.0d0-kkk*sl(3))*dm(3) + (1.0d0+kkk*sl(3))*du(3) )
!     pl   = dv(1,i)         + 0.25d0*sl(4)*( (1.0d0-kkk*sl(4))*dm(4) + (1.0d0+kkk*sl(4))*du(4) )
!     ql = ul*aij(j,i) + vl*bij(j,i)   ! Convariant Velocity
!
!     gam1 = dv(4,i) - 1.0d0
!     ggm1 = dv(4,i)/gam1
!     hl  = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)
!     rhl = hl*rl
!     
!!    RIGHT STATE DEFINITION
!     rr   = cv(1,p)           - 0.25d0*sr(1)*( (1.0d0-kkk*sr(1))*dp(1) + (1.0d0+kkk*sr(1))*du(1))
!     ur   = cv(2,p)/cv(1,p)   - 0.25d0*sr(2)*( (1.0d0-kkk*sr(2))*dp(2) + (1.0d0+kkk*sr(2))*du(2))
!     vr   = cv(3,p)/cv(1,p)   - 0.25d0*sr(3)*( (1.0d0-kkk*sr(3))*dp(3) + (1.0d0+kkk*sr(3))*du(3))
!     pr   = dv(1,p)           - 0.25d0*sr(4)*( (1.0d0-kkk*sr(4))*dp(4) + (1.0d0+kkk*sr(4))*du(4))
!     qr =  ur*aij(j,i)+ vr*bij(j,i)
!
!     gam1 = dv(4,p) - 1.0d0
!     ggm1 = dv(4,p)/gam1
!     hr  = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
!     rhr = rr*hr
!     
!!    roe averaging
!     rav  = sqrt(rl*rr)
!     gam1 = 0.5d0*(dv(4,i)+dv(4,p)) - 1.0d0
!     dd   = rav/rl
!     dd1  = 1.0d0/(1.0d0+dd)
!     uav  = (ul+dd*ur)*dd1
!     vav  = (vl+dd*vr)*dd1
!     hav  = (hl+dd*hr)*dd1
!     q2a  = 0.5d0*(uav*uav+vav*vav)
!     c2a  = gam1*(hav-q2a)
!     
!     cav  = sqrt(c2a)
!     uv   = uav*nx + vav*ny
!     dup   = (ur-ul)*nx + (vr-vl)*ny
!     
!!    eigenvalues
!     h1 = abs(uv - cav)
!     h2 = abs(uv)
!     h4 = abs(uv + cav)
!     delta =  epsentr*h4
!
!!    entropy correction
!     eabs1 = EntropyCorr2( h1,delta )
!     eabs2 = EntropyCorr2( h2,delta )
!     eabs4 = EntropyCorr2( h4,delta )
!     
!     h1 = rav*cav*dup
!     h2 = eabs1*(pr-pl - h1)/(2.D0*c2a)
!     h3 = eabs2*(rr-rl - (pr-pl)/c2a)
!     h4 = eabs2*rav
!     h5 = eabs4*(pr-pl + h1)/(2.D0*c2a)
!     
!     fd(1) = h2 + h3 + h5
!     fd(2) = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-dup*nx) + h5*(uav+cav*nx)
!     fd(3) = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-dup*ny) + h5*(vav+cav*ny)
!     fd(4) = h2*(hav-cav*uv) + h3*q2a + h4*(uav*(ur-ul)+vav*(vr-vl)-uv*dup) + h5*(hav+cav*uv)
!
!!   individual flux components of cloud components
!    pav = 0.5d0*(pl + pr)
!    fc(1) = 0.5d0*( ql*rl     + qr*rr     ) 
!    fc(2) = 0.5d0*( ql*rl*ul  + qr*rr*ur  ) + aij(j,i)*pav 
!    fc(3) = 0.5d0*( ql*rl*vl  + qr*rr*vr  ) + bij(j,i)*pav 
!    fc(4) = 0.5d0*( ql*rhl    + qr*rhr    )
!    
!   !individual flux components of cloud components
!    fi(1) = ql*rl   
!    fi(2) = ql*rl*ul + aij(j,i)*pl 
!    fi(3) = ql*rl*vl + bij(j,i)*pl 
!    fi(4) = ql*rhl  
!    
!    fij(1) = fc(1) - 0.5d0*fd(1)*ds
!    fij(2) = fc(2) - 0.5d0*fd(2)*ds
!    fij(3) = fc(3) - 0.5d0*fd(3)*ds
!    fij(4) = fc(4) - 0.5d0*fd(4)*ds
!
!   !residual
!    rhs(1,i) = rhs(1,i) + 2.0d0*( fij(1) - fi(1) )
!    rhs(2,i) = rhs(2,i) + 2.0d0*( fij(2) - fi(2) )
!    rhs(3,i) = rhs(3,i) + 2.0d0*( fij(3) - fi(3) )
!    rhs(4,i) = rhs(4,i) + 2.0d0*( fij(4) - fi(4) )
!    !
!    !do ii=1,4
!    !if( isnan(rhs(ii,i)) )then
!    !write(*,*) i
!    !pause
!    !exit
!    !endif
!    !enddo
!    
!    enddo
!    endif
!    enddo
!    
!contains
!
!    real(rtype) function valimit(a, b)
!    use ModDataTypes
!    use prmflow
!    use ModInterfaces
!    implicit none
!    real(rtype) :: a, b
!    double precision num, den, epsi
!    parameter(epsi=1.0d-15)
!
!    if(iorder.eq.1)then
!    valimit = 0.0d0
!    elseif(iorder.eq.2)then
!    valimit = 1.0d0
!    elseif(iorder.eq.3)then
!    num     = 2.0d0*a*b + epsi
!    den     = a**2 + b**2 + epsi
!    valimit = dmax1(0.0d0, num/den)
!    endif
!
!    return
!    end function valimit
!    
!    real(rtype) function valimit2(a, b, ds)
!    use ModDataTypes
!    use prmflow
!    use ModInterfaces
!    implicit none
!    real(rtype) :: a, b
!    double precision num, den, epsi,ds
!    parameter(epsi=1.0d-15)
!
!    valimit2  = 1 - abs( ( a-b )/( abs(a) + abs(b) + epsi*ds ) )**3.0d0
!
!    return
!    end function valimit2
!    
!
!  real(rtype) function EntropyCorr2( z,d )
!    implicit none
!    real(rtype) :: z, d
!
!    if (z > d) then
!      EntropyCorr2 = z
!    else
!      EntropyCorr2 = 0.5D0*(z*z+d*d)/d
!    endif
!  end function EntropyCorr2
!    
!    end subroutine residual_roe2b
!