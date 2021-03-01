
    Subroutine residual_cusp
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    
    integer :: i,p,j,jj
    real(rtype) :: ds,nx,ny,rrho,rl,ul,vl,pl,hl,ccl,qsl,mml,rhl,rr,ur,vr,pr,hr,ccr,qsr
    real(rtype) :: mmr,rhr,c1p2,vn1p,mn1p,rav,gam1,dd,dd1,uav,vav,hav,q2a,c2a,cav,uv
    real(rtype) :: eigp,eigm,termbeta1,betaf,alphaf,fd(4)
    real(rtype) :: ql,qr,pav,fc(4),fi(4),fij(4)
    
    do i=1,phynod 
    rhs(:,i) = 0.0d0
    if( ptype(i) .eq. 1 )then
    do j=1,nbers(i)
    
     p = conn(j,i)  
     ds = sqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
     nx = aij(j,i)/ds
     ny = bij(j,i)/ds
    
!   left state
    rrho = 1.D0/cv(1,i)
    rl   = cv(1,i)
    ul   = cv(2,i)*rrho
    vl   = cv(3,i)*rrho
    pl   = dv(1,i)
    hl   = (pl+cv(4,i))*rrho
    ccl  = dv(3,i) ! C_l
    qsl  = ( cv(2,i)*nx + cv(3,i)*ny )/cv(1,i) ! V_l
    rhl = dv(1,i) + cv(4,i)
    
!   right state
    rrho = 1.D0/cv(1,p)
    rr   = cv(1,p)
    ur   = cv(2,p)*rrho
    vr   = cv(3,p)*rrho
    pr   = dv(1,p)
    hr   = (pr+cv(4,p))*rrho
    ccr  = dv(3,p)
    qsr  = ( cv(2,p)*nx + cv(3,p)*ny )/cv(1,p)
    rhr = dv(1,p) + cv(4,p)
    
!   Interface Mach, convariant velocity
    c1p2 = 0.5d0*( ccl + ccr )
    vn1p = 0.5d0*( (ul+ur)*nx + (vl+vr)*ny )
    mn1p = vn1p/c1p2
    
!   roe averaging
    rav   = Sqrt(rl*rr)
    gam1  = 0.5D0*(dv(4,i)+dv(4,p)) - 1.D0
    dd    = rav/rl
    dd1   = 1.D0/(1.D0+dd)
    uav   = (ul+dd*ur)*dd1
    vav   = (vl+dd*vr)*dd1
    hav   = (hl+dd*hr)*dd1
    q2a   = 0.5D0*(uav*uav+vav*vav)
    c2a   = gam1*(hav-q2a)
    cav   = Sqrt(c2a)
    uv    = uav*nx + vav*ny           ! convariant Roe 
    
!   eigenvalues
    eigp = 0.5d0*uv*( gamma + 1.0d0 )/gamma + sqrt( (0.5d0*uv*(gamma + 1.0d0)/gamma)**2 + c2a/gamma )
    eigm = 0.5d0*uv*( gamma + 1.0d0 )/gamma - sqrt( (0.5d0*uv*(gamma + 1.0d0)/gamma)**2 + c2a/gamma )
    
!    eigp = ((gamma+1.0d0)/(2.0d0*gamma))*uv + sqrt( (uv*(gamma+1.0d0)/(2.0d0*gamma))**2 + (cav**2-uv**2)/gamma )
!    eigm = ((gamma+1.0d0)/(2.0d0*gamma))*uv - sqrt( (uv*(gamma+1.0d0)/(2.0d0*gamma))**2 + (cav**2-uv**2)/gamma )
    
!   factor beta
    if( mn1p >= 0.0d0 .and. mn1p < 1.0d0 )then  ! 0<M<1
     betaf = max1( 0.0d0, (vn1p + eigm)/(vn1p - eigm) ) 
    elseif(  mn1p > -1.0d0 .and. mn1p < 0.0d0 )then
     betaf = - max1( 0.0d0, (vn1p + eigp)/(vn1p - eigp) )
    elseif( abs(mn1p) >= 1.0d0 )then
     betaf = sign(1.0d0,mn1p) 
    endif 
    
!   factor alpha
    if( betaf == 0.0d0 )then
     alphaf = abs(vn1p)
    elseif( betaf > 0.0d0 .and. mn1p > 0.0d0 .and. mn1p < 1.0d0 )then
     alphaf = - ( 1.0d0 + betaf )*eigm
    elseif( betaf < 0.0d0 .and. mn1p > -1.0d0 .and. mn1p < 0.0d0 )then
     alphaf = ( 1.0d0 - betaf )*eigp
    elseif( abs( mn1p ) > 1.0d0 )then
     alphaf = 0.0d0
    endif 
    
!   CUSP dissipation term
    fd(1) = alphaf*(cv(1,p)-cv(1,i)) + betaf*( qsr*cv(1,p) - qsl*cv(1,i) )
    fd(2) = alphaf*(cv(2,p)-cv(2,i)) + betaf*( (qsr*cv(2,p)+nx*pr) - (qsl*cv(2,i)+nx*pl) )
    fd(3) = alphaf*(cv(3,p)-cv(3,i)) + betaf*( (qsr*cv(3,p)+ny*pr) - (qsl*cv(3,i)+ny*pl) )
    fd(4) = alphaf*(  rhr - rhl    ) + betaf*( qsr*rhr - qsl*rhl )
    
!    LEFT STATE
     rhl   = dv(1,i) + cv(4,i)                           ! rH 
     ql = (cv(2,i)*aij(j,i)+cv(3,i)*bij(j,i))/cv(1,i)   ! Convariant Velocity
!    RIGHT STATE
     rhr   = dv(1,p) + cv(4,p)
     qr =  (cv(2,p)*aij(j,i)+cv(3,p)*bij(j,i))/cv(1,p)
     
!    pressure term
     pav = 0.5d0*(dv(1,i) + dv(1,p))

!   individual flux components of cloud components
    fc(1) = 0.5d0*( ql*cv(1,i)  + qr*cv(1,p)  ) 
    fc(2) = 0.5d0*( ql*cv(2,i)  + qr*cv(2,p)  ) + aij(j,i)*pav 
    fc(3) = 0.5d0*( ql*cv(3,i)  + qr*cv(3,p)  ) + bij(j,i)*pav 
    fc(4) = 0.5d0*(  ql*rhl     + qr*rhr      )
    
   !individual flux components of cloud components
    fi(1) = ql*cv(1,i)   
    fi(2) = ql*cv(2,i) + aij(j,i)*dv(1,i) 
    fi(3) = ql*cv(3,i) + bij(j,i)*dv(1,i) 
    fi(4) = ql*rhl  
    
    fij(1) = fc(1) - 0.5d0*fd(1)*ds
    fij(2) = fc(2) - 0.5d0*fd(2)*ds
    fij(3) = fc(3) - 0.5d0*fd(3)*ds
    fij(4) = fc(4) - 0.5d0*fd(4)*ds
    
   !residual
    rhs(1,i) = rhs(1,i) + 2.0d0*( fij(1) - fi(1) )
    rhs(2,i) = rhs(2,i) + 2.0d0*( fij(2) - fi(2) )
    rhs(3,i) = rhs(3,i) + 2.0d0*( fij(3) - fi(3) )
    rhs(4,i) = rhs(4,i) + 2.0d0*( fij(4) - fi(4) )
    
    enddo
    endif
    enddo
    
    end subroutine residual_cusp
    

    Subroutine residual_cusp2
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    
    integer :: i,p,j,jj
    real(rtype) :: ds,nx,ny,rrho,rl,ul,vl,pl,hl,ccl,qsl,mml,rhl,rr,ur,vr,pr,hr,ccr,qsr
    real(rtype) :: mmr,rhr,c1p2,vn1p,mn1p,rav,gam1,dd,dd1,uav,vav,hav,q2a,c2a,cav,uv
    real(rtype) :: eigp,eigm,termbeta1,betaf,alphaf,fd(4)
    real(rtype) :: ql,qr,pav,fc(4),fi(4),fij(4)
    real(rtype) :: eps,dQ(4),dWirj(4),dWjrj(4),sl(4),sr(4)
    
    real(rtype) :: Fp(4),Fr(4),Fl(4)
    real(rtype) :: va,unr,unl,una,ua,sx,sy,st,rr12,rl12,rer,rel,rd,qr2,ql2,qa2,nl,mach,m1
    real(rtype) :: m0,lpos,lneg,l2,l1,ha,gg,gamma1,ct,b3,b1,ar2,al2,aa2,aa,a3,a2,a1
    real(rtype) :: dqh(4),aqh(4),b2,ggm1
    
    do i=1,phynod 
    rhs(:,i) = 0.0d0
    if( ptype(i) .eq. 1 )then
    do j=1,nbers(i)
    
     p = conn(j,i)  
     nl = sqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
     ct = aij(j,i)/nl
     st = bij(j,i)/nl
    
     eps = 1.0d-12
    
!     dQ(1) = cv(1,p) - cv(1,i)
!     dQ(2) = cv(2,p) - cv(2,i)
!     dQ(3) = cv(3,p) - cv(3,i)
!     dQ(4) = cv(4,p) - cv(4,i)
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
    
!    left state
     rl   = cv(1,i)  
     ul   = cv(2,i)/rl
     vl   = cv(3,i)/rl
     ql2  = ul**2 + vl**2
     pl   = dv(1,i)
     gam1 = dv(4,i) - 1.0d0
     ggm1 = dv(4,i)/gam1
     hl  = ggm1*pl/rl + 0.5d0*ql2

!    right state
     rr   = cv(1,p)  
     ur   = cv(2,p)/rr
     vr   = cv(3,p)/rr
     qr2  = ur**2 + vr**2
     pr   = dv(1,p)
     gam1 = dv(4,p) - 1.0d0
     ggm1 = dv(4,p)/gam1
     hr  = ggm1*pr/rr + 0.5d0*qr2
    
!    rotated velocity
     unl = ul*ct + vl*st
     unr = ur*ct + vr*st
     
!    Centred Flux
     Fl(1) = rl*unl
     Fl(2) = pl*ct + rl*ul*unl
     Fl(3) = pl*st + rr*vl*unl
     Fl(4) = rl*hl*unl
     
     Fr(1) = rr*unr
     Fr(2) = pr*ct + rr*ur*unr
     Fr(3) = pr*st + rr*vr*unr
     Fr(4) = rr*hr*unr
     
     Fp(1) = 0.0d0
     Fp(2) = (pr+pl)*ct
     Fp(3) = (pr+pl)*st
     Fp(4) = 0.0d0
     
     dqh(1) = cv(1,p) - cv(1,i)
     dqh(2) = cv(2,p) - cv(2,i)
     dqh(3) = cv(3,p) - cv(3,i)
     dqh(4) = rr*hr - rl*hl
     
     aqh(1) = 0.5d0*( cv(1,p) - cv(1,i) )
     aqh(2) = 0.5d0*( cv(2,p) - cv(2,i) )
     aqh(3) = 0.5d0*( cv(3,p) - cv(3,i) )
     aqh(4) = 0.5d0*( rr*hr + rl*hl )
     
!    Roe Averaging
     rl12 = dsqrt(rl)
     rr12 = dsqrt(rr)
     rd = 1.0d0/(rl12+rr12)
     
     ua   = (ul*rl12 + ur*rr12)*rd
     va   = (vl*rl12 + vr*rr12)*rd
     ha   = (hl*rl12 + hr*rr12)*rd
     qa2  = ua**2 + va**2
     aa2  = GAMMA1*(ha - 0.5d0*qa2)
     
     aa = dsqrt(aa2)
     una = ua*ct + va*st
     mach = una/aa/nl
    
!    Eigenvalues with entropy fix
     gg = 0.5d0*(gamma + 1.0d0 )/gamma
     l1 = (gg*una)**2 + ((aa*nl)**2 - una**2)/gamma
     l2 = dsqrt(l1)
     lpos = gg*una + l2
     lneg = gg*una - l2
     
     m0 = 0.01d0
     m1 = dabs(mach)
     if( m1 .gt. m0 )then
      alphaf = m1
     else
      alphaf = 0.5d0*( m0 + mach**2.0d0/m0 ) 
     endif
     
     if( mach .le. -1.0d0 )then
      betaf = -1.0d0
     elseif( mach .gt. -1.0d0 .and. mach .le. 0.0d0 )then
      a1 = una + lpos
      a2 = una - lpos
      a3 = a1/a2
      betaf = -dmax1(0.0d0,a3)
     else
      betaf = 1.0d0
     endif

!    Total Flux
     b1 = 0.5d0*alphaf*aa*nl
     b2 = 0.5d0*betaf*(unr-unl)
     b3 = 0.5d0*betaf
     
     do jj=1,4
      fij(jj) = 0.5d0*( Fl(jj) + Fr(jj)) - b1*dqh(jj) - b2*aqh(jj) - b3*Fp(jj)
      fi(jj)  = Fl(jj) - b3*Fp(jj)
     enddo 
    
   !residual
    rhs(1,i) = rhs(1,i) + 2.0d0*( fij(1) - fi(1) )
    rhs(2,i) = rhs(2,i) + 2.0d0*( fij(2) - fi(2) )
    rhs(3,i) = rhs(3,i) + 2.0d0*( fij(3) - fi(3) )
    rhs(4,i) = rhs(4,i) + 2.0d0*( fij(4) - fi(4) )
    
    enddo
    endif
    enddo
    
    end subroutine residual_cusp2
    
    