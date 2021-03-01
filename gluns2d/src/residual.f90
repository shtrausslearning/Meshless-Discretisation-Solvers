! 
!    Subroutine residual_roe2
!    use ModDataTypes
!    use prmflow
!    use ModInterfaces
!    implicit none
!    
!    integer :: i,j,p,m,ic,mm,jjj
!    real(rtype) :: ds,nx,ny,rrho,rl,ul,vl,pl,hl,gaml,rr,ur,vr,pr,hr,gamr
!    real(rtype) :: rav,gam1,dd,dd1,uav,hav,vav,q2a,c2a,cav,uv,encorf
!    real(rtype) :: h1,h2,h4,delta,eabs1,eabs2,eabs4,fd(4),h3,h5
!    real(rtype) :: rhoul,rhovl,rhl,ql,rhour,rhovr,rhr,qr,pav,fi(4),fc(4),fij(4)
!    real(rtype) :: fx(4),fy(4)
!    real(rtype) :: sx,sy,du(4),dm(4),dp(4),sl(4),sr(4),eps,kkk,ggm1,qsrl,qsrr,dub,qsl,qsr
!
!    do i=1,phynod 
!    rhs(:,i) = 0.0d0
!    if( ptype(i) .eq. 1 )then
!    do j=1,nbers(i)
!    
!    p = conn(j,i)  
!    ds = sqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
!    nx = aij(j,i)/ds
!    ny = bij(j,i)/ds
!    sx = x(p) - x(i)
!    sy = y(p) - y(i)
!    
!    du(1) = cv(1,p) - cv(1,i)
!    du(2) = cv(2,p)/cv(1,p) - cv(2,i)/cv(1,i)
!    du(3) = cv(3,p)/cv(1,p) - cv(3,i)/cv(1,i)
!    du(4) = dv(1,p) - dv(1,i)
!!    
!    if( p .le. phynod )then
!    do m=1,4
!    dm(m) = 2.0d0*(sx*gradx(m,i) + sy*grady(m,i)) - du(m)
!    dp(m) = 2.0d0*(sx*gradx(m,j) + sy*grady(m,p)) - du(m)
!    sl(m) = dmax1( 0.0d0, (2.0d0*dm(m)*du(m) + eps)/(dm(m)**2 + du(m)**2 + eps) )
!    sr(m) = dmax1( 0.0d0, (2.0d0*dp(m)*du(m) + eps)/(dp(m)**2 + du(m)**2 + eps) )
!    enddo
!    else
!    dm(:) = 0.0d0 ; dp(:) = 0.0d0 ; sl(:) = 0.0d0 ; sr(:) = 0.0d0 
!    endif
!
!!   LEFT STATE DEFINITION
!!   #########################################################
!    !rl   = cv(1,i)      
!    !ul   = cv(2,i)/cv(1,i)
!    !vl   = cv(3,i)/cv(1,i)
!    !pl   = dv(1,i)
!    !hl   = (pl+cv(4,i))/cv(1,i)
!    !gaml = dv(4,i)
!!   #########################################################
!    rl   = cv(1,i)         
!    ul   = cv(2,i)/cv(1,i) 
!    vl   = cv(3,i)/cv(1,i)
!    pl   = dv(1,i)        
!    gam1 = dv(4,i) - 1.D0
!    ggm1 = dv(4,i)/gam1
!    hl   = ggm1*pl/rl + 0.5D0*(ul*ul+vl*vl)
!    qsl = (ul*aij(j,i)+vl*bij(j,i))
!
!!   RIGHT STATE DEFINITION
!!   ######################################################### 
!    !rr   = cv(1,p)        
!    !ur   = cv(2,p)/cv(1,p) 
!    !vr   = cv(3,p)/cv(1,p) 
!    !pr   = dv(1,p)      
!    !hr   = (pr+cv(4,p))/cv(1,p)
!    !gamr = dv(4,p)
!!   #########################################################
!    rr   = cv(1,p)         
!    ur   = cv(2,p)/cv(1,p) 
!    vr   = cv(3,p)/cv(1,p) 
!    pr   = dv(1,p)         
!    gam1 = dv(4,p) - 1.D0
!    ggm1 = dv(4,p)/gam1
!    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
!    qsr = (ur*aij(j,i)+vr*bij(j,i))
!     
!!   roe averaging
!    rav  = sqrt(rl*rr)
!    gam1 = 0.5d0*(dv(4,i)+dv(4,j)) - 1.0d0
!    dd   = rav/rl
!    dd1  = 1.0d0/(1.0d0+dd)
!    uav  = (ul+dd*ur)*dd1
!    vav  = (vl+dd*vr)*dd1
!    hav  = (hl+dd*hr)*dd1
!    q2a  = 0.5d0*(uav*uav+vav*vav)
!    c2a  = gam1*(hav-q2a)
!     
!    cav  = sqrt(c2a)
!    uv   = uav*nx + vav*ny
!    dub   = (ur-ul)*nx + (vr-vl)*ny
!     
!!   eigenvalues
!    h1 = abs(uv - cav)
!    h2 = abs(uv)
!    h4 = abs(uv + cav)
!    delta =  0.001d0*h4
!     
!    if(h1>delta) then
!    eabs1 = h1
!    else
!    eabs1 = 0.5D0*(h1*h1+delta*delta)/delta
!    endif
!    if(h2>delta)then
!    eabs2 = h2 
!    else
!    eabs2 = 0.5D0*(h2*h2+delta*delta)/delta
!    endif
!    if(h4>delta)then
!    eabs4 = h4
!    else
!    eabs4 = 0.5D0*(h4*h4+delta*delta)/delta
!    endif
!     
!    h1 = rav*cav*dub
!    h2 = eabs1*(pr-pl - h1)/(2.D0*c2a)
!    h3 = eabs2*(rr-rl - (pr-pl)/c2a)
!    h4 = eabs2*rav
!    h5 = eabs4*(pr-pl + h1)/(2.D0*c2a)
!     
!    fd(1) = h2 + h3 + h5
!    fd(2) = h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-dub*nx) + h5*(uav+cav*nx)
!    fd(3) = h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-dub*ny) + h5*(vav+cav*ny)
!    fd(4) = h2*(hav-cav*uv) + h3*q2a + h4*(uav*(ur-ul)+vav*(vr-vl)-uv*dub) + &
!            h5*(hav+cav*uv)
!            
!!   Convective Part
!    pav   = 0.5D0*(pl+pr)
!    fc(1) = 0.5D0*(qsl*cv(1,i)+qsr*cv(1,p))
!    fc(2) = 0.5D0*(qsl*cv(2,i)+qsr*cv(2,p)) + aij(j,i)*pav
!    fc(3) = 0.5D0*(qsl*cv(3,i)+qsr*cv(3,p)) + bij(j,i)*pav
!    fc(4) = 0.5D0*(qsl*hl*rl    +qsr*hr*rr    )
!    
!   !individual flux components of cloud components
!    fi(1) = qsrl   
!    fi(2) = qsrl*ul + aij(j,i)*pav
!    fi(3) = qsrl*vl + bij(j,i)*pav
!    fi(4) = qsrl*hl
!    
!    fij(1) = fc(1)/ds - 0.50d0*fd(1)
!    fij(2) = fc(2)/ds - 0.50d0*fd(2)
!    fij(3) = fc(3)/ds - 0.50d0*fd(3)
!    fij(4) = fc(4)/ds - 0.50d0*fd(4)
!
!   !residual
!    rhs(1,i) = rhs(1,i) + 2.0d0*( fij(1) - fi(1)/ds )
!    rhs(2,i) = rhs(2,i) + 2.0d0*( fij(2) - fi(2)/ds )
!    rhs(3,i) = rhs(3,i) + 2.0d0*( fij(3) - fi(3)/ds )
!    rhs(4,i) = rhs(4,i) + 2.0d0*( fij(4) - fi(4)/ds )
!    
!    enddo
!    endif
!    enddo
!
!    
!    do i=1,phynod
!    if( rhs(1,i) .eq. -777.0d0 )then
!    write(*,*) i
!    write(*,*) ntype(i)
!    write(*,*) nbers(i)
!    write(*,*) ( conn(j,i),j=1,nbers(i) )
!    pause '                   '
!    endif
!    enddo
!    
!    end subroutine residual_roe2
!    
!!    Subroutine residual2
!!    use ModDataTypes
!!    use ModGeometry
!!    use ModNumerics
!!    use ModPhysics
!!    use ModInterfaces
!!    implicit none
!!    
!!    integer :: i,j,p
!!    real(rtype) :: ds,nx,ny,shigma,Ku,Kp,rrho,gam1,ggm1,rl,ul,vl,pl,ql2,hl,cl2
!!    real(rtype) :: rr,ur,vr,pr,qr2,hr,cr2,cbar,Coff
!!    real(rtype) :: VnLFT,VnRHT,sMLFT,sMRHT,Mbar2,Mbar,Mo2,Mo,faMo,alphaco,fppLa
!!    real(rtype) :: fpmRa,fMpLb,fMmRb,R1over2,Mp,M1over2,Pu,PT,sm,fij(4),fx(4),fy(4)
!!    
!!    do i=1,phynod
!!    rhs(1,i) = -777.0d0
!!    enddo
!!
!!    do i=1,phynod 
!!    rhs(:,i) = 0.0d0
!!    
!!!   Fx,Fy flux at i
!!    fx(1) = cv(2,i)
!!    fx(2) = cv(2,i)*(cv(2,i)/cv(1,i)) + dv(1,i)
!!    fx(3) = cv(2,i)*(cv(3,i)/cv(1,i)) 
!!    fx(4) = (dv(1,i) + cv(4,i))*(cv(2,i)/cv(1,i))
!!    
!!    fy(1) = cv(3,i)
!!    fy(2) = cv(3,i)*(cv(2,i)/cv(1,i))
!!    fy(3) = cv(3,i)*(cv(3,i)/cv(1,i)) + dv(1,i)
!!    fy(4) = (dv(1,i) + cv(4,i))*(cv(3,i)/cv(1,i))
!!    
!!    if( ptype(i) .eq. 1 )then
!!    do j=1,nbers(i)
!!    
!!     p = conn(j,i)                                        
!!     ds = sqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
!!     nx = aij(j,i)/ds
!!     ny = bij(j,i)/ds
!!     
!!!   AUSM CONSTANTS [ VARIABLE ]
!!    shigma = 1.0d0
!!    Ku = 0.75d0
!!    Kp = 0.25d0
!!    
!!!   left state 
!!    rrho = 1.D0/cv(1,i)
!!    gam1 = dv(4,i) - 1.0d0
!!    ggm1 = dv(4,i)/gam1
!!    
!!    rl   = cv(1,i)
!!    ul   = cv(2,i)*rrho
!!    vl   = cv(3,i)*rrho
!!    pl   = dv(1,i)
!!    ql2  = ul*ul + vl*vl
!!    hl   = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)
!!    cl2  = dv(3,i)**2
!!    
!!!   rights state
!!    rrho = 1.D0/cv(1,p)
!!    gam1 = dv(4,p) - 1.0d0
!!    ggm1 = dv(4,p)/gam1
!!    
!!    rr   = cv(1,p)
!!    ur   = cv(2,p)*rrho
!!    vr   = cv(3,p)*rrho
!!    pr   = dv(1,p)
!!    qr2   = ur*ur + vr*vr
!!    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
!!    cr2  = dv(3,p)**2
!!    
!!!   c_bar, average speed of sound
!!    cbar = 0.5d0 * ( sqrt(cl2) + sqrt(cr2) ) !SI
!!    
!!!   %%%%%%%%%%%%%%%%%%%%%%
!!!   SETTING AUSM VARIABLES
!!!   %%%%%%%%%%%%%%%%%%%%%%
!!    
!!!   Cutoff Mach ( need to set )
!!    Coff = 1.0d0 ! decide
!!    
!!!   lhs,rhs convariant velocity 
!!    VnLFT = nx*ul + ny*vl !SI
!!    VnRHT = nx*ur + ny*vr 
!!!   lhs,rhs Mach states
!!    sMLFT = VnLFT/cbar !Ml
!!    sMRHT = VnRHT/cbar !Mr
!!    
!!!   added
!!    Mbar2 = 0.5d0/(cbar**2) * (VnLFT**2 + VnRHT**2)
!!    Mbar = DSQRT(Mbar2)
!!    Mo2 = dmin1(1.0d0,dmax1(Mbar2,(Coff)**2))
!!     
!!    Mo = DSQRT(Mo2)
!!    faMo = Mo*(2.0d0-Mo)
!!    alphaco = 0.1875d0*(-4.0d0 + 5.0d0*(faMo**2))
!!    
!!!    F +/-_p | ALPHA FUNCTION
!!     IF( DABS(sMLFT) >= 1.0d0 )THEN
!!       fppLa = 0.5d0 * ( 1.0d0 + dmax1(-1.0d0,dmin1(1.0d0,sMLFT) ))
!!     ELSE
!!      fppLa = 0.25d0*(sMLFT+1.0d0)**2*(2.0d0-sMLFT) + alphaco*sMLFT*(sMLFT**2 - 1.0d0)**2
!!     ENDIF 
!!     IF( DABS(sMRHT) >= 1.0d0 )THEN
!!      fpmRa = 0.5d0 * ( 1.0d0 - dmax1(-1.0d0,dmin1(1.0d0,sMRHT)) )
!!     ELSE
!!      fpmRa = 0.25d0*(sMRHT-1.0d0)**2*(2.0d0+sMRHT) - alphaco*sMRHT*(sMRHT**2 - 1.0d0)**2
!!     ENDIF  
!!    
!!!    F+/-_M | BETA FUNCTION
!!     IF( DABS(sMLFT) >= 1.0d0 )THEN
!!       fMpLb = 0.5d0*(sMLFT+dabs(sMLFT))
!!     ELSE
!!       fMpLb = 0.25d0*(sMLFT+1.0d0)**2 + 0.125d0*(sMLFT**2 - 1.0d0)**2
!!     ENDIF
!!     IF( DABS(sMRHT) >= 1.0d0 )THEN
!!       fMmRb = 0.5d0*(sMRHT-dabs(sMRHT))
!!     ELSE
!!       fMmRb = - 0.25d0*(sMRHT-1.0d0)**2 - 0.125d0*(sMRHT**2 - 1.0d0)**2
!!     ENDIF  
!!     
!!     R1over2 = 0.5d0*(rl+rr)
!!     Mp = - Kp/faMo*dmax1((1.0d0-shigma*Mbar2),0.0d0)*(pr-pl)/R1over2/(cbar**2)
!!     M1over2 = fMpLb + fMmRb + Mp
!!
!!     Pu = - Ku * fppLa * fpmRa * (rl+rr)*(faMo*cbar)*(VnRHT-VnLFT)
!!    
!!!    PRESSURE FLUX FUNCTION
!!     PT = fppLa*pl + fpmRa*pr + Pu
!!
!!!    MASS FLUX FUNCTION
!!     IF( M1over2 > 0.0d0 )THEN
!!      sm = M1over2 * cbar * rl
!!     ELSE
!!      sm = M1over2 * cbar * rr
!!     ENDIF 
!!     
!!!   AUSM+up flux
!!    fij(1) = (0.5d0*(sm+abs(sm))    +0.5d0*(sm-abs(sm))           )
!!    fij(2) = (0.5d0*(sm+abs(sm))*ul +0.5d0*(sm-abs(sm))*ur+ PT*nx )
!!    fij(3) = (0.5d0*(sm+abs(sm))*vl +0.5d0*(sm-abs(sm))*vr+ PT*ny )
!!    fij(4) = (0.5d0*(sm+abs(sm))*hl +0.5d0*(sm-abs(sm))*hr        )
!!    
!!!   residual
!!    rhs(1,i) = rhs(1,i) + 2.0d0*( fij(1) - aij(j,i)*fx(1) - bij(j,i)*fy(1) )
!!    rhs(2,i) = rhs(2,i) + 2.0d0*( fij(2) - aij(j,i)*fx(2) - bij(j,i)*fy(2) )
!!    rhs(3,i) = rhs(3,i) + 2.0d0*( fij(3) - aij(j,i)*fx(3) - bij(j,i)*fy(3) )
!!    rhs(4,i) = rhs(4,i) + 2.0d0*( fij(4) - aij(j,i)*fx(4) - bij(j,i)*fy(4) )    
!!    enddo
!!    endif
!!    enddo
!!    
!!    do i=1,phynod
!!    if( rhs(1,i) .eq. -777.0d0 )then
!!    write(*,*) i
!!    write(*,*) ntype(i)
!!    write(*,*) nbers(i)
!!    write(*,*) ( conn(j,i),j=1,nbers(i) )
!!    pause '                   '
!!    endif
!!    enddo
!!    
!!    return
!!    end subroutine residual2
!!    
!    Subroutine residual3
!    use ModDataTypes
!    use prmflow
!    use ModInterfaces
!    implicit none
!    
!    integer :: i,j,p,m,ic,mm,jjj
!    real(rtype) :: ds,nx,ny,rrho,rl,ul,vl,pl,hl,gaml,rr,ur,vr,pr,hr,gamr
!    real(rtype) :: rav,gam1,dd,dd1,uav,hav,vav,q2a,c2a,cav,uv,du,encorf
!    real(rtype) :: h1,h2,h4,delta,eabs1,eabs2,eabs4,fd(4),h3,h5
!    real(rtype) :: rhoul,rhovl,rhl,ql,rhour,rhovr,rhr,qr,pav,fi(4),fc(4),fij(4)
!    real(rtype) :: fx(4),fy(4)
!        
!    !open(23,file='sort_residual.dat',form='formatted')
!    !do i=1,phynod
!    !write(23,*) ( conn(j,i),j=1,nbers(i) )
!    !enddo
!    !close(23)
!    
!    !open(23,file='sort_residual.dat',form='formatted')
!    !do j=1,nbers(148)
!    !ic=conn(j,148)
!    !write(23,*) x(ic),y(ic)
!    !enddo
!    !close(22)
!    
!    do i=1,phynod
!    rhs(1,i) = -777.0d0
!    enddo
!
!    do i=1,phynod 
!    rhs(:,i) = 0.0d0
!    
!    if( ptype(i) .eq. 1 )then
!    do j=1,nbers(i)
!    
!     p = conn(j,i)  
!     ds = sqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
!     nx = aij(j,i)/ds
!     ny = bij(j,i)/ds
!
!!    LEFT STATE DEFINITION
!     rl   = cv(1,i)
!     ul   = cv(2,i)/cv(1,i)
!     vl   = cv(3,i)/cv(1,i)
!     pl   = dv(1,i)
!     hl   = (pl+cv(4,i))/cv(1,i)
!     gaml = dv(4,i)
!
!!    RIGHT STATE DEFINITION
!     rr   = cv(1,p)
!     ur   = cv(2,p)/cv(1,p)
!     vr   = cv(3,p)/cv(1,p)
!     pr   = dv(1,p)
!     hr   = (pr+cv(4,p))/cv(1,p)
!     gamr = dv(4,p)
!     
!     
!!    roe averaging
!     rav  = sqrt(rl*rr)
!     gam1 = 0.5d0*(gaml+gamr) - 1.0d0
!     dd   = rav/rl
!     dd1  = 1.0d0/(1.0d0+dd)
!     uav  = (ul+dd*ur)*dd1
!     vav  = (vl+dd*vr)*dd1
!     hav  = (hl+dd*hr)*dd1
!     q2a  = 0.5d0*(uav*uav+vav*vav)
!     c2a  = gaml*(hav-q2a)
!     !c2a  = gaml*abs(hav-q2a)
!     
!     if( c2a .gt. 0.0d0 )then
!     cav  = sqrt(c2a)
!     else
!     pause 'sonic speed is negative'
!     endif
!     
!     uv   = uav*nx + vav*ny
!     du   = (ur-ul)*nx + (vr-vl)*ny
!     
!!    eigenvalues
!     h1 = dmin1(uv - cav, 0.0d0)
!     h2 = dmin1(uv      , 0.0d0)
!     h4 = dmin1(uv + cav, 0.0d0)
!
!     eabs1 = h1
!     eabs2 = h2 
!     eabs4 = h4
!     
!!    diffusive flux
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
!     do m=1,4
!     if( isnan( fd(m) ))then
!     fd(m) = 0.0d0
!     endif
!     enddo
!             
!!    LEFT STATE
!     rhl   = dv(1,i) + cv(4,i)                           ! rH 
!     ql = (cv(2,i)*aij(j,i)+cv(3,i)*bij(j,i))/cv(1,i)   ! Convariant Velocity
!!    RIGHT STATE
!     rhr   = dv(1,p) + cv(4,p)
!     qr =  (cv(2,p)*aij(j,i)+cv(3,p)*bij(j,i))/cv(1,p)
!        
!!   individual flux components of cloud components
!    fi(1) = ql*cv(1,i)   
!    fi(2) = ql*cv(2,i) + aij(j,i)*dv(1,i) 
!    fi(3) = ql*cv(3,i) + bij(j,i)*dv(1,i) 
!    fi(4) = ql*rhl  
!    
!    fij(1) = fi(1) - fd(1)
!    fij(2) = fi(2) - fd(2)
!    fij(3) = fi(3) - fd(3)
!    fij(4) = fi(4) - fd(4)
!    !fij(1) = fc(1)
!    !fij(2) = fc(2)
!    !fij(3) = fc(3)
!    !fij(4) = fc(4)
!    
!!   residual
!    rhs(1,i) = rhs(1,i) + fij(1)
!    rhs(2,i) = rhs(2,i) + fij(2)
!    rhs(3,i) = rhs(3,i) + fij(3)
!    rhs(4,i) = rhs(4,i) + fij(4)
!    
!    do m=1,4
!    if( isnan( rhs(m,i) ))then
!    write(*,*) (fd(jjj),jjj=1,4)
!    write(*,*) ''
!    write(*,*) gaml,hav,q2a
!    pause
!    endif
!    enddo
!    
!    enddo
!    endif
!    enddo
!
!    
!    do i=1,phynod
!    if( rhs(1,i) .eq. -777.0d0 )then
!    write(*,*) i
!    write(*,*) ntype(i)
!    write(*,*) nbers(i)
!    write(*,*) ( conn(j,i),j=1,nbers(i) )
!    pause '                   '
!    endif
!    enddo
!    
!    end subroutine residual3
!
!     
!     
!    
!    
!    