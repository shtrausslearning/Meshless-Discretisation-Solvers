
    subroutine residual_slau
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    
    integer :: i,j,p
    real(rtype) :: ds,nx,ny
    real(rtype) :: rrho,rl,ul,vl,pl,hl,ccl,qsl,rhl
    real(rtype) :: rr,ur,vr,pr,hr,ccr,qsr,rhr
    real(rtype) :: fc(4),fij(4),fi(4)
    
    real(rtype) :: mhat,Xai,cbar,Ml,Mr,ppl,pml,ppr,pmr,ptau
    real(rtype) :: ABSVnB,gf,VnBP,mf,VnBm,VnB
    
    real(rtype) :: VnRHTBARa,VnRHTa,VnRHT,VnLFTBARa,VnLFTa,VnLFT,VnBARa,sMRHT2
    real(rtype) :: sMRHT,sMLFT2,sMLFT,sMhat,sm,skai,qr2,ql2,PT,ggm1,gam1,G,cr2,cl2
    real(rtype) :: brht,blft
    
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
    gam1 = dv(4,i) - 1.0d0
    ggm1 = dv(4,i)/gam1
    
    rl   = cv(1,i)
    ul   = cv(2,i)*rrho
    vl   = cv(3,i)*rrho
    pl   = dv(1,i)
    ql2  = ul*ul + vl*vl
    hl   = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)
    cl2  = dv(3,i)**2

!   rights state
    rrho = 1.D0/cv(1,p)
    gam1 = dv(4,p) - 1.0d0
    ggm1 = dv(4,p)/gam1
    
    rr   = cv(1,p)
    ur   = cv(2,p)*rrho
    vr   = cv(3,p)*rrho
    pr   = dv(1,p)
    qr2   = ur*ur + vr*vr
    hr   = ggm1*pr/rr + 0.5D0*(ur*ur+vr*vr)
    cr2  = dv(3,p)**2
    
!   c_bar, average speed of sound
    cbar = 0.5d0 * ( dsqrt(cl2) + dsqrt(cr2) ) !SI
    
    sMhat = dmin1(1.0d0,dsqrt(0.5d0 *( ( ul*ul + vl*vl ) + (ur*ur + vr*vr )))/cbar)

!   %%%%%%%%%%%%%%%%%%%%%%
!   SETTING SLAU VARIABLES
!   %%%%%%%%%%%%%%%%%%%%%%

!   lhs,rhs convariant velocity 
    VnLFT = nx*ul + ny*vl !SI
    VnRHT = nx*ur + ny*vr 
!   lhs,rhs Mach states
    sMLFT = VnLFT/cbar !Ml
    sMRHT = VnRHT/cbar !Mr
    
    sMLFT2 = (VnLFT)**0.5d0/cbar
    sMRHT2 = (VnRHT)**0.5d0/cbar 
    
    skai = (1.0d0-sMhat)**2.0d0
    
    VnLFTa = abs(VnLFT) 
    VnRHTa = abs(VnRHT)
!   Vn_bar
    VnBARa = ( rl*abs(VnLFT) + rr*abs(VnRHT) )/(rl+rr)
    
    G = -max(min(sMLFT,0.0d0),-1.0d0)*min(max(VnRHT,0.0d0),1.0d0)
!    G = min(max(0.0d0,G),1.0d0)
    
    VnLFTBARa = (1.0d0-G)*VnBARa + G*VnLFTa
    VnRHTBARa = (1.0d0-G)*VnBARa + G*VnRHTa
    
!   MASS FLUX FUNCTION
    sm = 0.5d0*( rl*(VnLFT + VnLFTBARa) + rr*(VnRHT - VnRHTBARa) - skai*(pl - pr)/cbar )
    
!   P+/P- w/ alpha = 0
    IF( abs( sMLFT ) >= 1.0d0 )THEN
     blft = 0.5d0 * ( 1.0d0 + max(-1.0d0,min(1.0d0,sMLFT) ))
    ELSE
     blft = 0.25d0 * ((sMLFT + 1.0d0)**2) * (2.0d0 - sMLFT)  ! P+|a=0
    ENDIF  
      
    IF( abs( sMRHT ) >= 1.0d0 )THEN
     brht = 0.5d0 * ( 1.0d0 - max(-1.0d0,min(1.0d0,sMRHT)) )
    ELSE
     brht = 0.25d0 * ((sMRHT - 1.0d0)**2) * (2.0d0 + sMRHT)  ! P-|a=0
    ENDIF    
    
!   SLAU1 PRESSURE FLUX FUNCTION
    PT = 0.5d0*(pl  + pr) &
     + 0.5d0*(blft  - brht)*(pl - pr) &
     + 0.5d0*(1.0d0 - skai)*(blft + brht - 1.0d0)*(pl + pr)
   !SLAU2 PRESSURE FLUX FUNCTION
!    PT = 0.5d0*(pl  + pr)          &                             
!      + 0.5d0*(blft  - brht)*(pl - pr)   &                      
!      + SQRT(0.5d0*(ql2+qr2))*(blft+brht-1.0d0)*(0.5d0*(rl+rr))*cbar
     
!   SLAU flux
    fij(1) = (0.5d0*(sm+abs(sm))    +0.5d0*(sm-abs(sm))           )
    fij(2) = (0.5d0*(sm+abs(sm))*ul +0.5d0*(sm-abs(sm))*ur+ PT*nx )
    fij(3) = (0.5d0*(sm+abs(sm))*vl +0.5d0*(sm-abs(sm))*vr+ PT*ny )
    fij(4) = (0.5d0*(sm+abs(sm))*hl +0.5d0*(sm-abs(sm))*hr        )
    
!   SLAU flux
    fi(1) = (0.5d0*(sm+abs(sm))            )
    fi(2) = (0.5d0*(sm+abs(sm))*ul + PT*nx )
    fi(3) = (0.5d0*(sm+abs(sm))*vl + PT*ny )
    fi(4) = (0.5d0*(sm+abs(sm))*hl         )
    
   !residual
    rhs(1,i) = rhs(1,i) + 2.0d0*( fij(1) - fi(1) )/ds
    rhs(2,i) = rhs(2,i) + 2.0d0*( fij(2) - fi(2) )/ds
    rhs(3,i) = rhs(3,i) + 2.0d0*( fij(3) - fi(3) )/ds
    rhs(4,i) = rhs(4,i) + 2.0d0*( fij(4) - fi(4) )/ds
    
    enddo
    endif
    enddo
    
    end subroutine residual_slau
    