
    Subroutine residual_ausm
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    
    integer :: i,j,p
    real(rtype) :: ds,nx,ny,shigma,Ku,Kp
	real(rtype) :: rho,gam1,ggm1,rl,ul,vl,pl,ql2,hl,cl2
	real(rtype) :: rr,ur,vr,pr,qr2,hr,cr2,cbar,Coff
	real(rtype) :: VnLFT,VnRHT,sMLFT,sMRHT
	real(rtype) :: Mbar2,Mbar,Mo2,Mo,faMo,alphaco,fppLa,fpmRa,fMpLb,fMmRb
	real(rtype) :: R1over2,Mp,M1over2,Pu,PT,sm,fc(4),fi(4),rrho
    
    do i=1,phynod 
    rhs(:,i) = 0.0d0
    if( ptype(i) .eq. 1 )then
    do j=1,nbers(i)
    
     p = conn(j,i)  
     ds = sqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
     nx = aij(j,i)/ds
     ny = bij(j,i)/ds
     
!   AUSM CONSTANTS [ VARIABLE ]
    shigma = 1.0d0
    Ku = 0.75d0
    Kp = 0.25d0
    
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
    cbar = 0.5d0 * ( sqrt(cl2) + sqrt(cr2) ) !SI

!   %%%%%%%%%%%%%%%%%%%%%%
!   SETTING AUSM VARIABLES
!   %%%%%%%%%%%%%%%%%%%%%%
    
!   Cutoff Mach ( need to set )
    Coff = 1.0d0 ! decide
    
!   lhs,rhs convariant velocity 
    VnLFT = nx*ul + ny*vl !SI
    VnRHT = nx*ur + ny*vr 
!   lhs,rhs Mach states
    sMLFT = VnLFT/cbar !Ml
    sMRHT = VnRHT/cbar !Mr
    
!   added
    Mbar2 = 0.5d0/(cbar**2) * (VnLFT**2 + VnRHT**2)
    Mbar = DSQRT(Mbar2)
    Mo2 = dmin1(1.0d0,dmax1(Mbar2,(Coff)**2))
     
    Mo = DSQRT(Mo2)
    faMo = Mo*(2.0d0-Mo)
    alphaco = 0.1875d0*(-4.0d0 + 5.0d0*(faMo**2))
    
!    F +/-_p | ALPHA FUNCTION
     IF( DABS(sMLFT) >= 1.0d0 )THEN
       fppLa = 0.5d0 * ( 1.0d0 + dmax1(-1.0d0,dmin1(1.0d0,sMLFT) ))
     ELSE
      fppLa = 0.25d0*(sMLFT+1.0d0)**2*(2.0d0-sMLFT) 
     ENDIF 
     IF( DABS(sMRHT) >= 1.0d0 )THEN
      fpmRa = 0.5d0 * ( 1.0d0 - dmax1(-1.0d0,dmin1(1.0d0,sMRHT)) )
     ELSE
      fpmRa = 0.25d0*(sMRHT-1.0d0)**2*(2.0d0+sMRHT) 
     ENDIF  
     
!    F+/-_M | BETA FUNCTION
     IF( DABS(sMLFT) >= 1.0d0 )THEN
       fMpLb = 0.5d0*(sMLFT+dabs(sMLFT))
     ELSE
       fMpLb = 0.25d0*(sMLFT+1.0d0)**2 + 0.125d0*(sMLFT**2 - 1.0d0)**2
     ENDIF
     IF( DABS(sMRHT) >= 1.0d0 )THEN
       fMmRb = 0.5d0*(sMRHT-dabs(sMRHT))
     ELSE
       fMmRb = - 0.25d0*(sMRHT-1.0d0)**2 - 0.125d0*(sMRHT**2 - 1.0d0)**2
     ENDIF  
     
     R1over2 = 0.5d0*(rl+rr)
     Mp = - Kp/faMo*dmax1((1.0d0-shigma*Mbar2),0.0d0)*(pr-pl)/R1over2/(cbar**2)
     M1over2 = fMpLb + fMmRb + Mp
!     M1over2 = ( sMLFT + sMRHT )

     Pu = - Ku * fppLa * fpmRa * (rl+rr)*(faMo*cbar)*(VnRHT-VnLFT)
    
!    PRESSURE FLUX FUNCTION
!     PT = fppLa*pl + fpmRa*pr + Pu
     PT = 0.5d0*(pl+pr)
     
!    MASS FLUX FUNCTION
     IF( M1over2 > 0.0d0 )THEN
      sm = M1over2 * cbar * rl
     ELSE
      sm = M1over2 * cbar * rr
     ENDIF 
     
!   AUSM+up flux
    fc(1) = (0.5d0*(sm+abs(sm))    +0.5d0*(sm-abs(sm))           )*ds
    fc(2) = (0.5d0*(sm+abs(sm))*ul +0.5d0*(sm-abs(sm))*ur+ PT*nx )*ds
    fc(3) = (0.5d0*(sm+abs(sm))*vl +0.5d0*(sm-abs(sm))*vr+ PT*ny )*ds
    fc(4) = (0.5d0*(sm+abs(sm))*hl +0.5d0*(sm-abs(sm))*hr        )*ds
    
!   AUSM+up flux
    fi(1) = 0.5d0*((sm+abs(sm))               )*ds
    fi(2) = 0.5d0*((sm+abs(sm))*ul + pl*nx    )*ds
    fi(3) = 0.5d0*((sm+abs(sm))*vl + pl*ny    )*ds
    fi(4) = 0.5d0*((sm+abs(sm))*hl            )*ds
    
   !residual
    rhs(1,i) = rhs(1,i) + 2.0d0*( fc(1) - fi(1) )
    rhs(2,i) = rhs(2,i) + 2.0d0*( fc(2) - fi(2) )
    rhs(3,i) = rhs(3,i) + 2.0d0*( fc(3) - fi(3) )
    rhs(4,i) = rhs(4,i) + 2.0d0*( fc(4) - fi(4) )

     
     
    
    enddo
    endif
    enddo
    
    End subroutine residual_ausm