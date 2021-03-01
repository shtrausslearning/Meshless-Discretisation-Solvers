
    Subroutine residual_roeB
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
    real(rtype) :: dup(4),dm(4),dp(4),sl(4),sr(4),kkk,rx,ry
    

    do i=1,phynod 
    rhs(:,i) = 0.0d0
    if( ptype(i) .eq. 1 )then
    do j=1,nbers(i)
    
     p = conn(j,i)  
     ds = dsqrt(aij(j,i)*aij(j,i)+bij(j,i)*bij(j,i))
     nx = aij(j,i)/ds
     ny = bij(j,i)/ds
     rx = 0.5d0*(x(p)-x(i))
     ry = 0.5d0*(y(p)-y(i))
     
!    LEFT STATE DEFINITION
     rl   = cv(1,i)         + lim(1,i)*(gradx(1,i)*rx+grady(1,i)*ry)
     ul   = cv(2,i)/cv(1,i) + lim(2,i)*(gradx(2,i)*rx+grady(2,i)*ry)
     vl   = cv(3,i)/cv(1,i) + lim(3,i)*(gradx(3,i)*rx+grady(3,i)*ry)
     pl   = dv(1,i)         + lim(4,i)*(gradx(4,i)*rx+grady(4,i)*ry)
     ql = (ul*nx + vl*ny)   ! Convariant Velocity
     gam1 = dv(4,i) - 1.0d0
     ggm1 = dv(4,i)/gam1
     hl  = ggm1*pl/rl + 0.5d0*(ul*ul+vl*vl)

!    RIGHT STATE DEFINITION
     rr   = cv(1,p)         - lim(1,p)*(gradx(1,p)*rx+grady(1,p)*ry)
     ur   = cv(2,p)/cv(1,p) - lim(2,p)*(gradx(2,p)*rx+grady(2,p)*ry)
     vr   = cv(3,p)/cv(1,p) - lim(3,p)*(gradx(3,p)*rx+grady(3,p)*ry)
     pr   = dv(1,p)         - lim(4,p)*(gradx(4,p)*rx+grady(4,p)*ry)
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
    
    fij(1) = fc(1) - 0.25d0*fd(1)
    fij(2) = fc(2) - 0.25d0*fd(2)
    fij(3) = fc(3) - 0.25d0*fd(3)
    fij(4) = fc(4) - 0.25d0*fd(4)

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

    end subroutine residual_roeB