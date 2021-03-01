module covars
  implicit none
  integer covarflag 
  real(8), parameter :: root3=sqrt(3.d0) 
  real(8), parameter :: root5=sqrt(5.d0) 
  real(8), parameter :: third=1.D0/3.D0
contains

!   One dimensional covariance function
    function covar1d(X,Xt,theta)
    implicit none
    real(8), intent(in) :: X,Xt,theta
    
    real(8) covar1d
    real(8) r
    
    r=abs((X-Xt)/theta)
    
    select case(covarflag)
    case(0)
    covar1d=exp(-r)
    case(1)
    covar1d=(1.D0+root3*r)*exp(-root3*r)
    case(2)
    covar1d=(1.D0+root5*r+5.D0*third*r**2)*exp(-root5*r)
    end select

    end function covar1d

    function dcovar1d(X,Xt,theta)
    implicit none
    real(8), intent(in) :: X,Xt,theta
    real(8) dcovar1d
    real(8) r
    real(8) dr

    r=abs((X-Xt)/theta)
    dr=sign(1.D0,(X-Xt)/theta)*(1.D0/theta)
    
    select case(covarflag)
    case(0)
    dcovar1d=-exp(-r)*dr
    case(1)
    dcovar1d=-3.D0*r*exp(-root3*r)*dr
    case(2)
    dcovar1d=-5.D0*third*r*(1.D0+root5*r)*exp(-root5*r)*dr
    end select
    
    end function dcovar1d

    function d2covar1d(X,Xt,theta)
    implicit none
    real(8), intent(in) :: X,Xt,theta
    real(8) d2covar1d
    
    real(8) r
    real(8) dr1,dr2
    real(8) d2r

    r=abs((X-Xt)/theta)
    dr1=sign(1.D0,(X-Xt)/theta)*(1.D0/theta)
    dr2=sign(1.D0,(X-Xt)/theta)*(-1.D0/theta)
    d2r=0.D0
    
    select case(covarflag)
    case(0)
    WRITE(*,*) 'Covariance Function is not twice differentiable, try another'
    stop
    case(1)
    d2covar1d=-3.D0*((1.D0-root3*r)*dr1*dr2+r*d2r)*exp(-root3*r)
    case(2)
    d2covar1d=-5.D0*third*((1.D0+root5*r-5.D0*r**2)*dr2*dr1+r*(1.D0+root5*r)*d2r)*exp(-root5*r)
    end select
       
    end function d2covar1d

!   multidimensional covariance function
    subroutine covarfunc(ndim,X,Xt,theta,k)
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: X(ndim),Xt(ndim),theta(ndim)

    real(8), intent(out) :: k

    real(8) K1
    integer i

    if (ndim==1) then
    K=covar1d(X(1),Xt(1),theta(1))
    else    
    k=1.D0
    do i=1,ndim
    K1=covar1d(X(i),Xt(i),theta(i))
    K=K*K1
    enddo
    endif

    return
    end subroutine covarfunc

    subroutine dcovarfunc(ndim,X,Xt,theta,dk)
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: X(ndim),Xt(ndim),theta(ndim)

    real(8), intent(out) :: dk(ndim)

    integer i

    real(8) L(ndim),U(ndim),Fs(ndim)

    real(8) dkr

    if (ndim==1) then
    dK(1)=dcovar1d(X(1),Xt(1),theta(1))
    else
    do i=1,ndim
    Fs(i)=covar1d(X(i),Xt(i),theta(i))
    enddo

    L(1)=1.D0
    do i=1,ndim-1
    L(i+1)=L(i)*Fs(i)
    enddo

    U(ndim)=1.d0
    do i=ndim,2,-1
    U(i-1)=U(i)*Fs(i)
    enddo

    do i=1,ndim
    dKr=dcovar1d(X(i),Xt(i),theta(i))
    dK(i)=L(i)*dKr*U(i)
    enddo
    endif

    return
    end subroutine dcovarfunc

    subroutine d2covarfunc(ndim,X,Xt,theta,d2k)
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: X(ndim),Xt(ndim),theta(ndim)

    real(8), intent(out) :: d2k(ndim,ndim)

    integer i,j

    real(8) Fs(ndim),Gs(ndim),Hs(ndim)
    real(8) L(ndim),U(ndim),M(ndim,ndim)
    real(8) dkr,d2kr

    if (ndim==1) then
    d2k=d2covar1d(X(1),Xt(1),theta(1))
    else
    do i=1,ndim
    Fs(i)=covar1d(X(i),Xt(i),theta(i))
    enddo
    do i=1,ndim
    Gs(i)=dcovar1d(X(i),Xt(i),theta(i))
    enddo
    do i=1,ndim
    Hs(i)=d2covar1d(X(i),Xt(i),theta(i))
    enddo

    L(1)=1.D0
    do i=1,ndim-1
    L(i+1)=L(i)*Fs(i)
    enddo

    U(ndim)=1.d0
    do i=ndim,2,-1
    U(i-1)=U(i)*Fs(i)
    enddo

    do i=1,ndim-1
    M(i,i+1)=1.D0
    do j=i+2,ndim
    M(i,j)=M(i,j-1)*Fs(j-1)
    enddo
    enddo

    do i=1,ndim
    d2k(i,i)=L(i)*Hs(i)*U(i)
    enddo

    do i=1,ndim
    do j=i+1,ndim
    d2k(i,j)=-L(i)*Gs(i)*M(i,j)*Gs(j)*U(j)
    d2k(j,i)=d2k(i,j)
    enddo
    enddo
    endif

    return
    end subroutine d2covarfunc
    
    end module covars

!   Create the covariance matrix between training points (function kriging)
    subroutine CovarT(ndim,ntot,X,theta,Kmat)
    use covars
    implicit none
    integer, intent(in) :: ndim, ntot
    real(8), intent(in) :: X(ndim,ntot)
    real(8), intent(in) :: theta(ndim)

    real(8), intent(out) :: Kmat(ntot,ntot)

    integer i,j
    real(8) kr

!   Output just lower part of entire matrix as it's symmetric
    do i=1,ntot
    do j=1,i
    call covarfunc(ndim,X(:,i),X(:,j),theta,kr)
    Kmat(i,j)=kr
    enddo
    enddo
    
    return
    end subroutine CovarT

!   Build a general covariance matrix between two sets of points (function kriging)
!   ie. construct the covariance between training points & test points
    subroutine CovarG(ndim,ntot,X,mtot,Xm,theta,Kmat)
    use covars
    implicit none
    integer, intent(in) :: ndim, ntot
    real(8), intent(in) :: X(ndim,ntot)
    integer, intent(in) :: mtot
    real(8), intent(in) :: Xm(ndim,mtot)
    real(8), intent(in) :: theta(ndim)

    real(8), intent(out) :: Kmat(ntot,mtot)

    integer i,j
    real(8) kr

    do i=1,ntot
    do j=1,mtot
    call covarfunc(ndim,X(:,i),Xm(:,j),theta,kr)
    Kmat(i,j)=kr
    enddo
    enddo

    return
    end subroutine CovarG

    subroutine GCovar(ndim,ntot,X,mtot,htot,ptsm,dimsm,Xm,theta,Kmat)
    use covars
    implicit none
    integer, intent(in) :: ndim, ntot
    real(8), intent(in) :: X(ndim,ntot)

    integer, intent(in) :: mtot,htot
    integer, intent(in) :: ptsm(htot),dimsm(htot)
    real(8), intent(in) :: Xm(ndim,mtot)

    real(8), intent(in) :: theta(ndim)

    real(8), intent(out) :: Kmat(ntot,htot)

    integer i,j,k,l
    real(8) dkr(ndim)

    integer nderm(mtot)
    integer npointm(mtot+1)

    nderm=0
    do i=1,htot
       nderm(ptsm(i))=nderm(ptsm(i))+1
    enddo

    npointm(1)=1
    do i=1,mtot
       npointm(i+1)=npointm(i)+nderm(i)
    enddo

    do i=1,ntot
       do j=1,mtot
          call dcovarfunc(ndim,X(:,i),Xm(:,j),theta,dkr)
          do k=npointm(j),npointm(j+1)-1
             l=dimsm(k)
             Kmat(i,k)=-dkr(l)
          enddo
       enddo
    enddo

    return
  end subroutine GCovar
  


module covarmatrix_grad_mod
  implicit none
  interface covarmatrix_grad
     module procedure covarmatrix_grad_gek, covarmatrix_grad_func
  end interface covarmatrix_grad
contains
!> \brief Subroutine to create the covariance matrix between a set of function values and a set of derivative values (typically located at the test points).
  subroutine covarmatrix_grad_func(ndim,ntot,X,mtot,htot,ptsm,dimsm,Xm,theta,Kmat)
    use covars
    implicit none
    integer, intent(in) :: ndim, ntot
    real(8), intent(in) :: X(ndim,ntot)

    integer, intent(in) :: mtot,htot
    integer, intent(in) :: ptsm(htot),dimsm(htot)
    real(8), intent(in) :: Xm(ndim,mtot)

    real(8), intent(in) :: theta(ndim)

    real(8), intent(out) :: Kmat(ntot,htot)

    integer i,j,k,l
    real(8) dkr(ndim)

    integer nderm(mtot)
    integer npointm(mtot+1)

    nderm=0
    do i=1,htot
       nderm(ptsm(i))=nderm(ptsm(i))+1
    enddo

    npointm(1)=1
    do i=1,mtot
       npointm(i+1)=npointm(i)+nderm(i)
    enddo

    do i=1,ntot
       do j=1,mtot
          call dcovarfunc(ndim,X(:,i),Xm(:,j),theta,dkr)
          do k=npointm(j),npointm(j+1)-1
             l=dimsm(k)
             Kmat(i,k)=-dkr(l)
          enddo
       enddo
    enddo

    return
  end subroutine covarmatrix_grad_func
!> \brief Subroutine to compute the covariance between a set of training data and a set of derivative values for a gradient-enhanced model. For a gradient-enhanced model, the set of training data includes function and derivative values evaluated at the training points. 
  subroutine covarmatrix_grad_gek(ndim,ntot,gtot,pts,dims,X,mtot,htot,ptsm,dimsm,Xm,theta,Kmat)
    use covars
    implicit none
    integer, intent(in) :: ndim, ntot,gtot
    integer, intent(in) :: pts(gtot),dims(gtot)
    real(8), intent(in) :: X(ndim,ntot)

    integer, intent(in) :: mtot,htot
    integer, intent(in) :: ptsm(htot),dimsm(htot)
    real(8), intent(in) :: Xm(ndim,mtot)

    real(8), intent(in) :: theta(ndim)

    real(8), intent(out) :: Kmat(ntot+gtot,htot)

    integer i,j,k,l
    real(8) kr,dkr(ndim),d2kr(ndim,ndim)

    integer nder(ntot), nderm(mtot)
    integer npoint(ntot+1),npointm(mtot+1)

    nder=0
    do i=1,gtot
       nder(pts(i))=nder(pts(i))+1
    enddo
    nderm=0
    do i=1,htot
       nderm(ptsm(i))=nderm(ptsm(i))+1
    enddo

    npoint(1)=1
    do i=1,ntot
       npoint(i+1)=npoint(i)+nder(i)
    enddo
    npointm(1)=1
    do i=1,mtot
       npointm(i+1)=npointm(i)+nderm(i)
    enddo

    do i=1,ntot
       do j=1,mtot
          call dcovarfunc(ndim,X(:,i),Xm(:,j),theta,dkr)
          do k=npointm(j),npointm(j+1)-1
             l=dimsm(k)
             Kmat(i,k)=-dkr(l)
          enddo
       enddo
    enddo

    do i=1,ntot
       do j=1,mtot
          call d2covarfunc(ndim,X(:,i),Xm(:,j),theta,d2kr)
          do k=npoint(i),npoint(i+1)-1
             do l=npointm(j),npointm(j+1)-1
                Kmat(ntot+k,l)=d2kr(dims(k),dimsm(l))
             enddo
          enddo
       enddo
    enddo

    return
  end subroutine covarmatrix_grad_gek
end module covarmatrix_grad_mod
