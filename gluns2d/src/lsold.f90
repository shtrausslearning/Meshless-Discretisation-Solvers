    subroutine LeastSquares
!   ###################################################################
    use prmflow
    implicit none

    integer          :: i, j, nn,ierr,jj,inn
    double precision :: x0, y0, x1(12), y1(12)
    
    ierr=0;allocate( aij(12,phynod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for bij()" )
    ierr=0;allocate( bij(12,phynod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for bij()" )
    
    aij(:,:) = 0.0d0
    bij(:,:) = 0.0d0
    
    do i=1,phynod
    
     x0 = x(i);y0 = y(i)
     do j=1,nbers(i)                   ! for all neighbours of main cube
     nn    = conn(j,i)                 ! neighbour cube index
     x1(j) = x(nn)
     y1(j) = y(nn)
     enddo
     
     call LSCoeff(nbers(i),x0,y0,x1,y1,aij(1,i),bij(1,i))
  !   call LSCoeff2(nbers(i),x0,y0,x1,y1,aij(1,i),bij(1,i))
     
    enddo
      
    return
!   ###################################################################
    end subroutine LeastSquares
      
!    subroutine LSCoeff(n,x0,y0,x1,y1,cx,cy)
!!   ###################################################################
!    use prmflow
!!   Find coefficients for first derivative using least squares
!!   n           = Number of neighbours
!!   x0,y0       = coordinates of point at which grad is required
!!   x1(n),y1(n) = coordinates of neighbouring points
!!   cx(n),cy(n) = coefficients
!!   ###################################################################
!    implicit none
!    integer          :: n
!    double precision :: x0, y0, x1(*), y1(*), cx(*), cy(*)
!
!    integer          :: i, j
!    double precision :: sdx2, sdy2, sdxdy, dx, dy, dr, w(n), a(n,2), det
!    double precision :: ainv(2,2)
!
!    sdx2  = 0.0d0
!    sdy2  = 0.0d0
!    sdxdy = 0.0d0
!
!    do i=1,n
!    dx     = x1(i) - x0
!    dy     = y1(i) - y0
!    dr     = dsqrt(dx**2 + dy**2)
!    w(i)   = 1.0d0/dr
!    a(i,1) = w(i)*dx
!    a(i,2) = w(i)*dy
!    sdx2   = sdx2  + a(i,1)**2
!    sdy2   = sdy2  + a(i,2)**2
!    sdxdy  = sdxdy + a(i,1)*a(i,2)
!    enddo
!
!!   invert 2x2 symmetric least squares matrix
!    det       = sdx2*sdy2 - sdxdy**2
!    ainv(1,1) = sdy2/det
!    ainv(1,2) =-sdxdy/det
!    ainv(2,1) = ainv(1,2)
!    ainv(2,2) = sdx2/det
!
!!     finally the coefficients,
!!     [ cx ]
!!     [    ] = ainv * a^T * diag(w)
!!     [ cy ]
!    do i=1,n
!     cx(i) = 0.0d0
!     cy(i) = 0.0d0
!     do j=1,2
!     cx(i) = cx(i) + ainv(1,j)*a(i,j)
!     cy(i) = cy(i) + ainv(2,j)*a(i,j)
!     enddo
!    enddo
!
!    return
!!   ###################################################################
!    End subroutine LSCoeff
    
    subroutine LSCoeff(n,x0,y0,x1,y1,cx,cy)
    use prmflow
!     ======================================================================
!     Find coefficients for first derivative using least squares
!     n           = Number of neighbours
!     x0,y0       = coordinates of point at which grad is required
!     x1(n),y1(n) = coordinates of neighbouring points
!     cx(n),cy(n) = coefficients
!     ======================================================================
    implicit none
    integer          :: n
    double precision :: x0, y0, x1(*), y1(*), cx(*), cy(*)
    integer          :: flag
    integer          :: i, j
    double precision :: dx, dy, dr, w(n), det
    double precision :: ainv(2,2)
    double precision :: c(n,2), a(2,2)
      
    a(:,:) = 0.0d0 ; c(:,:) = 0.0d0
      
    do i=1,n ! for all neighbours
      
    dx     = x1(i) - x0            ! check
    dy     = y1(i) - y0            ! check
    dr     = dsqrt(dx**2 + dy**2)  
    w(i)   = 1.0d0/dr
         
    !       MATRIX C individual components i=1,n
    c(i,1) = w(i)*dx
    c(i,2) = w(i)*dy
         
    !       LEAST SQUARES MATRIX COMPONENTS 
    a(1,1) = a(1,1) +       w(i)*dx**2
    a(1,2) = a(1,2) +       w(i)*dx*dy
    a(2,1) = a(1,2)
    a(2,2) = a(2,2) +       w(i)*(dy**2)
      
    enddo   

    CALL m22inv(a,ainv,flag)

    do i=1,n ! for all neighbours
    cx(i) = 0.0d0
    cy(i) = 0.0d0
    do j=1,2
    cx(i) = cx(i) + ainv(1,j)*c(i,j)
    cy(i) = cy(i) + ainv(2,j)*c(i,j)
    enddo
    enddo

    return
    End Subroutine LSCoeff
      
    subroutine m22inv(a,ainv,flag)
    implicit none

    double precision, dimension(2,2), intent(in)  :: a
    double precision, dimension(2,2), intent(out) :: ainv
    integer, intent(out) :: flag

    double precision, parameter :: eps = 1.0d-10
    double precision :: det
    double precision, dimension(2,2) :: cofactor

    det =   a(1,1)*a(2,2) - a(1,2)*a(2,1)

!   Singular Matrix Condition Check
    if( abs(det) .le. eps)then
    ainv = 0.0d0 ! set default
    flag = 0  ! non invertable matrix [ singular matrix ]
    return
    end if

    cofactor(1,1) = +a(2,2)
    cofactor(1,2) = -a(2,1)
    cofactor(2,1) = -a(1,2)
    cofactor(2,2) = +a(1,1)

    ainv = transpose(cofactor)/det

    flag = 1 ! invertable matrix

    return
    end subroutine m22inv
    
    subroutine LSCoeff2(n,x0,y0,x1,y1,cx,cy)
    use prmflow
!     ======================================================================
!     Find coefficients for first derivative using least squares
!     n           = Number of neighbours
!     x0,y0       = coordinates of point at which grad is required
!     x1(n),y1(n) = coordinates of neighbouring points
!     cx(n),cy(n) = coefficients
!     ======================================================================
    implicit none
    integer          :: n
    double precision :: x0, y0, x1(*), y1(*), cx(*), cy(*)
    integer          :: flag
    integer          :: i, j
    double precision :: dx, dy, dr, w(n), det
    double precision :: ainv(5,5), c(n,5), a(5,5)

    a(:,:) = 0.0d0 ; c(:,:) = 0.0d0

    do i=1,n ! for all neighbours
      
        dx     = x1(i) - x0            ! check
        dy     = y1(i) - y0            ! check
        dr     = dsqrt(dx**2 + dy**2)  
        w(i)   = 1.0d0/dr
         
    !        MATRIX C individual components
        c(i,1) = w(i)*dx
        c(i,2) = w(i)*dy
        c(i,3) = w(i)*0.5d0*dx**2
        c(i,4) = w(i)*0.5d0*dy**2
        c(i,5) = w(i)*dx*dy
         
    !        LEAST SQUARES MATRIX COMPONENTS
        a(1,1) = a(1,1) +        w(i)*dx**2
        a(1,2) = a(1,2) +        w(i)*dx*dy
        a(1,3) = a(1,3) + 0.50d0*w(i)*dx**3
        a(1,4) = a(1,4) + 0.50d0*w(i)*(dx**2)*dy
        a(1,5) = a(1,5) +        w(i)*(dx**2)*dy
         
        a(2,1) = a(1,2)
        a(2,2) = a(2,2) +        w(i)*(dy**2)
        a(2,3) = a(1,4) 
        a(2,4) = a(2,4) + 0.50d0*w(i)*(dy**3)
        a(2,5) = a(2,5) +        w(i)*dx*(dy**2)
         
        a(3,1) = a(3,1) + 0.50d0*w(i)*(dx**3)
        a(3,2) = a(1,4)
        a(3,3) = a(3,3) + 0.25d0*w(i)*(dx**4)
        a(3,4) = a(3,4) + 0.25d0*w(i)*(dx**2)*(dy**2)
        a(3,5) = a(3,5) + 0.50d0*w(i)*(dx**3)*dy
         
        a(4,1) = a(1,4)
        a(4,2) = a(2,4)
        a(4,3) = a(3,4)
        a(4,4) = a(3,3)
        a(4,5) = a(4,5) + 0.50d0*w(i)*dx*(dy**3)
         
        a(5,1) = a(1,5)
        a(5,2) = a(2,5)
        a(5,3) = a(3,5)
        a(5,4) = a(4,5)
        a(5,5) = a(5,5) +        w(i)*(dx**2)*(dy**2)
         

    enddo
      
    !call inverse(a,ainv,5)
    call m55inv(a,ainv,flag)
      
    !     Check
    !if (flag .eq. 0) then
    ! write (*,'(/a)') ' singular matrix.'
    ! pause
    !end if

    do i=1,n ! for all neighbours
    cx(i) = 0.0d0
    cy(i) = 0.0d0
    do j=1,5
    cx(i) = cx(i) + ainv(1,j)*c(i,j) ! Fx0 component
    cy(i) = cy(i) + ainv(2,j)*c(i,j) ! Fy0 component
    enddo
    enddo

    return
    End Subroutine LSCoeff2
      
      
    subroutine inverse(a,c,n)
    !============================================================
    ! Inverse matrix
    ! Method: Based on Doolittle LU factorization for Ax=b
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! n      - dimension
    ! output ...
    ! c(n,n) - inverse matrix of A
    ! comments ...
    ! the original matrix a(n,n) will be destroyed 
    ! during the calculation
    !===========================================================
    implicit none 
    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
    do i=k+1,n
        coeff=a(i,k)/a(k,k)
        L(i,k) = coeff
        do j=k+1,n
            a(i,j) = a(i,j)-coeff*a(k,j)
        end do
    end do
    end do

    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
    L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
    do i=1,j
    U(i,j) = a(i,j)
    end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
    b(k)=1.0
    d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
    do i=2,n
    d(i)=b(i)
    do j=1,i-1
        d(i) = d(i) - L(i,j)*d(j)
    end do
    end do
    ! Step 3b: Solve Ux=d using the back substitution
    x(n)=d(n)/U(n,n)
    do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
        x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
    end do
    ! Step 3c: fill the solutions x(n) into column k of C
    do i=1,n
    c(i,k) = x(i)
    end do
    b(k)=0.0
    end do
    
    end subroutine inverse
    
    subroutine m55inv (a,ainv,flag)
    implicit none

    double precision, dimension(5,5), intent(in)  :: a
    double precision, dimension(5,5), intent(out) :: ainv
    integer, intent(out) :: flag

    double precision, parameter :: eps = 1.0d-10
    double precision :: det, a11, a12, a13, a14, a15, a21, a22, a23, a24, &
                 a25, a31, a32, a33, a34, a35, a41, a42, a43, a44, a45,   &
                 a51, a52, a53, a54, a55
    double precision, dimension(5,5) :: cofactor

    a11=a(1,1); a12=a(1,2); a13=a(1,3); a14=a(1,4); a15=a(1,5)
    a21=a(2,1); a22=a(2,2); a23=a(2,3); a24=a(2,4); a25=a(2,5)
    a31=a(3,1); a32=a(3,2); a33=a(3,3); a34=a(3,4); a35=a(3,5)
    a41=a(4,1); a42=a(4,2); a43=a(4,3); a44=a(4,4); a45=a(4,5)
    a51=a(5,1); a52=a(5,2); a53=a(5,3); a54=a(5,4); a55=a(5,5)

    det = a15*a24*a33*a42*a51-a14*a25*a33*a42*a51-a15*a23*a34*a42*a51+    &
        a13*a25*a34*a42*a51+a14*a23*a35*a42*a51-a13*a24*a35*a42*a51-       &
        a15*a24*a32*a43*a51+a14*a25*a32*a43*a51+a15*a22*a34*a43*a51-       &
        a12*a25*a34*a43*a51-a14*a22*a35*a43*a51+a12*a24*a35*a43*a51+       &
        a15*a23*a32*a44*a51-a13*a25*a32*a44*a51-a15*a22*a33*a44*a51+       &
        a12*a25*a33*a44*a51+a13*a22*a35*a44*a51-a12*a23*a35*a44*a51-       &
        a14*a23*a32*a45*a51+a13*a24*a32*a45*a51+a14*a22*a33*a45*a51-       &
        a12*a24*a33*a45*a51-a13*a22*a34*a45*a51+a12*a23*a34*a45*a51-       &
        a15*a24*a33*a41*a52+a14*a25*a33*a41*a52+a15*a23*a34*a41*a52-       &
        a13*a25*a34*a41*a52-a14*a23*a35*a41*a52+a13*a24*a35*a41*a52+       &
        a15*a24*a31*a43*a52-a14*a25*a31*a43*a52-a15*a21*a34*a43*a52+       &
        a11*a25*a34*a43*a52+a14*a21*a35*a43*a52-a11*a24*a35*a43*a52-       &
        a15*a23*a31*a44*a52+a13*a25*a31*a44*a52+a15*a21*a33*a44*a52-       &
        a11*a25*a33*a44*a52-a13*a21*a35*a44*a52+a11*a23*a35*a44*a52+       &
        a14*a23*a31*a45*a52-a13*a24*a31*a45*a52-a14*a21*a33*a45*a52+       &
        a11*a24*a33*a45*a52+a13*a21*a34*a45*a52-a11*a23*a34*a45*a52+       &
        a15*a24*a32*a41*a53-a14*a25*a32*a41*a53-a15*a22*a34*a41*a53+       &
        a12*a25*a34*a41*a53+a14*a22*a35*a41*a53-a12*a24*a35*a41*a53-       &
        a15*a24*a31*a42*a53+a14*a25*a31*a42*a53+a15*a21*a34*a42*a53-       &
        a11*a25*a34*a42*a53-a14*a21*a35*a42*a53+a11*a24*a35*a42*a53+       &
        a15*a22*a31*a44*a53-a12*a25*a31*a44*a53-a15*a21*a32*a44*a53+       &
        a11*a25*a32*a44*a53+a12*a21*a35*a44*a53-a11*a22*a35*a44*a53-       &
        a14*a22*a31*a45*a53+a12*a24*a31*a45*a53+a14*a21*a32*a45*a53-       &
        a11*a24*a32*a45*a53-a12*a21*a34*a45*a53+a11*a22*a34*a45*a53-       &
        a15*a23*a32*a41*a54+a13*a25*a32*a41*a54+a15*a22*a33*a41*a54-       &
        a12*a25*a33*a41*a54-a13*a22*a35*a41*a54+a12*a23*a35*a41*a54+       &
        a15*a23*a31*a42*a54-a13*a25*a31*a42*a54-a15*a21*a33*a42*a54+       &
        a11*a25*a33*a42*a54+a13*a21*a35*a42*a54-a11*a23*a35*a42*a54-       &
        a15*a22*a31*a43*a54+a12*a25*a31*a43*a54+a15*a21*a32*a43*a54-       &
        a11*a25*a32*a43*a54-a12*a21*a35*a43*a54+a11*a22*a35*a43*a54+       &
        a13*a22*a31*a45*a54-a12*a23*a31*a45*a54-a13*a21*a32*a45*a54+       &
        a11*a23*a32*a45*a54+a12*a21*a33*a45*a54-a11*a22*a33*a45*a54+       &
        a14*a23*a32*a41*a55-a13*a24*a32*a41*a55-a14*a22*a33*a41*a55+       &
        a12*a24*a33*a41*a55+a13*a22*a34*a41*a55-a12*a23*a34*a41*a55-       &
        a14*a23*a31*a42*a55+a13*a24*a31*a42*a55+a14*a21*a33*a42*a55-       &
        a11*a24*a33*a42*a55-a13*a21*a34*a42*a55+a11*a23*a34*a42*a55+       &
        a14*a22*a31*a43*a55-a12*a24*a31*a43*a55-a14*a21*a32*a43*a55+       &
        a11*a24*a32*a43*a55+a12*a21*a34*a43*a55-a11*a22*a34*a43*a55-       &
        a13*a22*a31*a44*a55+a12*a23*a31*a44*a55+a13*a21*a32*a44*a55-       &
        a11*a23*a32*a44*a55-a12*a21*a33*a44*a55+a11*a22*a33*a44*a55

    !if( abs(det) .le. eps)then
    if( abs(det) .eq. 0.0d0 )then
     ainv = 0.0d0
     flag = 0
     return
    endif

    cofactor(1,1) = a25*a34*a43*a52-a24*a35*a43*a52-a25*a33*a44*a52+      &
        a23*a35*a44*a52+a24*a33*a45*a52-a23*a34*a45*a52-a25*a34*a42*a53+   &
        a24*a35*a42*a53+a25*a32*a44*a53-a22*a35*a44*a53-a24*a32*a45*a53+   &
        a22*a34*a45*a53+a25*a33*a42*a54-a23*a35*a42*a54-a25*a32*a43*a54+   &
        a22*a35*a43*a54+a23*a32*a45*a54-a22*a33*a45*a54-a24*a33*a42*a55+   &
        a23*a34*a42*a55+a24*a32*a43*a55-a22*a34*a43*a55-a23*a32*a44*a55+   &
        a22*a33*a44*a55

    cofactor(2,1) = -a15*a34*a43*a52+a14*a35*a43*a52+a15*a33*a44*a52-     &
        a13*a35*a44*a52-a14*a33*a45*a52+a13*a34*a45*a52+a15*a34*a42*a53-   &
        a14*a35*a42*a53-a15*a32*a44*a53+a12*a35*a44*a53+a14*a32*a45*a53-   &
        a12*a34*a45*a53-a15*a33*a42*a54+a13*a35*a42*a54+a15*a32*a43*a54-   &
        a12*a35*a43*a54-a13*a32*a45*a54+a12*a33*a45*a54+a14*a33*a42*a55-   &
        a13*a34*a42*a55-a14*a32*a43*a55+a12*a34*a43*a55+a13*a32*a44*a55-   &
        a12*a33*a44*a55

    cofactor(3,1) = a15*a24*a43*a52-a14*a25*a43*a52-a15*a23*a44*a52+      &
        a13*a25*a44*a52+a14*a23*a45*a52-a13*a24*a45*a52-a15*a24*a42*a53+   &
        a14*a25*a42*a53+a15*a22*a44*a53-a12*a25*a44*a53-a14*a22*a45*a53+   &
        a12*a24*a45*a53+a15*a23*a42*a54-a13*a25*a42*a54-a15*a22*a43*a54+   &
        a12*a25*a43*a54+a13*a22*a45*a54-a12*a23*a45*a54-a14*a23*a42*a55+   &
        a13*a24*a42*a55+a14*a22*a43*a55-a12*a24*a43*a55-a13*a22*a44*a55+   &
        a12*a23*a44*a55

    cofactor(4,1) = -a15*a24*a33*a52+a14*a25*a33*a52+a15*a23*a34*a52-     &
        a13*a25*a34*a52-a14*a23*a35*a52+a13*a24*a35*a52+a15*a24*a32*a53-   &
        a14*a25*a32*a53-a15*a22*a34*a53+a12*a25*a34*a53+a14*a22*a35*a53-   &
        a12*a24*a35*a53-a15*a23*a32*a54+a13*a25*a32*a54+a15*a22*a33*a54-   &
        a12*a25*a33*a54-a13*a22*a35*a54+a12*a23*a35*a54+a14*a23*a32*a55-   &
        a13*a24*a32*a55-a14*a22*a33*a55+a12*a24*a33*a55+a13*a22*a34*a55-   &
        a12*a23*a34*a55

    cofactor(5,1) = a15*a24*a33*a42-a14*a25*a33*a42-a15*a23*a34*a42+      &
        a13*a25*a34*a42+a14*a23*a35*a42-a13*a24*a35*a42-a15*a24*a32*a43+   &
        a14*a25*a32*a43+a15*a22*a34*a43-a12*a25*a34*a43-a14*a22*a35*a43+   &
        a12*a24*a35*a43+a15*a23*a32*a44-a13*a25*a32*a44-a15*a22*a33*a44+   &
        a12*a25*a33*a44+a13*a22*a35*a44-a12*a23*a35*a44-a14*a23*a32*a45+   &
        a13*a24*a32*a45+a14*a22*a33*a45-a12*a24*a33*a45-a13*a22*a34*a45+   &
        a12*a23*a34*a45

    cofactor(1,2) = -a25*a34*a43*a51+a24*a35*a43*a51+a25*a33*a44*a51-     &
        a23*a35*a44*a51-a24*a33*a45*a51+a23*a34*a45*a51+a25*a34*a41*a53-   &
        a24*a35*a41*a53-a25*a31*a44*a53+a21*a35*a44*a53+a24*a31*a45*a53-   &
        a21*a34*a45*a53-a25*a33*a41*a54+a23*a35*a41*a54+a25*a31*a43*a54-   &
        a21*a35*a43*a54-a23*a31*a45*a54+a21*a33*a45*a54+a24*a33*a41*a55-   &
        a23*a34*a41*a55-a24*a31*a43*a55+a21*a34*a43*a55+a23*a31*a44*a55-   &
        a21*a33*a44*a55

    cofactor(2,2) = a15*a34*a43*a51-a14*a35*a43*a51-a15*a33*a44*a51+      &
        a13*a35*a44*a51+a14*a33*a45*a51-a13*a34*a45*a51-a15*a34*a41*a53+   &
        a14*a35*a41*a53+a15*a31*a44*a53-a11*a35*a44*a53-a14*a31*a45*a53+   &
        a11*a34*a45*a53+a15*a33*a41*a54-a13*a35*a41*a54-a15*a31*a43*a54+   &
        a11*a35*a43*a54+a13*a31*a45*a54-a11*a33*a45*a54-a14*a33*a41*a55+   &
        a13*a34*a41*a55+a14*a31*a43*a55-a11*a34*a43*a55-a13*a31*a44*a55+   &
        a11*a33*a44*a55

    cofactor(3,2) = -a15*a24*a43*a51+a14*a25*a43*a51+a15*a23*a44*a51-     &
        a13*a25*a44*a51-a14*a23*a45*a51+a13*a24*a45*a51+a15*a24*a41*a53-   &
        a14*a25*a41*a53-a15*a21*a44*a53+a11*a25*a44*a53+a14*a21*a45*a53-   &
        a11*a24*a45*a53-a15*a23*a41*a54+a13*a25*a41*a54+a15*a21*a43*a54-   &
        a11*a25*a43*a54-a13*a21*a45*a54+a11*a23*a45*a54+a14*a23*a41*a55-   &
        a13*a24*a41*a55-a14*a21*a43*a55+a11*a24*a43*a55+a13*a21*a44*a55-   &
        a11*a23*a44*a55

    cofactor(4,2) = a15*a24*a33*a51-a14*a25*a33*a51-a15*a23*a34*a51+      &
        a13*a25*a34*a51+a14*a23*a35*a51-a13*a24*a35*a51-a15*a24*a31*a53+   &
        a14*a25*a31*a53+a15*a21*a34*a53-a11*a25*a34*a53-a14*a21*a35*a53+   &
        a11*a24*a35*a53+a15*a23*a31*a54-a13*a25*a31*a54-a15*a21*a33*a54+   &
        a11*a25*a33*a54+a13*a21*a35*a54-a11*a23*a35*a54-a14*a23*a31*a55+   &
        a13*a24*a31*a55+a14*a21*a33*a55-a11*a24*a33*a55-a13*a21*a34*a55+   &
        a11*a23*a34*a55

    cofactor(5,2) = -a15*a24*a33*a41+a14*a25*a33*a41+a15*a23*a34*a41-     &
        a13*a25*a34*a41-a14*a23*a35*a41+a13*a24*a35*a41+a15*a24*a31*a43-   &
        a14*a25*a31*a43-a15*a21*a34*a43+a11*a25*a34*a43+a14*a21*a35*a43-   &
        a11*a24*a35*a43-a15*a23*a31*a44+a13*a25*a31*a44+a15*a21*a33*a44-   &
        a11*a25*a33*a44-a13*a21*a35*a44+a11*a23*a35*a44+a14*a23*a31*a45-   &
        a13*a24*a31*a45-a14*a21*a33*a45+a11*a24*a33*a45+a13*a21*a34*a45-   &
        a11*a23*a34*a45

    cofactor(1,3) = a25*a34*a42*a51-a24*a35*a42*a51-a25*a32*a44*a51+      &
        a22*a35*a44*a51+a24*a32*a45*a51-a22*a34*a45*a51-a25*a34*a41*a52+   &
        a24*a35*a41*a52+a25*a31*a44*a52-a21*a35*a44*a52-a24*a31*a45*a52+   &
        a21*a34*a45*a52+a25*a32*a41*a54-a22*a35*a41*a54-a25*a31*a42*a54+   &
        a21*a35*a42*a54+a22*a31*a45*a54-a21*a32*a45*a54-a24*a32*a41*a55+   &
        a22*a34*a41*a55+a24*a31*a42*a55-a21*a34*a42*a55-a22*a31*a44*a55+   &
        a21*a32*a44*a55

    cofactor(2,3) = -a15*a34*a42*a51+a14*a35*a42*a51+a15*a32*a44*a51-     &
        a12*a35*a44*a51-a14*a32*a45*a51+a12*a34*a45*a51+a15*a34*a41*a52-   &
        a14*a35*a41*a52-a15*a31*a44*a52+a11*a35*a44*a52+a14*a31*a45*a52-   &
        a11*a34*a45*a52-a15*a32*a41*a54+a12*a35*a41*a54+a15*a31*a42*a54-   &
        a11*a35*a42*a54-a12*a31*a45*a54+a11*a32*a45*a54+a14*a32*a41*a55-   &
        a12*a34*a41*a55-a14*a31*a42*a55+a11*a34*a42*a55+a12*a31*a44*a55-   &
        a11*a32*a44*a55

    cofactor(3,3) = a15*a24*a42*a51-a14*a25*a42*a51-a15*a22*a44*a51+      &
        a12*a25*a44*a51+a14*a22*a45*a51-a12*a24*a45*a51-a15*a24*a41*a52+   &
        a14*a25*a41*a52+a15*a21*a44*a52-a11*a25*a44*a52-a14*a21*a45*a52+   &
        a11*a24*a45*a52+a15*a22*a41*a54-a12*a25*a41*a54-a15*a21*a42*a54+   &
        a11*a25*a42*a54+a12*a21*a45*a54-a11*a22*a45*a54-a14*a22*a41*a55+   &
        a12*a24*a41*a55+a14*a21*a42*a55-a11*a24*a42*a55-a12*a21*a44*a55+   &
        a11*a22*a44*a55

    cofactor(4,3) = -a15*a24*a32*a51+a14*a25*a32*a51+a15*a22*a34*a51-     &
        a12*a25*a34*a51-a14*a22*a35*a51+a12*a24*a35*a51+a15*a24*a31*a52-   &
        a14*a25*a31*a52-a15*a21*a34*a52+a11*a25*a34*a52+a14*a21*a35*a52-   &
        a11*a24*a35*a52-a15*a22*a31*a54+a12*a25*a31*a54+a15*a21*a32*a54-   &
        a11*a25*a32*a54-a12*a21*a35*a54+a11*a22*a35*a54+a14*a22*a31*a55-   &
        a12*a24*a31*a55-a14*a21*a32*a55+a11*a24*a32*a55+a12*a21*a34*a55-   &
        a11*a22*a34*a55

    cofactor(5,3) = a15*a24*a32*a41-a14*a25*a32*a41-a15*a22*a34*a41+      &
        a12*a25*a34*a41+a14*a22*a35*a41-a12*a24*a35*a41-a15*a24*a31*a42+   &
        a14*a25*a31*a42+a15*a21*a34*a42-a11*a25*a34*a42-a14*a21*a35*a42+   &
        a11*a24*a35*a42+a15*a22*a31*a44-a12*a25*a31*a44-a15*a21*a32*a44+   &
        a11*a25*a32*a44+a12*a21*a35*a44-a11*a22*a35*a44-a14*a22*a31*a45+   &
        a12*a24*a31*a45+a14*a21*a32*a45-a11*a24*a32*a45-a12*a21*a34*a45+   &
        a11*a22*a34*a45

    cofactor(1,4) = -a25*a33*a42*a51+a23*a35*a42*a51+a25*a32*a43*a51-     &
        a22*a35*a43*a51-a23*a32*a45*a51+a22*a33*a45*a51+a25*a33*a41*a52-   &
        a23*a35*a41*a52-a25*a31*a43*a52+a21*a35*a43*a52+a23*a31*a45*a52-   &
        a21*a33*a45*a52-a25*a32*a41*a53+a22*a35*a41*a53+a25*a31*a42*a53-   &
        a21*a35*a42*a53-a22*a31*a45*a53+a21*a32*a45*a53+a23*a32*a41*a55-   &
        a22*a33*a41*a55-a23*a31*a42*a55+a21*a33*a42*a55+a22*a31*a43*a55-   &
        a21*a32*a43*a55

    cofactor(2,4) = a15*a33*a42*a51-a13*a35*a42*a51-a15*a32*a43*a51+      &
        a12*a35*a43*a51+a13*a32*a45*a51-a12*a33*a45*a51-a15*a33*a41*a52+   &
        a13*a35*a41*a52+a15*a31*a43*a52-a11*a35*a43*a52-a13*a31*a45*a52+   &
        a11*a33*a45*a52+a15*a32*a41*a53-a12*a35*a41*a53-a15*a31*a42*a53+   &
        a11*a35*a42*a53+a12*a31*a45*a53-a11*a32*a45*a53-a13*a32*a41*a55+   &
        a12*a33*a41*a55+a13*a31*a42*a55-a11*a33*a42*a55-a12*a31*a43*a55+   &
        a11*a32*a43*a55

    cofactor(3,4) = -a15*a23*a42*a51+a13*a25*a42*a51+a15*a22*a43*a51-     &
        a12*a25*a43*a51-a13*a22*a45*a51+a12*a23*a45*a51+a15*a23*a41*a52-   &
        a13*a25*a41*a52-a15*a21*a43*a52+a11*a25*a43*a52+a13*a21*a45*a52-   &
        a11*a23*a45*a52-a15*a22*a41*a53+a12*a25*a41*a53+a15*a21*a42*a53-   &
        a11*a25*a42*a53-a12*a21*a45*a53+a11*a22*a45*a53+a13*a22*a41*a55-   &
        a12*a23*a41*a55-a13*a21*a42*a55+a11*a23*a42*a55+a12*a21*a43*a55-   &
        a11*a22*a43*a55

    cofactor(4,4) = a15*a23*a32*a51-a13*a25*a32*a51-a15*a22*a33*a51+      &
        a12*a25*a33*a51+a13*a22*a35*a51-a12*a23*a35*a51-a15*a23*a31*a52+   &
        a13*a25*a31*a52+a15*a21*a33*a52-a11*a25*a33*a52-a13*a21*a35*a52+   &
        a11*a23*a35*a52+a15*a22*a31*a53-a12*a25*a31*a53-a15*a21*a32*a53+   &
        a11*a25*a32*a53+a12*a21*a35*a53-a11*a22*a35*a53-a13*a22*a31*a55+   &
        a12*a23*a31*a55+a13*a21*a32*a55-a11*a23*a32*a55-a12*a21*a33*a55+   &
        a11*a22*a33*a55

    cofactor(5,4) = -a15*a23*a32*a41+a13*a25*a32*a41+a15*a22*a33*a41-     &
        a12*a25*a33*a41-a13*a22*a35*a41+a12*a23*a35*a41+a15*a23*a31*a42-   &
        a13*a25*a31*a42-a15*a21*a33*a42+a11*a25*a33*a42+a13*a21*a35*a42-   &
        a11*a23*a35*a42-a15*a22*a31*a43+a12*a25*a31*a43+a15*a21*a32*a43-   &
        a11*a25*a32*a43-a12*a21*a35*a43+a11*a22*a35*a43+a13*a22*a31*a45-   &
        a12*a23*a31*a45-a13*a21*a32*a45+a11*a23*a32*a45+a12*a21*a33*a45-   &
        a11*a22*a33*a45

    cofactor(1,5) = a24*a33*a42*a51-a23*a34*a42*a51-a24*a32*a43*a51+      &
        a22*a34*a43*a51+a23*a32*a44*a51-a22*a33*a44*a51-a24*a33*a41*a52+   &
        a23*a34*a41*a52+a24*a31*a43*a52-a21*a34*a43*a52-a23*a31*a44*a52+   &
        a21*a33*a44*a52+a24*a32*a41*a53-a22*a34*a41*a53-a24*a31*a42*a53+   &
        a21*a34*a42*a53+a22*a31*a44*a53-a21*a32*a44*a53-a23*a32*a41*a54+   &
        a22*a33*a41*a54+a23*a31*a42*a54-a21*a33*a42*a54-a22*a31*a43*a54+   &
        a21*a32*a43*a54

    cofactor(2,5) = -a14*a33*a42*a51+a13*a34*a42*a51+a14*a32*a43*a51-     &
        a12*a34*a43*a51-a13*a32*a44*a51+a12*a33*a44*a51+a14*a33*a41*a52-   &
        a13*a34*a41*a52-a14*a31*a43*a52+a11*a34*a43*a52+a13*a31*a44*a52-   &
        a11*a33*a44*a52-a14*a32*a41*a53+a12*a34*a41*a53+a14*a31*a42*a53-   &
        a11*a34*a42*a53-a12*a31*a44*a53+a11*a32*a44*a53+a13*a32*a41*a54-   &
        a12*a33*a41*a54-a13*a31*a42*a54+a11*a33*a42*a54+a12*a31*a43*a54-   &
        a11*a32*a43*a54

    cofactor(3,5) = a14*a23*a42*a51-a13*a24*a42*a51-a14*a22*a43*a51+      &
        a12*a24*a43*a51+a13*a22*a44*a51-a12*a23*a44*a51-a14*a23*a41*a52+   &
        a13*a24*a41*a52+a14*a21*a43*a52-a11*a24*a43*a52-a13*a21*a44*a52+   &
        a11*a23*a44*a52+a14*a22*a41*a53-a12*a24*a41*a53-a14*a21*a42*a53+   &
        a11*a24*a42*a53+a12*a21*a44*a53-a11*a22*a44*a53-a13*a22*a41*a54+   &
        a12*a23*a41*a54+a13*a21*a42*a54-a11*a23*a42*a54-a12*a21*a43*a54+   &
        a11*a22*a43*a54

    cofactor(4,5) = -a14*a23*a32*a51+a13*a24*a32*a51+a14*a22*a33*a51-     &
        a12*a24*a33*a51-a13*a22*a34*a51+a12*a23*a34*a51+a14*a23*a31*a52-   &
        a13*a24*a31*a52-a14*a21*a33*a52+a11*a24*a33*a52+a13*a21*a34*a52-   &
        a11*a23*a34*a52-a14*a22*a31*a53+a12*a24*a31*a53+a14*a21*a32*a53-   &
        a11*a24*a32*a53-a12*a21*a34*a53+a11*a22*a34*a53+a13*a22*a31*a54-   &
        a12*a23*a31*a54-a13*a21*a32*a54+a11*a23*a32*a54+a12*a21*a33*a54-   &
        a11*a22*a33*a54

    cofactor(5,5) = a14*a23*a32*a41-a13*a24*a32*a41-a14*a22*a33*a41+      &
        a12*a24*a33*a41+a13*a22*a34*a41-a12*a23*a34*a41-a14*a23*a31*a42+   &
        a13*a24*a31*a42+a14*a21*a33*a42-a11*a24*a33*a42-a13*a21*a34*a42+   &
        a11*a23*a34*a42+a14*a22*a31*a43-a12*a24*a31*a43-a14*a21*a32*a43+   &
        a11*a24*a32*a43+a12*a21*a34*a43-a11*a22*a34*a43-a13*a22*a31*a44+   &
        a12*a23*a31*a44+a13*a21*a32*a44-a11*a23*a32*a44-a12*a21*a33*a44+   &
        a11*a22*a33*a44

    ainv = transpose(cofactor)/det
    flag = 1

    return
    end subroutine m55inv