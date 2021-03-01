!
!      subroutine LSCoeff(n,x0,y0,x1,y1,cx,cy)
!      use prm_flow
!!     ======================================================================
!!     Find coefficients for first derivative using least squares
!!     n           = Number of neighbours
!!     x0,y0       = coordinates of point at which grad is required
!!     x1(n),y1(n) = coordinates of neighbouring points
!!     cx(n),cy(n) = coefficients
!!     ======================================================================
!      implicit none
!      integer          :: n
!      double precision :: x0, y0, x1(*), y1(*), cx(*), cy(*)
!      integer          :: flag
!      integer          :: i, j
!      double precision :: dx, dy, dr, w(n), det
!      double precision :: ainv(2,2)
!      double precision :: c(n,2), a(2,2)
!      
!      a(:,:) = 0.0d0 ; c(:,:) = 0.0d0
!      
!      do i=1,n ! for all neighbours
!      
!        dx     = x1(i) - x0            ! check
!        dy     = y1(i) - y0            ! check
!        dr     = dsqrt(dx**2 + dy**2)  
!        w(i)   = 1.0d0/dr
!         
!!       MATRIX C individual components i=1,n
!        c(i,1) = w(i)*dx
!        c(i,2) = w(i)*dy
!         
!!       LEAST SQUARES MATRIX COMPONENTS 
!        a(1,1) = a(1,1) +       w(i)*dx**2
!        a(1,2) = a(1,2) +       w(i)*dx*dy
!        a(2,1) = a(1,2)
!        a(2,2) = a(2,2) +       w(i)*(dy**2)
!      
!      enddo   
!
!      CALL m22inv(a,ainv,flag)
!
!      do i=1,n ! for all neighbours
!      cx(i) = 0.0d0
!      cy(i) = 0.0d0
!      do j=1,2
!       cx(i) = cx(i) + ainv(1,j)*c(i,j)
!       cy(i) = cy(i) + ainv(2,j)*c(i,j)
!      enddo
!      enddo
!
!      return
!      End Subroutine LSCoeff
!      
!      subroutine LSCoeff2(n,x0,y0,x1,y1,cx,cy)
!      use prm_flow
!!     ======================================================================
!!     Find coefficients for first derivative using least squares
!!     n           = Number of neighbours
!!     x0,y0       = coordinates of point at which grad is required
!!     x1(n),y1(n) = coordinates of neighbouring points
!!     cx(n),cy(n) = coefficients
!!     ======================================================================
!      implicit none
!      integer          :: n
!      double precision :: x0, y0, x1(*), y1(*), cx(*), cy(*)
!      integer          :: flag
!      integer          :: i, j
!      double precision :: dx, dy, dr, w(n), det
!      double precision :: ainv(5,5), c(n,5), a(5,5)
!
!      a(:,:) = 0.0d0 ; c(:,:) = 0.0d0
!
!      do i=1,n ! for all neighbours
!      
!         dx     = x1(i) - x0            ! check
!         dy     = y1(i) - y0            ! check
!         dr     = dsqrt(dx**2 + dy**2)  
!         w(i)   = 1.0d0/dr
!         
!!        MATRIX C individual components
!         c(i,1) = w(i)*dx
!         c(i,2) = w(i)*dy
!         c(i,3) = w(i)*0.5d0*dx**2
!         c(i,4) = w(i)*0.5d0*dy**2
!         c(i,5) = w(i)*dx*dy
!         
!!        LEAST SQUARES MATRIX COMPONENTS
!         a(1,1) = a(1,1) +        w(i)**2*dx**2
!         a(1,2) = a(1,2) +        w(i)**2*dx*dy
!         a(1,3) = a(1,3) + 0.50d0*w(i)**2*dx**3
!         a(1,4) = a(1,4) + 0.50d0*w(i)**2*(dx**2)*dy
!         a(1,5) = a(1,5) +        w(i)**2*(dx**2)*dy
!         
!         a(2,1) = a(1,2)
!         a(2,2) = a(2,2) +        w(i)**2*(dy**2)
!         a(2,3) = a(1,4) 
!         a(2,4) = a(2,4) + 0.50d0*w(i)**2*(dy**3)
!         a(2,5) = a(2,5) +        w(i)**2*dx*(dy**2)
!         
!         a(3,1) = a(3,1) + 0.50d0*w(i)**2*(dx**3)
!         a(3,2) = a(1,4)
!         a(3,3) = a(3,3) + 0.25d0*w(i)**2*(dx**4)
!         a(3,4) = a(3,4) + 0.25d0*w(i)**2*(dx**2)*(dy**2)
!         a(3,5) = a(3,5) + 0.50d0*w(i)**2*(dx**3)*dy
!         
!         a(4,1) = a(1,4)
!         a(4,2) = a(2,4)
!         a(4,3) = a(3,4)
!         a(4,4) = a(3,3)
!         a(4,5) = a(4,5) + 0.50d0*w(i)**2*dx*(dy**3)
!         
!         a(5,1) = a(1,5)
!         a(5,2) = a(2,5)
!         a(5,3) = a(3,5)
!         a(5,4) = a(4,5)
!         a(5,5) = a(5,5) +        w(i)**2*(dx**2)*(dy**2)
!         
!
!      enddo
!      
!      call inverse(a,ainv,5)
!      !call m55inv(a,ainv,flag)
!      
!!     Check
!      !if (flag .eq. 0) then
!      ! write (*,'(/a)') ' singular matrix.'
!      ! pause
!      !end if
!
!      do i=1,n ! for all neighbours
!       cx(i) = 0.0d0
!       cy(i) = 0.0d0
!       do j=1,5
!        cx(i) = cx(i) + ainv(1,j)*c(i,j)*w(i) ! Fx0 component
!        cy(i) = cy(i) + ainv(2,j)*c(i,j)*w(i) ! Fy0 component
!       enddo
!      enddo
!
!      return
!      End Subroutine LSCoeff2
!      
!      
!    subroutine inverse(a,c,n)
!    !============================================================
!    ! Inverse matrix
!    ! Method: Based on Doolittle LU factorization for Ax=b
!    ! Alex G. December 2009
!    !-----------------------------------------------------------
!    ! input ...
!    ! a(n,n) - array of coefficients for matrix A
!    ! n      - dimension
!    ! output ...
!    ! c(n,n) - inverse matrix of A
!    ! comments ...
!    ! the original matrix a(n,n) will be destroyed 
!    ! during the calculation
!    !===========================================================
!    implicit none 
!    integer n
!    double precision a(n,n), c(n,n)
!    double precision L(n,n), U(n,n), b(n), d(n), x(n)
!    double precision coeff
!    integer i, j, k
!
!    ! step 0: initialization for matrices L and U and b
!    ! Fortran 90/95 aloows such operations on matrices
!    L=0.0
!    U=0.0
!    b=0.0
!
!    ! step 1: forward elimination
!    do k=1, n-1
!    do i=k+1,n
!        coeff=a(i,k)/a(k,k)
!        L(i,k) = coeff
!        do j=k+1,n
!            a(i,j) = a(i,j)-coeff*a(k,j)
!        end do
!    end do
!    end do
!
!    ! Step 2: prepare L and U matrices 
!    ! L matrix is a matrix of the elimination coefficient
!    ! + the diagonal elements are 1.0
!    do i=1,n
!    L(i,i) = 1.0
!    end do
!    ! U matrix is the upper triangular part of A
!    do j=1,n
!    do i=1,j
!    U(i,j) = a(i,j)
!    end do
!    end do
!
!    ! Step 3: compute columns of the inverse matrix C
!    do k=1,n
!    b(k)=1.0
!    d(1) = b(1)
!    ! Step 3a: Solve Ld=b using the forward substitution
!    do i=2,n
!    d(i)=b(i)
!    do j=1,i-1
!        d(i) = d(i) - L(i,j)*d(j)
!    end do
!    end do
!    ! Step 3b: Solve Ux=d using the back substitution
!    x(n)=d(n)/U(n,n)
!    do i = n-1,1,-1
!    x(i) = d(i)
!    do j=n,i+1,-1
!        x(i)=x(i)-U(i,j)*x(j)
!    end do
!    x(i) = x(i)/u(i,i)
!    end do
!    ! Step 3c: fill the solutions x(n) into column k of C
!    do i=1,n
!    c(i,k) = x(i)
!    end do
!    b(k)=0.0
!    end do
!    
!    
!    end subroutine inverse
!          

      subroutine LSCoeff(n,x0,y0,x1,y1,cx,cy)
      use prm_flow
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
        a(1,1) = a(1,1) +       w(i)**2*dx**2
        a(1,2) = a(1,2) +       w(i)**2*dx*dy
        a(2,1) = a(1,2)
        a(2,2) = a(2,2) +       w(i)**2*(dy**2)
      
      enddo   

      CALL m22inv(a,ainv,flag)

      do i=1,n ! for all neighbours
      cx(i) = 0.0d0
      cy(i) = 0.0d0
      do j=1,2
       cx(i) = cx(i) + ainv(1,j)*c(i,j)*w(i)
       cy(i) = cy(i) + ainv(2,j)*c(i,j)*w(i)
      enddo
      enddo

      return
      End Subroutine LSCoeff
      
!      subroutine LSCoeff2(n,x0,y0,x1,y1,cx,cy)
!      use prm_flow
!!     ======================================================================
!!     Find coefficients for first derivative using least squares
!!     n           = Number of neighbours
!!     x0,y0       = coordinates of point at which grad is required
!!     x1(n),y1(n) = coordinates of neighbouring points
!!     cx(n),cy(n) = coefficients
!!     ======================================================================
!      implicit none
!      integer          :: n
!      double precision :: x0, y0, x1(*), y1(*), cx(*), cy(*)
!      integer          :: flag
!      integer          :: i, j
!      double precision :: dx, dy, dr, w(n), det
!      double precision :: ainv(5,5), c(n,5), a(5,5)
!
!      a(:,:) = 0.0d0 ; c(:,:) = 0.0d0
!
!      do i=1,n ! for all neighbours
!      
!         dx     = x1(i) - x0            ! check
!         dy     = y1(i) - y0            ! check
!         dr     = dsqrt(dx**2 + dy**2)  
!         w(i)   = 1.0d0/dr
!         
!!        MATRIX C individual components
!         c(i,1) = w(i)*dx
!         c(i,2) = w(i)*dy
!         c(i,3) = w(i)*0.5d0*dx**2
!         c(i,4) = w(i)*0.5d0*dy**2
!         c(i,5) = w(i)*dx*dy
!         
!!        LEAST SQUARES MATRIX COMPONENTS
!         a(1,1) = a(1,1) +        w(i)**2*dx**2
!         a(1,2) = a(1,2) +        w(i)**2*dx*dy
!         a(1,3) = a(1,3) + 0.50d0*w(i)**2*dx**3
!         a(1,4) = a(1,4) + 0.50d0*w(i)**2*(dx**2)*dy
!         a(1,5) = a(1,5) +        w(i)**2*(dx**2)*dy
!         
!         a(2,1) = a(1,2)
!         a(2,2) = a(2,2) +        w(i)**2*(dy**2)
!         a(2,3) = a(1,4) 
!         a(2,4) = a(2,4) + 0.50d0*w(i)**2*(dy**3)
!         a(2,5) = a(2,5) +        w(i)**2*dx*(dy**2)
!         
!         a(3,1) = a(3,1) + 0.50d0*w(i)**2*(dx**3)
!         a(3,2) = a(1,4)
!         a(3,3) = a(3,3) + 0.25d0*w(i)**2*(dx**4)
!         a(3,4) = a(3,4) + 0.25d0*w(i)**2*(dx**2)*(dy**2)
!         a(3,5) = a(3,5) + 0.50d0*w(i)**2*(dx**3)*dy
!         
!         a(4,1) = a(1,4)
!         a(4,2) = a(2,4)
!         a(4,3) = a(3,4)
!         a(4,4) = a(3,3)
!         a(4,5) = a(4,5) + 0.50d0*w(i)**2*dx*(dy**3)
!         
!         a(5,1) = a(1,5)
!         a(5,2) = a(2,5)
!         a(5,3) = a(3,5)
!         a(5,4) = a(4,5)
!         a(5,5) = a(5,5) +        w(i)**2*(dx**2)*(dy**2)
!         
!
!      enddo
!      
!      call inverse(a,ainv,5)
!      !call m55inv(a,ainv,flag)
!      
!!     Check
!      !if (flag .eq. 0) then
!      ! write (*,'(/a)') ' singular matrix.'
!      ! pause
!      !end if
!
!      do i=1,n ! for all neighbours
!       cx(i) = 0.0d0
!       cy(i) = 0.0d0
!       do j=1,5
!        cx(i) = cx(i) + ainv(1,j)*c(i,j)*w(i) ! Fx0 component
!        cy(i) = cy(i) + ainv(2,j)*c(i,j)*w(i) ! Fy0 component
!       enddo
!      enddo
!
!      return
!      End Subroutine LSCoeff2
      
      
!    subroutine inverse(a,c,n)
!    !============================================================
!    ! Inverse matrix
!    ! Method: Based on Doolittle LU factorization for Ax=b
!    ! Alex G. December 2009
!    !-----------------------------------------------------------
!    ! input ...
!    ! a(n,n) - array of coefficients for matrix A
!    ! n      - dimension
!    ! output ...
!    ! c(n,n) - inverse matrix of A
!    ! comments ...
!    ! the original matrix a(n,n) will be destroyed 
!    ! during the calculation
!    !===========================================================
!    implicit none 
!    integer n
!    double precision a(n,n), c(n,n)
!    double precision L(n,n), U(n,n), b(n), d(n), x(n)
!    double precision coeff
!    integer i, j, k
!
!    ! step 0: initialization for matrices L and U and b
!    ! Fortran 90/95 aloows such operations on matrices
!    L=0.0
!    U=0.0
!    b=0.0
!
!    ! step 1: forward elimination
!    do k=1, n-1
!    do i=k+1,n
!        coeff=a(i,k)/a(k,k)
!        L(i,k) = coeff
!        do j=k+1,n
!            a(i,j) = a(i,j)-coeff*a(k,j)
!        end do
!    end do
!    end do
!
!    ! Step 2: prepare L and U matrices 
!    ! L matrix is a matrix of the elimination coefficient
!    ! + the diagonal elements are 1.0
!    do i=1,n
!    L(i,i) = 1.0
!    end do
!    ! U matrix is the upper triangular part of A
!    do j=1,n
!    do i=1,j
!    U(i,j) = a(i,j)
!    end do
!    end do
!
!    ! Step 3: compute columns of the inverse matrix C
!    do k=1,n
!    b(k)=1.0
!    d(1) = b(1)
!    ! Step 3a: Solve Ld=b using the forward substitution
!    do i=2,n
!    d(i)=b(i)
!    do j=1,i-1
!        d(i) = d(i) - L(i,j)*d(j)
!    end do
!    end do
!    ! Step 3b: Solve Ux=d using the back substitution
!    x(n)=d(n)/U(n,n)
!    do i = n-1,1,-1
!    x(i) = d(i)
!    do j=n,i+1,-1
!        x(i)=x(i)-U(i,j)*x(j)
!    end do
!    x(i) = x(i)/u(i,i)
!    end do
!    ! Step 3c: fill the solutions x(n) into column k of C
!    do i=1,n
!    c(i,k) = x(i)
!    end do
!    b(k)=0.0
!    end do
!    
!    
!    end subroutine inverse
          
!