  subroutine cholesky(n,A)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in out) :: A(n,n)

    integer i,j,k

    real(8) Ljj, Lij

    do j=1,n
       Ljj=A(j,j)
       do k=1,j-1
          Ljj=Ljj-A(j,k)**2
       enddo
       Ljj=sqrt(Ljj)
       do i=j+1,n
          Lij=A(i,j)
          do k=1,j-1
             Lij=Lij-A(i,k)*A(j,k)
          enddo
          A(i,j)=Lij/Ljj
       enddo
       A(j,j)=Ljj
    enddo

    do i=1,n
       do j=i+1,n
          A(i,j)=0.D0
       enddo
    enddo
    
    return
  end subroutine cholesky

subroutine choleskysolve(n,A,b,x)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: A(n,n),b(n)
  real(8), intent(out) :: x(n)

  integer i,j

  do i=1,n
     X(i)=B(i)
     do j=1,i-1
        X(i)=X(i)-A(i,j)*X(j)
     enddo
     X(i)=X(i)/A(i,i)
  enddo
  
  do i=n,1,-1
     do j=i+1,n
        X(i)=X(i)-A(j,i)*X(j)
     enddo
     X(i)=X(i)/A(i,i)
  enddo
  
  return
end subroutine choleskysolve

subroutine invertcholesky(n,A,Ai)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: A(n,n)
  
  real(8), intent(out) :: Ai(n,n)
  
  integer i,j
  real(8) X(n), b(n)
  
  do i=1,n
     B(:)=0.d0
     B(i)=1.D0
     call choleskysolve(n,A,b,x)
     Ai(:,i)=x(:)
  enddo
  
  do i=1,n
     do j=i+1,n
        Ai(i,j)=Ai(j,i)
     enddo
  enddo
  
  return
end subroutine invertcholesky
  
subroutine invertsymmetric(n,A)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in out) :: A(n,n)

  real(8) Ai(n,n)

  call cholesky(n,A)

  call invertcholesky(n,A,Ai)

  A=Ai

  return
end subroutine invertsymmetric

subroutine symmetricsolve(n,A,x)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in out) :: A(n,n)
  real(8), intent(in out) :: X(n)

  real(8) B(n)

  B=X

  call cholesky(n,A)
  call choleskysolve(n,A,B,x)
  
  return
end subroutine symmetricsolve

subroutine matrixmult(n,m,A,p,B,C)
  implicit none
  integer, intent(in) :: n,m,p
  real(8), intent(in) :: A(n,m), B(m,p)
  real(8), intent(out) :: C(n,p)

  integer i,j,k

  do i=1,n
     do j=1,p
        C(i,j)=0.D0
        do k=1,m
           C(i,j)=C(i,j)+A(i,k)*B(k,j)
        enddo
     enddo
  enddo

  return
end subroutine matrixmult

subroutine symmatrixmult(n,A,p,B,C)
  implicit none
  integer, intent(in) :: n,p
  real(8), intent(in) :: A(n,n), B(n,p)
  real(8), intent(out) :: C(n,p)

  integer i,j,k

  do i=1,n
     do j=1,p
        C(i,j)=0.D0
        do k=1,i
           C(i,j)=C(i,j)+A(i,k)*B(k,j)
        enddo
        do k=i+1,n
           C(i,j)=C(i,j)+A(k,i)*B(k,j)
        enddo
     enddo
  enddo

  return
end subroutine symmatrixmult

subroutine matrixmulttrans(n,m,A,p,B,C)
  implicit none
  integer, intent(in) :: n,m,p
  real(8), intent(in) :: A(n,m), B(p,m)
  real(8), intent(out) :: C(p,n)

  integer i,j,k

  do j=1,n
     do i=1,p
        C(i,j)=0.D0
        do k=1,m
           C(j,i)=C(j,i)+A(i,k)*B(j,k)
        enddo
     enddo
  enddo

  return
end subroutine matrixmulttrans

subroutine symmatrixmulttrans(n,A,p,B,C)
  implicit none
  integer, intent(in) :: n,p
  real(8), intent(in) :: A(n,n), B(p,n)
  real(8), intent(out) :: C(p,n)

  integer i,j,k

  do i=1,n
     do j=1,p
        C(i,j)=0.D0
        do k=1,i
           C(j,i)=C(j,i)+A(i,k)*B(j,k)
        enddo
        do k=i+1,n
           C(j,i)=C(j,i)+A(k,i)*B(j,k)
        enddo
     enddo
  enddo

  return
end subroutine symmatrixmulttrans

subroutine matvec(n,m,A,x,B)
  implicit none
  integer, intent(in) :: n,m
  real(8), intent(in) :: A(n,m),x(m)
  real(8), intent(out) :: B(n)

  integer i,j

  do i=1,n
     B(i)=0.D0
     do j=1,m
        B(i)=B(i)+A(i,j)*X(j)
     enddo
  enddo

  return
end subroutine matvec

subroutine symmatvec(n,A,x,B)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: A(n,n),x(n)
  real(8), intent(out) :: B(n)

  integer i,j

  do i=1,n
     B(i)=0.D0
     do j=1,i
        B(i)=B(i)+A(i,j)*X(j)
     enddo
     do j=i+1,n
        B(i)=B(i)+A(j,i)*X(j)
     enddo
  enddo

  return
end subroutine symmatvec

subroutine matvectrans(n,m,A,x,B)
  implicit none
  integer, intent(in) :: n,m
  real(8), intent(in) :: A(n,m),x(n)
  real(8), intent(out) :: B(m)

  integer i,j

  do i=1,m
     B(i)=0.D0
     do j=1,n
        B(i)=B(i)+A(j,i)*X(j)
     enddo
  enddo

  return
end subroutine matvectrans