
    !program m22inv_main
    !
    !implicit none
    !
    !integer :: i, j
    !double precision, dimension(2,2) :: mat, matinv
    !integer :: flag
    !
    !logical :: m22inv
    !
    !  
    !write (unit=*, fmt='(/a/)') ' enter matrix:'
    !
    !do i = 1, 2
    !do j = 1, 2
    !write (unit=*, fmt='(a,i1,1h,,i1,a)', advance='no') ' a(', i, j, ') = '
    !read (unit=*, fmt=*) mat(i,j)
    !end do
    !end do
    !
    !call m22inv (mat,matinv,flag)
    !
    !if( flag .eq. 1 )then
    !!write (*,'(/a/)') ' inverse:'
    !!write (*,'(2es25.15)') ((matinv(i,j), j=1,2), i=1,2)
    !elseif( flag .eq. 0 )then
    !write (unit=*, fmt='(/a)') ' singular matrix.'
    !end if
    !
    !stop
    !end program m22inv_main


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

