

    !program m33inv_main
    !
    !implicit none
    !
    !integer :: i, j
    !double precision, dimension(3,3) :: mat, matinv
    !integer :: flag
    !
    !write (*,'(/a/)') ' enter matrix:'
    !
    !do i = 1, 3
    !do j = 1, 3
    !write (unit=*, fmt='(a,i1,1h,,i1,a)', advance='no') ' a(', i, j, ') = '
    !read (unit=*, fmt=*) mat(i,j)
    !end do
    !end do
    !
    !call m33inv(mat, matinv,flag)
    !
    !if(flag .eq. 1)then
    ! !write (*,'(/a/)') ' inverse:'
    ! !write (*,'(3es25.15)') ((matinv(i,j), j=1,3), i=1,3)
    !elseif( flag .eq. 0 )then
    ! write (*,'(/a)') ' singular matrix.'
    !end if
    !
    !stop
    !
    !end program m33inv_main

    subroutine m33inv (a,ainv,flag)
    implicit none

    double precision, dimension(3,3), intent(in)  :: a
    double precision, dimension(3,3), intent(out) :: ainv
    integer, intent(out) :: flag

    double precision, parameter :: eps = 1.0d-10
    double precision :: det
    double precision, dimension(3,3) :: cofactor


    det =   a(1,1)*a(2,2)*a(3,3)  &
        - a(1,1)*a(2,3)*a(3,2)  &
        - a(1,2)*a(2,1)*a(3,3)  &
        + a(1,2)*a(2,3)*a(3,1)  &
        + a(1,3)*a(2,1)*a(3,2)  &
        - a(1,3)*a(2,2)*a(3,1)

!   Singular Matrix Condition Check
    if(abs(det) .le. eps)then
    ainv = 0.0d0
    flag = 0 ! non invertible matrix
    return
    end if

    cofactor(1,1) = +(a(2,2)*a(3,3)-a(2,3)*a(3,2))
    cofactor(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
    cofactor(1,3) = +(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    cofactor(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
    cofactor(2,2) = +(a(1,1)*a(3,3)-a(1,3)*a(3,1))
    cofactor(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
    cofactor(3,1) = +(a(1,2)*a(2,3)-a(1,3)*a(2,2))
    cofactor(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
    cofactor(3,3) = +(a(1,1)*a(2,2)-a(1,2)*a(2,1))

    ainv = transpose(cofactor)/det

    flag = 1

    return
    end subroutine m33inv