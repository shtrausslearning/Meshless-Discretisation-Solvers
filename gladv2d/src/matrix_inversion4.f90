
    !program m44inv_main
    !implicit none
    !
    !integer :: i, j
    !double precision, dimension(4,4) :: mat, matinv
    !integer :: flag
    !
    !write (*,'(/a/)') ' enter matrix:'
    !
    !do i = 1, 4
    !do j = 1, 4
    !write (*,'(a,i1,1h,,i1,a)', advance='no') ' a(', i, j, ') = '
    !read(*,*) mat(i,j)
    !end do
    !end do
    !
    !call m44inv (mat,matinv,flag)
    !
    !if (flag .eq. 1) then
    ! !write (*,'(/a/)') ' inverse:'
    ! !write (*,'(4es25.15)') ((matinv(i,j), j=1,4), i=1,4)
    !else
    ! write (*,'(/a)') ' singular matrix.'
    !end if
    !
    !stop
    !
    !end program m44inv_main

    subroutine m44inv (a,ainv,flag)
    implicit none

    double precision, dimension(4,4), intent(in)  :: a
    double precision, dimension(4,4), intent(out) :: ainv
    integer, intent(out) :: flag

    double precision, parameter :: eps = 1.0d-10
    double precision :: det
    double precision, dimension(4,4) :: cofactor

    det =  a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)*(a(3,4)*a(4,2)-a(3,2)*a(4,4))+a(2,4)*(a(3,2)*a(4,3)- &
            a(3,3)*a(4,2)))-a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)*(a(3,4)*a(4,1)-a(3,1)*a(4,4))+ &
            a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))+a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,2)*(a(3,4)*a(4,1)- &
            a(3,1)*a(4,4))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))-a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))+ &
            a(2,2)*(a(3,3)*a(4,1)-a(3,1)*a(4,3))+a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))

!   Singular Matrix Condition Check
    if( abs(det) .le. eps )then
     ainv = 0.0d0
     flag = 0 ! Singular matrix
     return
    endif

    cofactor(1,1) = a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)*(a(3,4)*a(4,2)-a(3,2)*a(4,4))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))
    cofactor(1,2) = a(2,1)*(a(3,4)*a(4,3)-a(3,3)*a(4,4))+a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,3)*a(4,1)-a(3,1)*a(4,3))
    cofactor(1,3) = a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,2)*(a(3,4)*a(4,1)-a(3,1)*a(4,4))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))
    cofactor(1,4) = a(2,1)*(a(3,3)*a(4,2)-a(3,2)*a(4,3))+a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(2,3)*(a(3,2)*a(4,1)-a(3,1)*a(4,2))
    cofactor(2,1) = a(1,2)*(a(3,4)*a(4,3)-a(3,3)*a(4,4))+a(1,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(1,4)*(a(3,3)*a(4,2)-a(3,2)*a(4,3))
    cofactor(2,2) = a(1,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(1,3)*(a(3,4)*a(4,1)-a(3,1)*a(4,4))+a(1,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))
    cofactor(2,3) = a(1,1)*(a(3,4)*a(4,2)-a(3,2)*a(4,4))+a(1,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(1,4)*(a(3,2)*a(4,1)-a(3,1)*a(4,2))
    cofactor(2,4) = a(1,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))+a(1,2)*(a(3,3)*a(4,1)-a(3,1)*a(4,3))+a(1,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))
    cofactor(3,1) = a(1,2)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))+a(1,3)*(a(2,4)*a(4,2)-a(2,2)*a(4,4))+a(1,4)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))
    cofactor(3,2) = a(1,1)*(a(2,4)*a(4,3)-a(2,3)*a(4,4))+a(1,3)*(a(2,1)*a(4,4)-a(2,4)*a(4,1))+a(1,4)*(a(2,3)*a(4,1)-a(2,1)*a(4,3))
    cofactor(3,3) = a(1,1)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))+a(1,2)*(a(2,4)*a(4,1)-a(2,1)*a(4,4))+a(1,4)*(a(2,1)*a(4,2)-a(2,2)*a(4,1))
    cofactor(3,4) = a(1,1)*(a(2,3)*a(4,2)-a(2,2)*a(4,3))+a(1,2)*(a(2,1)*a(4,3)-a(2,3)*a(4,1))+a(1,3)*(a(2,2)*a(4,1)-a(2,1)*a(4,2))
    cofactor(4,1) = a(1,2)*(a(2,4)*a(3,3)-a(2,3)*a(3,4))+a(1,3)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))+a(1,4)*(a(2,3)*a(3,2)-a(2,2)*a(3,3))
    cofactor(4,2) = a(1,1)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))+a(1,3)*(a(2,4)*a(3,1)-a(2,1)*a(3,4))+a(1,4)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
    cofactor(4,3) = a(1,1)*(a(2,4)*a(3,2)-a(2,2)*a(3,4))+a(1,2)*(a(2,1)*a(3,4)-a(2,4)*a(3,1))+a(1,4)*(a(2,2)*a(3,1)-a(2,1)*a(3,2))
    cofactor(4,4) = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))+a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

    ainv = transpose(cofactor)/det
    flag = 1

    return
    end subroutine m44inv

