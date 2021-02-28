
!* 4. module my_allocation
!*
!* This module defines some useful subroutines used for dynamic allocation.
!*
!*  - my_alloc_int_ptr       : Allocate/reallocate an integer 1D array
!*  - my_alloc_p2_ptr        : Allcoate/reallocate a real 1D array
!*  - my_alloc_p2_matrix_ptr : Allcoate/reallocate a real 2D array

 module edu2d_my_allocation

  implicit none

  private

  public :: my_alloc_int_ptr
  public :: my_alloc_p2_ptr
  public :: my_alloc_p2_matrix_ptr

  contains

!********************************************************************************
!* This subroutine is useful to expand or shrink integer arrays.
!*
!*  Array, x, will be allocated if the requested dimension is 1 (i.e., n=1)
!*  Array, x, will be expanded to the requested dimension, n, if (n > dim(x)).
!*  Array, x, will be shrinked to the requested dimension, n, if (n < dim(x)).
!*
!********************************************************************************
  subroutine my_alloc_int_ptr(x,n)

  implicit none
  integer, intent(in) :: n
  integer, dimension(:), pointer :: x
  integer, dimension(:), pointer :: temp
  integer :: i

  if (n <= 0) then
   write(101,*) "my_alloc_int_ptr received non-positive dimension. Stop."
   stop
  endif

! If not allocated, allocate and return
  if (.not.(associated(x))) then
   allocate(x(n))
   return
  endif

! If reallocation, create a pointer with a target of new dimension.
  allocate(temp(n))
   temp = 0

! (1) Expand the array dimension
  if ( n > size(x) ) then

   do i = 1, size(x)
    temp(i) = x(i)
   end do

! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
  else

   do i = 1, n
    temp(i) = x(i)
   end do

  endif

! Destroy the target of x
!  deallocate(x)

! Re-assign the pointer
   x => temp

  return

  end subroutine my_alloc_int_ptr
!********************************************************************************


!********************************************************************************
!* This subroutine is useful to expand or shrink real arrays.
!*
!*  Array, x, will be allocated if the requested dimension is 1 (i.e., n=1)
!*  Array, x, will be expanded to the requested dimension, n, if (n > dim(x)).
!*  Array, x, will be shrinked to the requested dimension, n, if (n < dim(x)).
!*
!********************************************************************************
  subroutine my_alloc_p2_ptr(x,n)

  use edu2d_constants   , only : p2

  implicit none
  integer, intent(in) :: n
  real(p2), dimension(:), pointer :: x
  real(p2), dimension(:), pointer :: temp
  integer :: i

  if (n <= 0) then
   write(101,*) "my_alloc_int_ptr received non-positive dimension. Stop."
   stop
  endif

! If not allocated, allocate and return
  if (.not.(associated(x))) then
   allocate(x(n))
   return
  endif

! If reallocation, create a pointer with a target of new dimension.
  allocate(temp(n))
   temp = 0

! (1) Expand the array dimension
  if ( n > size(x) ) then

   do i = 1, size(x)
    temp(i) = x(i)
   end do

! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
  else

   do i = 1, n
    temp(i) = x(i)
   end do

  endif

! Destroy the target of x
  deallocate(x)

! Re-assign the pointer
   x => temp

  return

  end subroutine my_alloc_p2_ptr


!********************************************************************************
!* This subroutine is useful to expand or shrink real arrays.
!*
!*  Array, x, will be allocated if the requested dimension is 1 (i.e., n=1)
!*  Array, x, will be expanded to the requested dimension, n, if (n > dim(x)).
!*  Array, x, will be shrinked to the requested dimension, n, if (n < dim(x)).
!*
!********************************************************************************
  subroutine my_alloc_p2_matrix_ptr(x,n,m)

  use edu2d_constants   , only : p2

  implicit none
  integer, intent(in) :: n, m
  real(p2), dimension(:,:), pointer :: x
  real(p2), dimension(:,:), pointer :: temp
  integer :: i, j

  if (n <= 0) then
   write(101,*) "my_alloc_int_ptr received non-positive dimension. Stop."
   stop
  endif

! If not allocated, allocate and return
  if (.not.(associated(x))) then
   allocate(x(n,m))
   return
  endif

! If reallocation, create a pointer with a target of new dimension.
  allocate(temp(n,m))
   temp = 0.0_p2

  do i = 1, min(n, size(x,1))
   do j = 1, min(m, size(x,2))
    temp(i,j) = x(i,j)
   end do
  end do

! Destroy the target of x
  deallocate(x)

! Re-assign the pointer
   x => temp

  return

  end subroutine my_alloc_p2_matrix_ptr

 end module edu2d_my_allocation
!********************************************************************************