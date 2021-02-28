
 module edu2d_constants

  implicit none

  private

  public :: p2
  public :: zero, one, two, three, four, five, six, seven, eight, nine
  public :: ten, eleven
  public :: half, third, fourth, fifth, sixth, two_third, four_third
  public :: three_fourth, twelfth, pi, one_twentyfourth

  integer , parameter :: sp = kind(1.0)
  integer , parameter :: p2 = selected_real_kind(2*precision(1.0_sp))

  real(p2), parameter :: zero = 0.0_p2, &
                          one = 1.0_p2, &
                          two = 2.0_p2, &
                        three = 3.0_p2, &
                         four = 4.0_p2, &
                         five = 5.0_p2, &
                          six = 6.0_p2, &
                        seven = 7.0_p2, &
                        eight = 8.0_p2, &
                         nine = 9.0_p2, &
                          ten = 10.0_p2, &
                       eleven = 11.0_p2, &
                         half = 0.5_p2, &
                        third = 1.0_p2/ 3.0_p2, &
                       fourth = 1.0_p2/ 4.0_p2, &
                        fifth = 1.0_p2/ 5.0_p2, &
                        sixth = 1.0_p2/ 6.0_p2, &
                    two_third = 2.0_p2/ 3.0_p2, &
                   four_third = 4.0_p2/ 3.0_p2, &
                 three_fourth = 3.0_p2/ 4.0_p2, &
                      twelfth = 1.0_p2/12.0_p2, &
             one_twentyfourth = 1.0_p2/24.0_p2

  real(p2), parameter :: pi = 3.141592653589793238_p2

 end module edu2d_constants