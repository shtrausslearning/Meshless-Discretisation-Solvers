
!!   INITIAL SOLUTION
!    subroutine set_initial_solution
!
!    use edu2d_constants   , only : zero, one, pi, p2
!    use prmflow
!    use edu2d_my_main_data, only : nnodes, node, rho_inf, u_inf, v_inf, p_inf
!    use convertUW
!
!    implicit none
!
!    !Local variables
!    integer  :: i
!    real(p2) :: rad
!
!    write(101,*)
!    write(101,*) "Setting up initial solution (uniform stream)..."
!    write(101,*)
!
!    ! Uniform stream values (nondimensionalized variables)
!    rad = 180.0d0/pi
!    alphIn = alphIn/rad
!    rho_inf = one
!    u_inf = M_inf*cos(alphIn)
!    v_inf = M_inf*sin(alphIn)
!    p_inf = one/gamma
!
!    ! Specify the uniform stream values at all nodes.
!    NODES: do i=1,nnodes+nGC+nIP
!    node(i)%w = (/ rho_inf, u_inf, v_inf, p_inf /) ! Primitive variables
!    node(i)%u = w2u(node(i)%w)                     ! Conservative variables
!    enddo NODES
!
!    end subroutine set_initial_solution