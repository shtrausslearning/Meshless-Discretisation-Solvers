
!   subroutine to calculate timestep for each physical node in the domain
!   required for time integration

    subroutine gtimestep
    use edu2d_constants, only : p2
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    use prmflow
    implicit none
    
    integer :: i,j
    real(p2) :: term,a,ds,u,v,p,sx,sy,ter1,ter2
    
    do i=1,nnodes ! all physical nodes only
    
    term = 0.0d0
    do j=1,node(i)%nnghbrs    ! for all neighbours
    
    sx = node(i)%aij(j);sy = node(i)%bij(j)
    ds = sqrt(sx**2+sy**2)
    v = node(i)%cv(2)/node(i)%cv(1)
    v = node(i)%cv(3)/node(i)%cv(1)
    ter1 = abs(sx*u+sy*v)
    ter2 = node(i)%dv(3)*sqrt(sx**2+sy**2)
    term = term + (ter1+ter2)

    enddo

    node(i)%dt = CFLexp/term
    
    enddo
    
    return
    end subroutine gtimestep