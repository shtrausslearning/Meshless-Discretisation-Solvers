
    subroutine initsolution
    use edu2d_constants    , only : p2, zero
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    use prmflow
    implicit none
    
    integer :: i
    real(p2) :: gam1,rgas,g1cp,rhoq
    
    gam1 = gamma - 1.0d0
    rgas = gam1*cpgas/gamma
    g1cp = gam1*cpgas
    
!   set initial solution for all physical nodes
    do i=1,nnodes+nGC+nIP
    node(i)%cv(1) = rhoinf
    node(i)%cv(2) = rhoinf*uinf
    node(i)%cv(3) = rhoinf*vinf
    node(i)%cv(4) = pinf/gam1 + 0.5d0*rhoinf*qinf*qinf
    enddo
    
!   set dependent variables for all physical nodes
    do i=1,nnodes+ngc+nIP
    rhoq = node(i)%cv(2)*node(i)%cv(2) + node(i)%cv(3)*node(i)%cv(3)
    node(i)%dv(1) = gam1*(node(i)%cv(4)-0.5d0*rhoq/node(i)%cv(1))
    node(i)%dv(2) = node(i)%dv(1)/(rgas*node(i)%cv(1))
    node(i)%dv(3) = sqrt(g1cp*node(i)%dv(2))
    node(i)%dv(4) = gamma
    node(i)%dv(5) = cpgas
    enddo
    
    return
    end subroutine initsolution
    

    