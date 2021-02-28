
    subroutine initconstants
    use edu2d_constants    , only : p2, zero, pi
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    use prmflow
    implicit none
    
    real(p2) :: gam1,rgas
    
    gamma = 1.4d0
    cpgas = 1004.5d0
    pinf = 1.0d5
    tinf = 288.0d0
    

    rad = 180.0d0/pi
    
    nconv = 4
    ndepv = 5
    gam1 = gamma - 1.0d0
    rgas = gam1*cpgas/gamma
    
    alpha = alpha/rad
    rhoinf = pinf/(rgas*tinf)
    qinf = M_inf*sqrt(gam1*cpgas*tinf)
    uinf = qinf*cos(alpha)
    vinf = qinf*sin(alpha)
    
    xref = 0.25d0
    yref = 0.0d0 
    cref = 1.0d0 
    
    return
    end subroutine initconstants