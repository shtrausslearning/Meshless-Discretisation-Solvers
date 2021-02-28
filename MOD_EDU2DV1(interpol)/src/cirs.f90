
!   subroutine intented to formulate CIRS for convergence acceleration

    subroutine cirs
    use edu2d_constants, only : p2
    use prmflow
    use edu2d_my_main_data , only : node,nnodes,edge,nedges
    use edu2d_my_allocation
    implicit none
    
    integer :: itirs,i,j,ie,ib,ibn,ibegn,iendn,jj,ierr
    real(p2) :: den
    
!   (1) store old RHS
    rhsold = 0.0d0
    do i=1,nnodes
    ncontr(i) = 0
    do j=1,4
    rhsold(j,i) = node(i)%rhs(j)
    enddo
    enddo
    
!   (2) JACOBI TERATIONS for CIRS
    do itirs=1,nitirs
    
    rhsit=0.0d0
    if( itirs == 1 )then
    EDGES1: do ie=1,nedges
    i = edge(ie)%n1
    j = edge(ie)%n2
    ncontr(i) = ncontr(i) + 1
    ncontr(j) = ncontr(j) + 1
    do jj=1,4
    rhsit(jj,i) = rhsit(jj,i) + node(j)%rhs(jj)
    rhsit(jj,j) = rhsit(jj,j) + node(i)%rhs(jj)
    enddo
    enddo EDGES1
    
    else
    EDGES2: do ie=1,nedges
    i = edge(ie)%n1
    j = edge(ie)%n2
    do jj=1,4
    rhsit(jj,i) = rhsit(jj,i) + node(j)%rhs(jj)
    rhsit(jj,j) = rhsit(jj,j) + node(i)%rhs(jj)
    enddo
    enddo EDGES2
    endif
    
!   (3) New residual
    do i=1,nnodes
    den = 1.0d0/(1.0d0+epsirs*real(ncontr(i)))
    do jj=1,4
    node(i)%rhs(jj) = (rhsit(jj,i)*epsirs+rhsold(jj,i))*den
    enddo
    enddo
    
    enddo ! loop over itirs
    
    return
    end subroutine CIRS
    