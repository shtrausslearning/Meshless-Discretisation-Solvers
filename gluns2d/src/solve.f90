
    Subroutine solve
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    
    integer :: irk,i,ierr
    real(rtype) :: lst1,ltt2
    character(chrlen) :: fname2
    
    if(lrest .ne. 1) iter = 0
 
    write(fname2,"(A,A)") Trim(title),"_timestamp.dat"
    open(99, file=fname2, status="unknown", action="write", iostat=ierr)
    if (ierr /= 0) call EM( "cannot open plot file (fname2)" )
    
    call CPU_TIME(lst1) 
    
    do  
    iter = iter + 1
    cvold(:,:) = cv(:,:)
    call TimeStep

    do irk=1,nrk

    call Gradients     ! Least Square Gradient
!    call GradientsG     ! Green Gauss Gradient
!    call LimiterInit    ! calculates min/max around main node
!    call Limiter
!    call residual_roe   ! flux calculation
!    call residual_roer   ! flux calculation
    call residual_roeMU
    
! - implicit residual smoothing
    if (epsirs > 0.D0) call irs
    
!   Conservative Updated ###########################
    do i=1,phynod
    cv(1,i) = cvold(1,i) - ark(irk)*tstep(i)*rhs(1,i)
    cv(2,i) = cvold(2,i) - ark(irk)*tstep(i)*rhs(2,i)
    cv(3,i) = cvold(3,i) - ark(irk)*tstep(i)*rhs(3,i)
    cv(4,i) = cvold(4,i) - ark(irk)*tstep(i)*rhs(4,i)
    enddo
!   Dependable Updated #############################
    call dvall          ! update dv for [1,allnod]
    
!    call BoundaryConditions
    call bc_ff  ! farfield dummy node update
    call bc_wgc         ! Update wall ghost cell cv(:)  
    enddo

    call Convergence

    if( mod(iter,outstep) == 0)then
    call CPU_TIME(ltt2) ! outstep total time
    write(99,1006) iter,ltt2
    endif
 
    if (Mod(iter,outstep) == 0) then
!!    call WriteSolution
    if( vtkout == 1 )then
    call output_vtkv3  ! overwrite type / not iterative
    endif
    call PlotSurfaces
    endif

    if(iter>=maxiter .or. drho<=convtol) exit
    enddo
    
    call output_vtkv3

1006  format(I0,1P,2X,F0.3)
    
    end subroutine solve
    
    
    
    
    
    