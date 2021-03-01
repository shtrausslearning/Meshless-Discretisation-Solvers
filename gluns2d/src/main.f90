
  program Unstruct2D
  use ModDataTypes
  use prmflow
  use ModInterfaces
  implicit none

! local variables
  integer           :: ierr
  integer, allocatable     :: niedge(:), & ! temporary edge structures (DO NOT overwrite
                              iedge(:,:)   ! between EdgesInitialize and InitMetrics)
    real(rtype) :: ltt2
! #######################################################################################
  
!   Ntype
!   1: Physical Wall
!   2: Physical Farfield Wall
!   3. Physical Fluid
!   4. Dummy Wall Ghost Cell
!   5. Dummy Farfield Ghost Cell
    
!   Ptype
!   1. All Physical Domain Nodes
!   2. All Ghost Cell Domain Nodes
  
    !call Getarg( 1,fname )
    !IF (Len_trim(fname) == 0) call Usage
    
    call readparams
    call InitConstants       ! initialize some constants
    call PrintParams         ! print input parameters for checking

    nconv = 4 ! set no. of equations (rho, rho*u, rho*v, rho*E, ...);
    ndepv = 5 ! set no. of dependent variables (p, T, c, gamma, cpgas, ...)
      !ndepv = ndepv + 2   ! laminar viscosity, heat conduction coeff.

    call ReadGrid         ! read grid dimensions, coordinates, and triangle nodes

!   generate edge list
    allocate( niedge(phynod),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate temporary edge pointer" )
    allocate( iedge(3,2*nt),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate temporary edge list" )

    call EdgesInitialize( niedge,iedge )
    call EdgesFinalize( niedge,iedge )  

    call createConnections ! gridless structure
    
    call InitMetrics( niedge,iedge )
    call InitMetricsBound( niedge,iedge )
    deallocate( iedge  ) ; deallocate( niedge )

!   Make Farfield Ghost Cells
    call make_ffgc          ! update x,y ghost cell location of farfield ghost cells
    call gcidarray            ! make gcid() array (1) fluid node (2) fluid node related ghost cell  
    call make_wgc          ! make wall ghost cells
    call alloc                    ! allocation of global calculation arrays
    cv(:,:) = -777.0d0;dv(:,:) = -777.0d0 ! reset

!   First Step Solution
    !call ReadSolution
    call InitialSolution   ! 1 to allnod initial solution 

!   pre calculation
    call ptypealloc         ! main masking 
 !   call LeastSquares       ! aij,bij gridless coefficients
    call lsWTLS
!    call lsWPLS
   ! call lsCWTLS
    
    !call confirms
    !pause 'confirms'
    
    write(*,*) aij(2,100),bij(2,100)
    pause

    write(fnConv,"(A,A,A)") "convall_",Trim(title),".dat"
    open(50, file=fnConv, form='formatted')
    write(fnConv,"(A,A,A)") "convrho_",Trim(title),".dat"
    open(51, file=fnConv, form='formatted')
    write(fnConv,"(A,A,A)") "convcl_",Trim(title),".dat"
    open(52,file=fnConv,form='formatted')
    write(fnConv,"(A,A,A)") "convcd_",Trim(title),".dat"
    open(53,file=fnConv,form='formatted')
    
    call LimiterRefvals ! reference values for limiter

!   Calculation
    write(*,*)
    write(*,*) 'Starting RK Time Integration...'
    write(*,*) ''
    call solve
    
!   Post Calculation
    call output_vtkv3         ! write final vtk
    call PlotSurfaces       ! write final surface points
    call WriteSolution      ! write input file
    close(50);close(51);close(52);close(53)
    
    call CPU_TIME(ltt2) ! outstep total time
    write(99,1006) iter,ltt2  ! write final timestamp
    close(99)

    write(*,*) 'Ending Program'
    
1006  format(I0,1P,2X,F0.3)
    
contains

    subroutine Usage
    write(*,"(/,A,/)") "Usage:"
    write(*,"(A,/)")   "Unstruct2D <input file>"
    stop
    end subroutine Usage
    
    end program Unstruct2D
    
    Subroutine Confirms
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none
    
    integer :: i,j,inn,nn,id
    
    call output_vtk3
    call output_vtk4
    call output_vtk5
    
    pause 'vtkoutputs'
    
    !open(33,file='xyall.dat',form='formatted')
    !do i=1,allnod+igcs
    !write(33,*) x(i),y(i)
    !enddo
    !close(33)
    
    !open(34,file='xyphynod.dat',form='formatted')
    !do i=1,phynod
    !write(34,*) x(i),y(i)
    !enddo
    !close(34)
    
    !open(35,file='xyallnod.dat',form='formatted')
    !do i=1,allnod
    !write(35,*) x(i),y(i)
    !enddo
    !close(35)
    
    !open(36,file='idneighbours.dat',form='formatted')
    !id=260
    !do j=1,nbers(id)
    !inn = conn(j,id)
    !write(36,*) x(inn),y(inn), inn
    !enddo
    !close(36)
    
    pause 'idneighbours'
    
    
    do i=1,phynod
    if( ntype(i) .eq. 1 )then  ! if wall node
    nn=0    
    do j=1,nbers(i)            ! do for all neighbours
    inn = conn(j,i)            ! corresponding neighbour id
    if( ptype(inn) .eq. 2 )then  
    nn=nn+1
    endif
    enddo
    
    if( nn .eq. 1 )then 
    write(*,*) 'i ',i
    write(*,*) x(i),y(i)
    pause
    endif
    
    if( nn .eq. 1 )then
    open(88,file='neighbours.dat',form='formatted')
    do j=1,nbers(i)
    inn=conn(j,i) 
    write(88,*) x(inn),y(inn),inn,ntype(inn)
    enddo
    close(88) 
    pause '88'
    endif
    
    write(*,*) i,nn
    endif
    enddo
    
    return
    end subroutine confirms
    
  subroutine EM( message )
  implicit none

! parameters
  character(*), intent(in) :: message

  write(*,"(A,A,/)") " Error: ",Trim( message )
  pause 'ddd'
  stop

  end subroutine EM
    
    
    