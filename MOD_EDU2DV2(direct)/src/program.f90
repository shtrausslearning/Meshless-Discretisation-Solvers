
    program EDU2D
    use edu2d_constants    , only : p2, zero
    use prmflow
    use edu2d_grid_data    , only : construct_grid_data, check_grid_data
    use edu2d_my_main_data , only : nnodes,node
    use modLSQ, only : compute_lsq_coeff_nc, check_lsq_coeff_nc
    implicit none

    integer       :: i,ierr, j
!   ###################################################################### 

!   PRECALCULATION #####################################################
    call inputs  ; write(*,*) '...inputs read'
    call initconstants  ; write(*,*) '...initial constants set'
    
    call preRgrid               ! preread grid to know roughly how many node() to alloc
    call read_grid            ; write(*,*) '...grid read'
    call construct_grid_data  ; write(*,*) '...grid constructed'
    call check_grid_data      ; write(*,*) '...grid checked' 
    call ptypePHYS ; write(*,*) '... allocated ptype for (physical) domain' 
!   ####################################################################
    call solveralloc ; write(*,*) '... allocated solver arrays'

!   GENERATE ADDITIONAL GRID -------------------------------------------
    call gcidRS  ! set boundary ghost x,y,ptype,gcid data
    call gcnxy   ! set gcid nx,ny
    
!   AFTER GRID ORIENTATION ( PHYSICAL / GHOST NODES ) ARE COMPLETELY SET

!   Meshless Coefficient Generation    
!    call compute_lsq_coeff_nc ; write(*,*) '...lq coefficients computed' 
!    call lsnWTLS ! node based WTLS (general)
!    call lsnPLS ! node based PLS ( needs weighting adjustment )
    call eOWTLS  ! edge based WTLS
!   --------------------------------------------------------------------
!   Meshless Coefficient Sort Only
    call lsqedge  ; write(*,*) '... grouped local edge LSQ coefficients'
    
!   initial CV and DV of all points in the domain
    call initsolution
    
    iter=0
    do 
        iter = iter + 1
        
        call gsolver
        call convergence
        
        if( mod(iter,outstep) == 0 )then
         call write_tecplot_file
         call plotsurfaces
        endif
        
        if( iter >= maxiter .or. Rdrho <= tolerance ) exit ! exit loop
        
    enddo
    
!   POST CALCULATION
    call write_tecplot_file ; write(*,*) '...fieldview file written'
    call plotsurfaces

    close(101) ! main summary file
    write(*,*) '... program finished.'
!   STOP Program
    STOP
    
    end program EDU2D
    
    subroutine del_spaces(s)
    character (*), intent (inout) :: s
    character (len=len(s)) tmp
    integer i,j
    j=1
    do i=1,len(s)
    if( s(i:i)==' ') cycle
    tmp(j:j) = s(i:i)
    j=j+1
    enddo
    s = tmp(1:j-1)
    end subroutine del_spaces