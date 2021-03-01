   
    subroutine Solve
    use prm_flow
    implicit none

    integer          :: i, nrks, iter
    double precision :: time, varmin, varmax
    double precision :: tc1,tc2,avt,mat,ti1,ti2
    double precision :: drho,dr

    nrks = 3
    time = 0.0d0
    iter = 0
    
    mat = 1d-10
        
    call cpu_time(tc1)

    open(20, file='./output/hist.dat')
    do while(time < Ttime)
        time = time + dt
        iter = iter + 1
        var_old(1:npts) = var(1:npts)
        
        call cpu_time(ti1)
        
        do i=1,nrks
        call Gradients
        call Residual
        call Update(i)
        enddo
        
        call CPU_TIME(ti2)
        
        mat = max(ti2-ti1,mat) ! cycled max dt
        
        varmin = 1.0d20
        varmax =-1.0d20

        do i=1,np0
        varmin = dmin1(varmin, var(i))
        varmax = dmax1(varmax, var(i))
        enddo
 
        write(*,10) iter,time,varmin,varmax
        write(20,10) iter,time,varmin,varmax
        
         vdrho = 0.0d0
         drho  = 0.0d0
         do i=1,npts
         if( ptype(i) .eq. 1 )then
         dr   = varini(i) - var(i)
         drho = drho + dr**2
         endif
         enddo
         vdrho  = Sqrt(drho)
         write(78,*) vdrho

        if(iter.eq. 1000 .or. iter .eq. 2000 .or. iter .eq. 3000 )then
!        if( mod(iter,100) .eq. 0 )then
         call output_vtk(iter)   ! variable
         call output_vtkdiff(iter) ! error accumulation
         
        endif

    enddo
    close(20)
    
    call cpu_time(tc2)
    
    avt = (tc2-tc1)/iter ! Averaged Iteration Cost
    
!   ######################################
    write(77,*) 'avt:     ',avt
    write(77,*) 'mat:     ',mat
    write(77,*) 'Totalt   ',tc2-tc1
    write(77,*) 'Initial  ',varmin0,varmax0
    write(77,*) 'Final    ',varmin,varmax
!   #######################################
    
10      format(i8,3e16.6)
    
    return
    end subroutine Solve