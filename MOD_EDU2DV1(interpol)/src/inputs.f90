
    subroutine inputs
    use edu2d_constants    , only : p2, zero
    use prmflow
    use edu2d_grid_data    , only : construct_grid_data, check_grid_data
    use edu2d_my_main_data , only : nq,gradient_type
    
    implicit none
    
    
    write(*,*) 'edu2d solver started...'
    write(*,*) '#######################################'

    open(1,file='input.dat',form='formatted')
    read(1,*) casename
    read(1,*) gridInF
    read(1,*) bcInF
    read(1,*) datafile_tec
    read(1,*) M_inf
    read(1,*) alphIn
    read(1,*) gamma
    read(1,*) CFLexp
    read(1,*) maxiter
    read(1,*) outstep
    read(1,*) tolerance
    read(1,*) nq
    read(1,*) inviscid_flux
    read(1,*) gradient_type
    read(1,*) limiter_type
    read(1,*) gradient_weight
    read(1,*) gradient_weight_p
    read(1,*) nitirs
    read(1,*) epsirs
    close(1)

    open(101,file='summary.dat',form='formatted') ! opem main console summary file
    write(101,*)
    write(101,*) " Summary of Input Parameteres:"
    write(101,*)
    write(101,*) "                  M_inf = ", M_inf
    write(101,*) "                  gamma = ", gamma
    write(101,*) "                 CFLexp = ", CFLexp !CFL for explicit method
    write(101,*) "              tolerance = ", tolerance
    write(101,*) "         max_iterations = ", maxiter
    write(101,*) "          inviscid_flux = ", trim(inviscid_flux)
    write(101,*) "          gradient_type = ", trim(gradient_type)
    write(101,*) "        gradient_weight = ", trim(gradient_weight)
    write(101,*) "      gradient_weight_p = ", gradient_weight_p
    write(101,*) ''
    
    return
    end subroutine inputs