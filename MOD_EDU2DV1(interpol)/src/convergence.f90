
!   subroutine to calculate the L2N convergence
    
    subroutine convergence
    use edu2d_constants 
    use edu2d_my_main_data, only : nnodes, node
    use prmflow
    implicit none
    
    real(p2) :: dr,drho1
    integer :: i
    logical :: exist
    character(256) :: linebuf
    
    Rdrho = 0.0d0
    do i=1,nnodes
    dr = node(i)%cv(1) - node(i)%cvold(1)
    Rdrho = Rdrho + dr*dr
    enddo
    
    if( iter == 1 )then
    drho1 = sqrt(Rdrho) + 1.0d-32
    Rdrho=1.0d0
    else
    Rdrho = sqrt(Rdrho)
    endif
    
    call forces
    
    inquire(file='convergence.dat',exist=exist)
    if(exist)then
    open(112,file='convergence.dat',status='old',position='append',action='write')
    else
    open(112,file='convergence.dat',status='new',action='write')
    endif
    
    write(*,100) iter,log10(Rdrho),clp,cdp
    if( ISNAN(Rdrho) ) pause 'rdrho NaN' 
    write(linebuf,*) iter,',',log10(Rdrho)
    call del_spaces(linebuf)
    write(112,*) trim(linebuf)
    close(112)
    
100 format(I5,1X,F7.4,1X,F7.4,1X,F7.4)
    
    return
    end subroutine convergence