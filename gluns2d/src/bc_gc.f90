
    Subroutine bc_wgc
!   ################################################################
    use ModDataTypes
    use prmflow
    use ModInterfaces
    implicit none

    integer :: ic,i, id1, id2, nn, ierr, inn, ii, iii,nnn, j, p1, p2, ps,ie,ie2,id,k, ps2
    integer, allocatable :: wnod(:)
    real(rtype) :: gam1,vnorm
    real(rtype) :: xt(5), yt(5) ! temporary x,y coordinates
    real(rtype) :: dxf,dyf,dsf,unx,uny,nx,ny,ui,vi,ut,un,ug,vg
    real(rtype) :: dist,x0,y0,x1,y1
    real(rtype) :: rgas,g1cp,rhoq    
    real(rtype) :: HSTFS,p_gc,rho_gc
!   ################################################################
    
1006  FORMAT("# vtk DataFile Version 2.0" )
1007  FORMAT("VTK FORMAT")
1008  FORMAT("ASCII")
1009  FORMAT("DATASET UNSTRUCTURED_GRID")
1004  FORMAT('POINTS',I7,' float ')
1011  FORMAT(1E14.6,1P,1E14.6,1P,1E14.6)
!
    !open(55,file='code2.vtk',form='formatted')
    !write(55,1006)
    !write(55,1007)
    !write(55,1008) 
    !write(55,*)""
    !write(55,1009)
    !write(55,1004) igcs
!    
!    do i=1,igcs
!    id1 = gcid(1,i) ! fluid node
!    id2 = gcid(2,i) ! dummy node
!    write(55,1011) x(id2),y(id2),0.0d0
!    enddo
!    close(55)
!    pause '55'

!    open(55,file='code2.vtk',form='formatted')
!    write(55,1006)
!    write(55,1007)
!    write(55,1008) 
!    write(55,*)""
!    write(55,1009)
!    
!    id=160;k=0
!    do j=1,nbers(id)
!    inn = conn(j,id)
!    k=k+1
!    enddo
!
!    write(55,1004) k
!!    
!    do j=1,nbers(id)
!    inn = conn(j,id)
!    write(55,1011) x(inn),y(inn),0.0d0
!    enddo
!    close(55)
!    pause '55'
    
!   Wall Ghost Cell
    do i=1,igcs  ! do for all ghost cell
   
    id1 = gcid(1,i) ! fluid node
    id2 = gcid(2,i) ! dummy node
    
    nn=0
!   count how many wall neighbours are connected to the ghost cell
    do j=1,phynod
    if( ntype(j) .eq. 1 )then  ! if wall node
     do k=1,nbers(j)         ! for all neighbours of wall node
     inn = conn(k,j)         ! neighbour of wall node
     if( inn .eq. id2 )then     ! if wall neighbour is the dummy node
     nn=nn+1                    ! how many ghost cells a wall neighbour has
     endif
     enddo
    endif
    enddo
    
    ierr=0;allocate( wnod(nn),stat=ierr )
    if (ierr /= 0) call EM( "cannot allocate memory for local wnod()" )
   
    nn=0
    do ii=1,phynod
    if( ntype(ii) .eq. 1 )then ! if wall node
     do iii=1,nbers(ii)        ! go through all neighbours
     inn = conn(iii,ii)        ! neighbour global number
     if( inn .eq. id2 )then    ! if neighbour is the ghost cell
     nn=nn+1
     wnod(nn) = ii             ! ghost cell's wall neighbour
     endif
     enddo
    endif
    enddo
    
!   select wall node for ghost cell
    if( nn .eq. 1 )then
    ps = wnod(1) 
    nx = nx_wc(ps)
    ny = ny_wc(ps)
    elseif( nn .eq. 2 )then
    ps = wnod(2) 
    ps2 = wnod(1)
    nx = (nx_wc(ps) + nx_wc(ps2))/2.0d0
    ny = (ny_wc(ps) + ny_wc(ps2))/2.0d0
!    nx = nx_wc(ps)
!    ny = ny_wc(ps)
    elseif( nn .eq. 3 )then
    ps = wnod(2) 
    ps2 = wnod(1)
    nx = (nx_wc(ps) + nx_wc(ps2))/2.0d0
    ny = (ny_wc(ps) + ny_wc(ps2))/2.0d0
!    nx = nx_wc(ps)
!    ny = ny_wc(ps)
    elseif( nn .gt. 3 )then
     write(*,*) 'nn>3'
!    ps = wnod(2) 
!    ps2 = wnod(3)
!    nx = (nx_wc(ps) + nx_wc(ps2))/2.0d0
!    ny = (ny_wc(ps) + ny_wc(ps2))/2.0d0
    endif
    
    ui= cv(2,id1)/cv(1,id1) ! u_ip
    vi= cv(3,id1)/cv(1,id1) ! v_ip  
!    if( ui .eq. -777.0d0 ) write(*,*) 'ui = -777'
!    if( nx .eq. -777.0d0 ) write(*,*) 'nx = -777'
    
    vnorm = ui*nx + vi*ny
    ug = ui - 2.d0*vnorm*nx
    vg = vi - 2.d0*vnorm*ny

!   CONSERVATIVE VARIABLES
    
!    cv(1,id2) = 2.0d0*cv(1,id1) - cv(1,ps) ! pushes shock back
!    cv(1,id2) = cv(1,id1)  ! more closesly resembles fvm
    cv(1,id2) = 2.0d0*cv(1,ps) - cv(1,id1) ! Wghost = Wsurf - Wfluidip
!    cv(1,id2) = abs(cv(1,id2))

    cv(2,id2) = ug*cv(1,id2)
    cv(3,id2) = vg*cv(1,id2)
    cv(4,id2) = 2.0d0*cv(4,id1) - cv(4,ps)
!    cv(4,id2) = 2.0d0*cv(4,ps) - cv(4,id1) ! completely reduces cp
!    cv(4,id2) = abs(cv(4,id2))
!    cv(4,id2) = cv(4,id1)  ! 

    gam1 = gamma - 1.D0
    rgas = gam1*cpgas/gamma
    g1cp = gam1*cpgas
    rhoq    = cv(2,id2)*cv(2,id2) + cv(3,id2)*cv(3,id2)
    
!   DEPENDABLE VARIABLES

!    dv(1,id2) = gam1*(cv(4,id2)-0.5D0*rhoq/cv(1,id2))
!    dv(1,id2) = dv(1,id1)
    dv(1,id2) = 2.0d0*dv(1,id1) - dv(1,ps)
!    dv(1,id2) = 2.0d0*dv(1,ps) - dv(1,id1)

    dv(2,id2) = dv(1,id2)/(rgas*cv(1,id2))
    dv(3,id2) = Sqrt(g1cp*dv(2,id2))
    dv(4,id2) = gamma
    dv(5,id2) = cpgas
    
    deallocate( wnod )
    enddo

    return
    end subroutine bc_wgc
    
    
