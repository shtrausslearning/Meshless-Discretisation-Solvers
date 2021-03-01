
    Subroutine make_ffgc
    use ModDataTypes
    use prmflow
    use ModInterfaces, only : EM
    implicit none
    
    integer :: i,ie,jj,p1,p2,bb,ic,bb2,itemp,j,ierr
    real(rtype) :: x0,y0,dxx,dyy,mid,cid,x_gc,y_gc
    integer, allocatable :: fnp(:) ! loop fluid node allocation
    
    
    ! adjust dummy node location of farfield bc
    do i=1,phynod
    if( ntype(i) .eq. 2 )then ! if ff boundary node
  
    !  find dummy node of ff boundary node
    ! each boundary node has a ghost dummy pair (from FVM formulation)
       do ie=nedint+1,nedges
       if( edge(1,ie) .eq. i )then
       itemp = edge(2,ie)  ! used at very end to adjust x,y position
       endif
       enddo
  
    !  find the edge face nodes associated with farfield bc node
       do jj=1,nbfaces
       if( bface(1,jj) .eq. i )then  
       p1 = i                                  ! corresponding wall neighbour
       p2 = bface(2,jj)                   ! bface edge partner
       endif
       enddo
  
       !write(*,*) (conn(j,i),j=1,nbers(i))
       !pause 'neighbours'
       !write(*,*) ( ntype(conn(j,i)),j=1,nbers(i)) 
       !pause 'point type of neighbours' 
       
       bb=0
    !  count how many fluid neighbours 
       do j=1,nbers(i) ! i is the ff boundary node
       ic = conn(j,i)
       if( ntype(ic) .eq. 3 )then ! if neighbour is fluid
       bb=bb+1
       endif
       enddo
       
       !write(*,*) 'number of fluid neighbours',bb

!      allocate/deallocate within loop
       allocate( fnp(bb),stat=ierr )
       if(ierr /= 0) call EM( "cannot allocate memory for local fnp()" )
  
       bb2=0
    !  store neighbour node number of fluid
       do j=1,nbers(i)
       ic = conn(j,i)
       if( ntype(ic) .eq. 3 )then
       bb2=bb2+1
       fnp(bb2) = conn(j,i)
       endif
       enddo
  
       !write(*,*) (fnp(j),j=1,bb)
       !pause 'fluid neighbour node values'
   
       !do j=1,bb
       !write(*,*) x( fnp(j) ),y( fnp(j) )
       !enddo
       !pause 'coordinates of these fluid neighbours'
   
       x0=0.0d0;y0=0.0d0
    !  find average fluid node location
       do jj=1,bb
       x0 = x0 + x( fnp(jj) )
       y0 = y0 + y( fnp(jj) )
       enddo
       x0 = x0/bb ; y0 = y0/bb
   
       dxx = x(p2) - x(p1)
       dyy = y(p2) - y(p1)
       if( dxx == 0.0d0 )then
       dxx = 1d-20
       endif
  
       mid = dyy/dxx
       cid = y(p1) - mid*x(p1)
   
      !write(*,*) itemp
     !pause 'itemp'
   
     !write(*,*) x(itemp),y(itemp)
     !pause 'eeee'
   
!    Reflection Point X,Y Coordinates 
     x_gc = ((1.0d0-mid**2)*x0 + 2.0d0*mid*(y0-cid))/(1.0d0+mid**2)
     y_gc = ( 2.0d0*mid*x0 - (1.0d0-mid**2)*y0 + 2.0d0*cid )/(1.0d0+mid**2)
   
     !write(*,*) x_gc,y_gc
     !pause 'x_gc,y_gc'
   
!    storage of gc value
     x(itemp) = x_gc
     y(itemp) = y_gc
  
     deallocate( fnp )
  
    endif
    enddo    
    
    !open(434,file='random.dat',form='formatted')
    !do j=1,nbers(126)
    !write(434,*) x( conn(j,126) ), y( conn(j,126) )
    !enddo
    !close(434)
  
    !open(100,file='777.dat',form='formatted')
    !do j=1,allnod
    !write(100,*) x(j),y(j)
    !enddo
    !close(100)
    
!   confirm ntype array fullness
    do i=1,allnod
    if( ntype(i) .eq.  -777 )then
    pause 'ntype array not fully defined'
    endif
    enddo
    
    return
    end subroutine make_ffgc