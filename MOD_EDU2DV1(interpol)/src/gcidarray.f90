
!   subroutine required to generate gcid(1:2,:) relation required for boundary
!   treatment ( cycle through all nodes required to set boundary value )

    subroutine gcidARRAY
    use edu2d_constants
    use prmflow
    use edu2d_my_main_data , only : node,bound,nbound,nnodes
    implicit none
    
    integer :: i,j,ii,ierr
    integer :: nGCid,nIPid
    character(256) :: linebuf
    
!   nIP and nGC must be the same to use gcid
    if( nIP .ne. nGC ) pause 'gcid array cant be generated'
    
!   if ngc and nIP are the same proceed
    igcs = nIP  ! specific to gcid array max ; its only valid if nIP = nGC
    
    ierr=0;allocate( gcid(2,igcs),stat=ierr )
    if( ierr /= 0 ) pause 'allocation error gcid'
    gcid=0 ! reset
    
!   now lets fill the array; the order of nGC generation via the boundary nodes
!   is identical to that for nIP, therefore their values will coincide in order
!   lets loop again over boundary nodes

!   what are the starts of both GC and IP
    nGCid = nnodes
    nIPid = nnodes+ngc

!   now fill in the gcid() array
    ii=0  ! need counter for storage
    do i=1,nbound                 ! cycle through all boundary segments
    do j=1,bound(i)%nbnodes-1     ! boundary node # exlude last as always
    
    ii=ii+1        ! counter for gcid storage
    nGCid=nGCid+1    ! counter for gc nodes
    nIPid=nIPid+1    ! counter for IP nodes
    
    gcid(1,ii) = nIPid  ! fluid node ( IP point )
    gcid(2,ii) = nGCid  ! ghost cell node 
    
    enddo
    enddo
    
!   [check] array is correctly filled in
    open(44,file='./info/gcid.dat',form='formatted')
        
    do i=1,igcs
    write(44,*) gcid(1,i),node(gcid(1,i))%ptype,gcid(2,i),node(gcid(2,i))%ptype
    enddo
    
    close(44)
    
!   [check]

    open(44,file='./info/gcidxyALL.dat',form='formatted')
    do j=1,2
    do i=1,igcs
    write(linebuf,*) node(gcid(j,i))%x,',',node(gcid(j,i))%y
    call del_spaces(linebuf)
    write(44,*) trim(linebuf)
    enddo
    enddo
    close(44)
!    
!    pause 'gcid'


    return
    end subroutine gcidARRAY

