
      program main
      use prm_flow
      implicit none
!     2D Meshless Linear Advection w/ Periodic BC [ X&Y ]
!     Cartesian point distribution w/ perturbation option
!     MUSCL reconstruction
!     Advection velocity can be set in AdvecVel.f by introducing a new veltype.
!     Initial condition can be set in the function InitCond in file SetInitCond.f90 by introducing a new ictype.
      
      integer :: ierr,i
      double precision :: lst1,lst2,drho,dr
      
      open(77,file='output/summary.dat',form='formatted')
      open(78,file='L2error.dat',form='formatted')
      
!     limiting array parameter
      nxmax = 1000 ; nymax = 1000 ; npmax = nxmax*nymax
      ntmax = 2*npmax ; nnbrmax = 20 

      call alloc    ! allocate arrays

!     Grid Settings
      xmin = 0.0d0 ! Global Xmin
      xmax = 1.0d0 ! Global Xmax
      ymin = 0.0d0 ! Global Ymin
      ymax = 1.0d0 ! Global Ymax
      nx   = 100   ! Points in X Direction
      ny   = 100   ! Points in Y Direction
!     Grid Perturbation Option 
      perturb = 0
      pertx   = 0.5d0
      perty   = 0.2d0

!     Solution Settings
      order= 3          ! 1=first order, 2=higher, 3=higher+limiter
      kkk = 1.0d0/3.0d0 ! MUSCL HO Constant
      ictype = 1        ! Initial Condition [ SetInitCond.f ]
      veltype = 1       ! Advection Velocity [ advecvel.f ]
      dt = 0.001d0      ! Timestep
      Ttime = 1.1d0       ! Final Time

      makeanim = 1
      animinterval = 2000

!     SOLUTION ###################################################
      call GenPtsConn    ! Generate grid and connectivity
      call cpu_time(lst1)
!      call LeastSquares  ! Compute First Order Derivatives Using Least Squares
      call lsTLS
!      call lsPLS
!      call lsCWTLS
      call cpu_time(lst2)
      
      write(77,*) 'lstime  : ',lst2-lst1
      write(77,*) 'nnbrmax : ',maxval(nnbr)
      
      ierr=0;allocate( varini(np0),stat=ierr )
      if( ierr /= 0 ) pause 'allocation error varini'
      varini(:) = 0.0d0
      
      call SetInitCond   ! Initialise Solution 
      call output_vtk(0)
!      call output_vtkdiff(0)

      call InitTimeStep  ! Global Timestep
      call Solve         ! Solution Subroutine

      write(*,*) 'L2 Error ',vdrho
      
      write(77,*) ' '
      write(77,*) 'xmin    : ',xmin
      write(77,*) 'xmax    : ',xmax
      write(77,*) 'ymin    : ',ymin
      write(77,*) 'ymax    : ',ymax
      write(77,*) 'nx      : ',nx
      write(77,*) 'ny      : ',ny
      write(77,*) 'perturb : ',perturb
      write(77,*) 'pertx   : ',pertx
      write(77,*) 'perty   : ',perty
      write(77,*) 'order   : ',order
      write(77,*) 'dt      : ',dt
      write(77,*) 'time    : ',ttime
      close(77)

      stop
      end program main
