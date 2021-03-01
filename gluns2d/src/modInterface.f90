
    module ModInterfaces

    implicit none
    interface
        
!    subroutine BcondFarfield( ibegn,iendn )
!    use ModDataTypes
!    integer, intent(in) :: ibegn, iendn
!!    end subroutine BcondFarfield
!    
!    subroutine BoundaryConditions( work )
!    use ModDataTypes
!    real(rtype) :: work(:)
!    end subroutine BoundaryConditions

    subroutine EdgesFinalize( niedge,iedge )
    integer :: niedge(:), iedge(:,:)
    end subroutine EdgesFinalize

    subroutine EdgesInitialize( niedge,iedge )
    integer, intent(out) :: niedge(:), iedge(:,:)
    end subroutine EdgesInitialize

    subroutine EM( message )
    character(*), intent(in) :: message
    end subroutine EM

    subroutine InitMetrics( niedge,iedge )
    integer, intent(in) :: niedge(:), iedge(:,:)
    end subroutine InitMetrics

    subroutine InitMetricsBound( marker,btria )
    integer :: marker(:), btria(:,:)
    end subroutine InitMetricsBound

    function ReadChar( iunit )
    integer, intent(in) :: iunit
    character(1) :: ReadChar
    end function ReadChar
    
    end interface

    end module ModInterfaces