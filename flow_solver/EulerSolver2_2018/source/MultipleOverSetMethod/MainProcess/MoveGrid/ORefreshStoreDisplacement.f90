!***********************************/
!	Name:累積させた移動量の初期化
!	Alias:ORefreshStoreDisplacement
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.01.10
!	Update:
!	Other:
!***********************************/
subroutine ORefreshStoreDisplacement(MW)
    use StructVar_Mod
    implicit none

    type(MotionWait), intent(inout) :: MW

    MW%Accumulate%Translation = 0.0d0
    MW%Accumulate%Rotation = 0.0d0

return
end subroutine ORefreshStoreDisplacement
