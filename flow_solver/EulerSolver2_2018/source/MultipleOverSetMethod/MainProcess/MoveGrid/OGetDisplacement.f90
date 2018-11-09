!***********************************/
!	Name:変位量の入力
!	Alias:OGetDisplacementMK2
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.10
!	Update:
!	Other:
!***********************************/
subroutine OGetDisplacement(MG,MW)
    use StructVar_Mod
    use ConstantVar_Mod
    implicit none

    type(MoveGrid), intent(inout) :: MG
    type(MotionWait), intent(in) :: MW

    MG%GCD%Translation = MW%Accumulate%Translation

    MG%GCD%Rotation = MW%Accumulate%Rotation
    !ここでの入力変位量は局所座標における変位量になっているため， 逆回転行列を掛けてグローバル座標系での表記に戻してやる必要がある
    MG%GCD%Translation(:) = matmul(MG%TM%Top2Next(1:3,1:3),MG%GCD%Translation)

return
end subroutine OGetDisplacement
