!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:OOverSetS2U
!	Description:
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.11.24
!	Update:2017.11.24
!	Other:
!***********************************/
subroutine OGetDisplacement(MG,iStep)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod
    implicit none

    integer, intent(in) :: iStep
    type(MoveGrid), intent(inout) :: MG

    call SampleRotation

    if(iStep /= 1) then
        call CoordinateTransformOfDisplacement
    end if

return
contains
!ここでの入力変位量は局所座標における変位量になっているため， 逆回転行列を掛けてグローバル座標系での表記に戻してやる必要がある

    subroutine SampleRotation
    implicit none
    !n段階からn+1段階への変位をグローバル座標系における表現で受け取る
        MG%GCD%Translation(1) = 0.001d0*cos(2.0d0*dPi*dble(iStep)/100.0d0)
        MG%GCD%Translation(2) = 0.001d0*sin(2.0d0*dPi*dble(iStep)/100.0d0)
        MG%GCD%Translation(3) = 0.0d0

        MG%GCD%Rotation(1) =  0.0d0
        MG%GCD%Rotation(2) =  0.0d0
        MG%GCD%Rotation(3) =  0.01d0*dPi/50.0d0 !1 degree?

    return
    end subroutine SampleRotation


    subroutine CoordinateTransformOfDisplacement
    implicit none

        MG%GCD%Translation(:) = matmul(MG%TM%Top2Current(1:3,1:3),MG%GCD%Translation)

    return
    end subroutine CoordinateTransformOfDisplacement

end subroutine OGetDisplacement
