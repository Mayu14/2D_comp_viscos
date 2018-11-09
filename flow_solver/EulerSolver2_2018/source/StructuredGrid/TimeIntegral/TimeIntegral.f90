!***********************************/
!	Name:時間積分制御用プログラム
!	Alias:TimeIntegral
!	Description:Configulationによって時間積分の精度を変更できる
!	Type:Conf,Geom,CellEdge,CellCenter,iStep(RRK2用)
!	Input:Geom,CC,CE,CC,iStep
!	Output:CC%ConservedQuantity
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.06
!	Other:局所時間刻み法はまだ未対応
!***********************************/
subroutine TimeIntegral(Conf,Geom,CE,CC,iStep)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod
implicit none
    type(Configulation),intent(in) :: Conf
    type(Geometry),intent(in) :: Geom
    type(CellEdge), intent(in) :: CE
    type(CellCenter), intent(inout) :: CC
    integer, intent(in) :: iStep
    integer :: iSwitch
    nullify(iCX,iCY,iCZ,iFX,iFY,iFZ,iTX,iTY,iTZ)
        iCX => iCenterX
        iCY => iCenterY
        iCZ => iCenterZ
        iFX => iCenterX
        iFY => iCenterY
        iFZ => iCenterZ
        iTX => iOne
        iTY => iOne
        iTZ => iOne
    if(Conf%UseVariableTime == 1) then
        call CalcTimeWidth(Geom,CC)
        if(Conf%UseLocalTimeStep == 1) then
            iTX => iCenterX
            iTY => iCenterY
            iTZ => iCenterZ
        else
            CC%TimeWidth(1,1,1) = minval(CC%TimeWidth)

        end if

    else
        CC%TimeWidth(1,1,1) = FixedTimeStep
    end if

    if(Conf%UseRRK2 == 1) then
        iSwitch = mod(iStep,2)
        call RationalRungeKutta2(Geom,CE,CC,iSwitch)
    else
        call ExplicitEulerMethod(Geom,CE,CC)
    end if

    nullify(iCX,iCY,iCZ,iFX,iFY,iFZ,iTX,iTY,iTZ)

return
end subroutine TimeIntegral
