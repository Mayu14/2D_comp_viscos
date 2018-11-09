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
subroutine UTimeIntegral(UConf,UG,UCE,UCC,iStep)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod
implicit none
    type(Configulation),intent(in) :: UConf
    type(UnstructuredGrid),intent(in) :: UG
    type(CellEdge), intent(in) :: UCE
    type(CellCenter), intent(inout) :: UCC
    integer, intent(in) :: iStep
    integer :: iSwitch

    if(UConf%UseVariableTime == 1) then
        call UCalcTimeWidth(UG,UCC)
        if(UConf%UseLocalTimeStep == 1) then
            iTC => iCell
        else
            iTC => iOne
            UCC%TimeWidth(1,1,1) = minval(UCC%TimeWidth(:,1,1),1)
        end if

    else
        UCC%TimeWidth(1,1,1) = FixedTimeStep
        iTC => iOne
    end if

    if(iStep /= 1 .and. UConf%UseSteadyCalc == 1) then
        UCC%iEndFlag = 1
    else
        UCC%iEndFlag = 0
    end if

    if(UConf%UseRRK2 == 1) then
        iSwitch = mod(iStep,2)
        !call RationalRungeKutta2(UG,UCE,UCC,iSwitch)
    else
        call UExplicitEulerMethod(UG,UCE,UCC)
    end if

return
end subroutine UTimeIntegral
