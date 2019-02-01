!***********************************/
!	Name:局所時間刻み幅用の時間刻み計算プログラム
!	Alias:UCalcTimeWidth
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.22
!	Update:
!	Other:
!***********************************/
subroutine JPUCalcTimeWidth(UG,CC)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, CFL => CourantFriedrichsLewyCondition
implicit none
    type(UnstructuredGrid),intent(in) :: UG
    type(CellCenter), intent(inout) :: CC
    type(CalcTimeWidthWithOMP) :: CTW

    do iCell = 1, UG%GI%RealCells
        CTW%AbsoluteVelocity2 = dot_product(CC%ConservedQuantity(2:4,iCell,1,1),CC%ConservedQuantity(2:4,iCell,1,1))/(CC%ConservedQuantity(1,iCell,1,1)**2) !速度L2ノルム

        if(Gamma*Gmin1*(CC%ConservedQuantity(5,iCell,1,1) / CC%ConservedQuantity(1,iCell,1,1) &
                            & - 0.5d0 * CTW%AbsoluteVelocity2) < 0.0d0) then
                write(6,*) "CalcTimeWidth, Error "
                write(6,*) iCell
                write(6,*) CC%ConservedQuantity(1,iCell,1,1),CC%ConservedQuantity(5,iCell,1,1)
                write(6,*) CTW%AbsoluteVelocity2
                !stop
        end if

        CTW%LocalSoundSpeed = sqrt(Gamma*Gmin1*(CC%ConservedQuantity(5,iCell,1,1) / CC%ConservedQuantity(1,iCell,1,1) &
                            & - 0.5d0 * CTW%AbsoluteVelocity2)) !局所音速

        CTW%AbsoluteVelocity = sqrt(CTW%AbsoluteVelocity2) !速度絶対値

        CC%TmpTimeWidth(iCell,1,1,1) = CFL*UG%InscribedCircle(iCell) / (CTW%LocalSoundSpeed + CTW%AbsoluteVelocity) !衝撃波が通過する速度(の半分くらい)
    end do

    do iCell = UG%GI%RealCells+1, UG%GI%AllCells
        CC%TmpTimeWidth(iCell,1,1,1) = CC%TmpTimeWidth(UG%VC%Cell(iCell,1),1,1,1)
    end do

    do iCell = 1, UG%GI%RealCells
        CC%TimeWidth(iCell,1,1) = 0.5d0*min(CC%TmpTimeWidth(iCell,1,1,1),CC%TmpTimeWidth(UG%Tri%Cell(iCell,1),1,1,1),&
                     &   CC%TmpTimeWidth(UG%Tri%Cell(iCell,2),1,1,1),CC%TmpTimeWidth(UG%Tri%Cell(iCell,3),1,1,1))
    end do

    return
end subroutine JPUCalcTimeWidth
