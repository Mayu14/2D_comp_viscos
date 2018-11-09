!***********************************/
!	Name:1段階目の計算を行うにあたって必要な情報を作成する
!	Alias:DataFabrication4FirstStep
!	Description:初回の計算をさもループ計算の途中であったかのように振る舞わせる → データ捏造っぽさみ → fabrication(偽造)
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.27
!	Update:
!	Other:
!***********************************/
subroutine DataFabrication4FirstStep(UG,MG,CC)
    use StructVar_Mod
    use LoopVar_Mod

    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(MoveGrid), intent(inout) :: MG
    type(CellCenter), intent(inout) :: CC

    !MG%NSR%C%NextCoordsInG(:,1:3) = UG%CD%Cell(:,1:3)
    MG%NSR%C%NextCoordsInG(:,4) = 1.0d0

    !MG%NSR%P%NextCoordsInG(:,1:3) = UG%CD%Point(:,1:3)
    MG%NSR%P%NextCoordsInG(:,4) = 1.0d0

    do iCell=UG%GI%RealCells+1, UG%GI%AllCells
        CC%ConservedQuantity(:,iCell,1,1) = CC%ConservedQuantity(:,UG%VC%Cell(iCell,1),1,1)
    end do

!    do iCell=1, UG%GI%RealCells
!        MG%IS(iCell)%InfluenceDepth = 0
!    end do

return
end subroutine DataFabrication4FirstStep
