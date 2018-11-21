!***********************************/
!	Name:外部境界のうち流入流速との内積が0以下になる境界を無反射境界に設定しなおす
!	Alias:UPrepareBoundary
!	Description:
!	Type:CellCenter
!	Input:Geometry%Dimension,Geom%CellNumber, CellCenter%ConservedQuantity
!	Output:CellCenter%ConservedQuantity
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.11.11
!	Update:2017.11.11
!	Other:
!***********************************/
subroutine UPrepareBoundary(UG,UCC)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
implicit none
    type(UnstructuredGrid), intent(inout) :: UG
    type(CellCenter), intent(inout) :: UCC

    iDim = UG%GM%Dimension

    do iCell = UG%GI%RealCells+1, UG%GI%AllCells
        if(UG%GM%CellType(iCell,1,1) == 1) then !Inflow and Outflow
            iEdge = UG%VC%Edge(iCell)
            if(dot_product(UG%GM%BC%InFlowVariable(2:4), UG%GM%Normal(iEdge, :)) <= 0.0d0) then    !　流入流束との内積が負になるとき流入境界→無反射境界に変更
                UG%GM%CellType(iCell, 1, 1) = 10   ! non reflect boundary
            end if
        end if
    end do

    return
end subroutine UPrepareBoundary
