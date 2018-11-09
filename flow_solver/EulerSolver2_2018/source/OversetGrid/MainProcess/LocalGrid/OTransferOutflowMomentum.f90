!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:OOverSetS2U
!	Description:
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.11.16
!	Update:2017.11.16
!	Other:
!***********************************/
subroutine OTransferOutflowMomentum(UG,UCC,MG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(MoveGrid), intent(inout) :: MG

    do iCell=1, UG%GI%AllCells
        UCC%ConservedQuantity(2:4,iCell,1,1) = matmul(MG%TM%Top2Next(1:3,1:3),UCC%ConservedQuantity(2:4,iCell,1,1))
    end do

return
end subroutine OTransferOutflowMomentum
