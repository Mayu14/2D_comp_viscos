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
subroutine OTransferInflowMomentum(UG,UCC,MG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(MoveGrid), intent(inout) :: MG

!    write(6,*) "Inflow rotation"

    do iCell=1, UG%GI%AllCells
        !UCC%ConservedQuantity(2:4,iCell,1,1) = matmul(MG%TM%Next2Top(1:3,1:3),UCC%ConservedQuantity(2:4,iCell,1,1))
        UCC%ConservedQuantity(2:4,iCell,1,1) = matmul(MG%TM%Top2Next(1:3,1:3),UCC%ConservedQuantity(2:4,iCell,1,1))
    end do

    !do iCell=UG%GI%RealCells+1, UG%GI%AllCells !全仮想セルに対し勾配なし条件を適用する !間違い
    !    UCC%ConservedQuantity(:,iCell,1,1) = UCC%ConservedQuantity(:,UG%VC%Cell(iCell,1),1,1)
    !end do

return
end subroutine OTransferInflowMomentum
