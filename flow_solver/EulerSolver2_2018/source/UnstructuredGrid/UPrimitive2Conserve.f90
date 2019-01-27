!***********************************/
!	Name:基礎変数から保存変数へ変換するプログラム
!	Alias:UPrimitive2Conserve
!	Description:定義通りに計算するだけ
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.06
!	Other:
!***********************************/
subroutine UPrimitive2Conserve(UG,UCC)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter),intent(inout) :: UCC

    iVar => iVariable
    iDim = UG%GM%Dimension

    do iCell=1,UG%GI%AllCells
        UCC%ConservedQuantity(1,iCell,1,1) = UCC%PrimitiveVariable(1,iCell,1,1)

        UCC%ConservedQuantity(2:iDim+1,iCell,1,1) = &
            &   UCC%PrimitiveVariable(2:iDim+1,iCell,1,1)*UCC%PrimitiveVariable(1,iCell,1,1)


        UCC%ConservedQuantity(iDim+2,iCell,1,1) = &
            &   InverseGmin1 * (UCC%PrimitiveVariable(iDim+2,iCell,1,1) + 0.5d0 * UCC%PrimitiveVariable(1,iCell,1,1) &
            &   * dot_product(UCC%PrimitiveVariable(2:iDim+1,iCell,1,1),UCC%PrimitiveVariable(2:iDim+1,iCell,1,1)))
    end do

return
end subroutine UPrimitive2Conserve
