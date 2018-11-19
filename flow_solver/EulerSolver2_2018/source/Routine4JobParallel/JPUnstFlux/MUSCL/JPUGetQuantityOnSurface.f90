!***********************************/
!	Name:流束制限関数を求いて，界面にTVDを満足する基礎変数を再構築するためのプログラム
!	Alias:GetQuantityOnSurface
!	Description:
!	Type:CellEdge
!	Input:Configulation,Geometory,CellCenter,CellEdge
!	Output:CE%RebuildQuantity
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2018.01.11
!	Other:仮想格子に関するif文等の修正
!***********************************/
subroutine JPUGetQuantityOnSurface(UG,UCC,UCE)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, only:CoreNumberOfCPU
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(in) :: UCC
    type(CellEdge), intent(inout) :: UCE

    do iEdge=1, UG%GI%Edges
        iFrontCell = UG%Line%Cell(iEdge,1,1)
        iBackCell = UG%Line%Cell(iEdge,2,1)

        UCE%RebuildQunatity(:,1,1,1,iEdge) = UCC%PrimitiveVariable(:,iFrontCell,1,1) &
            & + UCC%LimiterFunction(:,iFrontCell,1,1)*UCE%NormalGradient(:,1,1,1,iEdge)

        if(iBackCell <= UG%GI%RealCells) then
            UCE%RebuildQunatity(:,1,1,2,iEdge) = UCC%PrimitiveVariable(:,iBackCell,1,1) &
                & + UCC%LimiterFunction(:,iBackCell,1,1)*UCE%NormalGradient(:,1,1,2,iEdge)
        else !iBackCellが仮想格子
            UCE%RebuildQunatity(:,1,1,2,iEdge) = UCC%PrimitiveVariable(:,iBackCell,1,1)
        end if

    end do

    return
end subroutine JPUGetQuantityOnSurface
