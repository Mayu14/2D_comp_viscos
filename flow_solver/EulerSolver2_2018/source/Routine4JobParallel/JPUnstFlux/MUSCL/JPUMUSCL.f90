!***********************************/
!	Name:基礎変数MUSCL法を用いて2次精度流束を求めるためのプログラム
!	Alias:MUSCL
!	Description:
!	Type:CellCenter,CellEdge
!	Input:Configulation,Geometory,CellCenter,CellEdge
!	Output:CE%RebuildQuantity
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2018.01.12
!	Other:
!***********************************/
subroutine JPUMUSCL(UConf,UG,UCC,UCE)
    use StructVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(CellEdge), intent(inout) :: UCE

!基礎変数型式への変換
        call UConserve2Primitive(UG,UCC)
    !end if
!勾配の計算
        call UGetGradient(UG,UCC)

!勾配制限関数の計算
        call JPUGetLimiter(UConf,UG,UCC,UCE)

!セル境界での値の計算

        call JPUGetQuantityOnSurface(UG,UCC,UCE)

    return
end subroutine JPUMUSCL
