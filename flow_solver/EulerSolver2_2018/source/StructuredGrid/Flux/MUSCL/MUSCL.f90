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
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine MUSCL(Conf,Geom,CC,CE)
    use StructVar_Mod
    implicit none
    type(Configulation), intent(in) :: Conf
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CellEdge), intent(inout) :: CE

!基礎変数型式への変換
        call Conserve2Primitive(Geom,CC)

!勾配の計算
        call GetGradient(Geom,CC)
!勾配制限関数の計算
        call GetLimiter(Conf,Geom,CC,CE)
!セル境界での値の計算
        !CC%LimiterFunction = 0.0d0
        call GetQuantityOnSurface(Geom,CC,CE)
return
end subroutine MUSCL
