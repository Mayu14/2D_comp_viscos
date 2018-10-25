!***********************************/
!	Name:MayuGrid2形式の格子情報読み込みメソッド
!	Alias:GetCellAndEdgeCoordinate
!	Description:バージョン0.9に対応
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.??
!	Update:2018.10.25
!	Other:
!***********************************/
subroutine Set_Region(Conf, Geom, CC, CE)
    use StructVar_Mod_mk2
    implicit none

    type(Configulation), intent(in) :: Conf
    type(Geometry), intent(inout) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CellEdge), intent(inout) :: CE

    ! 格子点情報読み込み
    call read_mayugrid2(Conf%GridFileName, Geom)
    ! この段階でどの変数がどの大きさで必要かが分かるため一斉にallocateしておく
    call MallocAll(Conf, Geom, CC, CE)
    ! セル・界面の中心座標を算出
    call GetCellAndEdgeCoordinate(Geom)
    ! ヤコビ行列式を計算
    call GetGridJacobian(Geom, CC, CE)

    return
end subroutine Set_Region
