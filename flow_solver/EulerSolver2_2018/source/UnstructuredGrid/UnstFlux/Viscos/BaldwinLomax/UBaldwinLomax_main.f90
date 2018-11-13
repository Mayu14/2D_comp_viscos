!***********************************/
!	Name:粘性流束を計算するためのプログラム
!	Alias:UBaldwinLomax_main
!	Description:
!	Type:CellEdge
!	Input:Configulation,Geometory,CellCenter,CellEdge
!	Output:CE%NormalFluxDiff
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:1次元と2次元,1次精度と2次精度のみに対応
!***********************************/
subroutine UBaldwinLomax_main(UConf,UG,UCC,UCE)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(CellEdge), intent(inout) :: UCE

! RANS common
    ! Calc Laminar Viscosity from Sutherland's Law

    ! Calc Strain Rate Tensor

! Baldwin-Lomax
    ! get y_max, F_max, and u_dif on each wall boundary respectively
    ! 壁番号→壁に所属する要素の総数，近い順に整列済みでセル番号の検索が可能，高速巡回が可能なように内部では配列にしておく
    ! Calc y_cross of Baldwin-Lomax on each wall boundary respectively

    ! Calc Turbulance Viscosity of Baldwin-Lomax Model


return
end subroutine UBaldwinLomax_main
