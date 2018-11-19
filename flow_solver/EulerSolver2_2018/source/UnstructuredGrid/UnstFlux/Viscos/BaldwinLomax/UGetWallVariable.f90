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
subroutine GetWallVariable(UG, UCC, UCE, iWall, Wall_Density, Wall_Viscosity, Wall_dudy)
    use StructVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(in) :: UCC
    type(CellEdge), intent(in) :: UCE
    integer, intent(in) :: iWall
    double precision, intent(out) :: Wall_Density, Wall_Viscosity, Wall_dudy
    integer :: iWallEdge, iWallEdge_as_Local, iFrontCell, iBackCell
    double precision :: front_tangent_velocity, back_tangent_velocity, vertical_distance

    iWallEdge = UG%GM%BC%VW(iWall)%iGlobalEdge  ! 壁の大域界面番号
    iFrontCell = UG%GM%BC%VW(iWall)%iMemberCell(1)  ! 壁の物体側セル番号
    iBackCell = UG%Line%Cell(iWallEdge, 2, 1) ! 壁の仮想セル番号
    iWallEdge_as_Local = UG%Line%Cell(iWallEdge, 1, 2)  ! 物体外側セルから壁面セルを見たときの局所界面番号

    front_tangent_velocity = GetAbsTangentialVelocity_2D(UCC%PrimitiveVariable(2:3, iFrontCell, 1, 1), UG%GM%Normal(iWallEdge, 1:2))    ! 壁(物体外側)における接線方向速度の絶対値
    back_tangent_velocity = GetAbsTangentialVelocity_2D(UCC%PrimitiveVariable(2:3, iBackCell, 1, 1), UG%GM%Normal(iWallEdge, 1:2))  ! 壁(物体内部)における接線方向速度の絶対値
    vertical_distance = 2.0d0 * AbsVector(UG%GM%Width(iFrontCell, iWallEdge_as_Local, :))    ! 壁面を挟む2セルの中心間距離

    Wall_Density = 0.5d0 * (UCC%ConservedQuantity(1, iFrontCell, 1, 1) + UCC%ConservedQuantity(1, iBackCell, 1, 1)) ! 壁における密度
    Wall_Viscosity = UCE%LaminarViscosity(iWallEdge, 1, 1)   ! 壁における粘性
    Wall_dudy = abs((front_tangent_velocity - back_tangent_velocity)) / vertical_distance

    return

contains

    double precision function GetAbsTangentialVelocity_2D(velocity, normal) result (tangential_velocity)
        implicit none
        double precision, intent(in) :: velocity(:), normal(:)

        tangential_velocity = abs(velocity(1) * (-normal(2)) + velocity(2) * normal(1))
        return
    end function GetAbsTangentialVelocity_2D

end subroutine GetWallVariable
