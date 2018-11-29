!***********************************/
!	Name:流れ場中の物体の空力特性を計算するプログラム
!	Alias:UCalcAeroCharacteristics
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.11.17
!	Update:
!	Other:
!***********************************/
subroutine UCalcAeroCharacteristics(UConf, UCC, UG, iPlotStep, UAC)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod
    use FrequentOperation
    implicit none
    type(Configulation), intent(in) :: UConf
    type(CellCenter), intent(inout) :: UCC
    type(UnstructuredGrid), intent(in) :: UG
    integer, intent(in) :: iPlotStep
    type(AeroCharacteristics), intent(inout) :: UAC
    integer :: iWall
    double precision :: WallPressure
    double precision :: lift, drag
    double precision :: scaling_factor, ObjLength = 1.0d0
    double precision :: cosA, sinA
    call JPUConserve2Primitive(UG, UCC)

    lift = 0.0d0
    drag = 0.0d0

    cosA = cos(UConf%dAttackAngle)
    sinA = sin(UConf%dAttackAngle)

    do iWall = 1, UG%GM%BC%iWallTotal
        iEdge = UG%GM%BC%VW(iWall)%iGlobalEdge
        iFrontCell = UG%Line%Cell(iEdge, 1, 1)
        iBackCell = UG%Line%Cell(iEdge, 2, 1)
        WallPressure = 0.5d0 * (UCC%PrimitiveVariable(5, iFrontCell, 1, 1) + UCC%PrimitiveVariable(5, iBackCell, 1, 1))


        drag = drag + WallPressure * (-UG%GM%Normal(iEdge, 1) * cosA + UG%GM%Normal(iEdge, 2) * sinA) * UG%GM%Area(iEdge)

        lift = lift + WallPressure * (-UG%GM%Normal(iEdge, 1) * sinA - UG%GM%Normal(iEdge, 2) * cosA) * UG%GM%Area(iEdge)

    end do

    scaling_factor = 0.5d0 * UG%GM%BC%InFlowVariable(1) * dot_product(UG%GM%BC%InFlowVariable(2:4), UG%GM%BC%InFlowVariable(2:4)) * ObjLength
    UAC%coefficient(1, iPlotStep) = drag / scaling_factor
    UAC%coefficient(2, iPlotStep) = lift / scaling_factor

    return
end subroutine UCalcAeroCharacteristics

