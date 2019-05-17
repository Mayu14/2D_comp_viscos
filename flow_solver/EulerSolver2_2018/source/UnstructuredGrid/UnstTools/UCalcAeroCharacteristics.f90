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
    double precision :: scaling_factor
    double precision :: cosA, sinA
    double precision :: static_pressure, kinetic_pressure

    call JPUConserve2Primitive(UG, UCC)

    lift = 0.0d0
    drag = 0.0d0
    ! 迎角を+θ <=> 座標系を+θ <=> データを-θ
    cosA = cos(-UConf%dAttackAngle)
    sinA = sin(-UConf%dAttackAngle)

    static_pressure = UG%GM%BC%InFlowVariable(5)
    kinetic_pressure = 0.5d0 * UG%GM%BC%InFlowVariable(1) * dot_product(UG%GM%BC%InFlowVariable(2:4), UG%GM%BC%InFlowVariable(2:4))   ! for Cp

    do iWall = 1, UG%GM%BC%iWallTotal
        iEdge = UG%GM%BC%VW(iWall)%iGlobalEdge
        iFrontCell = UG%Line%Cell(iEdge, 1, 1)
        iBackCell = UG%Line%Cell(iEdge, 2, 1)
        WallPressure = 0.5d0 * (UCC%PrimitiveVariable(5, iFrontCell, 1, 1) + UCC%PrimitiveVariable(5, iBackCell, 1, 1))

        ! 法線ベクトルは物体内部方向へ向いているものが正になっているため-1倍する必要あり
        drag = drag + WallPressure * (UG%GM%Normal(iEdge, 1) * cosA + UG%GM%Normal(iEdge, 2) * sinA) * (-1.0d0) * UG%GM%Area(iEdge)

        lift = lift + WallPressure * (-UG%GM%Normal(iEdge, 1) * sinA + UG%GM%Normal(iEdge, 2) * cosA) * (-1.0d0) * UG%GM%Area(iEdge)

        !UAC%pressure_coefficient(iWall, iPlotStep) = (WallPressure - static_pressure) / kinetic_pressure
        UAC%pressure_coefficient(iWall, 1) = (WallPressure - static_pressure) / kinetic_pressure
    end do

    scaling_factor = kinetic_pressure * ObjLength

    UAC%coefficient(1, iPlotStep) = drag / scaling_factor
    UAC%coefficient(2, iPlotStep) = lift / scaling_factor

    if(iPlotStep > 1) then
        UAC%residual = sqrt(max((UAC%coefficient(1, iPlotStep) - UAC%coefficient(1, iPlotStep-1))**2, &
                         & (UAC%coefficient(2, iPlotStep) - UAC%coefficient(2, iPlotStep-1))**2))
    end if

    return
end subroutine UCalcAeroCharacteristics

