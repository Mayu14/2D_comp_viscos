!***********************************/
!	Name:RANSの共通操作を行うプログラム
!	Alias:URansCommon
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.11.18
!	Update:
!	Other:
!***********************************/
subroutine URansCommon(UConf, UG, UCC, UCE)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, gamma => SpecificOfHeatRatio, KC => KarmanConstant
    use FrequentOperation
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(CellEdge), intent(inout) :: UCE
    logical :: debug = .true.
    ! RANS common
    ! Calc Laminar Viscosity from Sutherland's Law
    call UGetLaminarViscosity_mk2(UConf, UG, UCC, UCE)

    ! Calc Strain Rate Tensor & AbsoluteVortisity
    call UGetStrainRateTensor_edge(UConf, UG, UCC, UCE)

    if(debug == .true.) then
        do iCell = 1, UG%GI%RealCells
            UCC%StrainRateTensor(1,1,iCell,1,1) = (UCE%StrainRateTensor(1,1,UG%Tri%Edge(iCell,1),1,1) + UCE%StrainRateTensor(1,1,UG%Tri%Edge(iCell,2),1,1) + UCE%StrainRateTensor(1,1,UG%Tri%Edge(iCell,3),1,1)) / 3.0d0
            UCC%StrainRateTensor(1,2,iCell,1,1) = (UCE%StrainRateTensor(1,2,UG%Tri%Edge(iCell,1),1,1) + UCE%StrainRateTensor(1,2,UG%Tri%Edge(iCell,2),1,1) + UCE%StrainRateTensor(1,2,UG%Tri%Edge(iCell,3),1,1)) / 3.0d0
            UCC%StrainRateTensor(2,1,iCell,1,1) = (UCE%StrainRateTensor(2,1,UG%Tri%Edge(iCell,1),1,1) + UCE%StrainRateTensor(2,1,UG%Tri%Edge(iCell,2),1,1) + UCE%StrainRateTensor(2,1,UG%Tri%Edge(iCell,3),1,1)) / 3.0d0
            UCC%StrainRateTensor(2,2,iCell,1,1) = (UCE%StrainRateTensor(2,2,UG%Tri%Edge(iCell,1),1,1) + UCE%StrainRateTensor(2,2,UG%Tri%Edge(iCell,2),1,1) + UCE%StrainRateTensor(2,2,UG%Tri%Edge(iCell,3),1,1)) / 3.0d0
        end do

        do iCell = 1, UG%GI%RealCells
            UCC%AbsoluteVortisity(iCell,1,1) = (UCE%AbsoluteVortisity(UG%Tri%Edge(iCell,1),1,1)+UCE%AbsoluteVortisity(UG%Tri%Edge(iCell,2),1,1)+UCE%AbsoluteVortisity(UG%Tri%Edge(iCell,3),1,1))/3.0d0
        end do
    end if

    if(UConf%TurbulenceModel == 1) then
        call UBaldwinLomax_main(UConf, UG, UCC, UCE)
    end if

    call URansFlux(UG, UCC, UCE)

    return

contains

    subroutine URansFlux(UG, UCC, UCE)
        implicit none
        type(UnstructuredGrid), intent(in) :: UG
        type(CellCenter), intent(in) :: UCC
        type(CellEdge), intent(inout) :: UCE

        double precision :: invRe, invGmin1Mach2Pr
        double precision :: tau_xx, tau_xy, tau_yy, beta_x, beta_y, Viscosity, ThermalConductivity
        double precision, allocatable :: average_velocity(:)

        integer :: iFrontLocalEdge, iBackLocalEdge
        double precision :: dc2dn, dc2dx, dc2dy, length

        invRe = 1.0d0 / ReynoldsNumber

        allocate(average_velocity(2))

        do iEdge = 1, UG%GI%Edges
            Viscosity = UCE%LaminarViscosity(iEdge, 1, 1) + UCE%EddyViscosity(iEdge, 1, 1)
            ThermalConductivity = UCE%LaminarViscosity(iEdge, 1, 1) / LaminarPrandtlNumber &
                              & + UCE%EddyViscosity(iEdge, 1, 1) / TurbulentPrandtlNumber
            ! 界面の平均流速(2次精度中心差分)
            average_velocity = 0.5d0 *(UCE%RebuildQunatity(2:3, 1, 1, 1, iEdge) + UCE%RebuildQunatity(2:3, 1, 1, 2, iEdge))

            call UCentralDifferencePrepareAroundFace(UG, iEdge, iFrontCell, iFrontLocalEdge, iBackCell, iBackLocalEdge, length)
            ! 界面に対し法線方向(eta方向)の温度微分
            dc2dn = gamma * (UCC%PrimitiveVariable(5,iFrontCell,1,1) / UCC%PrimitiveVariable(1,iFrontCell,1,1)  &
                           - UCC%PrimitiveVariable(5,iBackCell,1,1) / UCC%PrimitiveVariable(1,iBackCell,1,1))

            ! 温度微分を座標変換により各座標方向に分配
            dc2dx = dc2dn * UG%GM%Normal(iEdge, 1)
            dc2dy = dc2dn * UG%GM%Normal(iEdge, 2)

            ! せん断応力の計算
            tau_xx = 2.0d0 / 3.0d0 * Viscosity * (2.0d0 * UCE%StrainRateTensor(1, 1, iEdge, 1, 1) - UCE%StrainRateTensor(2, 2, iEdge, 1, 1))
            tau_xy = Viscosity * (UCE%StrainRateTensor(1, 2, iEdge, 1, 1) + UCE%StrainRateTensor(2, 1, iEdge, 1, 1))
            tau_yy = 2.0d0 / 3.0d0 * Viscosity * (2.0d0 * UCE%StrainRateTensor(2, 2, iEdge, 1, 1) - UCE%StrainRateTensor(1, 1, iEdge, 1, 1))

            beta_x = average_velocity(1) * tau_xx + average_velocity(2) * tau_xy + ThermalConductivity * gamma * InverseGmin1 * dc2dn
            beta_y = average_velocity(1) * tau_xy + average_velocity(2) * tau_yy + ThermalConductivity * gamma * InverseGmin1 * dc2dy

            UCE%NormalFluxDiff(2, 1, 1, 1, iEdge) = UCE%NormalFluxDiff(2, 1, 1, 1, iEdge) &
                                                   & - invRe * (tau_xx * UG%GM%Normal(iEdge, 1) + tau_xy * UG%GM%Normal(iEdge, 2))

            UCE%NormalFluxDiff(3, 1, 1, 1, iEdge) = UCE%NormalFluxDiff(3, 1, 1, 1, iEdge) &
                                                   & - invRe * (tau_xy * UG%GM%Normal(iEdge, 1) + tau_yy * UG%GM%Normal(iEdge, 2))

            UCE%NormalFluxDiff(5, 1, 1, 1, iEdge) = UCE%NormalFluxDiff(5, 1, 1, 1, iEdge) &
                                                   & - invRe * (beta_x * UG%GM%Normal(iEdge, 1) + beta_y * UG%GM%Normal(iEdge, 2))
        end do

        return
    end subroutine URansFlux

end subroutine URansCommon
