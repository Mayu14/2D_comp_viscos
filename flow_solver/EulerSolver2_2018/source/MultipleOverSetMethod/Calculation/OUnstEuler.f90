!***********************************/
!	Name:非構造オイラーソルバーの重合格子で使う部分だけ
!	Alias:OUnstEuler
!	Description:
!	Type:
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.11.28
!	Update:2017.11.28
!	Other:
!***********************************/
subroutine OUnstEuler(iStep,OConf,UG,UCC,UCE,SW)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod
!$ use omp_lib
    implicit none
    integer, intent(in) :: iStep
    type(Configulation), intent(in) :: OConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(CellEdge), intent(inout) :: UCE
    type(StopWatch), intent(inout) :: SW

            !if(OConf%UseSteadyCalc == 1) call TestVelocityCorrection(1)
            call USetBoundary(UG,UCC)


!$ SW%LapTime(1,2) = omp_get_wtime()
            call UCalcFlux(OConf,UG,UCC,UCE)
!$ SW%LapTime(2,2) = omp_get_wtime()

!$ SW%LapTime(1,3) = omp_get_wtime()
            call UTimeIntegral(OConf,UG,UCE,UCC,iStep)
!$ SW%LapTime(2,3) = omp_get_wtime()
!$ SW%LapTime(3,2) = SW%LapTime(3,2) + SW%LapTime(2,2) - SW%LapTime(1,2)
!$ SW%LapTime(3,3) = SW%LapTime(3,3) + SW%LapTime(2,3) - SW%LapTime(1,3)
            !if(OConf%UseSteadyCalc == 1) call TestVelocityCorrection(2)

return
contains
    subroutine TestVelocityCorrection(iSwitch)
    implicit none
        integer, intent(in) :: iSwitch
        double precision, allocatable :: Velocity(:)
        double precision :: Distance, InternalWindWeight,Gain,SigmoidX
        allocate(Velocity(3))

        Velocity(1) = -UG%GM%BC%InFlowVariable(2)
        Velocity(2) = -UG%GM%BC%InFlowVariable(3)
        Velocity(3) = 0.0d0

        if(iSwitch == 2) Velocity = - Velocity

        call UConserve2Primitive(UG,UCC)



        do iCell = 1, UG%GI%RealCells
            !Distance = sqrt(dot_product(UG%CD%Cell(iCell,1:3),UG%CD%Cell(iCell,1:3)))

            !if(iCell < UG%GI%RealCells+1) then
                !Gain = -2.0d0/(UG%GM%AverageWidth(iCell))*log((3.0d0-sqrt(3.0d0))/(3.0d0+sqrt(3.0d0)))
                !SigmoidX = Gain*(Distance - (0.2d0 + UG%GM%AverageWidth(iCell)) )
            !else
                !Gain = -2.0d0/(UG%GM%AverageWidth(UG%VC%Cell(iCell,1)))*log((3.0d0-sqrt(3.0d0))/(3.0d0+sqrt(3.0d0)))
                !SigmoidX = Gain*(Distance - (0.2d0 + UG%GM%AverageWidth(UG%VC%Cell(iCell,1))))
            !end if

            !Error procedure
                !if(SigmoidX < - SigmoidRange) then !Error1
                    !InternalWindWeight = 1.0d0 - 10.0d0**(-15)

                !else if(SigmoidX > SigmoidRange) then !Error2
                    !InternalWindWeight = 0.0d0

                !else !Sigmoid
                    !InternalWindWeight = -1.0d0 / (1.0d0 + exp(-SigmoidX)) + 1.0d0

                !end if

            UCC%ConservedQuantity(2:4,iCell,1,1) = UCC%ConservedQuantity(1,iCell,1,1) &
                &   * (UCC%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3))!*InternalWindWeight)

            UCC%ConservedQuantity(5,iCell,1,1) = &
                & InverseGmin1*UCC%PrimitiveVariable(5,iCell,1,1) + 0.5d0*UCC%ConservedQuantity(1,iCell,1,1) &
                & * (dot_product(UCC%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3),UCC%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3)))
                !& * (dot_product(UCC%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3)*InternalWindWeight,UCC%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3)*InternalWindWeight))
        end do

    return
    end subroutine TestVelocityCorrection

end subroutine OUnstEuler
