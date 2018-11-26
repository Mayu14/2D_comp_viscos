!***********************************/
!	Name:計算設定を読み込むプログラム
!	Alias:ReadConfigulation
!	Description:専用のCalcConfigファイルが必要
!	Type:Configulation
!	Input:CalcConfig(外部入力)
!	Output:Configulation
!	Note:CalcConfigはプログラム本体と同じディレクトリに置くこと
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.07
!	Other:
!***********************************/
subroutine JobParallelNS(UConf)
    use StructVar_Mod
    use ConstantVar_Mod
    implicit none
    type(Configulation), intent(inout) :: UConf
    type(CellCenter) :: UCC
    type(CellEdge) :: UCE
    type(UnstructuredGrid) :: UG
    integer :: iStep, iStartStep, iLoop, iSplit, iStep4Plot, iCalcStep !時間分割
    integer :: iAccel = 1
    double precision, allocatable :: obj_velocity(:)
    type(AeroCharacteristics) :: UAC

    call UReadUnStrGrid(UConf,UCC,UCE,UG)
    call UInitialize(UConf,UG,UCC) !ReadInitial,MakeInternalBoundary
    iStep = 0
    call UReadInflowOutflowCondition(UG, UConf)
    call UPrepareBoundary(UG, UCC)
    call JPUConserve2Primitive(UG, UCC)

    allocate(UAC%coefficient(2, int(IterationNumber)))
    allocate(obj_velocity(3))

    ! for accelaration area
    UConf%UseLocalTimeStep = 0
    UConf%UseMUSCL = 0
    UConf%TurbulenceModel = 0
    !
    iCalcStep = 0
    iStartStep = 1
    do iStep = iStartStep, IterationNumber
        if(UConf%UseLocalTimeStep /= 0 .and. iAccel == 0) then  ! steady area
            iCalcStep = iCalcStep + 1   ! tmp
            call RelativeCoordinateTransform(UG, UCC, obj_velocity)
            call UnstNS(iStep,UConf,UG,UCC,UCE)
            call RelativeCoordinateTransform(UG, UCC, -obj_velocity)

        else    ! accelaration area
            call JPUCheckCFL4FixedTime(UG,UCC,iSplit)
            FixedTimeStep = FixedTimeStep / dble(iSplit)

            do iLoop=1, iSplit
                iCalcStep = iCalcStep + 1
                obj_velocity(:) = - min(0.01d0 * dble(iCalcStep), 1.0d0) * UG%GM%BC%InFlowVariable(2:4)

                if(iLoop == 1) call RelativeCoordinateTransform(UG, UCC, obj_velocity)

                    call UnstNS(iStep,UConf,UG,UCC,UCE)

                if(iLoop == iSplit) call RelativeCoordinateTransform(UG, UCC, -obj_velocity)
            end do
            FixedTimeStep = DefaultTimeStep
        end if

        if(iCalcStep > 100) then
            iAccel = 0  ! 加速区間終了
            UConf%UseLocalTimeStep = 1
            UConf%UseMUSCL = 1
            UConf%TurbulenceModel = 1
            obj_velocity = - UG%GM%BC%InFlowVariable(2:4)
        end if

        if(mod(iStep,OutputInterval) == 0) then
            iStep4Plot = iStep / OutputInterval
            !call JPUOutput(UConf,UG,UCC,iStep4Plot)
        end if
        call UCalcAeroCharacteristics(UCC, UG, iStep4Plot, UAC)
    end do

    call JPUOutput(UConf,UG,UCC,iStep)
    call UOutput_Characteristics(UConf, UAC)

    UConf%UseLocalTimeStep = 1
    UConf%UseMUSCL = 1
    UConf%TurbulenceModel = 1

    return
contains

    subroutine UnstNS(iStep, OConf, UG, UCC, UCE)
        implicit none
        integer, intent(in) :: iStep
        type(Configulation), intent(in) :: OConf
        type(UnstructuredGrid), intent(in) :: UG
        type(CellCenter), intent(inout) :: UCC
        type(CellEdge), intent(inout) :: UCE

            call USetBoundary(UG,UCC)

            call UCalcFlux(OConf,UG,UCC,UCE)

            call UTimeIntegral(OConf,UG,UCE,UCC,iStep)

        return
    end subroutine UnstNS

    subroutine RelativeCoordinateTransform(UG, CC4MB, ObjectVelocity)
        implicit none
        type(UnstructuredGrid), intent(in) :: UG
        type(CellCenter), intent(inout) :: CC4MB
        double precision, intent(in) :: ObjectVelocity(:)
        integer :: iCell

        do iCell=1, UG%GI%AllCells

            CC4MB%ConservedQuantity(2:4,iCell,1,1) = CC4MB%ConservedQuantity(1,iCell,1,1)&
                &  * (CC4MB%PrimitiveVariable(2:4,iCell,1,1) - ObjectVelocity(1:3))

            CC4MB%ConservedQuantity(5,iCell,1,1) = &
                & InverseGmin1*CC4MB%PrimitiveVariable(5,iCell,1,1) + 0.5d0*CC4MB%ConservedQuantity(1,iCell,1,1) &
                & * (dot_product(CC4MB%PrimitiveVariable(2:4,iCell,1,1) - ObjectVelocity(1:3),CC4MB%PrimitiveVariable(2:4,iCell,1,1) - ObjectVelocity(1:3)))

        end do

        return
    end subroutine RelativeCoordinateTransform

    function accel_area_velocity(velocity_gradient, time, steady_velocity) result(velocity)
        implicit none
        double precision, intent(in) :: velocity_gradient, time, steady_velocity
        double precision :: velocity

        velocity = min(velocity_gradient * time, steady_velocity)

        return
    end function accel_area_velocity

end subroutine JobParallelNS

