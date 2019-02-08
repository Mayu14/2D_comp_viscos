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
    integer :: iStep, iStartStep, iLoop, iSplit, iStep4Plot, iCalcStep, iTry !時間分割
    integer :: iAccel = 1
    double precision, allocatable :: obj_velocity(:)
    type(AeroCharacteristics) :: UAC

    call UReadUnStrGrid(UConf,UCC,UCE,UG)

    if(UConf%SwitchProgram == 7) then
        UConf%UseResume = 1
        IterationNumber = 1
    end if

    call UInitialize(UConf,UG,UCC) !ReadInitial,MakeInternalBoundary
    iStep = 0
    call UReadInflowOutflowCondition(UG, UConf)
    call UPrepareBoundary(UG, UCC)
    call JPUConserve2Primitive(UG, UCC)

    allocate(obj_velocity(3))
    obj_velocity(:) = - UG%GM%BC%InFlowVariable(2:4)
    call RelativeCoordinateTransform(UG, UCC, obj_velocity)

    call CheckNaN(UConf, UCC)   !

    do iTry = 1, 100
        !call JPUOutput(UConf,UG,UCC,0)  ! debug
        iStartStep = 1
        do iStep = (iTry-1)*IterationNumber + iStartStep, iTry*IterationNumber
            if(DetailedReport > 2) then
                if(mod(IterationNumber, 10) == 0) then
                    write(6,*) iStep, "/", IterationNumber, " th iteration"
                end if
            end if

            call UnstNS(iStep,UConf,UG,UCC,UCE)

            if(mod(iStep,OutputInterval) == 0) then
                iStep4Plot = iStep / OutputInterval
                call JPUOutput(UConf,UG,UCC,iStep4Plot)
            end if

            if(mod(iStep, CheckNaNInterval) == 0)then
                if(DetailedReport > 2) write(6,*) "NaN Checking..."
                call CheckNaN(UConf, UCC)
            end if

            if(UCC%iEndFlag > 1) exit
        end do
    !end do

        UCC%iEndFlag = 3

        if(UConf%SwitchProgram /= 7) then
            call JPUOutput(UConf,UG,UCC,0)
        end if

        if (UCC%iEndFlag == 2) then
            iCalcStep = 100
        else if(UCC%iEndFlag == 3) then
            iCalcStep = 1
        end if

        allocate(UAC%coefficient(2, iCalcStep))
        allocate(UAC%pressure_coefficient(UG%GM%BC%iWallTotal, iCalcStep))

        !UConf%UseLocalTimeStep = 0
        !UConf%UseSteadyCalc = 0
        if(UCC%iEndFlag == 2) then
            do iStep4Plot = 1, iCalcStep
                call UnstNS(iStep4Plot, UConf, UG, UCC, UCE)
                call UCalcAeroCharacteristics(UConf, UCC, UG, iStep4Plot, UAC)
            end do

        else if(UCC%iEndFlag == 3) then
            call UCalcAeroCharacteristics(UConf, UCC, UG, iCalcStep, UAC)  ! debug
        end if

        !call JPUOutput(UConf,UG,UCC,iCalcStep)
        call UOutput_Characteristics(UConf, UG, UAC, iTry)   ! debug
        deallocate(UAC%coefficient, UAC%pressure_coefficient)   ! debug
        !UConf%UseLocalTimeStep = 1
        !UConf%UseSteadyCalc = 1
    end do

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
            !write(6,*) "Compute Flux..."
            call UCalcFlux(OConf,UG,UCC,UCE)
            !write(6,*) "Time integration..."
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

    subroutine CheckNaN(UConf, UCC)
        implicit none
        type(Configulation), intent(inout) :: UConf
        type(CellCenter), intent(inout) :: UCC
        integer :: iVariable, iCell

            do iCell = 1, UG%GI%RealCells
                do iVariable = 1, 5
                    if(isnan(UCC%ConservedQuantity(iVariable,iCell,1,1))) then
                        RetryFlag = 1
                    end if
                    if(RetryFlag == 1) exit
                end do
            end do

            if(RetryFlag == 0) then
                UCC%PastQuantity = UCC%ConservedQuantity
            else
                UCC%ConservedQuantity = UCC%PastQuantity
                if(UConf%UseJobParallel == 1) then
                    call JPUConserve2Primitive(UG, UCC)
                else
                    call UConserve2Primitive(UG, UCC)
                end if
                CourantFriedrichsLewyCondition = 0.5d0*CourantFriedrichsLewyCondition
                RetryFlag = 0
                CheckNaNInterval = max(ceiling(float(CheckNaNInterval) / 2.0), 1)
                if(DetailedReport > 1) then
                    write(6,*) "NaN detected."
                    write(6,*) "CFL Number ", CourantFriedrichsLewyCondition*2, " -> ", CourantFriedrichsLewyCondition
                    write(6,*) "Check interval ", CheckNaNInterval
                    write(6,*) ""
                end if
            end if

            if(CourantFriedrichsLewyCondition < 10.0**(-6)) then
                UCC%iEndFlag = 3
                if(DetailedReport > 0) then
                    write(6,*) "CFL Number ", CourantFriedrichsLewyCondition
                    write(6,*) "Finish Calculate"
                end if
            end if

        return
    end subroutine CheckNaN

end subroutine JobParallelNS

