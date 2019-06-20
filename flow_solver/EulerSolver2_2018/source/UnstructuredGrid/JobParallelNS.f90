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
    !$ use omp_lib
    implicit none
    type(Configulation), intent(inout) :: UConf
    type(CellCenter) :: UCC
    type(CellEdge) :: UCE
    type(UnstructuredGrid) :: UG
    integer :: iStep, iStartStep, iLoop, iSplit, iStep4Plot, iCalcStep, iTry !時間分割
    integer :: iAccel = 1
    double precision, allocatable :: obj_velocity(:)
    type(AeroCharacteristics) :: UAC
    double precision :: time_start, time_step, time_elapsed, time_specified, residual
    time_specified = 70200.0d0
    !$ time_start = omp_get_wtime()
! 初期化
    UCC%ConvergeCondition = Converge_Method
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

    if (UConf%UseResume /= 1) then
        write(6,*) "Velocity correction"
        allocate(obj_velocity(3))
        obj_velocity(:) = - UG%GM%BC%InFlowVariable(2:4)
        call RelativeCoordinateTransform(UG, UCC, obj_velocity)
    end if

    write(6,*) UG%GM%BC%InFlowVariable
    write(6,*) UG%GM%BC%OutFlowVariable

    call CheckNaN(UConf, UCC)   !
    allocate(UAC%coefficient(2, IterationNumber))

    if(DetailedReport > 0) then
        call show_setting(UConf)
    end if

    do iTry = 1, IterationNumber
        !$ time_elapsed = omp_get_wtime()
        iStartStep = 1
        ! 時間積分
        do iStep = (iTry-1)*OutputInterval + iStartStep, iTry*OutputInterval
            if(DetailedReport > 2) then
                if(mod(iStep, 100) == 0) then
                    write(6,*) iStep, "/", OutputInterval*IterationNumber, " th iteration"
                end if
            end if

            call UnstNS(iStep,UConf,UG,UCC,UCE)

        ! 画面出力＆NaN発生確認（発生時はデータをロールバックする）
            if(mod(iStep, CheckNaNInterval) == 0)then
                if(DetailedReport > 2) write(6,*) "NaN Checking..."
                call CheckNaN(UConf, UCC)
            end if
            !$ time_step = omp_get_wtime()
            !$ time_elapsed = time_step - time_start
            if(time_elapsed > time_specified) then
                UCC%iEndFlag = 4    ! resume
                UConf%OutputStatus = 2
                exit
            end if
            if(UCC%iEndFlag == 3) exit
        end do

        call set_residual(UCC, UAC, residual)
        write(6,*) UCC%ConvergeCondition, residual, Converge_tolerance

        if(residual < Converge_tolerance) then
            if(UConf%UseMUSCL == 0) then
                UConf%UseMUSCL = 1
                UCC%iEndFlag = 3    ! debug
            else
                UCC%iEndFlag = 3    ! finish
                UConf%OutputStatus = 1
            end if
        end if

        if(UConf%SwitchProgram /= 7) then
            call JPUOutput(UConf,UG,UCC,iTry)
            !call UOutput_Residuals(UConf, UCC%RH, -1)
        end if

        allocate(UAC%pressure_coefficient(UG%GM%BC%iWallTotal, 1))

        if(UCC%iEndFlag == 2) then
            call UnstNS(iTry, UConf, UG, UCC, UCE)
        end if

        call UCalcAeroCharacteristics(UConf, UCC, UG, iTry, UAC)  ! debug

        !call JPUOutput(UConf,UG,UCC,iCalcStep)
        call UOutput_Characteristics(UConf, UG, UAC, iTry)   ! debug

        deallocate(UAC%pressure_coefficient)   ! debug
        !UConf%UseLocalTimeStep = 1
        !UConf%UseSteadyCalc = 1
        !$ time_step = omp_get_wtime()
        !$ time_elapsed = time_step - time_start
        !$ write(6,*) time_elapsed, time_specified
        if(time_elapsed > time_specified) then
            UCC%iEndFlag = 3
        end if
        if(UCC%iEndFlag == 3) exit
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

    subroutine show_setting(UConf)
        use StructVar_Mod
        use ConstantVar_Mod
        implicit none
            type(Configulation), intent(in) :: UConf
            write(6,*) "--Setting--"
            write(6,*) "Detailed Report :", DetailedReport
            if(UConf%UseFluxMethod == 0) then
                write(6,*) "Inviscid Flux: Roe"
            else if(UConf%UseFluxMethod == 1) then
                write(6,*) "Inviscid Flux: SLAU2"
            end if
            if(UConf%UseMUSCL == 0) then
                write(6,*) "Accuracy : space:1st, time:1st"
            else
                write(6,*) "Accuracy : space:2nd(Venkatakrishnan), time:1st"
            end if
            if(UConf%UseVariableTime == 0) then
                write(6,*) "Uniform time step :", DefaultTimeStep
            else
                write(6,*) "variable time step"
            end if
            write(6,*) "CFL Number :", CourantFriedrichsLewyCondition
            write(6,*) "OutputInterval :", OutputInterval
            write(6,*) "CheckNaNInterval :", CheckNaNInterval
            write(6,*) "ConvergeMethod :", UCC%ConvergeCondition
            write(6,*) "Resume :", UConf%UseResume
        return
    end subroutine show_setting

    subroutine set_residual(UCC, UAC, residual)
        use StructVar_Mod
        implicit none
        type(CellCenter), intent(in) :: UCC
        type(AeroCharacteristics), intent(in) :: UAC
        double precision, intent(out) :: residual

        if(UCC%ConvergeCondition == 0) then
            residual = UCC%RH%AveResidual(6, UCC%RH%iTime - 1)
        else if(UCC%ConvergeCondition == 1) then
            residual = UCC%RH%MaxResidual(6, UCC%RH%iTime - 1)
        else if(UCC%ConvergeCondition == 2) then
            residual = UCC%RH%TotalResidual(6, UCC%RH%iTime - 1)
        else if(UCC%ConvergeCondition == 3) then
            residual = UAC%residual
        else if(UCC%ConvergeCondition == 4) then
            residual = max(UCC%RH%AveResidual(6, UCC%RH%iTime - 1), UCC%RH%MaxResidual(6, UCC%RH%iTime - 1), &
                        &   UCC%RH%TotalResidual(6, UCC%RH%iTime - 1), UAC%residual)
        end if

        return
    end subroutine set_residual

end subroutine JobParallelNS

