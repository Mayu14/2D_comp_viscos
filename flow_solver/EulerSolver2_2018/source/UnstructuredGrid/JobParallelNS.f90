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
    type(Configulation), intent(in) :: UConf
    type(CellCenter) :: UCC
    type(CellEdge) :: UCE
    type(UnstructuredGrid) :: UG
    integer :: iStep, iStartStep, iLoop, iSplit, iStep4Plot !時間分割
    type(AeroCharacteristics) :: UAC

    call UReadUnStrGrid(UConf,UCC,UCE,UG)
    call UInitialize(UConf,UG,UCC) !ReadInitial,MakeInternalBoundary
    iStep = 0
    call UReadInflowOutflowCondition(UG, UConf)
    call JPUOutput(UConf,UG,UCC,iStep)

    allocate(UAC%coefficient(2, int(IterationNumber/OutputInterval)))
    iStartStep = 1
    do iStep = iStartStep, IterationNumber
        if(UConf%UseLocalTimeStep == 0) then
            call JPUCheckCFL4FixedTime(UG,UCC,iSplit)
            FixedTimeStep = FixedTimeStep / dble(iSplit)
        else
            iSplit = 1
        end if

        !MW%StepDistance(1:3) = - UG%GM%BC%InFlowVariable(2:4)*FixedTimeStep
        !call VelocityCorrection(1,UG,MG,UCC,MW,FixedTimeStep,99.0d0)

        do iLoop=1, iSplit
            call UnstNS(iStep,UConf,UG,UCC,UCE)
        end do

        !call VelocityCorrection(2,UG,MG,UCC,MW,FixedTimeStep,99.0d0)
        FixedTimeStep = DefaultTimeStep

        if(mod(iStep,OutputInterval) == 0) then
            iStep4Plot = iStep / OutputInterval
            call JPUOutput(UConf,UG,UCC,iStep4Plot)
            call UCalcAeroCharacteristics(UCC, UG, iStep4Plot, UAC)
            write(6,*) UAC%coefficient(1, iStep4Plot), UAC%coefficient(2, iStep4Plot)
        end if

    end do

    call JPUOutput(UConf,UG,UCC,iStep)

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

end subroutine JobParallelNS

