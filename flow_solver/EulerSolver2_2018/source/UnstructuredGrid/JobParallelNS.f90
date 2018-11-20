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
    call UReadInflowOutflowCondition(UG)
    call UOutput(UConf,UG,UCC,iStep)

    allocate(UAC%coefficient(2, int(IterationNumber/OutputInterval)))
    iStartStep = 1
    do iStep = iStartStep, IterationNumber
write(6,*) UConf%my_rank, "a"
        if(UConf%UseLocalTimeStep == 0) then
            call JPUCheckCFL4FixedTime(UG,UCC,iSplit)
            FixedTimeStep = FixedTimeStep / dble(iSplit)
        else
            iSplit = 1
        end if
write(6,*) UConf%my_rank, "b"
        do iLoop=1, iSplit
            write(6,*) iLoop, iSplit
            call UnstNS(iStep,UConf,UG,UCC,UCE)
        end do
write(6,*) UConf%my_rank, "c"
        FixedTimeStep = DefaultTimeStep

        if(mod(iStep,OutputInterval) == 0) then
            iStep4Plot = iStep / OutputInterval
            call UOutput(UConf,UG,UCC,iStep4Plot)
            call UCalcAeroCharacteristics(UCC, UG, iStep4Plot, UAC)
            write(6,*) UAC%coefficient(1, iStep4Plot), UAC%coefficient(2, iStep4Plot)
        end if
write(6,*) UConf%my_rank, "d"
    end do
write(6,*) UConf%my_rank, "e"
    call UOutput(UConf,UG,UCC,iStep)

    return
contains

    subroutine UnstNS(iStep, OConf, UG, UCC, UCE)
        implicit none
        integer, intent(in) :: iStep
        type(Configulation), intent(in) :: OConf
        type(UnstructuredGrid), intent(in) :: UG
        type(CellCenter), intent(inout) :: UCC
        type(CellEdge), intent(inout) :: UCE
write(6,*) UConf%my_rank, "aa"
            call USetBoundary(UG,UCC)
write(6,*) UConf%my_rank, "ab"
            call UCalcFlux(OConf,UG,UCC,UCE)
write(6,*) UConf%my_rank, "ac"
            call UTimeIntegral(OConf,UG,UCE,UCC,iStep)
write(6,*) UConf%my_rank, "ad"
            call UOutput(OConf, UG, UCC, iStep)

        return
    end subroutine UnstNS

end subroutine JobParallelNS

