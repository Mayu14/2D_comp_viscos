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
subroutine ReadConfigulation(Conf)
    use StructVar_Mod
    use ConstantVar_Mod
    implicit none
    type(Configulation), intent(out) :: Conf
    character :: cAnnotate
    logical :: debug = .true.

    open(unit=1, file='CalcConfig',status='unknown')
        read(1,*) cAnnotate
        read(1,*) Conf%SwitchProgram, cAnnotate
        read(1,*) Conf%UseReadRegion,cAnnotate
        read(1,*) Conf%UseResume,cAnnotate
        read(1,*) Conf%ResumeFileName
        read(1,*) Conf%ResumeFileFormat
        read(1,*) Conf%ResumeNumber
        read(1,*) Conf%UseReadBoundaryCond,cAnnotate
        read(1,*) Conf%UseVariableTime,cAnnotate
        read(1,*) Conf%UseSteadyCalc,cAnnotate
        read(1,*) Conf%UseLocalTimeStep,cAnnotate
        read(1,*) Conf%UseRRK2,cAnnotate
        read(1,*) Conf%UseMUSCL,cAnnotate
        read(1,*) Conf%KindLimiter,cAnnotate
        read(1,*) Conf%UseOverSet, cAnnotate
        read(1,*) FixedTimeStep,cAnnotate
        read(1,*) OutputInterval, cAnnotate
        read(1,*) IterationNumber,cAnnotate
        read(1,*) GridNumber, cAnnotate
        read(1,*) SigmoidGain, cAnnotate
        read(1,*) EntropyCorrection, cAnnotate
        read(1,*) VenkatakrishnanParameterK, cAnnotate
        read(1,*) CourantFriedrichsLewyCondition, cAnnotate
        read(1,*) Conf%TurbulenceModel, cAnnotate
        read(1,*) ReynoldsNumber, cAnnotate
    close(1)
    if(Conf%TurbulenceModel /= 0) invicid = .false.
    if(Conf%UseRRK2 == 1) IterationNumber = 2*IterationNumber
    DefaultTimeStep = FixedTimeStep

    call Show_Configulation(debug)

    return

contains
    subroutine Show_Configulation(debug)
        logical, intent(in) :: debug

        if (debug == .false.) then
            ! no action
        else
            write(6,*) "show_configulation == TRUE"
            write(6,*) "Conf%SwitchProgram", Conf%SwitchProgram
            write(6,*) "Conf%UseReadRegion", Conf%UseReadRegion
            write(6,*) "Conf%UseResume", Conf%UseResume
            write(6,*) "Conf%ResumeFileName", Conf%ResumeFileName
            write(6,*) "Conf%ResumeFileFormat", Conf%ResumeFileFormat
            write(6,*) "Conf%ResumeNumber", Conf%ResumeNumber
            write(6,*) "Conf%UseReadBoundaryCond", Conf%UseReadBoundaryCond
            write(6,*) "Conf%UseVariableTime", Conf%UseVariableTime
            write(6,*) "Conf%UseSteadyCalc", Conf%UseSteadyCalc
            write(6,*) "Conf%UseLocalTimeStep", Conf%UseLocalTimeStep
            write(6,*) "Conf%UseRRK2", Conf%UseRRK2
            write(6,*) "Conf%UseMUSCL", Conf%UseMUSCL
            write(6,*) "Conf%KindLimiter", Conf%KindLimiter
            write(6,*) "Conf%UseOverSet", Conf%UseOverSet
            write(6,*) "FixedTimeStep", FixedTimeStep
            write(6,*) "OutputInterval", OutputInterval
            write(6,*) "IterationNumber", IterationNumber
            write(6,*) "GridNumber", GridNumber
            write(6,*) "SigmoidGain", SigmoidGain
            write(6,*) "EntropyCorrection", EntropyCorrection
            write(6,*) "VenkatakrishnanParameterK", VenkatakrishnanParameterK
            write(6,*) "CourantFriedrichsLewyCondition", CourantFriedrichsLewyCondition
            write(6,*) "Conf%TurbulenceModel", Conf%TurbulenceModel
            write(6,*) "ReynoldsNumber", ReynoldsNumber
        end if

        return
    end subroutine Show_Configulation

end subroutine ReadConfigulation

