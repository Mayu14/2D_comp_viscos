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
use StructVar_Mod_mk2
use ConstantVar_Mod_mk2
implicit none
type(Configulation), intent(out) :: Conf
character :: cAnnotate

    open(unit=1, file='CalcConfig',status='unknown')
        read(1,*) cAnnotate
        read(1,*) Conf%GridFileName, cAnnotate
        read(1,*) Conf%UseReadRegion,cAnnotate
        read(1,*) Conf%UseResume,cAnnotate
        read(1,*) Conf%ResumeFileName
        read(1,*) Conf%ResumeFileFormat
        read(1,*) Conf%ResumeNumber
        read(1,*) Conf%UseReadBoundaryCond,cAnnotate
        read(1,*) Conf%UseVariableTime,cAnnotate
        read(1,*) Conf%UseLocalTimeStep,cAnnotate
        read(1,*) Conf%UseRRK2,cAnnotate
        read(1,*) Conf%UseMUSCL,cAnnotate
        read(1,*) Conf%KindLimiter,cAnnotate
        read(1,*) Conf%UseOverSet, cAnnotate
        read(1,*) FixedTimeStep,cAnnotate
        read(1,*) OutputInterval, cAnnotate
        read(1,*) IterationNumber,cAnnotate
        read(1,*) GridNumber, cAnnotate
    close(1)

    if(Conf%UseRRK2 == 1) IterationNumber = 2*IterationNumber
    DefaultTimeStep = FixedTimeStep

    return
end subroutine ReadConfigulation

