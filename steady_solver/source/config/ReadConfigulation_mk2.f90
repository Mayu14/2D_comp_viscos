!***********************************/
!	Name:�v�Z�ݒ��ǂݍ��ރv���O����
!	Alias:ReadConfigulation
!	Description:��p��CalcConfig�t�@�C�����K�v
!	Type:Configulation
!	Input:CalcConfig(�O������)
!	Output:Configulation
!	Note:CalcConfig�̓v���O�����{�̂Ɠ����f�B���N�g���ɒu������
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.07
!	Other:
!***********************************/
subroutine ReadConfigulation_mk2(Conf)
use StructVar_Mod_mk2
use ConstantVar_Mod_mk2
implicit none
type(Configulation), intent(out) :: Conf
character :: cAnnotate

    open(unit=1, file='CalcConfig.dat',status='unknown')
        read(1,*) Conf%GridFileName, cAnnotate
        read(1,*) Conf%UseResume,cAnnotate
        read(1,*) Conf%ResumeFileName
        read(1,*) Conf%UseVariableTime,cAnnotate
        read(1,*) Conf%UseLocalTimeStep,cAnnotate
        read(1,*) Conf%UseMUSCL,cAnnotate
        read(1,*) Conf%KindLimiter,cAnnotate
        read(1,*) FixedTimeStep,cAnnotate
        read(1,*) OutputInterval, cAnnotate
        read(1,*) IterationNumber,cAnnotate
    close(1)

    DefaultTimeStep = FixedTimeStep
    return
end subroutine ReadConfigulation_mk2

