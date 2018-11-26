subroutine UOutput_Characteristics(UConf, UAC)
    use StructVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(AeroCharacteristics), intent(in) :: UAC
    character(len=64) :: cDirectory,cFileName, cCaseName
    integer :: iTime, iMaxTime

    !cDirectory = "ResultC/" !UConf%SaveDirectiry
    cDirectory = "../../../Case1/ResultC/" ! 東北大スパコン用
    cFileName = trim(adjustl(cDirectory))//"AC.dat"
    iMaxTime = ubound(UAC%coefficient, 2)

    open(unit = 1, file =trim(adjustl(cFileName)), status = 'unknown')
        do iTime = 1, iMaxTime
            write(1, *) UAC%coefficient(1,iTime), UAC%coefficient(2, iTime)
        end do
    close(1)

    return
end subroutine UOutput_Characteristics

