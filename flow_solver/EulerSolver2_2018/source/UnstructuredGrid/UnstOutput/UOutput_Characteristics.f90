subroutine UOutput_Characteristics(UConf, UAC)
    use StructVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(AeroCharacteristics), intent(in) :: UAC
    character(len=64) :: cDirectory,cFileName, cCaseName
    integer :: iTime, iMaxTime

    !cDirectory = "ResultC/" !UConf%SaveDirectiry
    cDirectory = "../../../Case1/ResultC/" ! 東北大スパコン用
    cFileName = trim(adjustl(cDirectory))//trim(adjustl(UConf%cFileName))//"_AC.dat"
    iMaxTime = ubound(UAC%coefficient, 2)

    open(unit = UConf%my_rank+100, file =trim(adjustl(cFileName)), status = 'unknown')
        do iTime = 1, iMaxTime
            write(UConf%my_rank+100, *) UAC%coefficient(1,iTime), UAC%coefficient(2, iTime)
        end do
    close(UConf%my_rank+100)

    return
end subroutine UOutput_Characteristics

