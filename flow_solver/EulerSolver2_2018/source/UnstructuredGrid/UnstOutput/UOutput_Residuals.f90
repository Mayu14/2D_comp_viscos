subroutine UOutput_Residuals(UConf, RH, iStep)
    use StructVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(ResidualHistory), intent(inout) :: RH
    integer, intent(in) :: iStep
    character(len=256) :: cDirectory,cFileName, cCaseName, cUFileName, cStep
    integer :: iTime, iMaxTime, iWall, debug = 0

    if(debug == 1) then
        cUFileName = "test_NACA0012_course"
    else
        cUFileName = UConf%cFileName
    end if

    write(cStep, *) iStep
    cDirectory = trim(adjustl(UConf%cDirectory))//trim(adjustl("ResultR/")) !UConf%SaveDirectiry    ! Œ¤‹†ŽºPC—p

    cFileName = trim(adjustl(cDirectory))//"RES_"//trim(adjustl(cUFileName))//"_"//trim(adjustl(UConf%cCaseName))//trim(adjustl(cStep))//"th.csv"

    if(RH%iLastOutput /= 1) then
        open(unit = UConf%my_rank+100, file =trim(adjustl(cFileName)), status = 'unknown', position ="append")
    else
        open(unit = UConf%my_rank+100, file =trim(adjustl(cFileName)), status = 'unknown')
    end if

        do iTime = RH%iLastOutput, RH%iTime - 1
            write(UConf%my_rank+100, "(i7, 12(a1,f22.17))") iTime,",",RH%AveResidual(1,iTime),",",RH%AveResidual(2,iTime),",",RH%AveResidual(3,iTime),",",    &
                                                                    & RH%AveResidual(4,iTime),",",RH%AveResidual(5,iTime),",",RH%AveResidual(6,iTime),",",    &
                                                                    & RH%MaxResidual(1,iTime),",",RH%MaxResidual(2,iTime),",",RH%MaxResidual(3,iTime),",",    &
                                                                    & RH%MaxResidual(4,iTime),",",RH%MaxResidual(5,iTime),",",RH%MaxResidual(6,iTime)
        end do
    close(UConf%my_rank+100)
    RH%iLastOutput = RH%iTime

    return
end subroutine UOutput_Residuals

