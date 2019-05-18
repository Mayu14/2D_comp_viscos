subroutine UOutput_Characteristics(UConf, UG, UAC, iStep)
    use StructVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(AeroCharacteristics), intent(in) :: UAC
    integer, intent(in) :: iStep
    character(len=256) :: cDirectory,cFileName, cCaseName, cUFileName, cStep
    integer :: iTime, iMaxTime, iWall, debug = 0


    if(debug == 1) then
        cUFileName = "test_NACA0012_course"
    else
        cUFileName = UConf%cFileName
    end if

    if(UConf%OutputStatus == 1) then    ! final step
        cDirectory = trim(adjustl(UConf%cDirectory))//trim(adjustl("Complete/dat/")) !UConf%SaveDirectiry    ! Œ¤‹†ŽºPC—p
        write(cStep, *) "Final"
    else
        cDirectory = trim(adjustl(UConf%cDirectory))//trim(adjustl("ResultC/")) !UConf%SaveDirectiry    ! Œ¤‹†ŽºPC—p
        write(cStep, *) iStep
        cStep = trim(adjustl(cStep))//"th"
    end if
    if(iStep == 0) then
        cFileName = trim(adjustl(cDirectory))//trim(adjustl(cUFileName))//"_"//trim(adjustl(UConf%cCaseName))//"_AC.dat"
    else
        cFileName = trim(adjustl(cDirectory))//trim(adjustl(cUFileName))//"_"//trim(adjustl(UConf%cCaseName))//trim(adjustl(cStep))//"_AC.dat"
    end if
    iMaxTime = ubound(UAC%coefficient, 2)

    open(unit = UConf%my_rank+100, file =trim(adjustl(cFileName)), status = 'unknown')
        if(iStep > 1) then
            write(UConf%my_rank+100, "(2(1x,f22.17))") 0.5d0 * (UAC%coefficient(1,iStep) + UAC%coefficient(1,iStep-1)), &
                                                     & 0.5d0 * (UAC%coefficient(2,iStep) + UAC%coefficient(2,iStep-1))
        end if
        do iTime = 1, iStep
            write(UConf%my_rank+100, "(2(1x,f22.17))") UAC%coefficient(1,iTime), UAC%coefficient(2, iTime)
        end do
    close(UConf%my_rank+100)

    if(iStep == 0) then
        cFileName = trim(adjustl(cDirectory))//"CP_"//trim(adjustl(cUFileName))//"_"//trim(adjustl(UConf%cCaseName))//".dat"
    else
        cFileName = trim(adjustl(cDirectory))//"CP_"//trim(adjustl(cUFileName))//"_"//trim(adjustl(UConf%cCaseName))//trim(adjustl(cStep))//"th.dat"
    end if

    open(unit = UConf%my_rank+100, file = trim(adjustl(cFileName)), status = 'unknown')
        !do iTime = 1, iMaxTime
        iTime = 1
            do iWall = 1, UG%GM%BC%iWallTotal
                write(UConf%my_rank+100, "(2(1x,f22.17))") UG%CD%Edge(UG%GM%BC%VW(iWall)%iGlobalEdge, 1), UAC%pressure_coefficient(iWall,iTime)
            end do
            write(UConf%my_rank+100, *) " "
        !end do
    close(UConf%my_rank+100)

    return
end subroutine UOutput_Characteristics

