subroutine UOutput_Characteristics(UConf, UG, UAC)
    use StructVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(AeroCharacteristics), intent(in) :: UAC
    character(len=256) :: cDirectory,cFileName, cCaseName, cUFileName
    integer :: iTime, iMaxTime, iWall, debug = 0

    if(debug == 1) then
        cUFileName = "test_NACA0012_course"
    else
        cUFileName = UConf%cFileName
    end if

    cDirectory = trim(adjustl(UConf%cDirectory))//trim(adjustl("ResultC/")) !UConf%SaveDirectiry    ! Œ¤‹†ŽºPC—p
    cFileName = trim(adjustl(cDirectory))//trim(adjustl(cUFileName))//"_AC.dat"
    iMaxTime = ubound(UAC%coefficient, 2)

    open(unit = UConf%my_rank+100, file =trim(adjustl(cFileName)), status = 'unknown')
        do iTime = 1, iMaxTime
            write(UConf%my_rank+100, "(2(1x,f22.17))") UAC%coefficient(1,iTime), UAC%coefficient(2, iTime)
        end do
    close(UConf%my_rank+100)

    cFileName = trim(adjustl(cDirectory))//"CP_"//trim(adjustl(cUFileName))//".dat"

    open(unit = UConf%my_rank+100, file = trim(adjustl(cFileName)), status = 'unknown')
        !do iTime = 1, iMaxTime
        iTime = iMaxTime
            do iWall = 1, UG%GM%BC%iWallTotal
                write(UConf%my_rank+100, "(2(1x,f22.17))") UG%CD%Edge(UG%GM%BC%VW(iWall)%iGlobalEdge, 1), UAC%pressure_coefficient(iWall,iTime)
            end do
            write(UConf%my_rank+100, *) " "
        !end do
    close(UConf%my_rank+100)

    return
end subroutine UOutput_Characteristics

