!***********************************/
!	Name:計算設定を読み込むプログラム
!	Alias:JPCalcCaseAutoFill
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
subroutine JPCalcCaseAutoFill(UConf, PETOT)
    use StructVar_Mod
    use ConstantVar_Mod
    use mpi
    use omp_lib
    implicit none
    type(Configulation), intent(inout) :: UConf
    integer, intent(in) :: PETOT
    integer :: i34digit, i2digit, i1digit, iAngleDeg, iLoop, i12digit
    double precision :: AttackAngleRad, dAngleDeg
    character(len=256) :: cGridName, cResultName, cLoop, cAngle
    integer :: debug = 0
    integer :: access
    character(len=256) :: cDirectory,cFileName, cCaseName, cTmpDir
    character(len=256) :: cStep
    integer :: naca4digit = 3

    write(cTmpDir, '("/work/A/FMa/FMa037/20_800_0010_0200/")')

    if(UConf%UseJobParallel == 1) then
    !PETET = 0 ~ 1619を仮定    ! NACA5の場合1~559
    ! do i1digit = 0, 9
        ! do i2digit = 0, 9
            ! do i34digit = 1, 40
            if(naca4digit == 1) then
                i12digit = int(float(UConf%my_rank) / 20.0d0) + 11 + int(float(UConf%my_rank) / 180.0d0) ! 11~19, 21~29, 31~..., 91~99
                i34digit = 4 * mod(UConf%my_rank, 20) + 12  ! 12~88, 4k+12
            else if(naca4digit == 2) then
                i12digit = 10 * (21 + (4 - mod(int(float(UConf%my_rank)/7.0d0), 5))) + int(float(UConf%my_rank) / 5.0d0)
                i34digit = mod(UConf%my_rank, 80) + 11
            end if

            if(UConf%CalcEnv == 0) then
                if(naca4digit == 1) then
                    write(UConf%cGridName, '("NACA", i2.2, i2.2, ".mayu")') i12digit, i34digit ! 研究室PC用
                else
                    write(UConf%cGridName, '("NACA", i3.3, i2.2, ".mayu")') i12digit, i34digit ! 研究室PC用
                end if
            else if(UConf%CalcEnv == 1) then
                if(naca4digit == 1) then
                    write(UConf%cGridName, '("/work/A/FMa/FMa037/mayu_grid/NACA", i2.2, i2.2, ".mayu")') i12digit, i34digit ! 東北大スパコン用
                else
                    write(UConf%cGridName, '("/work/A/FMa/FMa037/mayu_grid/NACA", i3.3, i2.2, ".mayu")') i12digit, i34digit ! 東北大スパコン用
                end if
            end if
            !write(6,*) UConf%my_rank ,UConf%cGridName, i34digit
                ! write(UConf%cGridName, '("NACA", i1, i1, i2.2, ".mayu")') i1digit, i2digit, i34digit
                !do iAngleDeg = 39, 0, -3
                iAngleDeg = 0
                    UConf%dAttackAngle = dPi * dble(iAngleDeg) / 180.0d0
                    if(naca4digit == 1) then
                        write(UConf%cFileName, '("NACA", i2.2, i2.2,  "_", i2.2)') i12digit, i34digit, iAngleDeg
                    else
                        write(UConf%cFileName, '("NACA", i3.3, i2.2,  "_", i2.2)') i12digit, i34digit, iAngleDeg
                    end if

                    if (naca4digit == 3) then
                        call read_re_cal_file(UConf, iOffset)
                    end if

                    if(UConf%CalcEnv == 0) then
                        write(UConf%cDirectory, '("")')
                    else if(UConf%CalcEnv == 1) then
                        UConf%cDirectory = cTmpDir
                    end if

                    if(debug == 1) then
                        write(UConf%cGridName, '("NACA0012_course.mayu")')
                        write(UConf%cFileName, '("NACA0012_course_mk2_t5000k", i2.2)') iAngleDeg
                    else if(debug == 2) then
                        write(UConf%cGridName, '("NACA0012_medium.mayu")')
                        write(UConf%cFileName, '("NACA0012_medium_", i2.2)') iAngleDeg
                    else if(debug == 3) then
                        write(UConf%cGridName, '("circle_HD.mayu")')
                        if(UConf%UseFluxMethod == 0) then
                            write(UConf%cFileName, '("circle_Roe", i2.2)') iAngleDeg
                        else if(UConf%UseFluxMethod == 1) then
                            write(UConf%cFileName, '("circle_SLAU2", i2.2)') iAngleDeg
                        end if
                    end if

                    CourantFriedrichsLewyCondition = CFL_default
                    CheckNaNInterval = CheckNaNInterval_default
                    call JobParallelNS(Uconf)
                !end do
            ! end do
        ! end do
    ! end do
    else
        !do i1digit = 9, 1, -1
        !do i1digit = 1, 2
            !do i2digit = 9, 1, -1
                !do i34digit = 88, 12, -4
                UConf%CalcEnv = 0
                i1digit = 0
                i2digit = 0
                i34digit = 12
                    if(UConf%CalcEnv == 0) then
                        !write(UConf%cGridName, '("NACA", i1, i1, i2.2, ".mayu")') i1digit, i2digit, i34digit ! 研究室PC用
                        write(UConf%cGridName, '("NACA0012_10_200_0100_0200.mayu")')! valid
                    else if(UConf%CalcEnv == 1) then
                        !write(UConf%cGridName, '("/work/A/FMa/FMa037/mayu_grid/NACA", i1, i1, i2.2, ".mayu")') i1digit, i2digit, i34digit ! 東北大スパコン用
                        write(UConf%cGridName, '("/work/A/FMa/FMa037/mayu_grid/NACA0012_20_800_0010_0200.mayu")')! 東北大スパコン用 valid
                    end if
                    !do iAngleDeg = 39, 0, -3
                        iAngleDeg = 1.25
                        dAngleDeg = 1.25d0

                        !UConf%dAttackAngle = dPi * dble(iAngleDeg) / 180.0d0
                        UConf%dAttackAngle = dPi * dAngleDeg / 180.0d0

                        if(UConf%CalcEnv == 0) then
                            write(UConf%cFileName, '("NACA", i1, i1, i2.2,  "_", i2.2)') i1digit, i2digit, i34digit, iAngleDeg
                            write(UConf%cDirectory, '("")')
                        else if(UConf%CalcEnv == 1) then
                            write(UConf%cFileName, '("NACA", i1, i1, i2.2,  "_", i2.2)') i1digit, i2digit, i34digit, iAngleDeg ! 東北大スパコン用
                            UConf%cDirectory = cTmpDir
                        end if

                        if(debug == 1) then
                            write(UConf%cGridName, '("NACA0012_course_rev2.mayu")')
                            write(UConf%cFileName, '("NACA0012_course_rev2_", i2.2)') iAngleDeg
                        else if(debug == 2) then
                            write(UConf%cGridName, '("NACA0012_fine_rev2.mayu")')
                            write(UConf%cFileName, '("NACA0012_fine_rev2_", i2.2)') iAngleDeg
                        else if(debug == 4) then
                            write(UConf%cGridName, '("mirror_Square_Half_fine_rev2.mayu")')
                            if(UConf%UseFluxMethod == 0) then
                                write(UConf%cFileName, '("ST_rev2_Roe")')
                            else if(UConf%UseFluxMethod == 1) then
                                write(UConf%cFileName, '("ST_rev2_SLAU2")')
                            end if

                        end if

                        CourantFriedrichsLewyCondition = CFL_default
                        CheckNaNInterval = CheckNaNInterval_default
                        ! 出力先ファイルがないときのみ実行
                        write(cStep,*) 0
                        if(UConf%CalcEnv == 0) then
                            cFileName = trim(adjustl("ResultU/"))//trim(adjustl(UConf%cFileName))//trim(adjustl(UConf%cCaseName))//trim(adjustl("_"))//trim(adjustl(cStep))//"th.vtk"
                        else if(UConf%CalcEnv == 1) then    ! 東北大スパコン用
                            cFileName = trim(adjustl(UConf%cDirectory))//trim(adjustl("ResultU/"))//trim(adjustl(UConf%cFileName))//trim(adjustl(UConf%cCaseName))//trim(adjustl("_"))//trim(adjustl(cStep))//"th.vtk"
                        end if

                        if(UConf%SwitchProgram /= 7) then
                            if(access(cFileName, " ") /= 0) then
                                write(6,*) trim(adjustl(UConf%cFileName))//"_"//trim(adjustl(UConf%cCaseName))
                                call JobParallelNS(Uconf)
                            end if
                        else
                            call JobParallelNS(UConf)
                        end if
                    !end do
                !end do
            !end do
        !end do

    end if

    return
contains
    subroutine read_re_cal_file(UConf, iOffset)
        implicit none
        type(Configulation), intent(inout) :: UConf
        integer, intent(in) :: iOffset
        integer :: iUnit, iDeg
        character(len=256) :: cDigit, cDeg, cMid

            iUnit = UConf%my_rank+100
            open(unit = iUnit, file="re_cal_namelist.dat", status = "unknown")
            do iLoop = 1 + iOffset, 2 + UConf%my_rank + iOffset
                read(iUnit, *) cDigit, cDeg, cMid
            end do

            UConf%cFileName = "NACA"//trim(adjustl(cDigit))//"_"//trim(adjustl(cDeg))
            UConf%cGridName = "NACA"//trim(adjustl(cDigit))//".mayu"
            if(UConf%CalcEnv == 1) then
                UConf%cGridName = "/work/A/FMa/FMa037/mayu_grid/"//trim(adjustl(UConf%cGridName))
            end if

            read(cDeg, *) iDeg
            UConf%dAttackAngle = dPi * dble(iDeg) / 180.0d0
            if(cMid /= "None") then
                UConf%UseResume = 1
                UConf%ResumeFileName = trim(adjustl(cMid))
            end if

        return
    end subroutine read_re_cal_file

    subroutine grid_change(Uconf)
        implicit none
        type(Configulation), intent(inout) :: UConf

            if(UConf%my_rank == 0) then
                UConf%cGridName = "NACA0012"
                UConf%dAttackAngle = 0.0d0 / 180.0d0 * dPi
            else if(UConf%my_rank == 1) then
                UConf%cGridName = "NACA0012"
                UConf%dAttackAngle = 10.0d0 / 180.0d0 * dPi
            else if(UConf%my_rank == 2) then
                UConf%cGridName = "NACA0223"
            else if(UConf%my_rank == 3) then
                UConf%cGridName = "NACA0115"
            else if(UConf%my_rank == 4) then
                UConf%cGridName = "NACA0117"
            else if(UConf%my_rank == 5) then
                UConf%cGridName = "NACA0119"
            else if(UConf%my_rank == 6) then
                UConf%cGridName = "NACA1111"
            else if(UConf%my_rank == 7) then
                UConf%cGridName = "NACA1113"
            else if(UConf%my_rank == 8) then
                UConf%cGridName = "NACA1115"
            else if(UConf%my_rank == 9) then
                UConf%cGridName = "NACA1117"
            else if(UConf%my_rank == 10) then
                UConf%cGridName = "NACA1119"
            else if(UConf%my_rank == 11) then
                UConf%cGridName = "NACA1517"
            else if(UConf%my_rank == 12) then
                UConf%cGridName = "NACA1519"
            else if(UConf%my_rank == 13) then
                UConf%cGridName = "NACA1611"
            else if(UConf%my_rank == 14) then
                UConf%cGridName = "NACA2013"
            else if(UConf%my_rank == 15) then
                UConf%cGridName = "NACA2015"
            else if(UConf%my_rank == 16) then
                UConf%cGridName = "NACA2017"
            end if

        return
    end subroutine grid_change

end subroutine JPCalcCaseAutoFill

