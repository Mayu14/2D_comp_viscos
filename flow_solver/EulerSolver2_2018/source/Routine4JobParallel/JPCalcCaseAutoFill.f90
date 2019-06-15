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
    integer :: naca4digit = 3   ! 1で4桁翼のループ，2で5桁翼のループ, 3でre_cal_dataから読み込んで再計算，4でgrid_chanteルーチンを回す

    write(cTmpDir, '("/work/A/FMa/FMa037/20_800_0010_0200/")')

    if(UConf%UseJobParallel == 1) then
        UConf%CalcEnv = 1
        call grid_change(UConf)

        iAngleDeg = int(UConf%dAttackAngle)
        CourantFriedrichsLewyCondition = CFL_default
        CheckNaNInterval = CheckNaNInterval_default
        call JobParallelNS(Uconf)

    else
        !do i1digit = 9, 1, -1
        !do i1digit = 1, 2
            !do i2digit = 9, 1, -1
                !do i34digit = 88, 12, -4
                UConf%CalcEnv = 0
                UConf%my_rank = 0
                i1digit = 0
                i2digit = 0
                i34digit = 12
                    if(UConf%CalcEnv == 0) then
                        !write(UConf%cGridName, '("NACA", i1, i1, i2.2, ".mayu")') i1digit, i2digit, i34digit ! 研究室PC用
                        write(UConf%cGridName, '("NACA0012_03_2000_0018_0100.mayu")')! valid
                        !write(UConf%cGridName, '("NACA0012_20_400_0050_0200.mayu")')! valid
                    else if(UConf%CalcEnv == 1) then
                        !write(UConf%cGridName, '("/work/A/FMa/FMa037/mayu_grid/NACA", i1, i1, i2.2, ".mayu")') i1digit, i2digit, i34digit ! 東北大スパコン用
                        write(UConf%cGridName, '("/work/A/FMa/FMa037/mayu_grid/NACA0012_20_800_0010_0200.mayu")')! 東北大スパコン用 valid
                    end if
                    !do iAngleDeg = 39, 0, -3
                        iAngleDeg = 1.25d0
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

                        !call grid_change(UConf) ! debug
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

    subroutine grid_change(UConf)
        implicit none
        type(Configulation), intent(inout) :: UConf
            UConf%UseResume = 1
            UConf%UseMUSCL = 1
            allocate(UConf%ResumeInFlowVars(5), UConf%ResumeOutFlowVars(5))
            UConf%ResumeInFlowVars = 0.0d0
            UConf%ResumeInFlowVars(1) = 1.0d0
            UConf%ResumeInFlowVars(5) = 1.0d0
            if(Uconf%my_rank == 0) then
                Uconf%cGridName = trim(adjustl("NACA0012_10_1600_0018_0100.mayu"))
                Uconf%cFileName = trim(adjustl("NACA0012_01"))
                UConf%cCaseName = trim(adjustl("case5"))
                UConf%dAttackAngle = 1.25d0
                !UConf%ResumeFileName = trim(adjustl("/mnt/d/Toyota/Data/validR/wide_grid_CN2/Resume/NACA0012_01__Resume.vtk"))
                UConf%ResumeFileName = trim(adjustl("/mnt/g/Toyota/validQ/resume/10_1600_0018_0100/Resume/NACA0012_01__Resume.vtk"))
                UConf%ResumeInFlowVars(2) = 0.80d0
                UConf%ResumeOutFlowVars = UConf%ResumeInFlowVars
                write(6,*) "resume_setting load"
                CheckNaNInterval = 1
            else if(Uconf%my_rank == 1) then
                Uconf%cGridName = trim(adjustl("/mnt/g/Toyota/Data/grid_vtk/valid/mayu/NACA21092.mayu"))
                Uconf%cFileName = trim(adjustl("NACA21092__AoA23.5__Ma0.30__Re4000__case5"))
                UConf%cCaseName = trim(adjustl("case5"))
                UConf%dAttackAngle = 23.5d0
                UConf%UseResume = 0
                UConf%ResumeInFlowVars(2) = 0.30d0
            end if
            UConf%ResumeOutFlowVars = UConf%ResumeOutFlowVars
            UConf%dAttackAngle = dPi * UConf%dAttackAngle / 180.0d0
        return
    end subroutine grid_change
end subroutine JPCalcCaseAutoFill
