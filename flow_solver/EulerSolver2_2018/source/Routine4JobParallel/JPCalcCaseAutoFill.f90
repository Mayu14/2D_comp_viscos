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
subroutine JPCalcCaseAutoFill(UConf, PETOT)
    use StructVar_Mod
    use ConstantVar_Mod
    use mpi
    use omp_lib
    implicit none
    type(Configulation), intent(inout) :: UConf
    integer, intent(in) :: PETOT
    integer :: i34digit, i2digit, i1digit, iAngleDeg, iLoop, i12digit
    double precision :: AttackAngleRad
    character(len=32) :: cGridName, cResultName, cLoop, cAngle

    !PETET = 0 ~ 1799を仮定
    ! do i1digit = 0, 9
        ! do i2digit = 0, 9
            ! do i34digit = 1, 40
            i12digit = mod(UConf%my_rank, 100)  ! 00~99
            i34digit = 5 * int((UConf%my_rank / 100)) + 11  ! 11~96, 5k+1
                write(UConf%cGridName, '("../../../../grid/NACA", i2.2, i2.2, ".mayu")') i12digit, i34digit ! 東北大スパコン用
                ! write(UConf%cGridName, '("NACA", i1, i1, i2.2, ".mayu")') i1digit, i2digit, i34digit
                do iAngleDeg = 0, 39, 3
                    UConf%dAttackAngle = dPi * dble(iAngleDeg) / 180.0d0
                    write(UConf%cFileName, '("NACA", i2.2, i2.2,  "_", i2.2)') i12digit, i34digit, iAngleDeg ! 東北大スパコン用
                    !write(UConf%cFileName, '("NACA", i1, i1, i2.2,  "_", i2.2)') i1digit, i2digit, i34digit, iAngleDeg
                    call JobParallelNS(Uconf)
                end do
            ! end do
        ! end do
    ! end do

    return
contains
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

