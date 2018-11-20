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
subroutine JPCalcCaseAutoFill(UConf, my_rank, PETOT, my_thread, Total_threads)
    use StructVar_Mod
    use ConstantVar_Mod
    use mpi
    use omp_lib
    implicit none
    type(Configulation), intent(inout) :: UConf
    integer, intent(in) :: my_rank, PETOT, my_thread, Total_threads
    integer :: iParallelTotal
    integer :: i4digit, iAngleDeg
    double precision :: AttackAngleRad
    character(len=32) :: cGridName, cResultName

    iParallelTotal = int(PETOT * Total_threads)
    write(6,*) iParallelTotal
    do i4digit = 1, 9999, 2
        write(cGridName, '("NACA", i4.4, ".mayu")') i4digit
        do iAngleDeg = 0, 39, 3
            write(cResultName, '("NACA", i4.4, "_", i2.2, ".mayu")') i4digit, iAngleDeg
            !access()
        end do
    end do

    return
end subroutine JPCalcCaseAutoFill

