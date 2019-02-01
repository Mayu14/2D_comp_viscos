!***********************************/
!	Name:計算設定を読み込むプログラム
!	Alias:JobParallelControler
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
subroutine JobParallelControler(UConf, iOffset)
    use StructVar_Mod
    use ConstantVar_Mod
    use mpi
    !!$use omp_lib
    implicit none
    type(Configulation), intent(inout) :: UConf
    integer :: ierr, my_rank, PETOT, iOffset
    integer :: Total_threads

    if(UConf%UseJobParallel == 1) then
        call MPI_INIT(ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, PETOT, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
            UConf%my_rank = my_rank + iOffset
            UConf%UseJobParallel = 1
            if (PETOT > 16) then
                UConf%CalcEnv = 1
            end if
            !!$omp parallel default(private)
            ! my_rankと，omp_get_thread_num()とユーザ入力または環境変数から迎え角と格子を変更
            !Total_threads = omp_get_max_threads()
            !UConf%my_thread = omp_get_thread_num()

            call ReadConfigulation(UConf, UConf%my_rank)
            call JPCalcCaseAutoFill(UConf, PETOT)

        !!$omp end parallel

        call MPI_FINALIZE(ierr)

    else
        UConf%my_rank = iOffset
        UConf%CalcEnv = 1   ! for 流体研スパコン
        call JPCalcCaseAutoFill(UConf, PETOT)
    end if
    return
end subroutine JobParallelControler

