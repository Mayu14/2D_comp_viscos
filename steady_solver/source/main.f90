!***********************************/
!	Name:圧縮性NS方程式の定常解をRANSで求めるプログラム
!	Alias:steady_NS_rans_solver
!	Description:
!	Type:
!	Input:vtk or MAYUGrid2
!	Output:vtk
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.10.25
!	Update:-
!	Other:
!***********************************/
program steady_NS_rans_solver
    use StructVar_Mod_mk2
    use LoopVar_Mod_mk2
    use ConstantVar_Mod_mk2
    use FrequentOperation_mk2
    implicit none

    integer :: iStep,iStartStep,iStep4Plot
    double precision :: StartTime,EndTime,TotalTime
    double precision, allocatable :: WatchTime(:,:) !(1Boundary,2Flux,3Integral:1Start,2End,3Total)
    type(Configulation) :: Conf
    allocate(WatchTime(3,3))
    TotalTime = 0.0d0
    WatchTime = 0.0d0
    call ReadConfigulation_mk2(Conf)
    !iSwitch = 3
    call body_fitted_and_baldwin_lomax(Conf)

    stop
contains
    subroutine body_fitted_and_baldwin_lomax(Conf)
    implicit none
    type(Configulation), intent(in) :: Conf
    type(CellCenter) :: CC
    type(CellEdge) :: CE
    type(Geometry) :: Geom

    ! 格子情報読み込み
        ! 格子点座標読み込み
        !　セル中心，セル界面座標生成
        ! セル中心、界面におけるヤコビ行列生成
    call Set_Region(Conf, Geom, CC, CE)

    ! 初期値設定

    ! Iteration開始
        ! 境界条件入力

        ! 流束計算

        ! 時間積分

    ! 結果出力

    return
    end subroutine body_fitted_and_baldwin_lomax
end program steady_NS_rans_solver
