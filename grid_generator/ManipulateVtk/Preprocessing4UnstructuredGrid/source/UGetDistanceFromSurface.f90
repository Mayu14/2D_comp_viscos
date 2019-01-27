!***********************************/
!	Name:各セルに対し最近傍物体表面を求め，そこからの垂直距離を求めるプログラム
!	Alias:UGetDistanceFromSurface
!	Description:
!	Type:UG
!	Input:
!	Output:
!	Note:効率を度外視してO(N^2)で書いたため，ボトルネックになるようなら修正するかも
!	Author:Akitaka Toyota
!	Date:2018.11.06
!	Update:2018.11.13
!	Other:
!***********************************/
subroutine UGetDistanceFromSurface(UG, ExistInnerBoundary)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG
    logical, intent(out)  :: ExistInnerBoundary

    integer :: iSurfaceCell ! 物体表面セルの総数 = 物体表面界面の数
    integer, allocatable :: iSurfaceEdge(:) ! surface edge number → edge number
    double precision :: tmp_de_num
    double precision, allocatable :: dDistance(:, :)    ! iCell, 1:所属界面，2:距離

    integer :: iTmpWall, iMember, iMem1, iMem2, iTermination ! 仮置き壁番号, メンバセル用ループ変数, メンバセル1, メンバセル2, バブルソートの終端
    integer, allocatable :: iCount(:)   ! 各壁ごとの入力済みカウント用配列
    integer :: iFlag = 0    ! 並び替え終了フラグ

    iSurfaceCell = (UG%GI%AllCells - UG%GI%RealCells) - UG%GI%OutlineCells  ! 物体表面界面数の特定
    UG%GM%BC%iWallTotal = iSurfaceCell

    if(iSurfaceCell == 0) then
        write(6,*) "No Internal Boundary"
        ExistInnerBoundary = .false.
    else
        ExistInnerBoundary = .true.
        allocate(iSurfaceEdge(iSurfaceCell))    ! 物体表面の局所界面番号→大域界面番号
        allocate(UG%GM%BC%VW(iSurfaceCell))     ! 粘性壁登録変数の動的確保

        iEdge = 1   ! 局所界面番号の初期化
        do iCell = UG%GI%RealCells+1, UG%GI%AllCells    ! すべての仮想セルについて巡回
            if (UG%VC%Type(iCell) == 2) then    ! 物体表面の仮想セルについて
                iSurfaceEdge(iEdge) = UG%VC%Edge(iCell) ! 大域辺番号をiSurfaceEdgeに格納
                UG%GM%BC%VW(iEdge)%iGlobalEdge = UG%VC%Edge(iCell)  ! 1行上のiSurfaceEdgeと同じやつ
                iEdge = iEdge + 1
            end if
        end do
        ! この時点で物体表面の辺を順不同に全格納完了

        ! すべてのセルについて、所属表面および所属表面からの距離を求める
        allocate(dDistance(UG%GI%RealCells, 2))
        ! UG%GM%BC%VW内部のメンバ総数の初期化
        call viscosity_wall_initialize(UG%GM%BC, iSurfaceCell)

        do iCell = 1, UG%GI%RealCells   ! すべてのセルについて
            call get_minimum_distance_bruteforce(UG, iCell, iSurfaceEdge, tmp_de_num, dDistance(iCell, 2))    ! 直線距離が最短となる界面の局所界面番号tmp_de_numと，垂直方向距離dDistance(iCell, 2)を格納
            dDistance(iCell, 1) = int(tmp_de_num)    ! 局所界面番号であることに注意されたし
            UG%GM%BC%VW(int(tmp_de_num))%iNumberOfMember = UG%GM%BC%VW(int(tmp_de_num))%iNumberOfMember + 1 ! 局所界面に所属するセルの総数を+1
        end do
        UG%Tri%Belongs2Wall = int(dDistance(:, 1))  ! 三角形要素番号→所属する壁の局所界面番号
        UG%Tri%Distance = dDistance(:, 2)   ! 三角形要素番号→所属する壁までの距離
        ! ここまでで各界面に所属するセルの総数が得られたため，壁→セルの検索用配列を用意する
        call viscosity_wall_malloc(UG%GM%BC, iSurfaceCell)
        ! 壁にセルを仮登録(順不同)
        allocate(iCount(iSurfaceCell))  ! 入力済み確認配列
        iCount = 1
        do iCell = 1, UG%GI%RealCells   ! すべてのセルについて
            iTmpWall = UG%Tri%Belongs2Wall(iCell)   ! 所属壁番号を仮置きして
            UG%GM%BC%VW(iTmpWall)%iMemberCell(iCount(iTmpWall)) = iCell ! 壁番号→所属要素
            iCount(iTmpWall) = iCount(iTmpWall) + 1 ! 入力した回数を反映する
        end do

        ! 距離準並び替えを実施(もうめんどいしバブルソートでええやろ(適当))(そもそも入力がほぼ整列済みである可能性が高いから...(震え声))
        do iTmpWall = 1, iSurfaceCell   ! すべての界面について
            iFlag = 0
            do while(iFlag == 0)
                iFlag = 1   ! 終了フラグを点灯
                iTermination = 1
                do iMember = 1, UG%GM%BC%VW(iTmpWall)%iNumberOfMember - iTermination
                    iMem1 = UG%GM%BC%VW(iTmpWall)%iMemberCell(iMember)  ! iMember番目のメンバセルのセル番号
                    iMem2 = UG%GM%BC%VW(iTmpWall)%iMemberCell(iMember + 1)  ! iMember + 1番目のメンバセルのセル番号
                    if (UG%Tri%Distance(iMem1) > UG%Tri%Distance(iMem2)) then   ! iMem1の距離の方がiMem2より遠いならば
                        iFlag = 0   ! 終了フラグを消灯し
                        call SwapInteger(UG%GM%BC%VW(iTmpWall)%iMemberCell(iMember), UG%GM%BC%VW(iTmpWall)%iMemberCell(iMember + 1))    ! メンバセルの登録順序を変更
                    end if
                end do
                iTermination = iTermination + 1 ! ソート済み領域をソート対象から外す
            end do
        end do
    end if

    return
contains

    double precision function get_distance_e2c(edge, cell) result(distance)
        double precision, intent(in) :: edge(:), cell(:)
        distance = sqrt(dot_product(edge - cell, edge - cell))
    end function get_distance_e2c

    subroutine get_edge_id(edge_num, em1, e00, ep1)
        integer, intent(in) :: edge_num
        integer, intent(out) :: em1, e00, ep1

        e00 = edge_num

        if(edge_num == 1) then
            em1 = iSurfaceCell
            ep1 = edge_num + 1
        else if (edge_num == iSurfaceCell) then
            em1 = edge_num - 1
            ep1 = 1
        else
            em1 = edge_num - 1
            ep1 = edge_num + 1
        end if
        return
    end subroutine get_edge_id

    subroutine get_minimum_distance_bruteforce(UG, iCellNum, iSurfaceEdge, de_num, minDistance)
        type(UnstructuredGrid), intent(in) :: UG
        integer, intent(in) :: iCellNum
        integer, intent(in) :: iSurfaceEdge(:)
        double precision, intent(out) :: de_num, minDistance
        double precision :: tmpDistance
        minDistance = 100000
        do iEdge = 1, iSurfaceCell
            tmpDistance = get_distance_e2c(UG%CD%Edge(iSurfaceEdge(iEdge), :), UG%CD%Cell(iCellNum, :))
            if(tmpDistance < minDistance) then
                minDistance = tmpDistance
                de_num = dble(iEdge)
            end if
        end do

        minDistance = AbsCrossProduct(UG%CD%Cell(iCellNum, :) - UG%CD%Edge(iSurfaceEdge(int(de_num)), :),&
                                    & UG%CD%Point(UG%Line%Point(iSurfaceEdge(int(de_num)), 1), :) - UG%CD%Point(UG%Line%Point(iSurfaceEdge(int(de_num)), 2), :)) &
                                    & / UG%GM%Area(iSurfaceEdge(int(de_num)))    ! 物体表面の辺と，界面中心からセル中心を結ぶ辺とがなす平行四辺形の高さを求める

        return
    end subroutine get_minimum_distance_bruteforce

    ! メンバ総数の初期化
    subroutine viscosity_wall_initialize(BC, iWallNumber)
        type(BoundaryCondition), intent(inout) :: BC
        integer, intent(in) :: iWallNumber

        do iLoop = 1, iWallNumber
            BC%VW(iLoop)%iNumberOfMember = 0
        end do

        return
    end subroutine viscosity_wall_initialize

    ! メンバセル登録用配列の動的確保
    subroutine viscosity_wall_malloc(BC, iWallNumber)
        type(BoundaryCondition), intent(inout) :: BC
        integer, intent(in) :: iWallNumber

        do iLoop = 1, iWallNumber
            allocate(BC%VW(iLoop)%iMemberCell(BC%VW(iLoop)%iNumberOfMember))
        end do

        return
    end subroutine viscosity_wall_malloc

    subroutine check_sort
        do iTmpWall = 1, iSurfaceCell
            write(6, *) UG%GM%BC%VW(iTmpWall)%iNumberOfMember
            do iMember = 1, UG%GM%BC%VW(iTmpWall)%iNumberOfMember
                write(6, *) UG%GM%BC%VW(iTmpWall)%iMemberCell(iMember), UG%Tri%Distance(UG%GM%BC%VW(iTmpWall)%iMemberCell(iMember))
            end do
            write(6,*) ""
            write(6,*) ""
            write(6,*) ""
        end do
    return
    end subroutine check_sort

end subroutine UGetDistanceFromSurface
