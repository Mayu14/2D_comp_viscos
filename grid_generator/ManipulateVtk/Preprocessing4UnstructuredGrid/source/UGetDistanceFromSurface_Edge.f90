!***********************************/
!	Name:各辺に対し最近傍物体表面を求め，そこからの垂直距離を求めるプログラム
!	Alias:UGetDistanceFromSurface_Edge
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
subroutine UGetDistanceFromSurface_Edge(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    integer :: iSurfaceCell ! 物体表面セルの総数 = 物体表面界面の数
    integer, allocatable :: iSurfaceEdge(:) ! surface edge number → edge number
    double precision :: tmp_de_num
    double precision, allocatable :: dDistance(:, :)    ! iCell, 1:所属界面，2:距離

    integer :: iTmpWall, iMember, iMem1, iMem2, iTermination ! 仮置き壁番号, メンバセル用ループ変数, メンバセル1, メンバセル2, バブルソートの終端
    integer, allocatable :: iCount(:)   ! 各壁ごとの入力済みカウント用配列
    integer :: iFlag = 0    ! 並び替え終了フラグ

    iSurfaceCell = UG%GM%BC%iWallTotal  ! 物体表面界面数の特定

    ! 物体表面の辺は格納済みとする

    ! すべての界面について、所属表面および所属表面からの距離を求める
    allocate(dDistance(UG%GI%Edges, 2))
    ! UG%GM%BC%VW内部のメンバ総数の初期化
    call viscosity_wall_initialize(UG%GM%BC, iSurfaceCell)

    do iEdge = 1, UG%GI%Edges   ! すべての界面について
        call get_minimum_distance_bruteforce_edge(iEdge, tmp_de_num, dDistance(iEdge, 2))    ! 直線距離が最短となる界面の局所界面番号tmp_de_numと，垂直方向距離dDistance(iCell, 2)を格納
        dDistance(iEdge, 1) = int(tmp_de_num)    ! 局所界面番号であることに注意されたし
        UG%GM%BC%VW(int(tmp_de_num))%iNumberOfMemberEdge = UG%GM%BC%VW(int(tmp_de_num))%iNumberOfMemberEdge + 1 ! 局所界面に所属するセルの総数を+1
    end do

    UG%Line%Belongs2Wall = int(dDistance(:, 1))  ! 界面番号→所属する壁の局所界面番号
    UG%Line%Distance = dDistance(:, 2)   ! 界面番号→所属する壁までの距離
    ! ここまでで各界面に所属する界面の総数が得られたため，壁→セルの検索用配列を用意する
    call viscosity_wall_malloc(UG%GM%BC, iSurfaceCell)
    ! 壁に界面を仮登録(順不同)
    allocate(iCount(iSurfaceCell))  ! 入力済み確認配列
    iCount = 1
    do iEdge = 1, UG%GI%Edges   ! すべての界面について
        iTmpWall = UG%Line%Belongs2Wall(iEdge)   ! 所属壁番号を仮置きして
        UG%GM%BC%VW(iTmpWall)%iMemberEdge(iCount(iTmpWall)) = iEdge ! 壁番号→所属している界面の番号
        iCount(iTmpWall) = iCount(iTmpWall) + 1 ! 入力した回数を反映する
    end do

    ! 距離準並び替えを実施(もうめんどいしバブルソートでええやろ(適当))(そもそも入力がほぼ整列済みである可能性が高いから...(震え声))
    do iTmpWall = 1, iSurfaceCell   ! すべての界面について
        iFlag = 0
        do while(iFlag == 0)
            iFlag = 1   ! 終了フラグを点灯
            iTermination = 1
            do iMember = 1, UG%GM%BC%VW(iTmpWall)%iNumberOfMemberEdge - iTermination
                iMem1 = UG%GM%BC%VW(iTmpWall)%iMemberEdge(iMember)  ! iMember番目のメンバセルのセル番号
                iMem2 = UG%GM%BC%VW(iTmpWall)%iMemberEdge(iMember + 1)  ! iMember + 1番目のメンバセルのセル番号
                if (UG%Line%Distance(iMem1) > UG%Line%Distance(iMem2)) then   ! iMem1の距離の方がiMem2より遠いならば
                    iFlag = 0   ! 終了フラグを消灯し
                    call SwapInteger(UG%GM%BC%VW(iTmpWall)%iMemberEdge(iMember), UG%GM%BC%VW(iTmpWall)%iMemberEdge(iMember + 1))    ! メンバセルの登録順序を変更
                end if
            end do
            iTermination = iTermination + 1 ! ソート済み領域をソート対象から外す
        end do
    end do

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

    subroutine get_minimum_distance_bruteforce_edge(iEdgeNum, de_num, minDistance)
        integer, intent(in) :: iEdgeNum
        double precision, intent(out) :: de_num, minDistance
        double precision :: tmpDistance
        integer :: iTmpEdge
        minDistance = 100000
        do iTmpEdge = 1, iSurfaceCell   ! 物体表面の界面すべてについて
            tmpDistance = get_distance_e2c(UG%CD%Edge(iSurfaceEdge(iTmpEdge), :), UG%CD%Edge(iEdgeNum, :))  ! 与えられた辺との距離を計算
            if(tmpDistance < minDistance) then  ! 最小値が更新されたときのみ
                minDistance = tmpDistance
                de_num = dble(iTmpEdge)
            end if
        end do

        minDistance = AbsCrossProduct(UG%CD%Edge(iEdgeNum, :) - UG%CD%Edge(iSurfaceEdge(int(de_num)), :),&
                                    & UG%CD%Point(UG%Line%Point(iSurfaceEdge(int(de_num)), 1), :) - UG%CD%Point(UG%Line%Point(iSurfaceEdge(int(de_num)), 2), :)) &
                                    & / UG%GM%Area(iSurfaceEdge(int(de_num)))    ! 物体表面の辺と，界面中心からセル中心を結ぶ辺とがなす平行四辺形の高さを求める

        return
    end subroutine get_minimum_distance_bruteforce_edge

    ! メンバ総数の初期化
    subroutine viscosity_wall_initialize(BC, iWallNumber)
        type(BoundaryCondition), intent(inout) :: BC
        integer, intent(in) :: iWallNumber

        do iLoop = 1, iWallNumber
            BC%VW(iLoop)%iNumberOfMemberEdge = 0
        end do

        return
    end subroutine viscosity_wall_initialize

    ! メンバ界面登録用配列の動的確保
    subroutine viscosity_wall_malloc(BC, iWallNumber)
        type(BoundaryCondition), intent(inout) :: BC
        integer, intent(in) :: iWallNumber

        do iLoop = 1, iWallNumber
            allocate(BC%VW(iLoop)%iMemberEdge(BC%VW(iLoop)%iNumberOfMemberEdge))
        end do

        return
    end subroutine viscosity_wall_malloc

    subroutine check_sort
        do iTmpWall = 1, iSurfaceCell
            write(6, *) UG%GM%BC%VW(iTmpWall)%iNumberOfMemberEdge
            do iMember = 1, UG%GM%BC%VW(iTmpWall)%iNumberOfMemberEdge
                write(6, *) UG%GM%BC%VW(iTmpWall)%iMemberEdge(iMember), UG%Line%Distance(UG%GM%BC%VW(iTmpWall)%iMemberEdge(iMember))
            end do
            write(6,*) ""
            write(6,*) ""
            write(6,*) ""
        end do
    return
    end subroutine check_sort

end subroutine UGetDistanceFromSurface_Edge
