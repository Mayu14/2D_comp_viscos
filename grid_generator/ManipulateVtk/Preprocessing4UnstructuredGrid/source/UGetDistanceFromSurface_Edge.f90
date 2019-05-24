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
subroutine UGetDistanceFromSurface_Edge(UG, ExistInnerBoundary)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG
    logical, intent(out) :: ExistInnerBoundary

    integer :: iSurfaceCell ! 物体表面セルの総数 = 物体表面界面の数

    double precision :: tmp_de_num
    double precision, allocatable :: dDistance(:, :)    ! iCell, 1:所属界面，2:距離

    integer :: iTmpWall, iMember, iMem1, iMem2, iTermination ! 仮置き壁番号, メンバセル用ループ変数, メンバセル1, メンバセル2, バブルソートの終端
    integer, allocatable :: iSurfaceEdge(:), iCount(:)   ! 各壁ごとの入力済みカウント用配列
    integer :: iFlag = 0    ! 並び替え終了フラグ
    type(DummySE), allocatable :: DSE
    integer :: iLeftEdge, iRightEdge, iGLEdge, iGCEdge, iGREdge

    iSurfaceCell = (UG%GI%AllCells - UG%GI%RealCells) - UG%GI%OutlineCells  ! 物体表面界面数の特定
    UG%GM%BC%iWallTotal = iSurfaceCell

    if(iSurfaceCell == 0) then
        write(6,*) "No Internal Boundary"
        ExistInnerBoundary = .false.
    else
        ExistInnerBoundary = .true.
        allocate(DSE%SE(iSurfaceCell))
        allocate(iSurfaceEdge(iSurfaceCell))    ! 物体表面の局所界面番号→大域界面番号
        allocate(UG%GM%BC%VW(iSurfaceCell))     ! 粘性壁登録変数の動的確保

        ! 物体表面界面を反時計周りに並べ替え
        iLocalEdge = 1   ! 物体表面(=壁)上の局所界面番号の初期化
        do iCell = UG%GI%RealCells+1, UG%GI%AllCells    ! すべての仮想セルについて巡回
            if (UG%VC%Type(iCell) == 2) then    ! 物体表面の仮想セルについて
                DSE%SE(iLocalEdge)%iGEdgeNum = UG%VC%Edge(iCell) ! 大域界面番号をiSurfaceEdgeに格納
                DSE%SE(iLocalEdge)%iLEdgeNum = iLocalEdge   ! 物体表面のみの局所界面番号
                DSE%SE(iLocalEdge)%iyNormalSign = int(sign(1.0d0, UG%GM%Normal(UG%VC%Edge(iCell), 2)))  !法線ベクトルy成分の符号
                DSE%SE(iLocalEdge)%dxCoord = UG%CD%Edge(iCell, 1)   ! 界面中心x座標
                iLocalEdge = iLocalEdge + 1
            end if
        end do
        ! この時点で物体表面の辺を順不同に全格納完了
        call sort_edge_ccw(DSE, iSurfaceCell)
        ! すべての辺を反時計周りに格納

        iLeftEdge = iSurfaceCell
        iRightEdge = 2
        do iLocalEdge = 1, iSurfaceCell    ! すべての仮想セルについて巡回
            iSurfaceEdge(iLocalEdge) = DSE%SE(iLocalEdge)%iGEdgeNum ! 大域辺番号をiSurfaceEdgeに格納
            UG%GM%BC%VW(iLocalEdge)%iGlobalEdge = DSE%SE(iLocalEdge)%iGEdgeNum  ! 1行上のiSurfaceEdgeと同じやつ
            iGLEdge = DSE%SE(iLeftEdge)%iGEdgeNum
            iGCEdge = DSE%SE(iLocalEdge)%iGEdgeNum
            iGREdge = DSE%SE(iRightEdge)%iGEdgeNum

            UG%GM%BC%VW(iLocalEdge)%curvature = get_curvature_of_edge(UG%CD%Edge(iGCEdge,1:2), &
                                            & UG%CD%Edge(iGLEdge,1:2), UG%CD%Edge(iGREdge,1:2), &
                                            & UG%GM%Normal(iGLEdge,1:2), UG%GM%Normal(iGREdge,1:2))

            iLeftEdge = iLeftEdge + 1
            iRightEdge = iRightEdge + 1
            if(iLeftEdge == iSurfaceCell + 1) then
                iLeftEdge = 1
            end if
            if(iRightEdge == iSurfaceCell + 1) then
                iRightEdge = 1
            end if
        end do
        ! この時点で物体表面の辺を順不同に全格納完了

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
        allocate(iCount(iSurfaceCell))  ! 入力済み確認配列（局所壁番号）
        iCount = 1
        do iEdge = 1, UG%GI%Edges   ! すべての界面について
            iTmpWall = UG%Line%Belongs2Wall(iEdge)   ! 所属壁の局所壁番号を仮置きして
            UG%GM%BC%VW(iTmpWall)%iMemberEdge(iCount(iTmpWall)) = iEdge ! 壁番号→所属している界面の番号
            iCount(iTmpWall) = iCount(iTmpWall) + 1 ! 入力した回数を反映する
        end do

        ! 距離準並び替えを実施(もうめんどいしバブルソートでええやろ(適当))(ボトルネックになるようなら書き換える感じで)
        do iTmpWall = 1, iSurfaceCell   ! すべての界面について
            iFlag = 0
            iTermination = 1
            do while(iFlag == 0)
                iFlag = 1   ! 終了フラグを点灯
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
    end if

    return
contains

    double precision function get_distance_e2e(edge1, edge2) result(distance)
        double precision, intent(in) :: edge1(:), edge2(:)
        distance = sqrt(dot_product(edge1 - edge2, edge1 - edge2))
    end function get_distance_e2e

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
        double precision :: tmpDistance, minD_direct, minD_projection
        integer :: iTmpEdge
        minD_direct = 100000

        do iTmpEdge = 1, iSurfaceCell   ! 物体表面の界面すべてについて
            tmpDistance = get_distance_e2e(UG%CD%Edge(UG%GM%BC%VW(iTmpEdge)%iGlobalEdge, :), UG%CD%Edge(iEdgeNum, :))  ! 与えられた辺との直線距離を計算
            if(tmpDistance < minD_direct) then  ! 最小値が更新されたときのみ
                minD_direct = tmpDistance
                de_num = dble(iTmpEdge)
            end if
        end do

        minD_projection = AbsCrossProduct(UG%CD%Edge(iEdgeNum, :) - UG%CD%Edge(UG%GM%BC%VW(int(de_num))%iGlobalEdge, :),&
                                    & UG%CD%Point(UG%Line%Point(UG%GM%BC%VW(int(de_num))%iGlobalEdge, 1), :) - UG%CD%Point(UG%Line%Point(UG%GM%BC%VW(int(de_num))%iGlobalEdge, 2), :)) &
                                    & / UG%GM%Area(UG%GM%BC%VW(int(de_num))%iGlobalEdge)    ! 物体表面の辺と，界面中心からセル中心を結ぶ辺とがなす平行四辺形の高さを求める

        if(sqrt(2.0d0) * minD_projection > minD_direct) then
            minDistance = minD_projection
        else
            minDistance = minD_direct
        end if
        minDistance = minD_direct
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

    subroutine sort_edge_ccw(DSE, iEdgeTotal)
        implicit none
        type(DummySE), intent(inout) :: DSE
        type(DummySE) :: tmpDSE
        integer :: inner, outer, iEdgeTotal
        integer :: iTop, iBottom

        ! バブルソートにて辺のx座標規準で並べ替え
        do outer = 1, iEdgeTotal - 1
            do inner = 2, iEdgeTotal
                call swap_edge_SE(DSE, inner, inner + 1)
            end do
        end do

        ! 界面法線ベクトルのy方向成分が，正：上から，負：下から配列に登録することでccwにならべる
        allocate(tmpDSE%SE(iEdgeTotal))
        iTop = 1
        iBottom = iEdgeTotal
        do iLocalEdge = 1, iEdgeTotal
            if(DSE%SE(iLocalEdge)%iyNormalSign < 0) then
                tmpDSE%SE(iBottom) = DSE%SE(iLocalEdge)
                iBottom = iBottom - 1
            else
                tmpDSE%SE(iTop) = DSE%SE(iLocalEdge)
                iTop = iTop + 1
            end if
        end do
        DSE = tmpDSE
        deallocate(tmpDSE%SE)
        return
    end subroutine sort_edge_ccw

    subroutine swap_edge_SE(DSE,iNum1, iNum2)
        implicit none
        type(DummySE), intent(inout) :: DSE
        integer, intent(in) :: iNum1, iNum2
        type(SurfaceEdge) :: SE

            if DSE%SE(iNum1)%dxCoord > DSE%SE(iNum2)%dxCoord then
                SE = DSE%SE(iNum1)
                DSE%SE(iNum1) = DSE%SE(iNum2)
                DSE%SE(iNum2) = SE
            end if

        return
    end subroutine swap_edge_SE

    function get_curvature_of_edge(edge_center, center_i, center_j, normal_i, normal_j) result(curvature)
        use FrequentOperation
        implicit none
        double precision, intent(in) :: edge_center(:), center_i(:), center_j(:)
        double precision, intent(in) :: normal_i(:), normal_j(:)
        double precision :: curvature_radius, curvature
        double precision, allocatable :: coef(:,:), const(:), intersect(:)
        integer :: eps_exp = 13
        logical :: parallel

        allocate(coef(2,2), const(2), intersect(2))
        ! 係数行列
        coef(1,1) = normal_i(2)
        coef(1,2) = - normal_i(1)
        coef(2,1) = normal_j(2)
        coef(2,2) = - normal_j(1)
        ! 定数ベクトル
        const(1) = normal_i(1) * center_i(2) - normal_i(2) * center_i(1)
        const(2) = normal_j(1) * center_j(2) - normal_j(2) * center_j(1)
        ! 逆行列を作成する(反転不可能なときはparallel==.true.となる)
        call set_2d_inverse_matrix(coef, parallel)

        if (parallel == .false.) then
            intersect = matmul(coef, const) ! 交点
            curvature_radius = sqrt((intersect(1)**2 - edge_center(1)**2) + (intersect(2)**2 - edge_center(2)**2))  ! 曲率半径
        else
            curvature_radius = 10.0d0**eps_exp + 1
        end if
        if (curvature_radius > 10.0d0**eps_exp) then
            curvature = 0.0d0
        else
            curvature = 1.0d0 / curvature_radius
        end if
        return
    end function get_curvature_of_edge

end subroutine UGetDistanceFromSurface_Edge
