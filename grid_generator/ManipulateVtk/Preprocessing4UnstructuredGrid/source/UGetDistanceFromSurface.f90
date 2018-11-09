!***********************************/
!	Name:非構造格子用EulerSolverのための仮想格子属性割り当てプログラム
!	Alias:UGetDistanceFromSurface
!	Description:ユーザー指定基準で仮想セルに境界番号を割り当てる
!	Type:UGetDistanceFromSurface
!	Input:
!	Output:
!	Note:座標に基づいて色々やる，とりあえず，上界下界の座標基準と中心からの距離基準による境界属性割り当ては実装した，それ以外は今後の努力次第
!	Author:Akitaka Toyota
!	Date:2018.11.06
!	Update:-
!	Other:
!***********************************/
subroutine UGetDistanceFromSurface(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    integer :: iSurfaceCell

    integer :: iTmpEdge, iCount, iFlag, iSupport
    integer, allocatable :: iSurfaceEdge(:) ! surface edge number → edge number
    integer, allocatable :: iCWEdge(:) ! iSurfaceEdgeを時計周りor反時計周りに格納
    integer, allocatable :: iCheck(:)
    integer :: iNextCell, iYoungCell, iCandidateEdge
    double precision, allocatable :: dDistance(:, :)    ! iCell, 1:所属界面，2:距離

    iSurfaceCell = (UG%GI%AllCells - UG%GI%RealCells) - UG%GI%OutlineCells
    allocate(iSurfaceEdge(iSurfaceCell), iCWEdge(iSurfaceCell))

    iEdge = 1
    do iCell = UG%GI%RealCells+1, UG%GI%AllCells    ! すべての仮想セルについて巡回
        if (UG%VC%Type(iCell) == 2) then    ! 物体表面の仮想セルについて
            iSurfaceEdge(iEdge) = UG%VC%Edge(iCell) ! 大域辺番号をiSurfaceEdgeに格納
            iEdge = iEdge + 1
        end if
    end do

    ! この時点で物体表面の辺を順不同に全格納完了
    ! 物体表面の辺を時計回りor反時計周りに並べ替える
    iCount = 1
    iFlag = 1
    iEdge = iSurfaceEdge(15)
    iCWEdge(iCount) = iEdge

    do while(iCount < iSurfaceCell)
        iCount = iCount + 1
        iPoint = UG%Line%Point(iEdge, iFlag)    ! 探索対象辺の端点
        iFlag = 0   ! 隣接辺が見つかり次第flag = 1

        do iTmpEdge = 1, iSurfaceCell   ! 物体表面の辺についてループ
            iAdjacentEdge = iSurfaceEdge(iTmpEdge)  ! 隣接辺候補を格納
            if(iEdge /= iAdjacentEdge) then   ! 隣接辺候補が探索対象と異なる辺であり、かつ
                if(iPoint == UG%Line%Point(iAdjacentEdge, 1)) then   ! 端点を共有する辺なら隣接辺であるから
                    iCWEdge(iCount) = iAdjacentEdge ! iCWEdgeに隣接辺を登録し
                    iFlag = 2   ! 検索終了フラグを立てる
                else if (iPoint == UG%Line%Point(iAdjacentEdge, 2)) then
                    iCWEdge(iCount) = iAdjacentEdge ! iCWEdgeに隣接辺を登録し
                    iFlag = 1   ! 検索終了フラグを立てる
                end if
            end if
            if(iFlag /= 0) exit
        end do
        if(iFlag == 0) then
            write(6, *) "error", iCount, iEdge, iAdjacentEdge, iPoint
        end if
        iEdge = iAdjacentEdge
    end do

    ! すべてのセルについて、所属表面および所属表面からの距離を求める
    allocate(iCheck(UG%GI%RealCells), dDistance(UG%GI%RealCells, 2))
    iCount = 0  ! 初期化
    iCheck = 0  ! 初期化
    iSupport = 1

    do iCell = 1, UG%GI%RealCells
        if(iCheck(iCell) /= 1) then
            iCheck(iCell) = 1
            iCount = iCount + 1
            call get_minimum_distance_bruteforce(iCell, dDistance(iCell, 1), dDistance(iCell, 2))
        end if
    end do

    UG%Tri%Belongs2Wall = int(dDistance(:, 1))
    UG%Tri%Distance = dDistance(:, 2)
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

    subroutine get_minimum_distance_e2c(edge_num, de_num, minDistance)
        integer, intent(in) :: edge_num
        double precision, intent(out) :: de_num, minDistance
        integer, allocatable :: e_num(:)
        integer :: iTmp, iNum
        double precision, allocatable :: Distance(:)
        allocate(e_num(3), Distance(3))
        ! 調べる辺をリストアップ
        call get_edge_id(edge_num, e_num(1), e_num(2), e_num(3))
        ! 調べるべき3辺とセル中心の距離を求める
        do iTmp = 1, 3
            iCell = UG%Line%Cell(iCWEdge(e_num(iTmp)), 1, 1)  ! 最初のセル
            Distance(iTmp) = get_distance_e2c(UG%CD%Edge(iCWEdge(e_num(iTmp)), :), UG%CD%Cell(iCell, :))
        end do

        ! 最小の距離と，辺の所属を返す
        if (Distance(1) < Distance(2)) then
            if (Distance(1) < Distance(3)) then
                iNum = 1
            else
                iNum = 3
            end if
        else
            if (Distance(2) < Distance(3)) then
                iNum = 2
            else
                iNum = 3
            end if
        end if

        de_num = dble(iNum) !dble(e_num(iNum))
        minDistance = Distance(iNum)
        return
    end subroutine

    subroutine get_minimum_distance_bruteforce(iCellNum, de_num, minDistance)
        integer, intent(in) :: iCellNum
        double precision, intent(out) :: de_num, minDistance
        double precision :: tmpDistance
        minDistance = 100000
        do iEdge = 1, iSurfaceCell
            tmpDistance = get_distance_e2c(UG%CD%Edge(iSurfaceEdge(iEdge), :), UG%CD%Cell(iCellNum, :))
            if(tmpDistance < minDistance) then
                minDistance = tmpDistance
                de_num = iEdge
            end if
        end do

        ! call CrossProduct()
        minDistance = abs(dot_product(UG%CD%Edge(iSurfaceEdge(de_num), :) - UG%CD%Cell(iCellNum, :), -UG%GM%Normal(de_num, :)))    !法線ベクトルは物体外→物体内部に向いているため負号で反転
        !write(6, *) UG%CD%Edge(iSurfaceEdge(de_num), :) - UG%CD%Cell(iCellNum, :), -UG%GM%Normal(de_num, :)
        !write(6, *) minDistance

        return
    end subroutine get_minimum_distance_bruteforce

end subroutine UGetDistanceFromSurface
