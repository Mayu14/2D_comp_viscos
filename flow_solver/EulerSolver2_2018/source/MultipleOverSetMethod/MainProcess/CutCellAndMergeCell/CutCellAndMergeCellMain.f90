!***********************************/
!	Name:カットセル・マージセルの制御プログラム
!	Alias:OCutCellAndMergeCellMain
!	Description:カットセル・マージセルの制御プログラム
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.23
!	Update:
!	Other:not parallelize
!***********************************/
    subroutine OCutCellAndMergeCellMain!(UG,MG,RCN,EOM)
    use StructVar_Mod
    implicit none

    type(LocalGrids) :: LG
    type(Geometry) :: Geom
    integer :: iGrid, iEdge, iCell, iAdjacentRealCell
    integer :: GetGlobalNumber2D
    double precision :: MinimulArea
    !0.Initialize
        !call Initilize4CutCell

    !1.界面の接触判定
    do iGrid = 1, LG%TotalGridsNumber
        do iCell = 1, LG%UG(iGrid)%GI%RealCells+1, LG%UG(iGrid)%GI%OutlineCells
            LG%UG(iGrid)%VC%CutEdge(iCell)%InternalPointNumber = 0 !initialize
            iAdjacentRealCell = LG%UG(iGrid)%VC%Cell(iCell,1)

            call GetEndpointWithoutInterference(LG,iGrid,iCell,iAdjacentRealCell) !L格子の境界界面のうち実際にG格子と接触しているもののみを取り出す

            if(LG%UG(iGrid)%VC%CutEdge(iCell)%CrossPattern /= 2) then !G格子との接触またはG/L両方の格子と接触しているときカットセル処理対象

                LG%UG(iGrid)%VC%CutEdge(iCell)%GlobalNumber1 = GetGlobalNumber2D(LG%UG(iGrid)%VC%CutEdge(iCell)%EndPoint(1,1:2),&
                                                                &   LG%GG%GASG%Width(1:2),LG%GG%GASG%Bound(1,1:2),LG%GG%GASG%CellNumber(1))

                LG%UG(iGrid)%VC%CutEdge(iCell)%GlobalNumber2 = GetGlobalNumber2D(LG%UG(iGrid)%VC%CutEdge(iCell)%EndPoint(2,1:2),&
                                                                &   LG%GG%GASG%Width(1:2),LG%GG%GASG%Bound(1,1:2),LG%GG%GASG%CellNumber(1))


                if(LG%UG(iGrid)%VC%CutEdge(iCell)%GlobalNumber1 == LG%UG(iGrid)%VC%CutEdge(iCell)%GlobalNumber2) then !2つの端点が同じGセルに属する時

                    call RegisterCutCondition(0,LG%UG(iGrid)%VC%CutEdge(iCell),LG%GG%CutC(LG%UG(iGrid)%VC%CutEdge(iCell)%GlobalNumber1)) !Edge側には分割した点，Gセル側には新しい内点を登録する


                else !2つの端点が異なるGセルに属するとき
                    call PointSwap(LG%UG(iGrid)%VC%CutEdge(iCell)) !GlobalNumberの大小で点座標の入れ替え
                    call CheckContactEdge2Edge(Geom,LG%UG(iGrid)%VC%CutEdge(iCell),LG%GG) !接触判定&界面のリストを作成
                    !call RegisterCutCondition !引数で中身の処理制御してあげる，これはCheckContactに内包されている
                    call RegisterEdgeLengthRatio(LG%UG(iGrid)%VC%CutEdge(iCell),LG%GG)!界面に登録した点の番号を並び変えて各線分の長さ比率をGセルに登録する
                end if

            end if
        end do
    end do

    call MergeIEC2withIEC(LG%GG)

    do iCell=1, LG%GG%GASG%TotalCell !G格子の全セルに関するループ
        if(LG%GG%CutC(iCell)%UnsolvedFlag == 1) then !未解決フラグが1であるセルに対して
!            call GetNewAreaOfCuttedCell(LG%GG%CutC(iCell)) !カット後のセル面積を計算する

            !未解決フラグ = 0
            LG%GG%CutC(iCell)%UnsolvedFlag = 0

            if(LG%GG%CutC(iCell)%CutVolume < MinimulArea) then
                !call MergeBesideCell !もっとも長い界面で隣接する要素に接続する
                !隣接させた要素の未解決フラグ = 1
            end if

        end if


    end do

    return
    end subroutine OCutCellAndMergeCellMain
