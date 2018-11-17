!***********************************/
!	Name:非構造格子用EulerSolverのための前処理プログラム(計算格子データ生成)
!	Alias:Preprocessing4UnstructuredGrid
!	Description:非構造格子用のvtkファイルから計算用格子：UnStrGridを作成する・格子点番号の振り直し，領域の分割には対応していない
!	Type:Main
!	Input:vtkファイル
!	Output:UnStrGridファイル
!	Note:変換したいvtkファイルはプログラムと同じor下のディレクトリに置くこと．
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:2018.01.20
!	Other:現時点でハイブリッド格子・3次元格子への対応は未定
!***********************************/
program Preprocessing4UnstructuredGrid
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid) :: UG
    character(len=64) :: cFileName

    write(6,*) "Please input the name of VTK file that defined region.(***.vtk's ***)"
    read(5,*) cFileName
    !cFileName = 'Square2DMeshRough'
    !write(6,*) "Please input Dimension to calculate."
    !read(5,*) UG%GM%Dimension
    UG%GM%Dimension = 2

!vtkの読み込み
    call UReadRegionVTK(UG,cFileName)

!共有辺の抽出
    call UMakeEdgeNumber(UG)

!ここでとりあえず必要なデータが全部揃うのでデータ格納用配列を全部割り当て
    call AllocVariables_P4U

!各種幾何学的量を求める
!座標の中心
    call UCalcElementCenter(UG)

!セル体積
    call UCalcCellVolume(UG)

!界面面積
    call UCalcEdgeArea(UG)

!要素中心から界面までの距離ベクトル
    call UCalcWidthCell2Edge(UG)

!要素界面における法線ベクトルの計算
    call UCalcNormalVector(UG)

!Inscribed Circle Raduis !need to Volume and Area
    call UCalcInscribedCircleOfCell(UG)

!並べ替えの基準に用いる仮の仮想セル中心を求める
    !call UDataPointOfVirtualCell(1) !Preprocess

!仮想セルの並べ替え(外周部から内周部へ)
    !call UReSortVirtualCell(UG)

!実際の仮想セル中心を計算する
    call UDataPointOfVirtualCell(2) !MainProcess

!仮想セルの属性付け
    call UMarkingVirtualCell(UG)

!格子外周の凸包の点列を作成
    call UMakeConvexHull(UG)

    if(UG%VC%Total /= UG%GI%OutlineCells) then
        !内部に物体を持つ格子のみ，"物体"っぽく見えるものを出力する
        call UMakeInternalObject(UG)
    end if

!物体表面からの距離を与える
    call UGetDistanceFromSurface(UG)

!テスト出力
    !call UCheckGrid(UG)

!中間ファイルの出力
    call UOutputUnStrGrid(UG,cFileName)

!do iCell=UG%GI%RealCells+1, UG%GI%AllCells
!    write(6,*) "Length2",dot_product(UG%CD%Cell(iCell,:), UG%CD%Cell(iCell,:)), UG%GM%CellType(iCell,1,1)
!    write(6,*) "AdjReCell/Edge",UG%VC%Cell(iCell,1),UG%VC%Cell(iCell,2)
!    write(6,*) "EdgeNum,CellNum",UG%VC%Edge(iCell),UG%Tri%Edge(UG%VC%Cell(iCell,1),UG%VC%Cell(iCell,2))
!    write(6,*) "VCx,y",UG%CD%Cell(iCell,1),UG%CD%Cell(iCell,2)
!    write(6,*) "RCx,y",UG%CD%Cell(UG%VC%Cell(iCell,1),1),UG%CD%Cell(UG%VC%Cell(iCell,1),2)
!    write(6,*) "REx,y",UG%CD%Edge(UG%VC%Edge(iCell),1),UG%CD%Edge(UG%VC%Edge(iCell),2)
!    write(6,*) ""
!end do
!再構成済みファイル

    print *, "Optimized Grid Data Generated!"
    stop

contains
    subroutine AllocVariables_P4U
    implicit none

        allocate(UG%CD%Cell(UG%GI%AllCells,3)) !実格子の数でいいのか，仮想格子も数えあげるべきなのか
        allocate(UG%CD%Edge(UG%GI%Edges,3))

        allocate(UG%GM%Area(UG%GI%Edges))
        allocate(UG%GM%Volume(UG%GI%RealCells))
        allocate(UG%GM%Width(UG%GI%RealCells,3,3)) !三角形格子のみを仮定

        allocate(UG%GM%Normal(UG%GI%Edges,3))
        allocate(UG%GM%AverageWidth(UG%GI%RealCells))

        allocate(UG%GM%CellType(UG%GI%AllCells,1,1))
        allocate(UG%VC%Type(UG%GI%RealCells+1:UG%GI%AllCells))

        allocate(UG%InscribedCircle(UG%GI%RealCells)) !Radius of Cell's Inscribed Circle
        ! for test
        allocate(UG%Tri%Belongs2Wall(UG%GI%RealCells))
        allocate(UG%Tri%Distance(UG%GI%RealCells))

        allocate(UG%Line%Belongs2Wall(UG%GI%Edges))
        allocate(UG%Line%Distance(UG%GI%RealCells))
    return
    end subroutine AllocVariables_P4U


    subroutine UDataPointOfVirtualCell(iStep)
    implicit none
    integer :: iStep !1:Pre, 2:After"ReSort"
    !仮想格子の中心は，隣接セル中心と界面中心について点対称である位置に設定する
    if(iStep == 1) then
        do iElement = UG%GI%RealCells+1, UG%GI%AllCells
            iShareEdge = UG%VC%Edge(iElement)
            iAdjacentCell = UG%VC%Cell(iElement,1)
            iAdjacentEdge = UG%VC%Cell(iElement,2)
            UG%CD%Cell(iElement,:) = UG%CD%Edge(iShareEdge,:) + UG%GM%Width(iAdjacentCell,iAdjacentEdge,:)
        end do

    else if(iStep == 2) then
        do iElement = UG%GI%RealCells+1, UG%GI%AllCells
            iShareEdge = UG%VC%Edge(iElement) !GlobalEdgeNumber
            iAdjacentCell = UG%VC%Cell(iElement,1)
            iAdjacentEdge = UG%VC%Cell(iElement,2) !LocalEdgeNumber

            if(dot_product(UG%CD%Edge(iShareEdge,1:3)-UG%CD%Cell(iAdjacentCell,1:3),UG%GM%Normal(iShareEdge,1:3)) < 0.0d0) then !面の法線は常に格子の外を向いている．隣接実セルの中心から境界界面へ向かうベクトルと，界面法線ベクトルの内積の符号で内外判定できる
                UG%CD%Cell(iElement,:) = UG%CD%Edge(iShareEdge,:) + UG%GM%Width(iAdjacentCell,iAdjacentEdge,:) !outer boundary
            else
                UG%CD%Cell(iElement,:) = UG%CD%Edge(iShareEdge,:) - UG%GM%Width(iAdjacentCell,iAdjacentEdge,:) !internal boundary
            end if

        end do

    end if

    return
    end subroutine UDataPointOfVirtualCell

end program

