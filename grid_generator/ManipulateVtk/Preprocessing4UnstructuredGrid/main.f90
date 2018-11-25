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
    character(len=64) :: cFileName, c1, c2, c34, cWing, cInPath="", cOutPath=""
    logical :: ExistBound
    integer :: MultipleConvertMode = 0, i1, i2, i34

    MultipleConvertMode = 1    ! このオプションを使う場合は以下のelse内を書き換えて使う

    if(MultipleConvertMode /= 1) then
        write(6,*) "Please input the name of VTK file that defined region.(***.vtk's ***)"
        read(5,*) cFileName
        UG%GM%Dimension = 2
        call main_process(UG, cInPath, cOutPath, cFileName)

    else
        cInPath = "G:/Toyota/Data/grid_vtk/NACA4_vtk/"
        cOutPath = "G:/Toyota/Data/grid_vtk/NACA4_mayu/"
        do i1 = 0, 9
            do i2 =0, 9
                do i34 = 1, 40
                    write(cWing, '("NACA", i1, i1, i2, ".vtk")') i1, i2, i34
                    cFileName = adjustl(trim(cPath))//adjustl(trim(cWing))
                    call main_process(UG, cFileName)
                    call deallocate_UG
                end do
            end do
        end do
    end if
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
        allocate(UG%Line%Distance(UG%GI%Edges))
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

    subroutine main_process(UG, cInPath, cOutPath, cFileName)
        implicit none
        type(UnstructuredGrid), intent(out) :: UG
        character(len=64), intent(inout) :: cFileName, cInPath, cOutPath

    !vtkの読み込み
        call UReadRegionVTK(UG, cInPath, cFileName)

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
        call UGetDistanceFromSurface(UG, ExistBound)
        if(ExistBound == .True.) then
            call UGetDistanceFromSurface_Edge(UG)
        end if
    !テスト出力
        !call UCheckGrid(UG)

    !中間ファイルの出力
        call UOutputUnStrGrid(UG, cOutPath, cFileName, ExistBound)

    !再構成済みファイル
        print *, "Optimized Grid Data Generated!"

        return
    end subroutine main_process

    subroutine deallocate_UG
        implicit none

        if(allocated(UG%CD%Cell)) deallocate(UG%CD%Cell)
        if(allocated(UG%CD%Edge)) deallocate(UG%CD%Edge)

        if(allocated(UG%GM%Area)) deallocate(UG%GM%Area)
        if(allocated(UG%GM%Volume)) deallocate(UG%GM%Volume)
        if(allocated(UG%GM%Width)) deallocate(UG%GM%Width)

        if(allocated(UG%GM%Normal)) deallocate(UG%GM%Normal)
        if(allocated(UG%GM%AverageWidth)) deallocate(UG%GM%AverageWidth)

        if(allocated(UG%GM%CellType)) deallocate(UG%GM%CellType)
        if(allocated(UG%GM%CellType)) deallocate(UG%VC%Type)

        if(allocated(UG%GM%CellType)) deallocate(UG%InscribedCircle)
        ! for test
        if(allocated(UG%Tri%Belongs2Wall)) deallocate(UG%Tri%Belongs2Wall)
        if(allocated(UG%Tri%Distance)) deallocate(UG%Tri%Distance)

        if(allocated(UG%Line%Belongs2Wall)) deallocate(UG%Line%Belongs2Wall)
        if(allocated(UG%Line%Distance)) deallocate(UG%Line%Distance)

        if(allocated(UG%GM%BC%VW)) deallocate(UG%GM%BC%VW)
        if(allocated(UG%CH%PointNum)) deallocate(UG%CH%PointNum)
        if(allocated(UG%Tri%Cell)) deallocate(UG%Tri%Cell)
        if(allocated(UG%Tri%Edge)) deallocate(UG%Tri%Edge)
        if(allocated(UG%Line%Cell)) deallocate(UG%Line%Cell)
        if(allocated(UG%Line%Point)) deallocate(UG%Line%Point)
        if(allocated(UG%VC%Cell)) deallocate(UG%VC%Cell)
        if(allocated(UG%IO%PointNum)) deallocate(UG%IO%PointNum)
        if(allocated(UG%CD%Point)) deallocate(UG%CD%Point)
        if(allocated(UG%Tri%Point)) deallocate(UG%Tri%Point)

        if(allocated(UG%GM%BC%VW)) then
            do iLoop = 1, UG%GM%BC%iWallTotal
                deallocate(UG%GM%BC%VW(iLoop)%iMemberCell)
                deallocate(UG%GM%BC%VW(iLoop)%iMemberEdge)
            end do
            deallocate(UG%GM%BC%VW)
        end if

        return
    end subroutine deallocate_UG

end program

