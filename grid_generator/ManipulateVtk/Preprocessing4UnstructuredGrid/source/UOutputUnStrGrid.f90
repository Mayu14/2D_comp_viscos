!***********************************/
!	Name:非構造格子を計算用に整理したUnStrGridを出力するプログラム
!	Alias:UOutputUnStrGrid
!	Description:UnStrGrid(通称USG:ウサギ)を出力する．一応テキストで出力してるから読めるけど専用のreadプログラムが必要
!	Type:UnstructuredGrid
!	Input:UG%CD,UG%GI,UG%Tri,UG%VC,UG%Line,UG%GM (xxxは使用する種類のセルデータ)
!	Output:UnStrGrid
!	Note:拡張子は4文字でmayuにしたので，読み取りプログラムは*.mayuを探せばいいんじゃない？(適当)
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:2018.02.03
!	Other:GJKアルゴリズムによる接触判定を実行するため凸包出力を追加/複数の境界条件へ対応
!***********************************/
subroutine UOutputUnStrGrid(UG, cPath, cFileName, ExistBound)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    character(len=256), intent(inout) :: cFileName, cPath
    logical, intent(in) :: ExistBound
    !write(6,*) "Please input filename of the Computed Grid Data"
    !read(5,*) cFileName
    cFileName = trim(adjustl(cPath))//trim(adjustl(cFileName))//".mayu" !調べたら.usgも.ccdもなんか既に使われてるみたいだったのでキレた
    !cFileName = "UnStrGrid"

!点番号は1始まり
    open(unit=1,file=trim(adjustl(cFileName)),status='unknown')
!格子の基本情報
        write(1,*) "TotalPoint", UG%GI%Points
        write(1,*) "TotalEdge", UG%GI%Edges
        write(1,*) "TotalRealCell", UG%GI%RealCells
        write(1,*) "TotalVirtualCell", UG%VC%Total
        write(1,*) "TotalOutlineCell", UG%GI%OutlineCells

        write(1,*) "TotalTriangle", UG%Tri%Total
        write(1,*) "TotalSquare", UG%Quad%Total
        write(1,*) "TotalLine", UG%Line%Total
        write(1,*) "TotalConvexHullPoint", 0    !UG%CH%iTotal
        write(1,*) "InternalObjectPoint", UG%IO%iTotal

!三角形要素の幾何的関係
        write(1,*) "TriPoint ", UG%Tri%Total * 3
        do iCell=1, UG%Tri%Total
            write(1,"(3(1x,i7))") (UG%Tri%Point(iCell,iLoop),iLoop=1,3)
        end do

        write(1,*) "TriEdge ", UG%Tri%Total * 3
        do iCell=1, UG%Tri%Total
            write(1,"(3(1x,i7))") (UG%Tri%Edge(iCell,iLoop),iLoop=1,3)
        end do

        write(1,*) "TriCell ", UG%Tri%Total * 3
        do iCell=1, UG%Tri%Total
            write(1,"(3(1x,i7))") (UG%Tri%Cell(iCell,iLoop),iLoop=1,3)
        end do

!四辺形要素の幾何的関係
        if(UG%Quad%Total /= 0) then
            write(6,*) "coming soon..."
            stop
        end if

!線要素の幾何的関係
        write(1,*) "LinePoint ", UG%Line%Total * 2
        do iEdge=1, UG%Line%Total
            write(1,"(2(1x,i7))") (UG%Line%Point(iEdge,iLoop),iLoop=1,2)
        end do

        write(1,*) "LineCell ", UG%Line%Total * 4
        do iEdge=1, UG%Line%Total
            write(1,"(4(1x,i7))") ((UG%Line%Cell(iEdge,iSide,iLocalEdge),iLocalEdge=1,2),iSide=1,2)
        end do

!仮想格子の幾何的関係
        write(1,*) "VCCell ", (UG%GI%AllCells - UG%GI%RealCells) * 2
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            write(1,"(2(1x,i7))") (UG%VC%Cell(iCell,iLoop),iLoop=1,2)
        end do

        write(1,*) "VCEdge ", UG%GI%AllCells - UG%GI%RealCells
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            write(1,"(1x,i7)") UG%VC%Edge(iCell)
        end do

!格子点の座標情報
        write(1,*) "PointC ", UG%GI%Points * 3
        do iPoint=1, UG%GI%Points
            write(1,"(3(1x,f22.17))") (UG%CD%Point(iPoint,iLoop),iLoop=1,3)
        end do

!格子境界中心の座標情報
        write(1,*) "EdgeC ", UG%GI%Edges * 3
        do iEdge=1,UG%GI%Edges
            write(1,"(3(1x,f22.17))") (UG%CD%Edge(iEdge,iLoop),iLoop=1,3)
        end do

!格子要素中心の座標情報
        write(1,*) "CellC ", UG%GI%AllCells * 3
        do iCell=1,UG%GI%AllCells
            write(1,"(3(1x,f22.17))") (UG%CD%Cell(iCell,iLoop),iLoop=1,3)
        end do

!格子境界の面積
        write(1,*) "EdgeS ", UG%GI%Edges
        do iEdge=1,UG%GI%Edges
            write(1,"((2x,f22.14))") UG%GM%Area(iEdge)
        end do

!格子境界の体積
        write(1,*) "CellV ", UG%GI%RealCells
        do iCell=1,UG%GI%RealCells
            write(1,"((2x,f22.14))") UG%GM%Volume(iCell)
        end do

!セル界面における法線ベクトル
        write(1,*) "EdgeNormal ", UG%GI%Edges * 3
        do iEdge=1,UG%GI%Edges
            write(1,"(3(1x,f22.17))") (UG%GM%Normal(iEdge,iLoop),iLoop=1,3)
        end do

!セル中心から界面中心までの距離ベクトル(中心基準)
        write(1,*) "Width ", UG%GI%RealCells
        do iCell=1,UG%GI%RealCells
            write(1,"(3(1x,f22.17))") ((UG%GM%Width(iCell,iLocalEdge,iLoop),iLoop=1,3),iLocalEdge=1,3)
        end do

!仮想セルに適用する境界条件
        write(1,*) "BoudnaryCondition_of_VirtualCell ", UG%GI%AllCells - UG%GI%RealCells
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            write(1,"(1x,i7)") UG%GM%CellType(iCell,1,1)
        end do

!Radius of Cell's Inscribed-Circle for compute CFL Condition
        write(1,*) "InscribedCircleRadius ", UG%GI%RealCells
        do iCell=1,UG%GI%RealCells
            write(1,*) UG%InscribedCircle(iCell)
        end do

!内部に含む物体を構成する凸包のデータ
        write(1,*) "InternalObject ", UG%IO%iTotal
        do iPoint=1, UG%IO%iTotal
            write(1,"(1x,i7)") UG%IO%PointNum(iPoint)
        end do

!MUSCLのLimiterで使用する各セルにおける界面までの平均距離
        write(1,*) "AverageWid ", UG%GI%RealCells
        do iCell=1,UG%GI%RealCells
            write(1,"((2x,f22.14))") UG%GM%AverageWidth(iCell)
        end do

!格子外周を構成する凸包のデータ
        !write(1,*) "ConvexHull ", UG%CH%iTotal
        !do iPoint=1, UG%CH%iTotal
            !write(1,"(1x,i7)") UG%CH%PointNum(iPoint)
        !end do

        if(ExistBound .eqv. .true.) then
    ! セルから最も近い物体表面セル
            !write(1, *) "NearestSurfaceBoundaryEdgeNum ", UG%Tri%Total
            !do iCell = 1, UG%Tri%Total
                !write(1, *) UG%Tri%Belongs2Wall(iCell)
            !end do

    ! 物体表面からの距離
            !write(1, *) "DistanceFromObjectSurface ", UG%Tri%Total
            !do iCell = 1, UG%Tri%Total
                !write(1, *) UG%Tri%Distance(iCell)
            !end do
            !write(1,*) ""

    ! 界面から最も近い物体表面
            write(1, *) "NearestSurfaceBoundaryEdgeNum4Edge ", UG%Line%Total
            do iEdge = 1, UG%Line%Total
                write(1, "(1x,i7)") UG%Line%Belongs2Wall(iEdge)
            end do

    ! 物体表面からの距離
            write(1, *) "DistanceFromObjectSurface4Edge ", UG%Line%Total
            do iEdge = 1, UG%Line%Total
                write(1, "((2x,f22.14))") UG%Line%Distance(iEdge)
            end do
            write(1,*) ""

    ! 各壁に対応する
            !write(1,*) "Wall2Cell_data ", UG%GM%BC%iWallTotal
            !do iEdge = 1, UG%GM%BC%iWallTotal
                !write(1, *) UG%GM%BC%VW(iEdge)%iGlobalEdge, UG%GM%BC%VW(iEdge)%iNumberOfMember
                !do iCell = 1, UG%GM%BC%VW(iEdge)%iNumberOfMember
                    !write(1,*) UG%GM%BC%VW(iEdge)%iMemberCell(iCell)
                !end do
            !end do

            write(1,*) "Wall2Edge_data ", UG%GM%BC%iWallTotal
            do iLoop = 1, UG%GM%BC%iWallTotal
                write(1, *) UG%GM%BC%VW(iLoop)%iGlobalEdge, UG%GM%BC%VW(iLoop)%iNumberOfMemberEdge
                do iEdge = 1,  UG%GM%BC%VW(iLoop)%iNumberOfMemberEdge
                    write(1,"(1x,i7)") UG%GM%BC%VW(iLoop)%iMemberEdge(iEdge)
                end do
            end do
        end if

    close(1)

    return
end subroutine UOutputUnStrGrid
