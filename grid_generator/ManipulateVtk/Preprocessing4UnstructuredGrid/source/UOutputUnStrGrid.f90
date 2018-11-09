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
subroutine UOutputUnStrGrid(UG,cFileName)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    character(len=64), intent(inout) :: cFileName

    !write(6,*) "Please input filename of the Computed Grid Data"
    !read(5,*) cFileName
    cFileName = trim(adjustl(cFileName))//".mayu" !調べたら.usgも.ccdもなんか既に使われてるみたいだったのでキレた
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
        write(1,*) "TotalCeonvexHullPoint", UG%CH%iTotal
        write(1,*) "InternalObjectPoint", UG%IO%iTotal

!三角形要素の幾何的関係
        write(1,*) "TriPoint"
        do iCell=1, UG%Tri%Total
            write(1,*) (UG%Tri%Point(iCell,iLoop),iLoop=1,3)
        end do

        write(1,*) "TriEdge"
        do iCell=1, UG%Tri%Total
            write(1,*) (UG%Tri%Edge(iCell,iLoop),iLoop=1,3)
        end do

        write(1,*) "TriCell"
        do iCell=1, UG%Tri%Total
            write(1,*) (UG%Tri%Cell(iCell,iLoop),iLoop=1,3)
        end do

!四辺形要素の幾何的関係
        if(UG%Quad%Total /= 0) then
            write(1,*) "coming soon..."
            stop
        end if

!線要素の幾何的関係
        write(1,*) "LinePoint"
        do iEdge=1, UG%Line%Total
            write(1,*) (UG%Line%Point(iEdge,iLoop),iLoop=1,2)
        end do

        write(1,*) "LineCell"
        do iEdge=1, UG%Line%Total
            write(1,*) ((UG%Line%Cell(iEdge,iSide,iLocalEdge),iLocalEdge=1,2),iSide=1,2)
        end do

!仮想格子の幾何的関係
        write(1,*) "VCCell"
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            write(1,*) (UG%VC%Cell(iCell,iLoop),iLoop=1,2)
        end do

        write(1,*) "VCEdge"
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            write(1,*) UG%VC%Edge(iCell)
        end do

!格子点の座標情報
        write(1,*) "PointC"
        do iPoint=1, UG%GI%Points
            write(1,*) (UG%CD%Point(iPoint,iLoop),iLoop=1,3)
        end do

!格子境界中心の座標情報
        write(1,*) "EdgeC"
        do iEdge=1,UG%GI%Edges
            write(1,*) (UG%CD%Edge(iEdge,iLoop),iLoop=1,3)
        end do

!格子要素中心の座標情報
        write(1,*) "CellC"
        do iCell=1,UG%GI%AllCells
            write(1,*) (UG%CD%Cell(iCell,iLoop),iLoop=1,3)
        end do

!格子境界の面積
        write(1,*) "EdgeS"
        do iEdge=1,UG%GI%Edges
            write(1,*) UG%GM%Area(iEdge)
        end do

!格子境界の体積
        write(1,*) "CellV"
        do iCell=1,UG%GI%RealCells
            write(1,*) UG%GM%Volume(iCell)
        end do

!セル界面における法線ベクトル
        write(1,*) "EdgeNormal"
        do iEdge=1,UG%GI%Edges
            write(1,*) (UG%GM%Normal(iEdge,iLoop),iLoop=1,3)
        end do

!セル中心から界面中心までの距離ベクトル(中心基準)
        write(1,*) "Width"
        do iCell=1,UG%GI%RealCells
            write(1,*) ((UG%GM%Width(iCell,iLocalEdge,iLoop),iLoop=1,3),iLocalEdge=1,3)
        end do

!仮想セルに適用する境界条件
        write(1,*) "BoudnaryCondition of VirtualCell"
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            write(1,*) UG%GM%CellType(iCell,1,1)
        end do

!Radius of Cell's Inscribed-Circle for compute CFL Condition
        write(1,*) "InscribedCircleRadius"
        do iCell=1,UG%GI%RealCells
            write(1,*) UG%InscribedCircle(iCell)
        end do

!内部に含む物体を構成する凸包のデータ
        write(1,*) "InternalObject"
        do iPoint=1, UG%IO%iTotal
            write(1,*) UG%IO%PointNum(iPoint)
        end do

!MUSCLのLimiterで使用する各セルにおける界面までの平均距離
        write(1,*) "AverageWid"
        do iCell=1,UG%GI%RealCells
            write(1,*) UG%GM%AverageWidth(iCell)
        end do

!格子外周を構成する凸包のデータ
        write(1,*) "ConvexHull"
        do iPoint=1, UG%CH%iTotal
            write(1,*) UG%CH%PointNum(iPoint)
        end do



    close(1)
return
end subroutine UOutputUnStrGrid
