!***********************************/
!	Name:非構造格子用EulerSolverのための界面面積計算プログラム
!	Alias:UCalcEdgeArea
!	Description:界面を構成する点番号と各点の座標情報を元に面積を計算する
!	Type:UnstructuredGrid
!	Input:UG%GM%Dimension,UG%xxx%Point (xxxは使用する種類のセルデータ)
!	Output:UG%GM%Area
!	Note:2次元系における界面面積は辺の長さに相当する
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:2018.02.03 !Add AverageWidth
!	Other:現時点でハイブリッド格子・3次元格子への対応は未定(準備はある)
!***********************************/
subroutine UCalcEdgeArea(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    if(UG%GM%Dimension == 2) then
        call AreaOfLine
        call GetCellAverageArea

    else
        write(6,*) "coming soon..."
        stop
    end if
return

contains
    subroutine AreaOfLine
    implicit none
        do iEdge=1, UG%Line%Total
            iEndPoint1 = UG%Line%Point(iEdge,1)
            iEndPoint2 = UG%Line%Point(iEdge,2)
            UG%GM%Area(iEdge) = sqrt(sum((UG%CD%Point(iEndPoint2,:)-UG%CD%Point(iEndPoint1,:))**2))
        end do
    end subroutine AreaOfLine

    subroutine GetCellAverageArea
    implicit none

    do iCell=1, UG%GI%RealCells
        UG%GM%AverageWidth(iCell) = (UG%GM%Area(UG%Tri%Edge(iCell,1)) + UG%GM%Area(UG%Tri%Edge(iCell,2)) + UG%GM%Area(UG%Tri%Edge(iCell,3)))/3.0d0
    end do

    return
    end subroutine GetCellAverageArea

end subroutine UCalcEdgeArea
