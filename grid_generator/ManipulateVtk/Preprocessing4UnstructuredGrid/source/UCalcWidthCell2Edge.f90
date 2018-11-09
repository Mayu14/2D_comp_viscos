!***********************************/
!	Name:非構造格子用EulerSolverのための，要素中心界から面中心までの距離ベクトルを計算するプログラム
!	Alias:UCalcWidthCell2Edge
!	Description:セルごとに界面中心から面中心までの距離を求めると同時に，セルごとの平均距離を求める(Venkatakrishnanの流束制限関数用)
!	Type:UnstructuredGrid
!	Input:UG%GM%Dimension,UG%CD%Edge,UG%CD%Cell (xxxは使用する種類のセルデータ)
!	Output:UG%GM%Normal
!	Note:UG%GM%Widthは引数が3つもあってややこしいので気をつけること
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:2018.02.03 !delete AverageWidth
!	Other:現時点でハイブリッド格子・3次元格子への対応は未定(準備はある)
!***********************************/
subroutine UCalcWidthCell2Edge(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    if(UG%GM%Dimension == 2) then
        call WidthTriangle2Line
        if(UG%Quad%Total /= 0) then
            write(6,*) "coming soon...?"
        end if
    else
        write(6,*) "coming soon..."
        stop
    end if
return

contains
    subroutine WidthTriangle2Line
    implicit none
        !UG%GM%AverageWidth = 0.0d0
        do iCell=1, UG%Tri%Total
            do iLocalEdge=1,3
                iEdge = UG%Tri%Edge(iCell,iLocalEdge)
                UG%GM%Width(iCell,iLocalEdge,:) = UG%CD%Edge(iEdge,:) - UG%CD%Cell(iCell,:)
                !UG%GM%AverageWidth(iCell) = UG%GM%AverageWidth(iCell) + sqrt(sum(UG%GM%Width(iCell,iLocalEdge,:)**2)) !第2項はベクトルの絶対値
            end do
            !UG%GM%AverageWidth(iCell) = UG%GM%AverageWidth(iCell)/3.0d0
        end do
        return
    end subroutine WidthTriangle2Line

end subroutine UCalcWidthCell2Edge
