!***********************************/
!	Name:非構造格子用EulerSolverのためのセル体積計算プログラム
!	Alias:UCalcCellVolume
!	Description:セルを構成する点番号と各点の座標情報を元に体積を計算する
!	Type:Unstructured Grid
!	Input:UG%GM%Dimension,UG%xxx%Point
!	Output:UG%GM%Volume
!	Note:2次元系におけるセル体積はセルの面積に相当する
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:-
!	Other:現時点でハイブリッド格子・3次元格子への対応は未定(準備はある)
!***********************************/
subroutine UCalcCellVolume(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    if(UG%GM%Dimension == 2) then
        call VolumeOfTriangle
        if(UG%Quad%Total /= 0) then
            write(1,*) "coming soon..."
        end if
    end if


return

contains
    subroutine VolumeOfTriangle
    implicit none
        do iCell=1, UG%Tri%Total
            iVertex1 = UG%Tri%Point(iCell,1)
            iVertex2 = UG%Tri%Point(iCell,2)
            iVertex3 = UG%Tri%Point(iCell,3)

            UG%GM%Volume(iCell) = 0.5d0*AbsCrossProduct(UG%CD%Point(iVertex2,:)-UG%CD%Point(iVertex1,:), & !外積の1/2が面積
                &   UG%CD%Point(iVertex3,:)-UG%CD%Point(iVertex2,:))
        end do
    end subroutine VolumeOfTriangle

end subroutine UCalcCellVolume
