!***********************************/
!	Name:非構造格子用EulerSolverのための要素中心計算プログラム
!	Alias:UCalcElementCenter
!	Description:セル・界面を構成する点番号と各点の座標情報を元に，それぞれの中心座標を計算する
!	Type:UnstructuredGrid
!	Input:UG%GM%Dimension,UG%xxx%Point (xxxは使用する種類のセルデータ)
!	Output:UG%CD%Edge,UG%CD%Cell
!	Note:中心は重心を利用，次元の辺要素のみ中点を利用する
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:-
!	Other:現時点でハイブリッド格子・3次元格子への対応は未定(準備はある)
!***********************************/
subroutine UCalcElementCenter(UG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    if(UG%GM%Dimension == 2) then
        iELM => iElement
        call GravityCenterOfTriangle
        call MidPointOfLine
    else
        write(6,*) "coming soon..."
        stop
    end if

return
contains
!_________________________________________________
    subroutine GravityCenterOfTriangle
    implicit none
        do iElement=1,UG%Tri%Total
            iVertex1 = UG%Tri%Point(iELM,1) !三角形の頂点の点番号
            iVertex2 = UG%Tri%Point(iELM,2)
            iVertex3 = UG%Tri%Point(iELM,3)
            !重心
            UG%CD%Cell(iELM,:) = (UG%CD%Point(iVertex1,:) + UG%CD%Point(iVertex2,:) &
                                        &  + UG%CD%Point(iVertex3,:))/3.0d0
        end do
    return
    end subroutine GravityCenterOfTriangle

!_________________________________________________

    subroutine MidPointOfLine
    implicit none
        do iElement=1, UG%Line%Total
            iEndPoint1 = UG%Line%Point(iELM,1)
            iEndPoint2 = UG%Line%Point(iELM,2)
            UG%CD%Edge(iELM,:) = 0.5d0*(UG%CD%Point(iEndPoint1,:)+UG%CD%Point(iEndPoint2,:))
        end do
    return
    end subroutine MidPointOfLine

end subroutine UCalcElementCenter



