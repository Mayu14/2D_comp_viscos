!***********************************/
!	Name:非構造格子用EulerSolverのための要素中心計算プログラム
!	Alias:UCalcNormalVector
!	Description:界面を構成する点番号と各点の座標情報を元に，面ごとの単位法線ベクトルを計算する，なお表裏の判定のためにセル中心から界面までの距離ベクトルを用いる
!	Type:UnstructuredGrid
!	Input:UG%GM%Dimension,UG%xxx%Point,UG%GM%Width (xxxは使用する種類のセルデータ)
!	Output:UG%GM%Normal
!	Note:隣接要素の要素番号が大きい側が面の表である．表側のセル中心→界面までの距離ベクトルとの内積が負になるとき法線は表を向いていることになる．
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:-
!	Other:現時点でハイブリッド格子・3次元格子への対応は未定(準備はある)
!***********************************/
subroutine UCalcNormalVector(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    if(UG%GM%Dimension == 2) then
        call NormalVectorOfLine
    else
        write(6,*) "coming soon..."
        stop
    end if
return

contains
    subroutine NormalVectorOfLine
    implicit none
    double precision :: InverseLength
        do iEdge=1, UG%Line%Total
            !界面表方向の単位法線ベクトル = 隣接要素のうち要素番号が大きい方を向いた単位法線ベクトル
            !表側のセル中心における，中心から界面への距離ベクトルとの内積が負になる方向
            !※内積の値が-1からずれるほど格子の歪みが大きいといえる

            iEndPoint1 = UG%Line%Point(iEdge,1)
            iEndPoint2 = UG%Line%Point(iEdge,2)
            iAdjacentCell = UG%Line%Cell(iEdge,1,1)
            iLocalEdge = UG%Line%Cell(iEdge,1,2)
            !UG%GM%Width
            UG%GM%Normal(iEdge,1) =  UG%CD%Point(iEndPoint2,2) - UG%CD%Point(iEndPoint1,2) !内積=0より
            UG%GM%Normal(iEdge,2) = -UG%CD%Point(iEndPoint2,1) + UG%CD%Point(iEndPoint1,1)
            UG%GM%Normal(iEdge,3) = 0.0d0

            InverseLength = 1.0d0/sqrt(dot_product(UG%GM%Normal(iEdge,:),UG%GM%Normal(iEdge,:)))
            UG%GM%Normal(iEdge,:) = InverseLength * UG%GM%Normal(iEdge,:)

            if(dot_product(UG%GM%Normal(iEdge,:),UG%GM%Width(iAdjacentCell,iLocalEdge,:)) >= 0) then
                UG%GM%Normal(iEdge,:) = - UG%GM%Normal(iEdge,:)
            end if

        end do
        return
    end subroutine NormalVectorOfLine

end subroutine UCalcNormalVector
