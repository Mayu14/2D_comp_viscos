!***********************************/
!	Name:非構造格子のデータが入ったvtkファイルを読み込むプログラム
!	Alias:UReadRegionVTK
!	Description:読み込んだあと，点番号を+1する(配列回りの煩雑さをなくすため)
!	Type:UnstructuredGrid
!	Input:vtk
!	Output:UnstructuredGrid
!	Note:※このプログラムの末尾で読み出した点番号を1始まりに変更している
!           よって点の数だけ必要な配列を確保するときは1~総点数でよい
!           また，vtkに出力するときは，点番号-1とすること
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:-
!	Other:一応3次元も行けるはず
!***********************************/

subroutine UReadRegionVTK(UG,cInPath, cFileName)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    integer :: iCellType, iPointNumber
    integer :: iTmp
    character(len=64), intent(in) :: cFileName, cInPath
    character(len=64) :: cVTK,cAnnotate
    !点の番号と座標の対応
    !点の総数
    !実要素の総数
    !三角形要素の総数
    !三角形要素を構成する点の番号

    cVTK = trim(adjustl(cInPath))//trim(adjustl(cFileName))//".vtk"
    open(unit=1, file=cVTK, status='unknown') !配列確保用
    open(unit=2, file=cVTK, status='unknown') !データ読み出し用
        do iLoop=1,4
            read(1,*) cAnnotate !ヘッダー読み飛ばし
        end do

        read(1,*) cAnnotate, UG%GI%Points !点の総数
        !write(6,*)UG%GI%Points

        allocate(UG%CD%Point(UG%GI%Points,3)) !点座標格納用配列の動的割当
        read(1,*) ((UG%CD%Point(iPoint,iDim),iDim=1,3),iPoint=1,UG%GI%Points) !点の座標読み込み
        !write(6,*) UG%PointCoordinate

        !ここから下のプログラムは冗長だが，いずれ3次元やハイブリッド格子に手を出すために必要な工夫...のはず
        read(1,*) cAnnotate, UG%GI%RealCells

        do iLoop = 1, 5+(UG%GI%Points+UG%GI%RealCells+2) !unit=2をCell_Types直下まで移動させる!Header,座標，CELLS,注記
            read(2,*) cAnnotate
        end do

!以下if分岐の工夫でハイブリッド格子に対応可能
        iTriangleCell = 0
        do iCell=1,UG%GI%RealCells
            read(1,*) iPointNumber !CELLSの内容のうち1列目だけ読み続ける
            read(2,*) iCellType     !CELL_TYPESの内容を読み取る
            if(iPointNumber == 3) then !3点から構成されるもののうち
                if(iCellType == 5) then !3角形要素であるとき
                    iTriangleCell = iTriangleCell + 1 !3角系要素の総数を1つ増やす
                end if
            end if
        end do

        UG%Tri%Total = iTriangleCell !三角形要素の総数を記録
        allocate(UG%Tri%Point(UG%Tri%Total,3)) !三角形要素中の局所点番号と大域点番号を結びつける

        read(1,*) cAnnotate !unit=1をCellTypes直下へ移動

        rewind(2)
        do iLoop = 1, 5+(UG%GI%Points+1) !unit=2をCells直下まで移動
            read(2,*) cAnnotate
        end do

        iTriangleCell = 0
        do iCell=1, UG%GI%RealCells
            read(1,*) iCellType
            if(iCellType == 5) then !ハイブリッド格子の場合ここに条件追加
                iTriangleCell = iTriangleCell + 1
                read(2,*) iTmp,(UG%Tri%Point(iTriangleCell,iLocalPoint),iLocalPoint=1,3)
            end if
        end do
        UG%Tri%Point = UG%Tri%Point +1 !点番号は+1される

    close(1)
    close(2)
    return
end subroutine UReadRegionVTK
