!***********************************/
!	Name:重合格子法の対象となる要素の登録
!	Alias:ORegisterOverSetBounds
!	Description:
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.11.16
!	Update:2017.11.22
!	Other:
!***********************************/
subroutine ORegisterOverSetBounds(UG,Geom,OSG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(Geometry), intent(inout) :: Geom
    type(OverSetGrid), intent(inout) :: OSG
    integer :: iGlobalCellNumber,iOverSetCell
    integer :: iStructCell, iTotalStructCells
 !各変数の初期化
        OSG%RelatedU2S = 0
        OSG%RelatedS2U = 0
        OSG%FrequencyDistribution = 0
        Geom%CellType = 0
        iOverSetCell = 0

    do iCell=UG%GI%RealCells+1, UG%GI%AllCells !すべての仮想セルについて
        iAdjacentCell = UG%VC%Cell(iCell,1)     !仮想セルの隣接実セル = 非構造格子の外部境界セル

        !構造直交格子の座標は配列の添字と対応しているため，
        !(非構造格子の座標r,構造格子の座標r0，構造格子負端からのセル数n，格子幅dxについて
        !r-r0 > n*dxなる関係があることから，nについて解くと以下のような式になる(※2次元用で書いたためzは常に固定値になる)
        OSG%RelatedU2S(iAdjacentCell,1:2) = int((UG%CD%Cell(iAdjacentCell,1:2) - Geom%Bound(1,1:2))/(2.0d0*Geom%Width(1,1,1:2)))+1
        OSG%RelatedU2S(iAdjacentCell,3) = 1

        !取扱いの利便性を上げるためグローバルセル番号(詳細は↓)に変換する．
        !(i,j,k) <=>!2D:(j-1)*(Lx)+i !3D:(k-1)*(Ly*Lx)+(j-1)*(jmax-1)+i
        iGlobalCellNumber = (OSG%RelatedU2S(iAdjacentCell,2)-1)*Geom%CellNumber(1) + OSG%RelatedU2S(iAdjacentCell,1)

        if(OSG%FrequencyDistribution(iGlobalCellNumber) == 0) then !初めて参照されたとき
            iOverSetCell = iOverSetCell + 1 !構造格子における重合格子セルの数を1増やす
        end if
        !構造格子ijkが非構造格子から参照された回数を記録
        OSG%FrequencyDistribution(iGlobalCellNumber) = OSG%FrequencyDistribution(iGlobalCellNumber) +1

        !最後に構造格子に所属する非構造格子への対応を記録する
        !OSG%RelatedS2U(iGlobalCellNumber,OSG%FrequencyDistribution(iGlobalCellNumber)) = iAdjacentCell !11/21
        OSG%RelatedS2U(iGlobalCellNumber,OSG%FrequencyDistribution(iGlobalCellNumber)) = iCell !11/22
    end do

    if(ubound(OSG%OverSetStructCell,1) /= iOverSetCell) then !前回と比べて重合格子セルの総数に変化があったときのみ
        if(allocated(OSG%OverSetStructCell)) deallocate(OSG%OverSetStructCell)
        allocate(OSG%OverSetStructCell(iOverSetCell,5)) !新な大きさで再確保
    end if

    iXmax = max(1,Geom%CellNumber(1))
    iYmax = max(1,Geom%CellNumber(2))
    iZmax = max(1,Geom%CellNumber(3))

    iTotalStructCells = iXmax * iYmax * iZmax

!構造格子側の重合格子境界の番号付けと，格子番号との対応付け
    OSG%OverSetStructCell = -1
    iOverSetCell = 0
    do iStructCell = 1, iTotalStructCells
        if(OSG%FrequencyDistribution(iStructCell) == 0) cycle !重合格子に関係ないセルはループの先頭へ
        !iOverSetCell is Local-OverSet-Cell-Number
        iOverSetCell = iOverSetCell + 1
        OSG%OverSetStructCell(iOverSetCell,1) = iStructCell-int(iStructCell/iXmax)*iXmax
        OSG%OverSetStructCell(iOverSetCell,2) = int(iStructCell/iXmax)
        OSG%OverSetStructCell(iOverSetCell,3) = 1
        OSG%OverSetStructCell(iOverSetCell,4) = OSG%FrequencyDistribution(iStructCell)
        OSG%OverSetStructCell(iOverSetCell,5) = iStructCell !internal to Global
    end do

    !call TestOverSetOutput

return
contains
    subroutine TestOverSetOutput
    implicit none

    character(len=64) :: cFileName
    integer,allocatable :: CellNumber(:)

    allocate(CellNumber(2))

    cFileName = "TestOverSetGrid.vtk"

    CellNumber(:) = max(1,Geom%CellNumber(2:3))

    open(unit=1,file=cFileName,status='unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,"('Struct_ShockTube')")
        write(1,"('ASCII')")
        write(1,"('DATASET RECTILINEAR_GRID')")
        write(1,"('DIMENSIONS ',3(1x,i5))") Geom%CellNumber(1)+1,CellNumber(1)+1,CellNumber(2)

!座標情報
        write(1,"('X_COORDINATES ',i5,' float')") Geom%CellNumber(1)+1
            do iCenterX=1, Geom%CellNumber(1)+1
                write(1,"(f17.14)") Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * (dble(iCenterX)-1.0d0)
            end do

        write(1,"('Y_COORDINATES ',i5,' float')") CellNumber(1)+1
        do iCenterY=1, CellNumber(1)+1
            write(1,"(f17.14)") Geom%Bound(1,2) + 2.0d0 * Geom%Width(1,1,2) * (dble(iCenterY)-1.0d0)
        end do

        write(1,"('Z_COORDINATES ',i5,' float')") CellNumber(2)
        do iCenterZ=1, CellNumber(2)
            write(1,"(f17.14)") 0.0
        end do

        write(1,"('CELL_DATA ',i7)") (Geom%CellNumber(1)) * (CellNumber(1))


!密度プロット
        write(1,"('SCALARS OverSetGrid float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCenterZ=1,CellNumber(2)
            do iCenterY=1,CellNumber(1)
                do iCenterX=1,Geom%CellNumber(1)
                        iGlobalCellNumber = (iCenterY-1)*Geom%CellNumber(1) + iCenterX
                        write(6,*) iGlobalCellNumber
                        write(1,"((2x,f22.14))") float(OSG%FrequencyDistribution(iGlobalCellNumber))
                end do
            end do
        end do
        close(1)
        return
    end subroutine TestOverSetOutput

end subroutine ORegisterOverSetBounds
