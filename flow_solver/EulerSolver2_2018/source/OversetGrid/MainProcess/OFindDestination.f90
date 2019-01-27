!***********************************/
!	Name:移動先の格子を補間するための補間データ提供元格子検索
!	Alias:OFindDestination
!	Description:
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.11.16
!	Update:2017.11.16
!	Other:
!***********************************/
subroutine OFindDestination(UG,Geom,OSG,MG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(Geometry), intent(in) :: Geom
    type(OverSetGrid), intent(inout) :: OSG
    type(MoveGrid), intent(inout) :: MG

    call LocalGrid2AuxiliaryGrid
    call LocalGrid2GlobalGrid

return
contains

    subroutine LocalGrid2AuxiliaryGrid
    implicit none
    double precision, allocatable :: AuxiliaryGridNum(:)
    integer, allocatable :: iGlobalCellNumber(:)
    integer :: iOverflowCell

    allocate(AuxiliaryGridNum(3))
    allocate(iGlobalCellNumber(9))

        MG%NSR%InterpolateFromGlobal = 0

        do iCell=1, UG%GI%AllCells !すべてのローカルセルについて
            iOverflowCell = 0
            !Auxiliary格子の座標は配列の添字と対応しているため，
            !(非構造格子の座標r,Auxiliary格子の座標r0，構造格子負端からのセル数n，格子幅dxについて
            !r-r0 > n*dxなる関係があることから，nについて解くと以下のような式になる(※2次元用で書いたためzは常に固定値になる)
            AuxiliaryGridNum(1:2) = int((MG%NSR%C%NextCoordsInN(iCell,1:2) - MG%ASG%Bound(1,1:2))/(MG%ASG%Width(1:2)))+1
            AuxiliaryGridNum(3) = 1

            if(AuxiliaryGridNum(1) > MG%ASG%CellNumber(1)) then
                iOverflowCell = 9
            else if(AuxiliaryGridNum(1) < 1) then
                iOverflowCell = 9
            else if(AuxiliaryGridNum(2) > MG%ASG%CellNumber(2)) then
                iOverflowCell = 9
            else if(AuxiliaryGridNum(2) < 1) then
                iOverflowCell = 9
            end if

            if(iOverflowCell == 0) then
                !補助格子の隣接8セルと自セルをグローバル番号で表記する
                iGlobalCellNumber(1) = (AuxiliaryGridNum(2)-2)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) - 1 !(x-1,y-1)
                iGlobalCellNumber(2) = (AuxiliaryGridNum(2)-2)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)     !(x,y-1)
                iGlobalCellNumber(3) = (AuxiliaryGridNum(2)-2)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) + 1 !(x+1,y-1)

                iGlobalCellNumber(4) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) - 1 !(x-1,y)
                iGlobalCellNumber(5) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)     !(x,y)
                iGlobalCellNumber(6) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) + 1 !(x+1,y)

                iGlobalCellNumber(7) = (AuxiliaryGridNum(2))*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) - 1   !(x-1,y+1)
                iGlobalCellNumber(8) = (AuxiliaryGridNum(2))*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)       !(x,y+1)
                iGlobalCellNumber(9) = (AuxiliaryGridNum(2))*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) + 1   !(x+1,y+1)


            !補間の対象となるはずの補助格子のグローバル番号を格子に登録する
                do iLoop=1,9
                    if(iGlobalCellNumber(iLoop) > MG%ASG%TotalCell) then !計算したグローバル番号が補助格子中に存在していない場合補間対象から除外
                        iOverflowCell = iOverflowCell + 1

                    else if(iGlobalCellNumber(iLoop) < 0) then
                        iOverflowCell = iOverflowCell + 1

                    else if(MG%ASG%ContentsAS2U(iGlobalCellNumber(iLoop)) == 0) then !計算したグローバル番号の補助格子中にローカル格子が存在しない場合補間対象から除外
                        iOverflowCell = iOverflowCell + 1
                    end if
                    if(iOverflowCell == 9) exit
                end do

            end if

            if(iOverflowCell < 2) then
                MG%NSR%Previous2Next(iCell,:) = iGlobalCellNumber(:)
            else
                MG%NSR%Previous2Next(iCell,:) = 0
                MG%NSR%InterpolateFromGlobal(iCell) = 1
            end if

        end do

        !call CheckOfInterpolateSource
    return
    end subroutine LocalGrid2AuxiliaryGrid

        subroutine CheckOfInterpolateSource
        implicit none
        iCenterX=0
        iCenterY=0
        do iLoop=1,UG%GI%AllCells
            if(MG%NSR%InterpolateFromGlobal(iLoop) == 0) then
                iCenterX = iCenterX+1
            else
                iCenterY = iCenterY+1
            end if
        end do
        write(6,*) UG%GI%AllCells,iCenterX,iCenterY
        stop
        return
        end subroutine CheckOfInterpolateSource


    subroutine LocalGrid2GlobalGrid
    implicit none

        integer :: iGlobalCellNumber

        OSG%FrequencyDistribution = 0

        do  iCell=1, UG%GI%AllCells

            OSG%RelatedU2S(iCell,1:2) = int((MG%NSR%C%NextCoordsInG(iCell,1:2) - Geom%Bound(1,1:2))/(2.0d0*Geom%Width(1,1,1:2)))+1
            OSG%RelatedU2S(iCell,3) = 1

            iGlobalCellNumber = (OSG%RelatedU2S(iCell,2)-1)*Geom%CellNumber(1) + OSG%RelatedU2S(iCell,1)

            !構造格子ijkが非構造格子から参照された回数を記録
            OSG%FrequencyDistribution(iGlobalCellNumber) = OSG%FrequencyDistribution(iGlobalCellNumber) +1

            !最後に構造格子に所属する非構造格子への対応を記録する
            !OSG%RelatedS2U(iGlobalCellNumber,OSG%FrequencyDistribution(iGlobalCellNumber)) = iAdjacentCell !11/21
            OSG%RelatedS2U(iGlobalCellNumber,OSG%FrequencyDistribution(iGlobalCellNumber)) = iCell !11/22
        end do

    return
    end subroutine LocalGrid2GlobalGrid


end subroutine OFindDestination
