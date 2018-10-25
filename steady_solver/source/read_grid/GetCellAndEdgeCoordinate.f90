!***********************************/
!	Name:MayuGrid2形式の格子情報読み込みメソッド
!	Alias:GetCellAndEdgeCoordinate
!	Description:バージョン0.9に対応
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.??
!	Update:2018.10.25
!	Other:
!***********************************/
subroutine GetCellAndEdgeCoordinate(Geom)
    use StructVar_Mod_mk2
    use LoopVar_Mod_mk2
    implicit none

    type(Geometry), intent(inout) :: Geom
    double precision, allocatable :: tmpEdgeCoords(:, :, :, :, :)

    allocate(Geom%EdgeCoords(2, 0:Geom%CellNumber(1)+1, 0:Geom%CellNumber(2)+1, 1, 2))   ! West, South, Bottomの順  !
    allocate(tmpEdgeCoords(2, 0:Geom%CellNumber(1)+1, 0:Geom%CellNumber(2)+1, 1, 3:4))   ! West, South, Bottomの順  !
    allocate(Geom%CellCoords(2, 0:Geom%CellNumber(1)+1, 0:Geom%CellNumber(2)+1, 1))      !仮想セル分も確保しておく
    ! とりあえず補間してみる(精度が悪ければ格子生成の段階で格子点数を4倍にして，中点ごと格子生成する)
    ! Edge_Center
    do iCenterY = 1, Geom%CellNumber(2)
        do iCenterX = 1, Geom%CellNumber(1) ! VertexCoordsは(CellNum(1)+1, CellNum(2)+1)で定義されているためoverflowせず
            Geom%EdgeCoords(:, iCenterX, iCenterY, 1, 1) = 0.5 * (Geom%VertexCoords(:, iCenterX, iCenterY, 1) + Geom%VertexCoords(:, iCenterX, iCenterY + 1, 1))    ! West
            Geom%EdgeCoords(:, iCenterX, iCenterY, 1, 2) = 0.5 * (Geom%VertexCoords(:, iCenterX, iCenterY, 1) + Geom%VertexCoords(:, iCenterX + 1, iCenterY, 1))    ! South
            tmpEdgeCoords(:, iCenterX, iCenterY, 1, 3) = 0.5 * (Geom%VertexCoords(:, iCenterX + 1, iCenterY, 1) + Geom%VertexCoords(:, iCenterX + 1, iCenterY + 1, 1))    ! East
            tmpEdgeCoords(:, iCenterX, iCenterY, 1, 4) = 0.5 * (Geom%VertexCoords(:, iCenterX, iCenterY + 1, 1) + Geom%VertexCoords(:, iCenterX + 1, iCenterY + 1, 1)) !North
        end do
    end do
    ! y:fixed boundary(線形外挿)
    Geom%EdgeCoords(:, :, 0, 1, :) = 2.0 * Geom%EdgeCoords(:, :, 1, 1, :) - Geom%EdgeCoords(:, :, 2, 1, :)
    Geom%EdgeCoords(:, :, Geom%CellNumber(2) + 1, 1, :) = 2.0 * Geom%EdgeCoords(:, :, Geom%CellNumber(2), 1, :) - Geom%EdgeCoords(:, :, Geom%CellNumber(2) - 1, 1, :)

    ! x:loop boundary(反対側の値を入れる)
    Geom%EdgeCoords(:, 0, :, 1,  :) = Geom%EdgeCoords(:, Geom%CellNumber(1), :, 1, :)
    Geom%EdgeCoords(:, Geom%CellNumber(1) + 1, :, 1,  :) = Geom%EdgeCoords(:, 1, :, 1, :)

    ! y:fixed boundary(線形外挿)
    tmpEdgeCoords(:, :, 0, 1, :) = 2.0 * tmpEdgeCoords(:, :, 1, 1, :) - tmpEdgeCoords(:, :, 2, 1, :)
    tmpEdgeCoords(:, :, Geom%CellNumber(2) + 1, 1, :) = 2.0 * tmpEdgeCoords(:, :, Geom%CellNumber(2), 1, :) - tmpEdgeCoords(:, :, Geom%CellNumber(2) - 1, 1, :)

    ! x:loop boundary(反対側の値を入れる)
    tmpEdgeCoords(:, 0, :, 1,  :) = tmpEdgeCoords(:, Geom%CellNumber(1), :, 1, :)
    tmpEdgeCoords(:, Geom%CellNumber(1) + 1, :, 1,  :) = tmpEdgeCoords(:, 1, :, 1, :)

    ! Cell　Center
    do iCenterY = 1, Geom%CellNumber(2)
        do iCenterX = 1, Geom%CellNumber(1) ! GridPointは(CellNum(1)+1, CellNum(2)+1)で定義されているためoverflowせず
            Geom%CellCoords(:, iCenterX, iCenterY, 1) = 0.25 * (Geom%EdgeCoords(:, iCenterX, iCenterY, 1, 1)&
                                                             &+ Geom%EdgeCoords(:, iCenterX, iCenterY, 1, 2)&
                                                             &+ tmpEdgeCoords(:, iCenterX, iCenterY, 1, 3)&
                                                             &+ tmpEdgeCoords(:, iCenterX, iCenterY, 1, 4))
        end do
    end do
    deallocate(tmpEdgeCoords)

    ! y:fixed boundary(線形外挿)
    Geom%CellCoords(:, :, 0, 1, :) = 2.0 * Geom%CellCoords(:, :, 1, 1, :) - Geom%CellCoords(:, :, 2, 1, :)
    Geom%CellCoords(:, :, Geom%CellNumber(2) + 1, 1, :) = 2.0 * Geom%CellCoords(:, :, Geom%CellNumber(2), 1, :) - Geom%CellCoords(:, :, Geom%CellNumber(2) - 1, 1, :)

    ! x:loop boundary(反対側の値を入れる)
    Geom%CellCoords(:, 0, :, 1,  :) = Geom%CellCoords(:, Geom%CellNumber(1), :, 1, :)
    Geom%CellCoords(:, Geom%CellNumber(1) + 1, :, 1,  :) = Geom%CellCoords(:, 1, :, 1, :)

    return
end subroutine GetCellAndEdgeCoordinate
