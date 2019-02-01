!***********************************/
!	Name:MayuGrid2�`���̊i�q���ǂݍ��݃��\�b�h
!	Alias:GetCellAndEdgeCoordinate
!	Description:�o�[�W����0.9�ɑΉ�
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

    allocate(Geom%EdgeCoords(2, 0:Geom%CellNumber(1)+1, 0:Geom%CellNumber(2)+1, 1, 2))   ! West, South, Bottom�̏�  !
    allocate(tmpEdgeCoords(2, 0:Geom%CellNumber(1)+1, 0:Geom%CellNumber(2)+1, 1, 3:4))   ! West, South, Bottom�̏�  !
    allocate(Geom%CellCoords(2, 0:Geom%CellNumber(1)+1, 0:Geom%CellNumber(2)+1, 1))      !���z�Z�������m�ۂ��Ă���
    ! �Ƃ肠������Ԃ��Ă݂�(���x��������Ίi�q�����̒i�K�Ŋi�q�_����4�{�ɂ��āC���_���Ɗi�q��������)
    ! Edge_Center
    do iCenterY = 1, Geom%CellNumber(2)
        do iCenterX = 1, Geom%CellNumber(1) ! VertexCoords��(CellNum(1)+1, CellNum(2)+1)�Œ�`����Ă��邽��overflow����
            Geom%EdgeCoords(:, iCenterX, iCenterY, 1, 1) = 0.5 * (Geom%VertexCoords(:, iCenterX, iCenterY, 1) + Geom%VertexCoords(:, iCenterX, iCenterY + 1, 1))    ! West
            Geom%EdgeCoords(:, iCenterX, iCenterY, 1, 2) = 0.5 * (Geom%VertexCoords(:, iCenterX, iCenterY, 1) + Geom%VertexCoords(:, iCenterX + 1, iCenterY, 1))    ! South
            tmpEdgeCoords(:, iCenterX, iCenterY, 1, 3) = 0.5 * (Geom%VertexCoords(:, iCenterX + 1, iCenterY, 1) + Geom%VertexCoords(:, iCenterX + 1, iCenterY + 1, 1))    ! East
            tmpEdgeCoords(:, iCenterX, iCenterY, 1, 4) = 0.5 * (Geom%VertexCoords(:, iCenterX, iCenterY + 1, 1) + Geom%VertexCoords(:, iCenterX + 1, iCenterY + 1, 1)) !North
        end do
    end do
    ! y:fixed boundary(���`�O�})
    Geom%EdgeCoords(:, :, 0, 1, :) = 2.0 * Geom%EdgeCoords(:, :, 1, 1, :) - Geom%EdgeCoords(:, :, 2, 1, :)
    Geom%EdgeCoords(:, :, Geom%CellNumber(2) + 1, 1, :) = 2.0 * Geom%EdgeCoords(:, :, Geom%CellNumber(2), 1, :) - Geom%EdgeCoords(:, :, Geom%CellNumber(2) - 1, 1, :)

    ! x:loop boundary(���Α��̒l������)
    Geom%EdgeCoords(:, 0, :, 1,  :) = Geom%EdgeCoords(:, Geom%CellNumber(1), :, 1, :)
    Geom%EdgeCoords(:, Geom%CellNumber(1) + 1, :, 1,  :) = Geom%EdgeCoords(:, 1, :, 1, :)

    ! y:fixed boundary(���`�O�})
    tmpEdgeCoords(:, :, 0, 1, :) = 2.0 * tmpEdgeCoords(:, :, 1, 1, :) - tmpEdgeCoords(:, :, 2, 1, :)
    tmpEdgeCoords(:, :, Geom%CellNumber(2) + 1, 1, :) = 2.0 * tmpEdgeCoords(:, :, Geom%CellNumber(2), 1, :) - tmpEdgeCoords(:, :, Geom%CellNumber(2) - 1, 1, :)

    ! x:loop boundary(���Α��̒l������)
    tmpEdgeCoords(:, 0, :, 1,  :) = tmpEdgeCoords(:, Geom%CellNumber(1), :, 1, :)
    tmpEdgeCoords(:, Geom%CellNumber(1) + 1, :, 1,  :) = tmpEdgeCoords(:, 1, :, 1, :)

    ! Cell�@Center
    do iCenterY = 1, Geom%CellNumber(2)
        do iCenterX = 1, Geom%CellNumber(1) ! GridPoint��(CellNum(1)+1, CellNum(2)+1)�Œ�`����Ă��邽��overflow����
            Geom%CellCoords(:, iCenterX, iCenterY, 1) = 0.25 * (Geom%EdgeCoords(:, iCenterX, iCenterY, 1, 1)&
                                                             &+ Geom%EdgeCoords(:, iCenterX, iCenterY, 1, 2)&
                                                             &+ tmpEdgeCoords(:, iCenterX, iCenterY, 1, 3)&
                                                             &+ tmpEdgeCoords(:, iCenterX, iCenterY, 1, 4))
        end do
    end do
    deallocate(tmpEdgeCoords)

    ! y:fixed boundary(���`�O�})
    Geom%CellCoords(:, :, 0, 1, :) = 2.0 * Geom%CellCoords(:, :, 1, 1, :) - Geom%CellCoords(:, :, 2, 1, :)
    Geom%CellCoords(:, :, Geom%CellNumber(2) + 1, 1, :) = 2.0 * Geom%CellCoords(:, :, Geom%CellNumber(2), 1, :) - Geom%CellCoords(:, :, Geom%CellNumber(2) - 1, 1, :)

    ! x:loop boundary(���Α��̒l������)
    Geom%CellCoords(:, 0, :, 1,  :) = Geom%CellCoords(:, Geom%CellNumber(1), :, 1, :)
    Geom%CellCoords(:, Geom%CellNumber(1) + 1, :, 1,  :) = Geom%CellCoords(:, 1, :, 1, :)

    return
end subroutine GetCellAndEdgeCoordinate
