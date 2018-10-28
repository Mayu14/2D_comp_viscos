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
subroutine GetGridJacobian(Geom, CC, CE)
    use StructVar_Mod_mk2
    use LoopVar_Mod_mk2
    use FrequentOperation_mk2
    implicit none

    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CellEdge), intent(inout) :: CE
    double precision :: detJ
    double precision, allocatable :: r_xi(:), r_eta(:)

    allocate(r_xi(2), r_eta(2))

    ! Cell Edge
    do iCenterY = 1, Geom%CellNumber(2) + 1
        do iCenterX = 1, Geom%CellNumber(1) + 1
            if ((iCenterY == Geom%CellNumber(2) + 1) .and. (iCenterX == Geom%CellNumber(1) + 1)) then

            else
                do iEdge = 1, 2
                    if(iEdge == 1) then
                        r_xi(:) = Geom%CellCoords(:, iCenterX, iCenterY, 1) - Geom%CellCoords(:, iCenterX - 1, iCenterY, 1)
                        if (iCenterY == Geom%CellNumber(2) + 1) then
                            r_eta(:) = 2.0d0 * Geom%VertexCoords(:, iCenterX, iCenterY, 1) &
                                            &- Geom%VertexCoords(:, iCenterX, iCenterY - 1, 1)
                        else
                            r_eta(:) = Geom%VertexCoords(:, iCenterX, iCenterY + 1, 1) - Geom%VertexCoords(:, iCenterX, iCenterY, 1)
                        end if

                    else
                        if (iCenterX == Geom%CellNumber(1) + 1) then
                            r_xi(:) = 2.0d0 * Geom%VertexCoords(:, iCenterX, iCenterY, 1)&
                                           &- Geom%VertexCoords(:, iCenterX - 1, iCenterY, 1)
                        else
                            r_xi(:) = Geom%VertexCoords(:, iCenterX + 1, iCenterY, 1) - Geom%VertexCoords(:, iCenterX, iCenterY, 1)
                        end if

                        r_eta(:) = Geom%CellCoords(:, iCenterX, iCenterY, 1) - Geom%CellCoords(:, iCenterX, iCenterY - 1, 1)
                    end if

                    detJ = 1.0 / DeterminantOf2DMatrix(r_xi(1), r_eta(1), r_xi(2), r_eta(2))
                    CE%GridJacobian(iCenterX, iCenterY, 1, iEdge) = detJ
                    CE%GridJacobiMatrix(1, 1, iCenterX, iCenterY, 1, iEdge) = detJ * r_eta(2)   ! xi_x
                    CE%GridJacobiMatrix(1, 2, iCenterX, iCenterY, 1, iEdge) = -detJ * r_eta(1)   ! xi_y
                    CE%GridJacobiMatrix(2, 1, iCenterX, iCenterY, 1, iEdge) = -detJ * r_xi(2)   ! eta_x
                    CE%GridJacobiMatrix(2, 2, iCenterX, iCenterY, 1, iEdge) = detJ * r_xi(1)   ! eta_y

                end do
            end if
        end do
    end do

    ! Cell Center
    do iCenterY = 1, Geom%CellNumber(2)
        do iCenterX = 1, Geom%CellNumber(1)
            r_xi(:) = Geom%EdgeCoords(:, iCenterX + 1, iCenterY, 1, 1) - Geom%EdgeCoords(:, iCenterX, iCenterY, 1, 1)   !west
            r_eta(:) = Geom%EdgeCoords(:, iCenterX, iCenterY + 1, 1, 2) - Geom%EdgeCoords(:, iCenterX, iCenterY, 1, 2)  !south

            detJ = 1.0 / DeterminantOf2DMatrix(r_xi(1), r_eta(1), r_xi(2), r_eta(2))
            CC%GridJacobian(iCenterX, iCenterY, 1) = detJ
            CC%GridJacobiMatrix(1, 1, iCenterX, iCenterY, 1) = detJ * r_eta(2)   ! xi_x
            CC%GridJacobiMatrix(1, 2, iCenterX, iCenterY, 1) = -detJ * r_eta(1)   ! xi_y
            CC%GridJacobiMatrix(2, 1, iCenterX, iCenterY, 1) = -detJ * r_xi(2)   ! eta_x
            CC%GridJacobiMatrix(2, 2, iCenterX, iCenterY, 1) = detJ * r_xi(1)   ! eta_y
        end do
    end do


    return
end subroutine GetGridJacobian
