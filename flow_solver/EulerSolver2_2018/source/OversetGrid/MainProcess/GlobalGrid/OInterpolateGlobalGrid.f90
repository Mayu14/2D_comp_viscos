!***********************************/
!	Name:移動先の格子を補間するための補間データ提供元格子検索
!	Alias:OInterpolatedLocalGrid
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
subroutine OInterpolateGlobalGrid(Geom,OSG,MG,UCC,CC)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none

    type(Geometry), intent(inout) :: Geom
    type(OverSetGrid), intent(in) :: OSG
    type(MoveGrid), intent(in) :: MG
    type(CellCenter), intent(in) :: UCC
    type(CellCenter), intent(inout) :: CC

    double precision, allocatable :: SCellCoordinate(:),Moment(:)
    double precision :: InverseDistanceSum, InverseDistance
    integer :: iContent
    integer :: iGCell
    allocate(SCellCoordinate(3),Moment(5))

    iXmax = max(1,Geom%CellNumber(1))
    iYmax = max(1,Geom%CellNumber(2))
    iZmax = max(1,Geom%CellNumber(3))

    Geom%Interpolated= 0

    do iGCell=1, iXmax*iYmax*iZmax !Global格子のすべてのセルについて
        if(OSG%FrequencyDistribution(iGCell) /= 0) then !if the cell need to interpolate
            InverseDistanceSum = 0.0d0
            Moment = 0.0d0

            !この構造格子の座標は
            iCenterX = iGCell - int(iGCell/iXmax)*iXmax
            iCenterY = int(iGCell/iXmax)
            SCellCoordinate(1) = Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * (dble(iCenterX)-1.0d0)
            SCellCoordinate(2) = Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * (dble(iCenterY)-1.0d0)
            SCellCoordinate(3) = 0.0d0

            do iContent=1, OSG%FrequencyDistribution(iGCell)
                iAdjacentCell = OSG%RelatedS2U(iGCell,iContent)

                InverseDistance = 1.0d0/sqrt(sum((MG%NSR%C%NextCoordsInG(iAdjacentCell,1:3) - SCellCoordinate(1:3))**2))

                InverseDistanceSum = InverseDistanceSum + InverseDistance
                !(距離の逆数)*(物理量)!!(距離) 近いほど影響が大きい
                Moment = Moment + InverseDistance * UCC%ConservedQuantity(:,iAdjacentCell,1,1)
            end do

            CC%ConservedQuantity(1:3,iCenterX,iCenterY,1) = Moment(1:3)/InverseDistanceSum
            CC%ConservedQuantity(4,iCenterX,iCenterY,1) = Moment(5)/InverseDistanceSum
            Geom%Interpolated(iCenterX,iCenterY,1) = 1
        end if
    end do

return
end subroutine OInterpolateGlobalGrid
