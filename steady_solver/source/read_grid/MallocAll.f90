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
subroutine MallocAll(Conf, Geom, CC, CE)
    use StructVar_Mod_mk2
    use LoopVar_Mod_mk2
    implicit none

    type(Configulation), intent(in) :: Conf
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CellEdge), intent(inout) :: CE

    integer, allocatable :: iMin(:)

    iDim =Geom%iDimension

    allocate(iMin(3))
    iMin = 0
    if(Geom%CellNumber(2) == 0) iMin(2) = 1
    if(Geom%CellNumber(3) == 0) iMin(3) = 1

    iXmax = Geom%CellNumber(1)
    iYmax = max(1, Geom%CellNumber(2))
    iZmax = max(1, Geom%CellNumber(3))

    ! CellCenter
    allocate(CC%ConservedQuantity(4, 0:iXmax, 0:iYmax, 1))
    allocate(CC%PrimitiveVariable(4, 0:iXmax, 0:iYmax, 1))

    allocate(CE%NormalFluxDiff(4, 0:iXmax+1, 0:iYmax+1, 1, 2))
    CC%ConservedQuantity = 0.0d0
    CC%PrimitiveVariable = 0.0d0
    CE%NormalFluxDiff = 0.0d0

    ! body fitted coordinate
    allocate(CC%GridJacobian(iXmax, iYmax, 1))
    allocate(CC%GridJacobiMatrix(2, 2, iXmax, iYmax, 1))

    allocate(CE%GridJacobian(iXmax+1, iYmax+1, 1, 2))
    allocate(CE%GridJacobiMatrix(2, 2, iXmax+1, iYmax+1, 1, 2))


    if(Conf%UseLocalTimeStep == 1) then
        allocate(CC%TimeWidth(iXmax,iYmax,iZmax))
        allocate(CC%TmpTimeWidth(iMin(1):iXmax+1,iMin(2):iYmax,1,iDim))
        CC%TmpTimeWidth = 0.0d0
    else
        allocate(CC%TimeWidth(1,1,1))
    end if
        CC%TimeWidth = 0.0d0

    return
end subroutine MallocAll
