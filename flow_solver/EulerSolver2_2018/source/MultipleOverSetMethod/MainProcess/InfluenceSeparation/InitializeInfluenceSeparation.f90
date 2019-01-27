!***********************************/
!	Name:影響域分離法のための配列確保
!	Alias:InitializeInfluenceSeparation
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.31
!	Update:
!	Other:
!***********************************/
    subroutine InitializeInfluenceSeparation(UG,MG,CC,CC4MB)
    use StructVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(MoveGrid), intent(inout) :: MG
    type(CellCenter), intent(in) :: CC
    type(CellCenter), intent(inout) :: CC4MB

    integer :: iCell
    CC4MB = CC
    CC4MB%InterpolatedQuantity = CC%ConservedQuantity

    do iCell = 1, UG%GI%RealCells
        MG%IS(iCell)%InfluenceDepth = 0
        MG%IS(iCell)%BoundaryAdjacent = 0
    end do

    MG%MaxInfluenceRange = 0.0d0

    return
    end subroutine InitializeInfluenceSeparation
