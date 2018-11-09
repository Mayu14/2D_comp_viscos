!***********************************/
!	Name:影響域分離法のための保存量修正
!	Alias:ModifyQuantity
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.02.26
!	Update:
!	Other:
!***********************************/
    subroutine ModifyQuantity(CC4MB,CC,MG,UG)
    use StructVar_Mod
    implicit none
    type(CellCenter), intent(inout) :: CC4MB
    type(CellCenter), intent(in) :: CC !reference
    type(MoveGrid), intent(in) :: MG
    type(UnstructuredGrid), intent(in) :: UG
    integer :: iCell

    do iCell = 1, UG%GI%RealCells
        if(MG%IS(iCell)%InfluenceDepth == 0) then
            CC4MB%ConservedQuantity(:,iCell,1,1) = CC%ConservedQuantity(:,iCell,1,1)
        end if
    end do

    return
    end subroutine ModifyQuantity
