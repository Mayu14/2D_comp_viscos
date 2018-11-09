!***********************************/
!	Name:
!	Alias:PointSwap
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.24
!	Update:
!	Other:not parallelize
!***********************************/
    subroutine MergeBesideCell(tGCC,dGCC)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    !integer, intent(in) :: iLocalEdge !of Global Cell (=1,2,3,4)
    !integer, intent(in) :: iCenterX, iCenterY
    !type(Geometry), intent(in) :: Geom
    type(CutCell), intent(inout) :: tGCC, dGCC
    !type(InternalPointOfCell)
    type(InternalEdgeOfCell), pointer :: TempIEC
    integer :: iCell


        TempIEC => tGCC%IEC

        do while(associated(TempIEC%Next))
            TempIEC => TempIEC%Next
        end do

            TempIEC%Next => dGCC%IEC



        return
    end subroutine MergeBesideCell
