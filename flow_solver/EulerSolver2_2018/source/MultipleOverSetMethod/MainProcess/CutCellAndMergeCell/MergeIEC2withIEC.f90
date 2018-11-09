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
    subroutine MergeIEC2withIEC(GG)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    !integer, intent(in) :: iLocalEdge !of Global Cell (=1,2,3,4)
    !integer, intent(in) :: iCenterX, iCenterY
    !type(Geometry), intent(in) :: Geom
    type(GlobalGrid), intent(inout) :: GG

    type(InternalEdgeOfCell), pointer :: TempIEC
    integer :: iCell

    do iCell = 1, GG%GASG%TotalCell

        TempIEC => GG%CutC(iCell)%IEC

        do while(associated(TempIEC%Next))
            TempIEC => TempIEC%Next
        end do

            TempIEC%Next => GG%CutC(iCell)%IEC2
    end do


        return
    end subroutine MergeIEC2withIEC
