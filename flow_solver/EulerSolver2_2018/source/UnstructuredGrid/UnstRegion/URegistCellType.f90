!***********************************/
!	Name:非構造格子Boundaryセルデータをregistするプログラム
!	Alias:URegistCellType
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.09
!	Update:2018.01.20
!	Other:
!***********************************/
subroutine URegistCellType(UG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    do iCell=UG%GI%RealCells+1, UG%GI%AllCells
        if(UG%GM%CellType(iCell,1,1) == 4) then !重合格子境界なら
            UG%GM%CellType(UG%VC%Cell(iCell,1),1,1) = 99 ! OVERSET BOUNDARY CELL

        else
            UG%GM%CellType(UG%VC%Cell(iCell,1),1,1) = 100  !INTERNAL BOUNDARY CELL (REAL)
        end if
    end do

return
end subroutine URegistCellType
