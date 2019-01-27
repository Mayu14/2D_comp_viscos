!***********************************/
!	Name:
!	Alias: UCalcInscribedCircleOfCell
!	Description:
!	Type:UnstructuredGrid
!	Input:Volume and Area
!	Output:Radius
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.02.03
!	Update:-
!	Other:現時点でハイブリッド格子・3次元格子への対応は未定(準備はある)
!***********************************/
subroutine UCalcInscribedCircleOfCell(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    do iCell=1, UG%GI%RealCells
        UG%InscribedCircle(iCell) = 2.0d0/3.0d0*UG%GM%Volume(iCell)/(UG%GM%AverageWidth(iCell)) !AverageWidth = a+b+c !Three Edge Length
    end do

    return
end subroutine  UCalcInscribedCircleOfCell
