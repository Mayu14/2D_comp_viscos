!***********************************/
!	Name:界面を挟んで中心差分を取るための準備をするためのプログラム
!	Alias:UCentralDifferencePrepareAroundFace
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.11.18
!	Update:
!	Other:
!***********************************/
subroutine UCentralDifferencePrepareAroundFace(UG, iFrontCell, iFrontLocalEdge, iBackCell, iBackLocalEdge, length)
    use StructVar_Mod
    use FrequentOperation
    type(UnstructuredGrid), intent(in) :: UG
    integer, intent(out) :: iFrontCell, iFrontLocalEdge, iBackCell, iBackLocalEdge
    double precision, intent(out) :: length

        iFrontCell = UG%Line%Cell(iEdge, 1, 1)
        iFrontLocalEdge = UG%Line%Cell(iEdge, 1, 2)
        iBackCell =  UG%Line%Cell(iEdge, 2, 1)
        iBackLocalEdge =UG%Line%Cell(iEdge, 2, 2)
        length = AbsVector(UG%GM%Width(iFrontCell, iFrontLocalEdge, :)) + AbsVector(UG%GM%Width(iBackCell, iBackLocalEdge))

    return
end subroutine UCentralDifferencePrepareAroundFace
