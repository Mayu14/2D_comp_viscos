!***********************************/
!	Name:
!	Alias:GetGlobalNumberOfEndPoint
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
    function GetGlobalNumber2D(Coords,Width,Bound,xCellNumber) result(iGlobalNum)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    double precision, intent(in) :: Coords(:), Width(:), Bound(:)
    integer, intent(in) :: xCellNumber
    integer :: iGlobalNum
    integer, allocatable :: GridNum(:)

    allocate(GridNum(2))

        GridNum(1:2) = int((Coords(1:2) - Bound(1:2))/(Width(1:2)))+1 !補助格子におけるセル中心位置

        iGlobalNum = (GridNum(2)-1)*xCellNumber + GridNum(1)

    return
    end function GetGlobalNumber2D
