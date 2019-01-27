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
    subroutine PointSwap(CutE)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    type(CuttedEdge), intent(inout) :: CutE
    double precision, allocatable :: Temp(:)
    integer :: iTemp
    allocate(Temp(3))

    if(CutE%GlobalNumber1 > CutE%GlobalNumber2) then !界面が単一のGセル内に入っている
        Temp(:) = CutE%EndPoint(1,:)
        CutE%EndPoint(1,:) = CutE%EndPoint(2,:)
        CutE%EndPoint(2,:) = CutE%EndPoint(1,:)

        iTemp = CutE%GlobalNumber1
        CutE%GlobalNumber1 = CutE%GlobalNumber2
        CutE%GlobalNumber2 = iTemp
    end if

    return
    end subroutine PointSwap
