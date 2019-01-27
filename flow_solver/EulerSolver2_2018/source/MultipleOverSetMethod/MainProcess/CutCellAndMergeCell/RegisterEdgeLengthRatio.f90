!***********************************/
!	Name:
!	Alias:RegisterEdgeLengthRatio
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.31
!	Update:
!	Other:not parallelize
!***********************************/
    subroutine RegisterEdgeLengthRatio(CutE,GG)
    use StructVar_Mod

    implicit none

    type(GlobalGrid), intent(inout) :: GG
    type(CuttedEdge), intent(in) :: CutE
    type(InternalEdgeOfCell), pointer :: TempIEC
    type(InternalPointOfEdge), pointer :: TempIPE
    double precision, allocatable :: Temp(:), MidPoint(:),TempCooods(:,:)

    integer, allocatable :: iTemp(:), iInvTemp(:)
    integer :: iLoop,iCount, iSortNum,iEndFlag, iGlobalNumber, GetGlobalNumber2D

    allocate(Temp(CutE%InternalPointNumber))
    allocate(iTemp(CutE%InternalPointNumber),iInvTemp(CutE%InternalPointNumber))
    allocate(MidPoint(3),TempCooods(CutE%InternalPointNumber,3))
!1st step :: data pool
    TempIPE => CutE%IPE
        iCount = 0
        do while(associated(TempIPE)) !次のレコードがなくなるまで続行
            iCount = iCount + 1
            !Temp(iCount) = TempIEC%CutRatio !Length from one of the end point
            TempCooods(iCount,1:3) = TempIPE%Coords(1:3)
            Temp(iCount) = TempIPE%Distance !Length from one of the end point
            iTemp(iCount) = iCount
            TempIPE => TempIPE%Next
        end do


!2nd step :: data sort by length
    iEndFlag = 1 !initialize 1
    iSortNum = 0 !Initialize
    do while(iEndFlag == 0)
        iEndFlag = 0 !initialize 2
        iSortNum = iSortNum + 1

        do iLoop=1, iCount-iSortNum !bubble sort is STABLE sort
            if(Temp(iLoop) > Temp(iLoop+1)) then
                call swap(Temp(iLoop),Temp(iLoop+1))
                call iswap(iTemp(iLoop),iTemp(iLoop+1))
                iEndFlag = 1
            end if
        end do
        if(iEndFlag == 1 .and. iSortNum == iCount-1) write(6,*) "Sorting Error"
    end do

    call MakeInverseOfiTemp(iTemp,iInvTemp, CutE%InternalPointNumber)

!3rd step :: Calculate Edge cut ratio
    TempIPE => CutE%IPE


    iCount = 0
    do while(associated(TempIPE))
        iCount = iCount + 1
        allocate(TempIEC)
        TempIEC%EdgeNum = CutE%EdgeNumber
        TempIEC%GridNum = CutE%GridNumber

        if(iTemp(iCount) == 1) then !shortest edge
            TempIEC%CutRatio = TempIPE%Distance/CutE%Length
            MidPoint(:) = 0.5d0*(TempCooods(iTemp(iCount),1:3)+CutE%EndPoint(1,1:3))

        else if(iTemp(iCount) == CutE%InternalPointNumber) then !longest edge
            TempIEC%CutRatio = (CutE%Length - TempIPE%Distance)/CutE%Length
            MidPoint(:) = 0.5d0*(CutE%EndPoint(2,1:3) + TempCooods(iTemp(iCount),1:3))

        else !not shortest and longest edge
            TempIEC%CutRatio = (Temp(iCount+1) - Temp(iCount))/CutE%Length
            MidPoint(:) = 0.5d0*(TempCooods(iInvTemp(iTemp(iCount+1)),1:3) + TempCooods(iInvTemp(iTemp(iCount)),1:3))

        end if


        !Regist
        iGlobalNumber = GetGlobalNumber2D(MidPoint(1:2),GG%GASG%Width(1:2),GG%GASG%Bound(1,1:2),GG%GASG%CellNumber(1))

        if(associated(GG%CutC(iGlobalNumber)%IEC)) then
            TempIEC%Next => GG%CutC(iGlobalNumber)%IEC
        else
            nullify(GG%CutC(iGlobalNumber)%IEC%Next)
        end if
        GG%CutC(iGlobalNumber)%IEC%Next => TempIEC

        TempIPE => TempIPE%Next !next node
    end do

    return
contains

    subroutine MakeInverseOfiTemp(iTemp,iInvTemp,TotalPoint)
    implicit none
        integer, intent(in) :: TotalPoint
        integer, intent(in) :: iTemp(:)
        integer, intent(inout) :: iInvTemp(:)
        integer :: iLoop

        do iLoop=1, TotalPoint
            iInvTemp(iTemp(iLoop)) = iLoop
        end do

    return
    end subroutine MakeInverseOfiTemp

    subroutine swap(A,B)
    implicit none
        double precision, intent(inout) :: A,B
        double precision :: Tmp
        Tmp = A
        A = B
        B = Tmp
        return
    end subroutine swap

    subroutine iswap(A,B)
    implicit none
        integer, intent(inout) :: A,B
        integer :: Tmp
        Tmp = A
        A = B
        B = Tmp
        return
    end subroutine iswap

    end subroutine RegisterEdgeLengthRatio

