!***********************************/
!	Name:物体まわりの点を並び変えて格納する
!	Alias:UMakeInternalObject
!	Description:
!	Type:UnstructuredGrid
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.02.08
!	Update:
!	Other:
!***********************************/
subroutine UMakeInternalObject(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG
    integer, allocatable :: Point2(:),Point(:)
    double precision, allocatable :: Angle(:),Norm(:)
    double precision, allocatable :: RelativeVector(:)
    double precision :: ErrorNorm, MinYCD, LastMinYCD
    integer :: iCount, iProcess, iEndFlag, iCheckMin

    iCount = 0
    do iCell = UG%GI%RealCells+1, UG%GI%AllCells
        if(UG%VC%Type(iCell) == 2) then
            iCount = iCount + 1
        end if
    end do

    allocate(Point(iCount))
    allocate(Point2(iCount*2))

    iCount = 0
    do iCell = UG%GI%RealCells+1, UG%GI%AllCells
        if(UG%VC%Type(iCell) == 2) then
            do iLoop=1,2
                iCount = iCount + 1
                Point2(iCount) = UG%Line%Point(UG%VC%Edge(iCell),iLoop)

            end do
        end if
    end do


    iEndFlag = 0
    iProcess = 0
    do while(iEndFlag == 0)
        iEndFlag = 1
        iProcess = iProcess + 1
        do iPoint=1, ubound(Point2,1)-iProcess

            if(Point2(iPoint) > Point2(iPoint+1)) then
                iEndFlag = 0
                call SwapInteger(Point2(iPoint),Point2(iPoint+1))
            end if
        end do
    end do

    MinYCD = 10**30

    do iPoint = 1, ubound(Point,1)
        Point(iPoint) = Point2(2*iPoint)

        LastMinYCD = MinYCD
        MinYCD = min(MinYCD,UG%CD%Point(Point(iPoint),2))
        if(LastMinYCD /= MinYCD) iCheckMin = Point(iPoint)
    end do
        if(LastMinYCD /= MinYCD) iCheckMin = Point(ubound(Point,1)) !Final Step

    deallocate(Point2)
    allocate(Angle(ubound(Point,1)),Norm(ubound(Point,1)))
    allocate(RelativeVector(3))




    do iPoint = 1, ubound(Point,1)
        if(iPoint /= iCheckMin) then
            RelativeVector = UG%CD%Point(Point(iPoint),1:3) - UG%CD%Point(iCheckMin,1:3)
            Angle(iPoint) = atan2(RelativeVector(2),RelativeVector(1))
            Norm(iPoint) = dot_product(RelativeVector(1:3),RelativeVector(1:3))

        else
            Angle(iPoint) = 0.0d0
            Norm(iPoint) = 0.0d0
        end if
    end do

    ErrorNorm = 10**(-10)

    iEndFlag = 0
    iProcess = 0
    do while(iEndFlag == 0)
        iEndFlag = 1
        iProcess = iProcess + 1
        do iPoint=1, ubound(Point,1)-iProcess
            if(Angle(iPoint)-ErrorNorm < Angle(iPoint) .and.  Angle(iPoint) < Angle(iPoint)+ErrorNorm) then
                if(Norm(iPoint) > Norm(iPoint + 1)) then
                    iEndFlag = 0
                    call swapAll

                end if

            else if(Angle(iPoint) < Angle(iPoint+1)) then
                iEndFlag = 0
                call swapAll
            end if
        end do
    end do

    allocate(UG%IO%PointNum(ubound(Point,1)))

    UG%IO%iTotal = ubound(Point,1)
    UG%IO%PointNum = Point

    return
contains
    subroutine swapAll
    implicit none
        call SwapInteger(Point(iPoint),Point(iPoint+1))
        call SwapDoublePrecision(Angle(iPoint),Angle(iPoint+1))
        call SwapDoublePrecision(Norm(iPoint),Norm(iPoint+1))
        return
    end subroutine swapAll

end subroutine UMakeInternalObject
