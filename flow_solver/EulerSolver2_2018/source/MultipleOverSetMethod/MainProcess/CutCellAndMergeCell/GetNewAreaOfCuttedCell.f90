!***********************************/
!	Name:
!	Alias:GetNewAreaOfCuttedCell
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
    subroutine GetNewAreaOfCuttedCell(GCC)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    type(CutCell), intent(inout) :: GCC

    double precision, allocatable :: PointCoords(:,:)
    integer, allocatable :: iOnOff(:)
    type(InternalPointOfCell), pointer :: TempIPC

    double precision :: Error,dTmp, Area
    integer :: iCount, iYminNumber
    integer :: iTmp, iEndFlag, iEnd1, iEnd2
    allocate(PointCoords(GCC%NumberOfPoint,1:5)) !1:3:xyz,0:LengthFrom1stPoint,4:LengthFromOnePoint, 5:AngleByOneAxis
    allocate(iOnOff(GCC%NumberOfPoint+1))

    !Read Record
    TempIPC => GCC%IPC
    iCount = 0
    do while(associated(TempIPC))
        iCount = iCount + 1
        PointCoords(iCount,1:3) = TempIPC%Coords
        TempIPC => TempIPC%Next

    end do

    !Find Minimum Value of Y
    iYminNumber = minloc(PointCoords(:,2),1)
    !Coordinate Transform to Relative Coordinate
    do iCount=1, GCC%NumberOfPoint
        if(iCount /= iYminNumber) then
            PointCoords(iCount,1:3) = PointCoords(iCount,1:3)-PointCoords(iYminNumber,1:3)
        end if
    end do

    PointCoords(iYminNumber,1:5) = 0

    !Compute Length and Angle
    do iCount = 1, GCC%NumberOfPoint
        if(iCount /= iYminNumber) then
            PointCoords(iCount,4) = sqrt(dot_product(PointCoords(iCount,1:3),PointCoords(iCount,1:3)))
            PointCoords(iCount,5) = atan2(PointCoords(iCount,2),PointCoords(iCount,1))
        end if
    end do

    !Sort By Length and Angle
    !Bubble Sort
    iEndFlag = 0
    do while(iEndFlag == 0)
        iEndFlag = 1
        do iCount=1, GCC%NumberOfPoint-1
            if(PointCoords(iCount,5) > PointCoords(iCount+1,5)) then
                call SwapDbleArray(PointCoords(iCount,1:5),PointCoords(iCount+1,1:5))
                iEndFlag = 0
            end if
        end do
    end do

    !Make OnOff Switch
    Error = 10**(-12) !should change
    iOnOff = 1 !Initialize

    do iCount=2, GCC%NumberOfPoint
        if(abs(PointCoords(iCount,4) - PointCoords(iCount-1,4)) < Error ) then
            iOnOff(iCount-1) = 0
        end if
    end do

    if(iOnOff(1) /= 0) then
        if(abs(PointCoords(1,4) - PointCoords(iOnOff(GCC%NumberOfPoint),4)) < Error ) then
            iOnOff(GCC%NumberOfPoint) = 0
        end if
    end if

!Calculate Area
    GCC%CutVolume = 0.0d0
    do iCount=1, GCC%NumberOfPoint
        if(iOnOff(iCount) == 0) cycle
        iEnd1 = iCount !1st End Point

        do iTmp=iCount+1, GCC%NumberOfPoint
            if(iOnOff(iTmp) == 0) cycle
            iEnd2 = iTmp !2nd End Point
            exit
        end do

        GCC%CutVolume = GCC%CutVolume + AbsCrossProduct(PointCoords(iEnd1,1:3),PointCoords(iEnd2,1:3))
    end do

    if(iOnOff(1) /= 0 .and. iOnOff(GCC%NumberOfPoint) /= 0) then
        GCC%CutVolume = GCC%CutVolume + AbsCrossProduct(PointCoords(GCC%NumberOfPoint,1:3),PointCoords(1,1:3))
    end if

        GCC%CutVolume = 0.5d0*GCC%CutVolume

        return
    end subroutine GetNewAreaOfCuttedCell
