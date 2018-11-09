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
    subroutine RegisterCutCondition(iCrossCondition,LCE,GCC)!(Coords,Width,Bound,xCellNumber) result(iGlobalNum)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    type(CuttedEdge), intent(inout) :: LCE
    type(CutCell), intent(inout) :: GCC

    integer, intent(in) :: iCrossCondition
    integer :: iLoop
    type(InternalPointOfEdge), pointer :: TempIPE
    type(InternalPointOfCell), pointer :: TempIPC
    type(InternalEdgeOfCell), pointer :: TempIEC

    if(iCrossCondition == 0) then !界面が単一のGセル内に入っている
        nullify(LCE%IPE)

        allocate(TempIEC)
        TempIEC%CutRatio = 1.0d0
        TempIEC%EdgeNum = LCE%EdgeNumber
        TempIEC%GridNum = LCE%GridNumber
        nullify(TempIEC%Next)
        GCC%IEC => TempIEC

        do iLoop=1, 2
            allocate(TempIPC)
            TempIPC%Coords(1:3) = LCE%EndPoint(iLoop,1:3)
            if(GCC%UnsolvedFlag == 0) then
                nullify(TempIPC%Next)
                GCC%UnsolvedFlag = 1
            else
                TempIPC%Next => GCC%IPC
            end if
            GCC%IPC => TempIPC
        end do
        GCC%NumberOfEdge = GCC%NumberOfEdge + 1
        GCC%NumberOfPoint = GCC%NumberOfPoint + 2 !allow overlap


    else !界面が複数のGセルに跨る
        allocate(TempIPE)
        allocate(TempIPC)

        LCE%IPE%Coords(:) = LCE%TmpInternalPoint(:)
        TempIPC%Coords(:) = LCE%TmpInternalPoint(:)

        TempIPE%Distance = sqrt(dot_product(LCE%TmpInternalPoint(1:3)-LCE%EndPoint(1,1:3),LCE%TmpInternalPoint(1:3)-LCE%EndPoint(1,1:3))) !This is still not ratio here.

        if(LCE%InternalPointNumber == 1) then
            nullify(TempIPE%Next)
        else
            TempIPE%Next => LCE%IPE
        end if
        LCE%IPE%Next => TempIPE

        if(GCC%UnsolvedFlag == 0) then
            nullify(TempIPC%Next)
            GCC%UnsolvedFlag = 1

        else
            TempIPC%Next => GCC%IPC
        end if

        GCC%IPC%Next => TempIPC

        GCC%NumberOfEdge = GCC%NumberOfEdge + 1
        GCC%NumberOfPoint = GCC%NumberOfPoint + 1 !allow overlap

        !call GetAlivePoint

    end if

    GCC%UseCutCell = 1

    return
    end subroutine RegisterCutCondition
