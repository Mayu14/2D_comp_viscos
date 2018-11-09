!***********************************/
!	Name:GJKアルゴリズムのための，凸包構成点列生成プログラム
!	Alias:UMakeConvexHull
!	Description:仮想格子が接する界面の構成点を列挙して重複を削除する
!	Type:UnstructuredGrid
!	Input:UG%GM%Dimension,UG%CD%Edge,UG%CD%Cell (xxxは使用する種類のセルデータ)
!	Output:UG%GM%Normal
!	Note:格子点数が少ないときは度数分布ソートを用いる．ある程度以上格子点数が多い場合はクイックソート
!	Author:Akitaka Toyota
!	Date:2017.12.19
!	Update:2018.01.20
!	Other:現時点でハイブリッド格子・3次元格子への対応は未定(準備はある)
!***********************************/
subroutine UMakeConvexHull(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG
    integer :: iKind !1:OuterBoundary, 2:Surface
    iKind = 1

    if(UG%GM%Dimension == 2) then
        if(UG%GI%Points < 10**5) then
            call CountingSort
        else
            call QuickSortProcess(UG%GI%OutlineCells)
        end if
    else
        write(6,*) "coming soon..."
        stop
    end if
return

contains
    subroutine CountingSort
    implicit none
    integer, allocatable :: iCount(:) !点番号：出現回数を返す
    integer :: iTotal = 0

        allocate(iCount(UG%GI%Points))
        iCount = 0

        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            if(UG%VC%Type(iCell) == iKind) then !ここで外周と内部物体切り替え可能
                do iLoop=1,2
                    iCount(UG%Line%Point(UG%VC%Edge(iCell),iLoop)) = iCount(UG%Line%Point(UG%VC%Edge(iCell),iLoop)) + 1
                    if(iCount(UG%Line%Point(UG%VC%Edge(iCell),iLoop)) == 1) then
                        iTotal = iTotal + 1
                    end if
                end do

            end if
        end do

            UG%CH%iTotal = iTotal
            allocate(UG%CH%PointNum(iTotal))

            iTotal = 0
        do iPoint=1, UG%GI%Points
            if(iCount(iPoint) /= 0) then
                iTotal = iTotal + 1
                UG%CH%PointNum(iTotal) = iPoint
            end if
        end do

        return
    end subroutine CountingSort

    subroutine QuickSortProcess(iBoundaryCellNumber)
    implicit none
    double precision, allocatable :: DuplicateConvexPoint(:)
    integer, intent(in) :: iBoundaryCellNumber
    integer :: iTmp,iTotal

        allocate(DuplicateConvexPoint(2*(iBoundaryCellNumber)))

        do iCell=1, UG%VC%Total
            if(UG%VC%Type(iCell) == iKind) then
                DuplicateConvexPoint(2*iCell-1) = UG%Line%Point(UG%VC%Edge(iCell+UG%GI%RealCells),1)
                DuplicateConvexPoint(2*iCell) = UG%Line%Point(UG%VC%Edge(iCell+UG%GI%RealCells),2)
            end if
        end do

        call QuickSortMain(DuplicateConvexPoint,1,2*(UG%GI%OutlineCells))

        UG%CH%iTotal = UG%GI%OutlineCells
        allocate(UG%CH%PointNum(UG%CH%iTotal))

        iTotal = 0
        do iPoint=1, 2*(UG%GI%OutlineCells)
            if(iPoint == 1) then
                iTmp = DuplicateConvexPoint(iPoint)
                iTotal = iTotal + 1
                UG%CH%PointNum(iTotal) = iTmp
            else
                if(iTmp /= DuplicateConvexPoint(iPoint)) then
                    iTmp = DuplicateConvexPoint(iPoint)
                    iTotal = iTotal + 1
                    UG%CH%PointNum(iTotal) = iTmp
                end if
            end if
        end do
        UG%CH%iTotal = iTotal

    return
    end subroutine


    recursive subroutine QuickSortMain(a,iFirst, iLast)
    implicit none
    double precision, intent(inout) :: a(:)
    integer, intent(in) :: iFirst, iLast
    double precision MiddlePointValue, Temp
    integer iLeft, iRight

        MiddlePointValue = a( (iFirst+iLast) / 2 )
        iLeft = iFirst
        iRight = iLast
        do
            do while (a(iLeft) < MiddlePointValue)
                iLeft = iLeft + 1
            end do

            do while (MiddlePointValue < a(iRight))
                iRight = iRight - 1
            end do

            if (iLeft >= iRight) exit

            Temp = a(iLeft);  a(iLeft) = a(iRight);  a(iRight) = Temp
            iLeft=iLeft+1
            iRight=iRight-1
        end do

        if (iFirst < iLeft-1) call QuickSortMain(a,iFirst, iLeft-1)
        if (iRight+1 < iLast)  call QuickSortMain(a,iRight+1, iLast)
    end subroutine QuickSortMain

end subroutine UMakeConvexHull
