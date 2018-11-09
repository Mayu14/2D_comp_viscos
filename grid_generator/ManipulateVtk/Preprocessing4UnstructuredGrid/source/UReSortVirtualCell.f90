!***********************************/
!	Name:非構造格子用EulerSolverのための仮想格子並び換えプログラム
!	Alias:UReSortVirtualCell
!	Description:仮想セルの中心座標を基準に仮想セルの番号を振り直す(外部→内部)
!	Type:UnstructuredGrid
!	Input:UG%GM%Dimension,UG%xxx%Point (xxxは使用する種類のセルデータ)
!	Output:UG
!	Note:中心からの距離を算出して，その値の大きい順に番号を振り直す
!	Author:Akitaka Toyota
!	Date:2018.01.20
!	Update:-
!	Other:
!***********************************/
subroutine UReSortVirtualCell(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    integer, allocatable :: iO2N(:),iN2O(:),iTmp(:) !Old Number to New Number, New Number to Old Number

    double precision, allocatable :: Length(:) !中心からの距離

    allocate(iO2N(UG%GI%RealCells+1:UG%GI%AllCells),iTmp(UG%GI%RealCells+1:UG%GI%AllCells),Length(UG%GI%RealCells+1:UG%GI%AllCells))

    call GetLength !比較基準用の距離取得
    call MakeReSortRelation !元の番号と並び換え後の番号を比較するための交換用配列作成
    !call MakeiN2O
    call ReSortAndDataUpdate


return
contains

    subroutine GetLength
    implicit none

    do iCell=UG%GI%RealCells+1, UG%GI%AllCells
        iO2N(iCell) = iCell
        Length(iCell) = dot_product(UG%CD%Cell(iCell,:),UG%CD%Cell(iCell,:)) !距離の2乗 !どうせ距離は正の数だから２乗で比較しても問題なし(処理も軽くなるし)
    end do

    return
    end subroutine GetLength

    subroutine MakeReSortRelation !奇偶転置ソートを用いる !複雑なソートを書く気力が残っていなかったため
    implicit none
    integer :: iDebug = -1
    integer :: iFlag = 1
    !距離が遠い順に並べ直す
    do while(iFlag == 1)
        iFlag = 0
        iDebug = iDebug + 1
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells-1, 2
            if(Length(iCell) < Length(iCell+1)) then
                call dSwap(Length(iCell),Length(iCell+1))
                call iSwap(iO2N(iCell),iO2N(iCell+1))
                iFlag = 1
            end if
        end do

        do iCell=UG%GI%RealCells+2, UG%GI%AllCells-1,2
            if(Length(iCell) < Length(iCell+1)) then
                call dSwap(Length(iCell),Length(iCell+1))
                call iSwap(iO2N(iCell),iO2N(iCell+1))
                iFlag = 1
            end if
        end do
    end do

    return
    end subroutine MakeReSortRelation


    subroutine iSwap(A,B)
    implicit none
    integer, intent(inout) :: A,B
    integer :: Tmp

        Tmp = B
        B = A
        A = Tmp

    return
    end subroutine iSwap

    subroutine dSwap(A,B)
    implicit none
    double precision, intent(inout) :: A,B
    double precision :: Tmp

        Tmp = B
        B = A
        A = Tmp

    return
    end subroutine dSwap

    subroutine MakeiN2O
    implicit none

    allocate(iN2O(UG%GI%RealCells+1:UG%GI%AllCells))

    do iCell=UG%GI%RealCells+1, UG%GI%AllCells
        iN2O(iO2N(iCell)) = iCell
    end do

    return
    end subroutine MakeiN2O


    subroutine ResortAndDataUpdate
    implicit none
    integer, allocatable :: ReferencedLocalEdge(:)
    integer, allocatable :: CopyOfTriCell(:,:), CopyOfLineCell(:,:,:)
    integer, allocatable :: CopyOfVCEdge(:),CopyOfVCCell(:,:)

    allocate(ReferencedLocalEdge(UG%GI%RealCells+1:UG%GI%AllCells))
    allocate(CopyOfTriCell(UG%GI%RealCells,3),CopyOfLineCell(UG%GI%Edges,2,2))
    allocate(CopyOfVCEdge(UG%GI%RealCells+1:UG%GI%AllCells),CopyOfVCCell(UG%GI%RealCells+1:UG%GI%AllCells,2))

!stop
!UG%Tri%Cellのupdate
    CopyOfTriCell = UG%Tri%Cell
    do iCell=1, UG%GI%RealCells
        do iLocalEdge=1,3
            if(CopyOfTriCell(iCell,iLocalEdge) > UG%GI%RealCells) then
                UG%Tri%Cell(iCell,iLocalEdge) = iO2N(CopyOfTriCell(iCell,iLocalEdge))
            end if
        end do
    end do

!UG%Line%Cellのupdate
    CopyOfLineCell = UG%Line%Cell
    do iEdge=1, UG%GI%Edges
        do iSide=1,2
            if(CopyOfLineCell(iEdge,iSide,1) > UG%GI%RealCells) then
                ReferencedLocalEdge(CopyOfLineCell(iEdge,iSide,1)) = CopyOfLineCell(iEdge,iSide,2)
            end if
        end do
    end do

    do iEdge=1, UG%GI%Edges
        do iSide=1,2
            if(CopyOfLineCell(iEdge,iSide,1) > UG%GI%RealCells) then
                UG%Line%Cell(iEdge,iSide,1) = iO2N(CopyOfLineCell(iEdge,iSide,1))
                UG%Line%Cell(iEdge,iSide,2) = ReferencedLocalEdge(iO2N(CopyOfLineCell(iEdge,iSide,1)))
            end if
        end do
    end do

!UG%VC%Edgeのupdate
    CopyOfVCEdge = UG%VC%Edge
    do iCell = UG%GI%RealCells+1, UG%GI%AllCells
        UG%VC%Edge(iCell) = CopyOfVCEdge(iO2N(iCell))
    end do

!UG%VC%Cellのupdate
    CopyOfVCCell = UG%VC%Cell
    do iCell = UG%GI%RealCells+1, UG%GI%AllCells
        UG%VC%Cell(iCell,:) = CopyOfVCCell(iO2N(iCell),:)
    end do

    return
    end subroutine ResortAndDataUpdate

end subroutine UReSortVirtualCell
