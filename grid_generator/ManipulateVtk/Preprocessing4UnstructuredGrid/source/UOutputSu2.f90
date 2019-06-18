!***********************************/
!	Name:非構造格子を計算用に整理したUnStrGridを出力するプログラム
!	Alias:UOutputUnStrGrid
!	Description:UnStrGrid(通称USG:ウサギ)を出力する．一応テキストで出力してるから読めるけど専用のreadプログラムが必要
!	Type:UnstructuredGrid
!	Input:UG%CD,UG%GI,UG%Tri,UG%VC,UG%Line,UG%GM (xxxは使用する種類のセルデータ)
!	Output:UnStrGrid
!	Note:拡張子は4文字でmayuにしたので，読み取りプログラムは*.mayuを探せばいいんじゃない？(適当)
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:2018.02.03
!	Other:GJKアルゴリズムによる接触判定を実行するため凸包出力を追加/複数の境界条件へ対応
!***********************************/
subroutine UOutputSu2(UG, cPath, cFileName)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    character(len=256), intent(inout) :: cFileName, cPath
    integer :: iInner, iOuter
    integer, allocatable :: vEdge_farfield(:), vEdge_obj(:)
    !write(6,*) "Please input filename of the Computed Grid Data"
    !read(5,*) cFileName
    cFileName = trim(adjustl(cPath))//trim(adjustl(cFileName))//".su2" !調べたら.usgも.ccdもなんか既に使われてるみたいだったのでキレた
    !cFileName = "UnStrGrid"

!点番号は1始まり
    open(unit=1,file=trim(adjustl(cFileName)),status='unknown')
!格子の基本情報
        write(1,"('%')")
        write(1,"('% Problem dimension')")
        write(1,"('%')")
        write(1,"('NDIME= ',(1x,i1))") 2
        write(1,"('%')")
        write(1,"('% Inner element connectivity')")
        write(1,"('%')")
        write(1,"('NELEM= ',(1x,i7))") UG%GI%RealCells

        do iCell=1, UG%GI%RealCells
            write(1,"(4(1x,i7))") 5,(UG%Tri%Point(iCell,iLoop)-1,iLoop=1,3) !Adjusted Point Number
        end do

        write(1,"('%')")
        write(1,"('% Node coordinates')")
        write(1,"('%')")
        write(1,"('NPOIN= ',(1x,i7))") UG%GI%Points

        do iPoint=1, UG%GI%Points
            write(1,*) UG%CD%Point(iPoint,1), UG%CD%Point(iPoint,2), iPoint-1
        end do

        write(1,"('%')")
        write(1,"('% Boundary elements')")
        write(1,"('%')")

        allocate(vEdge_farfield(UG%GI%OutlineCells))    !from UMarkingVirtualCell
        allocate(vEdge_obj((UG%GI%AllCells - UG%GI%RealCells) - UG%GI%OutlineCells))

        iInner = 1
        iOuter = 1
        do iCell = UG%GI%RealCells + 1, UG%GI%AllCells
            if(UG%VC%Type(iCell) == 1) then !from UMarkingVirtualCell
                vEdge_farfield(iInner) = UG%VC%Edge(iCell)  !from UMakeEdgeNumber
                iInner = iInner + 1
            else if(UG%VC%Type(iCell) == 2) then
                vEdge_obj(iOuter) = UG%VC%Edge(iCell)
                iOuter = iOuter +1
            end if
        end do


        write(1,"('NMARK = 2')")
        write(1,"('MARKER_TAG= airfoil')")
        write(1,"('MARKER_ELEMS= ',(1x,i7))") UG%GM%BC%iWallTotal
        do iInner = 1, UG%GM%BC%iWallTotal
            write(1,"(3(1x,i7))") 3, UG%Line%Point(vEdge_farfield(iInner), 1), UG%Line%Point(vEdge_farfield(iInner), 2)
        end do

        write(1,"('MARKER_TAG= farfield')")
        write(1,"('MARKER_ELEMS= ',(1x,i7))") UG%GI%OutlineCells
        do iOuter = 1, UG%GI%OutlineCells
            write(1,"(3(1x,i7))") 3, UG%Line%Point(vEdge_obj(iOuter), 1), UG%Line%Point(vEdge_obj(iOuter), 2)
        end do

    close(1)

    return
end subroutine UOutputSu2
