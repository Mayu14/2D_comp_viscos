!***********************************/
!	Name:非構造格子用EulerSolverのための仮想格子属性割り当てプログラム
!	Alias:UMarkingVirtualCell
!	Description:ユーザー指定基準で仮想セルに境界番号を割り当てる
!	Type:UnstructuredGrid
!	Input:
!	Output:
!	Note:座標に基づいて色々やる，とりあえず，上界下界の座標基準と中心からの距離基準による境界属性割り当ては実装した，それ以外は今後の努力次第
!	Author:Akitaka Toyota
!	Date:2018.01.20
!	Update:-
!	Other:
!***********************************/
subroutine UMarkingVirtualCell(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    integer :: iKind,iUseObject = 0 ,iInternalBound, OutlineCell
    double precision, allocatable :: ObjectCenterCorrection(:)
    double precision :: ObjectRadius
    character(len=256) :: cBoundaryDataFileName,cAnnotate

    UG%GM%CellType = 0 !実セル
    UG%GI%OutlineCells = 0 !外周セル数の初期化

    if(UG%GM%Dimension /= 2) then
        write(6,*) "3D grid cannot used. Please write new preprocess program."
        stop
    end if

    !格子形状
    write(6,*) "What is the criteria for deciding the cell attributes? ;"
    write(6,*) "1:Coordinate 2:Distance from Center, 3:Other"
    ! read(5,*) iKind
    iKind = 2
    if(iKind == 3) then
        write(6,*) "I can't do it! Let's write new Program for Preprocess!! YEAH!!!"
        stop
    end if


    if(iKind == 1) then !座標値を基準に外周の境界条件を入力
        call CoordinateBasis

    else if(iKind == 2) then !中心からの距離を基準に外周の境界条件を入力
        call DistanceBasis

    end if
    open(unit=1, file=trim(adjustl(cBoundaryDataFileName)), status='unknown')
        do iLoop = 1,3+iKind
            read(1,*) cAnnotate
        end do

        read(1,*) cAnnotate, iUseObject
        if(iUseObject == 1) then
            read(1,*) cAnnotate,iInternalBound
            allocate(ObjectCenterCorrection(3))
            do iLoop=1,3
                read(1,*) cAnnotate,ObjectCenterCorrection(iLoop)
            end do
            read(1,*) cAnnotate,ObjectRadius

            call SetObjectBoundary

        else !iUseObject == 0
            if(maxval(UG%GM%CellType(:,1,1)) >= 999) then
                write(6,*) "undefined virtual cells remain. process will terminated."
                stop
            end if
        end if
    close(1)


return
contains
    subroutine CoordinateBasis !座標基準での属性づけ
    implicit none
    character :: cAnnotate
    integer, allocatable :: iBoundType(:)

    allocate(iBoundType(4),UG%GM%Bound(2,2))

    UG%GM%Bound(1,1) = minval(UG%CD%Edge(:,1))
    UG%GM%Bound(2,1) = maxval(UG%CD%Edge(:,1))
    UG%GM%Bound(1,2) = minval(UG%CD%Edge(:,2))
    UG%GM%Bound(2,2) = maxval(UG%CD%Edge(:,2))

    cBoundaryDataFileName = "BoundaryData.CoordinateBasis"
    !write(6,*) "Please input ""Boundary-Data-File"" written in CoordinateBasis format."
    !read(5,*) cBoundaryDataFileName

    open(unit=1, file=trim(adjustl(cBoundaryDataFileName)), status='unknown')
        do iLoop=1,4
            read(1,*) cAnnotate, iBoundType(iLoop)
        end do
    close(1)

    OutlineCell = 0
    do iCell = UG%GI%RealCells+1, UG%GI%AllCells

        if(UG%CD%Cell(iCell,1) < UG%GM%Bound(1,1)) then !x方向下端側仮想セル
            UG%GM%CellType(iCell,1,1) = iBoundType(1)

        else if(UG%CD%Cell(iCell,1) > UG%GM%Bound(2,1)) then !x方向上端側仮想セル
            UG%GM%CellType(iCell,1,1) = iBoundType(2)

        else if(UG%CD%Cell(iCell,2) < UG%GM%Bound(1,2)) then !y方向下端側仮想セル
            UG%GM%CellType(iCell,1,1) = iBoundType(3)

        else if(UG%CD%Cell(iCell,2) > UG%GM%Bound(2,2)) then !y方向上端側仮想セル
            UG%GM%CellType(iCell,1,1) = iBoundType(4)
        else !条件を満たしていない要素があるとき(矩形ではなかった or 内部要素があるとき) エラーコードとして
            UG%GM%CellType(iCell,1,1) = 999
        end if

        if(UG%GM%CellType(iCell,1,1) /= 999) then
            UG%VC%Type(iCell) = 1 !Outline
            UG%GI%OutlineCells = UG%GI%OutlineCells + 1 !外周セル数の足し上げ
        else
            UG%VC%Type(iCell) = 2 !Material Surface
        end if

    end do

    return
    end subroutine



    subroutine DistanceBasis
    implicit none
    integer :: iBoundType
    double precision, allocatable :: CenterCorrection(:)
    double precision :: MinimumRadiusOfGrid
    allocate(CenterCorrection(3))
    cBoundaryDataFileName = "BoundaryData.DistanceBasis"
    !write(6,*) "Please input ""Boundary-Data-File"" written in DistanceBasis format."
    !read(5,*) cBoundaryDataFileName
    open(unit=1, file=trim(adjustl(cBoundaryDataFileName)), status='unknown')
        read(1,*) cAnnotate,iBoundType
        do iLoop=1,3
            read(1,*) cAnnotate, CenterCorrection(iLoop)
        end do
        read(1,*) cAnnotate,MinimumRadiusOfGrid
    close(1)

    OutlineCell = 0
    do iCell = UG%GI%RealCells+1, UG%GI%AllCells
        if(sqrt(dot_product(UG%CD%Cell(iCell,:)-CenterCorrection,UG%CD%Cell(iCell,:)-CenterCorrection)) > MinimumRadiusOfGrid) then !補正済みセル中心からの距離が最小半径より外側にある仮想セルは
            UG%GM%CellType(iCell,1,1) = iBoundType
        else
            UG%GM%CellType(iCell,1,1) = 999
        end if

        if(UG%GM%CellType(iCell,1,1) /= 999) then
            UG%VC%Type(iCell) = 1 !Outline
            UG%GI%OutlineCells = UG%GI%OutlineCells + 1 !外周セル数の足し上げ
        else
            UG%VC%Type(iCell) = 2 !Material Surface
        end if
    end do

    return
    end subroutine DistanceBasis



    subroutine SetObjectBoundary
    implicit none

        do iCell = UG%GI%RealCells+1, UG%GI%AllCells
            if(sqrt(dot_product(UG%CD%Cell(iCell,:)-ObjectCenterCorrection,UG%CD%Cell(iCell,:)-ObjectCenterCorrection)) < ObjectRadius) then !物体の最大半径以下に存在するセルについて物体用の境界条件を適用する
                UG%GM%CellType(iCell,1,1) = iInternalBound
            end if
        end do

    return
    end subroutine SetObjectBoundary



end subroutine UMarkingVirtualCell
