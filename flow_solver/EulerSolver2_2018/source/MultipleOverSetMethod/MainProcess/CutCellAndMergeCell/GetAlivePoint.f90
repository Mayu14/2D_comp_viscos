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
    subroutine GetAlivePoint(iLocalEdge,iCenterX,iCenterY,Geom,LCE,GCC)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    integer, intent(in) :: iLocalEdge !of Global Cell (=1,2,3,4)
    integer, intent(in) :: iCenterX, iCenterY
    type(Geometry), intent(in) :: Geom
    type(CutCell), intent(inout) :: GCC
    type(CuttedEdge), intent(in) :: LCE
    double precision :: TmpEdgeLength
    double precision :: dotP
    double precision, allocatable :: CenterCoords(:),PointCoords(:)
    integer :: iLoop, iCount
    type(InternalEdgeOfCell), pointer :: TempIEC
    allocate(CenterCoords(2),PointCoords(3))

    !1.L界面の法線ベクトルとG界面の右回りベクトルとの内積を計算（左回りの方が馴染む可能性）
    !G-Cell Edge Vector's are CCW
    !Edge 1: (0,-1,0)
    !Edge 2: (1,0,0)
    !Edge 3: (0,1,0)
    !Edge 4: (-1,0,0)

    if(iLocalEdge == 1) then
        dotP = LCE%Normal(2)*(-1)
    else if(iLocalEdge == 2) then
        dotP = LCE%Normal(1)
    else if(iLocalEdge == 3) then
        dotP = LCE%Normal(2)
    else if(iLocalEdge == 4) then
        dotP = LCE%Normal(1)*(-1)
    end if


	!2.i番目のG界面と交叉して
		!法線ベクトルとの内積が正→i番，割合，始点，カウンタ＋1
		!同負										→i-1番，割合，終点，カウンタ＋1
        !※割合＝（i番の点から交叉点までの距離）／G格子の格子幅
    if(dotP > 0.0d0) then
        GCC%iAlivePoint(1) = iLocalEdge
        GCC%AliveCoords(1,1:3) = LCE%TmpInternalPoint(1:3)
    else
        GCC%iAlivePoint(2) = iLocalEdge-1
        GCC%AliveCoords(2,1:3) =LCE%TmpInternalPoint(1:3)
    end if

        GCC%AlivePointCheck = GCC%AlivePointCheck + 1

	!3.カウンタ＝1のとき処理戻す
	!カウンタ＝2のとき処理続行
	if(GCC%AlivePointCheck == 2) then
	!4.始点〜終点までに含まれる点をすべてALIVEとしてマーキングする
	iCount = 0
	call GetCellCenterCoords(iCenterX,iCenterY,Geom%Bound(1,1:2),2.0d0*Geom%Width(1,1,1:2),CenterCoords(1:2))

    TmpEdgeLength = 0.0d0 !Find Longest Edge

        if(GCC%iAlivePoint(1) <= GCC%iAlivePoint(2) .and. GCC%iAlivePoint(2)+1 <= 4) then
            do iLoop = GCC%iAlivePoint(1), GCC%iAlivePoint(2)+1
                if(iLoop /= GCC%iAlivePoint(1)) then
                    call RegisterCuttedGEdge(1)
                else if(iLoop /= GCC%iAlivePoint(2)+1) then
                    call RegisterCuttedGEdge(2)
                else
                    call RegisterCuttedGEdge(0)
                end if
            end do

        else
            do iLoop = GCC%iAlivePoint(1), 4
                if(iLoop /= GCC%iAlivePoint(1)) then
                    call RegisterCuttedGEdge(1)
                else
                    call RegisterCuttedGEdge(0)
                end if
            end do

            do iLoop = 1, GCC%iAlivePoint(2)+1
                if(iLoop /= GCC%iAlivePoint(2)+1) then
                    call RegisterCuttedGEdge(2)
                else
                    call RegisterCuttedGEdge(0)
                end if
            end do

        end if

    end if

    return
contains
    subroutine RegisterCuttedGEdge(iSwitch)
    implicit none
        integer, intent(in) :: iSwitch

        iCount = iCount + 1

        allocate(TempIEC)
            TempIEC%EdgeNum = iLoop
            TempIEC%GridNum = 0

            if(iSwitch == 0) then !not cutted Edge
                TempIEC%CutRatio = 1.0d0

            else
                call GetPointCoordsFromCellCoords(iLoop,CenterCoords(1:2),Geom%Width(1,1,1:2),PointCoords(1:2))
                PointCoords(3) = 0.0d0
                TempIEC%CutRatio=sqrt(dot_product(GCC%AliveCoords(iSwitch,1:3)-PointCoords(1:3),GCC%AliveCoords(iSwitch,1:3)-PointCoords(1:3)))/Geom%Area(mod(iLoop,2)+1)

            end if

        if(iCount == 1) then
            nullify(TempIEC%Next)
        else
            TempIEC%Next => GCC%IEC2
        end if
            GCC%IEC2%Next => TempIEC

            GCC%NumberOfPoint = GCC%NumberOfPoint + 1
            GCC%NumberOfEdge = GCC%NumberOfEdge + 1

            if(TmpEdgeLength < TempIEC%CutRatio*Geom%Area(mod(iLoop,2)+1)) then
                TmpEdgeLength = TempIEC%CutRatio*Geom%Area(mod(iLoop,2)+1)
                GCC%LongestEdgeNumber = iLoop
            end if

        return
    end subroutine RegisterCuttedGEdge

    end subroutine GetAlivePoint
