!***********************************/
!	Name:
!	Alias:CheckContactEdge2Edge
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
    subroutine CheckContactEdge2Edge(Geom,CutE,GG)!(iCrossCondition)!(Coords,Width,Bound,xCellNumber) result(iGlobalNum)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    type(Geometry), intent(in) :: Geom
    type(CuttedEdge), intent(inout) :: CutE
    type(GlobalGrid), intent(inout) :: GG
    !integer :: ixCellNumber
    integer :: iEdgeX, iEdgeY, iToF, iGlobalNumber
    integer, allocatable :: ixyPosition(:,:) !(GlobalCell, xy)1or2:Global1, Global2 , 1or2:xy
    double precision, allocatable :: EdgeEndPoint(:,:) !界面の端点座標の格納
    double precision, allocatable :: IntersectCoords(:)
    allocate(ixyPosition(2,2),EdgeEndPoint(2,2),IntersectCoords(3))
    !2つのGセル番号，およびGセルのx方向のセル数がわかれば

    call Global2Local(CutE%GlobalNumber1,Geom%CellNumber(1),ixyPosition(1,1:2))
    call Global2Local(CutE%GlobalNumber2,Geom%CellNumber(1),ixyPosition(2,1:2))

!セル間に含まれる全界面に関するループ
!x軸と垂直な界面
    do iEdgeY = ixyPosition(1,2), ixyPosition(2,2)-1 !界面jはセルjの上側界面
        do iEdgeX = ixyPosition(1,1), ixyPosition(2,1)-1 !界面iはセルiの右側界面
            !界面両端の座標をEdgeEndPointに格納
            call GetEdgeEndPoint(iEdgeX,iEdgeY,Geom%Bound(1,1:2),2.0d0*Geom%Width(1,1,1:2),EdgeEndPoint(1,1:2))
            call GetEdgeEndPoint(iEdgeX+1,iEdgeY,Geom%Bound(1,1:2),2.0d0*Geom%Width(1,1,1:2),EdgeEndPoint(2,1:2))

            !界面の交叉を判定して...
            if(CheckEdgeIntersect(CutE%EndPoint(1,:),CutE%EndPoint(2,:),EdgeEndPoint(1,:),EdgeEndPoint(2,:)) == 1) then !Cross "True"
                !真ならば交叉座標を取得し，入力する
                CutE%InternalPointNumber = CutE%InternalPointNumber + 1
                call GetEdgeIntersectionCoords(CutE%EndPoint(1,:),CutE%EndPoint(2,:),EdgeEndPoint(1,:),EdgeEndPoint(2,:),CutE%TmpInternalPoint(:))
                iGlobalNumber = Local2Global(iEdgeX,iEdgeY,Geom%CellNumber(1))
                call GetAlivePoint(4,iEdgeX,iEdgeY,Geom,CutE,GG%CutC(iGlobalNumber))
                call GetAlivePoint(2,iEdgeX,iEdgeY-1,Geom,CutE,GG%CutC(iGlobalNumber))
                call RegisterCutCondition(1,CutE,GG%CutC(iGlobalNumber))

            end if

        end do
    end do

!Edge of vertical with Y-axis
    do iEdgeX = ixyPosition(1,1), ixyPosition(2,1)-1 !界面iはセルiの右側界面
        do iEdgeY = ixyPosition(1,2), ixyPosition(2,2)-1 !界面jはセルjの上側界面
            !界面両端の座標をEdgeEndPointに格納
            call GetEdgeEndPoint(iEdgeX,iEdgeY,Geom%Bound(1,1:2),2.0d0*Geom%Width(1,1,1:2),EdgeEndPoint(1,1:2))
            call GetEdgeEndPoint(iEdgeX,iEdgeY+1,Geom%Bound(1,1:2),2.0d0*Geom%Width(1,1,1:2),EdgeEndPoint(2,1:2))

            !界面の交叉を判定して...
            if(CheckEdgeIntersect(CutE%EndPoint(1,:),CutE%EndPoint(2,:),EdgeEndPoint(1,:),EdgeEndPoint(2,:)) == 1) then !Cross "True"
                !真ならば交叉座標を取得し，入力する
                CutE%InternalPointNumber = CutE%InternalPointNumber + 1
                call GetEdgeIntersectionCoords(CutE%EndPoint(1,:),CutE%EndPoint(2,:),EdgeEndPoint(1,:),EdgeEndPoint(2,:),CutE%TmpInternalPoint(:))
                iGlobalNumber = Local2Global(iEdgeX,iEdgeY,Geom%CellNumber(1))
                call GetAlivePoint(1,iEdgeX,iEdgeY,Geom,CutE,GG%CutC(iGlobalNumber))
                call GetAlivePoint(3,iEdgeX-1,iEdgeY,Geom,CutE,GG%CutC(iGlobalNumber))
                call RegisterCutCondition(1,CutE,GG%CutC(iGlobalNumber))

            end if

        end do
    end do


    return
    end subroutine CheckContactEdge2Edge
