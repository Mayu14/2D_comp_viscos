!***********************************/
!	Name:不完全なオーバーラップをしているセルを内挿対象から除外する
!	Alias:ExcludeImperfectOverlapped
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:To delete relationship of interpolation on the Outer-Boundary-Cell
!	Author:Akitaka Toyota
!	Date:2018.02.04
!	Update:
!	Other:
!***********************************/
    subroutine ExcludeImperfectOverlapped(UG,MG,GG,Geom) !(iGrid,GG,iGlobalCellNumber,Geom,CUV)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    type(UnstructuredGrid), intent(in) :: UG
    type(MoveGrid), intent(in) :: MG
    type(GlobalGrid), intent(inout) :: GG
    type(Geometry), intent(inout) :: Geom

    double precision, allocatable :: PointCoords(:,:)
    integer, allocatable :: LocalGridNum(:) !iCenterX, iCenterY
    integer :: iGlobalNumber
    integer :: iPoint, iCell, iEdge, iLoop
    type(DonorElement), pointer :: DE
    type(DonorElement), pointer :: TempDE



    allocate(LocalGridNum(2),PointCoords(2,2))
    !do iPoint = 1, UG%CH%iTotal
    do iCell = UG%GI%RealCells+1, UG%GI%AllCells
        if(UG%GM%CellType(iCell,1,1) == 4) then !Overset boundary only
        !仮想セル番号でループ回して界面番号を得る，得られた界面番号を元に両端の所属するセルを...(凸包を使う場合に比べて処理量は倍になるけどとりあえず)
            iEdge = UG%VC%Edge(iCell)
            do iLoop = 1, 3
                if(iLoop < 3) then
                    iPoint = UG%Line%Point(iEdge,iLoop)
                    PointCoords(iLoop,1:2) = MG%NSR%P%NextCoordsInG(iPoint,1:2)
                    call GetLocalNumber2D(PointCoords(iLoop,1:2),&
                            &   2.0d0*Geom%Width(1,1,1:2),Geom%Bound(1,1:2),LocalGridNum(1:2))

                else
                    call GetLocalNumber2D(0.5d0*(PointCoords(1,1:2)+PointCoords(2,1:2)),&
                            &   2.0d0*Geom%Width(1,1,1:2),Geom%Bound(1,1:2),LocalGridNum(1:2))

                end if


                if(Geom%Interpolated(LocalGridNum(1),LocalGridNum(2),1) /= 0) then

                    iGlobalNumber = Local2Global(LocalGridNum(1),LocalGridNum(2),Geom%CellNumber(1))

                    DE => GG%GDE(iGlobalNumber)%DE !Delete Relation
                    do while(associated(DE))
                        TempDE => DE%Next
                        nullify(DE%Next) !?
                        deallocate(DE)
                        DE => TempDE
                    end do
                    nullify(GG%GDE(iGlobalNumber)%DE) !?

                    Geom%Interpolated(LocalGridNum(1),LocalGridNum(2),1) = 0 !Don't Interpolate
                end if
            end do
        end if
    end do
!stop
    return
    end subroutine ExcludeImperfectOverlapped
