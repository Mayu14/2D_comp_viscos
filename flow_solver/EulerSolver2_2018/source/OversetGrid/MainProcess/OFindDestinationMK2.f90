!***********************************/
!	Name:移動先の格子を補間するための補間データ提供元格子検索
!	Alias:OFindDestination
!	Description:
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.11.16
!	Update:2017.11.16
!	Other:
!***********************************/
subroutine OFindDestinationMK2(UG,MG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    !type(Geometry), intent(in) :: Geom
    !type(OverSetGrid), intent(inout) :: OSG
    type(MoveGrid), intent(inout) :: MG

    call LocalGrid2AuxiliaryGrid
    call LocalGrid2GlobalGrid

return
contains

    subroutine LocalGrid2AuxiliaryGrid
    implicit none
    double precision, allocatable :: AuxiliaryGridNum(:)
    double precision, allocatable :: CornerPosition(:)
    double precision :: TmpSquare
    integer, allocatable :: iGlobalCellNumber(:)
    integer :: iOverflowCell,iContent
    integer :: iCheckPoint, iInterpolatePointNum
    integer :: iCheckEdge, iIntersectEdgeNum
    integer :: iOwnEdge, iPartnerEdge, iVertex4
    integer :: iTotalVertexNum
    double precision :: InterpolateVertex(:,:)

    type(OMPMoveGrid) :: OMG
    type(Used4Interpolation) :: U4I

    allocate(OMG%CUV%TargetVertex(6,3)) !1:6 is PointNumber
    allocate(AuxiliaryGridNum(3))
    allocate(iGlobalCellNumber(4))
    allocate(CornerPosition(2))

        MG%NSR%InterpolateFromGlobal = 0

        do iCell=1, UG%GI%AllCells !すべてのローカルセルについて
            iOverflowCell = 0
            !Auxiliary格子の座標は配列の添字と対応しているため，
            !(非構造格子の座標r,Auxiliary格子の座標r0，構造格子負端からのセル数n，格子幅dxについて
            !r-r0 > n*dxなる関係があることから，nについて解くと以下のような式になる(※2次元用で書いたためzは常に固定値になる)
            AuxiliaryGridNum(1:2) = int((MG%NSR%C%NextCoordsInN(iCell,1:2) - MG%ASG%Bound(1,1:2))/(MG%ASG%Width(1:2)))+1
            AuxiliaryGridNum(3) = 1

            if(AuxiliaryGridNum(1) > MG%ASG%CellNumber(1)) then
                iOverflowCell = 4
            else if(AuxiliaryGridNum(1) < 1) then
                iOverflowCell = 4
            else if(AuxiliaryGridNum(2) > MG%ASG%CellNumber(2)) then
                iOverflowCell = 4
            else if(AuxiliaryGridNum(2) < 1) then
                iOverflowCell = 4
            end if

            if(iOverflowCell == 0) then
                CornerPosition(1:2) = &
                    &   (MG%NSR%C%NextCoordsInN(iCell,1:2) - MG%ASG%Bound(1,1:2))/(MG%ASG%Width(1:2)) &
                    &-int((MG%NSR%C%NextCoordsInN(iCell,1:2) - MG%ASG%Bound(1,1:2))/(MG%ASG%Width(1:2)))
                !補助格子の隣接8セルと自セルをグローバル番号で表記する

                if(CornerPosition(1) <= 0.5d0) then
                    if(CornerPosition(2) <= 0.5d0) then !LeftiBottom
                        iGlobalCellNumber(3) = (AuxiliaryGridNum(2)-2)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) - 1 !(x-1,y-1)
                        iGlobalCellNumber(4) = (AuxiliaryGridNum(2)-2)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)     !(x,y-1)
                        iGlobalCellNumber(2) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) - 1 !(x-1,y)
                        iGlobalCellNumber(1) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)     !(x,y)
                    else !Left-Top
                        iGlobalCellNumber(4) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) - 1 !(x-1,y)
                        iGlobalCellNumber(1) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)     !(x,y)
                        iGlobalCellNumber(3) = (AuxiliaryGridNum(2))*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) - 1   !(x-1,y+1)
                        iGlobalCellNumber(2) = (AuxiliaryGridNum(2))*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)       !(x,y+1)

                    end if
                else !Right
                    if(CornerPosition(2) <= 0.5d0) then !Right-Bottom
                        iGlobalCellNumber(2) = (AuxiliaryGridNum(2)-2)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)     !(x,y-1)
                        iGlobalCellNumber(3) = (AuxiliaryGridNum(2)-2)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) + 1 !(x+1,y-1)
                        iGlobalCellNumber(1) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)     !(x,y)
                        iGlobalCellNumber(4) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) + 1 !(x+1,y)

                    else !Right-Top
                        iGlobalCellNumber(1) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)     !(x,y)
                        iGlobalCellNumber(2) = (AuxiliaryGridNum(2)-1)*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) + 1 !(x+1,y)
                        iGlobalCellNumber(4) = (AuxiliaryGridNum(2))*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1)       !(x,y+1)
                        iGlobalCellNumber(3) = (AuxiliaryGridNum(2))*MG%ASG%CellNumber(1) + AuxiliaryGridNum(1) + 1   !(x+1,y+1)

                    end if

                end if

                iVertex1 = UG%Tri%Point(iCell,1)
                iVertex2 = UG%Tri%Point(iCell,2)
                iVertex3 = UG%Tri%Point(iCell,3)

            do iLoop=1,4 !Loop for AuxiliryGrid
                if(iGlobalCellNumber(iLoop) > MG%ASG%TotalCell .and. iGlobalCellNumber(iLoop) < 1) then !if exist
                    do iContent=1, MG%ASG%ContentsAS2U(iGlobalCellNumber(iLoop)) !Content Loop

                        iAdjacentCell = MG%ASG%RelatedAS2U(iGlobalCellNumber(iLoop),iContent) !this cell will use for interpolation
                        TmpSquare = dot_product(MG%RC%Cell(iAdjacentCell,1:3)-MG%NSR%C%NextCoordsInN(iCell,1:3), &
                            &   MG%RC%Cell(iAdjacentCell,1:3)-MG%NSR%C%NextCoordsInN(iCell,1:3))

                        if(TmpSquare <= UG%GM%Area(iCell)+UG%GM%Area(iAdjacentCell)) then !too longer
                            iTotalVertexNum = 0
                            InterpolateVertex = 0.0d0
                            !call CheckPoint
                            do iPoint=1,3
                                iInterpolatePointNum = UG%Tri%Point(iAdjacentCell,iPoint)

                                iCheckPoint = CheckPointInternalTriangle(MG%RC%Point(iInterpolatePointNum,1:3), &
                                &   MG%NSR%PointCoords(iVertex1,1:3), MG%NSR%PointCoords(iVertex2,1:3),MG%NSR%PointCoords(iVertex3,1:3))

                                if(iCheckPoint == 1) then
                                    iTotalVertexNum = iTotalVertexNum + 1
                                    InterpolateVertex(iTotalVertexNum,:) = MG%RC%Point(iInterpolatePointNum,1:3)
                                end if
                            end do

                            do iOwnEdge=1,3
                                    iVertex1 = UG%Line%Point(UG%Tri%Edge(iCell,iOwnEdge),1)
                                    iVertex2 = UG%Line%Point(UG%Tri%Edge(iCell,iOwnEdge),2)
                                do iPartnerEdge=1,3
                                    !iIntersectEdgeNum = UG%Tri%Edge(iAdjacentCell,iPartnerEdge)
                                    iVertex3 = UG%Line%Point(UG%Tri%Edge(iAdjacentCell,iPartnerEdge),1)
                                    iVertex4 = UG%Line%Point(UG%Tri%Edge(iAdjacentCell,iPartnerEdge),2)

                                call CheckEdgeIntersectTriangle(MG%NSR%PointCoords(iVertex3,1:3),MG%NSR%PointCoords(iVertex4,1:3))
                                end do
                            end do

                            !call CheckEdge
                        end if

                    end do
                end if
            end do

            !補間の対象となるはずの補助格子のグローバル番号を格子に登録する
                do iLoop=1,9
                    if(iGlobalCellNumber(iLoop) > MG%ASG%TotalCell) then !計算したグローバル番号が補助格子中に存在していない場合補間対象から除外
                        iOverflowCell = iOverflowCell + 1

                    else if(iGlobalCellNumber(iLoop) < 0) then
                        iOverflowCell = iOverflowCell + 1

                    else if(MG%ASG%ContentsAS2U(iGlobalCellNumber(iLoop)) == 0) then !計算したグローバル番号の補助格子中にローカル格子が存在しない場合補間対象から除外
                        iOverflowCell = iOverflowCell + 1
                    end if
                    if(iOverflowCell == 9) exit
                end do

            end if

            if(iOverflowCell < 2) then
                MG%NSR%Previous2Next(iCell,:) = iGlobalCellNumber(:)
            else
                MG%NSR%Previous2Next(iCell,:) = 0
                MG%NSR%InterpolateFromGlobal(iCell) = 1
            end if

        end do

        !call CheckOfInterpolateSource
    return
    end subroutine LocalGrid2AuxiliaryGrid

    function CheckPointInternalTriangle(Point,Vertex1, Vertex2, Vertex3) result(TorF)
    implicit none
        double precision, intent(in) :: Point(:), Vertex1(:), Vertex2(:), Vertex3(:)
        integer :: TorF
        double precision, dimension(3,3) :: VertexOnRelatedCoords
        double precision :: AB, BC, CA, AA,BB

        VertexOnRelatedCoords(1,:) = Vertex1(:) - Point(:)
        VertexOnRelatedCoords(2,:) = Vertex2(:) - Point(:)
        VertexOnRelatedCoords(3,:) = Vertex3(:) - Point(:)


        AB = dot_product(VertexOnRelatedCoords(1,:),VertexOnRelatedCoords(2,:))
        BA = dot_product(VertexOnRelatedCoords(2,:),VertexOnRelatedCoords(3,:))
        CA = dot_product(VertexOnRelatedCoords(3,:),VertexOnRelatedCoords(1,:))
        AA = dot_product(VertexOnRelatedCoords(1,:),VertexOnRelatedCoords(1,:))

        if(AB*BC-AA*BC < 0.0d0) then
            TorF = 0
        else
            BB = dot_product(VertexOnRelatedCoords(2,:),VertexOnRelatedCoords(2,:))
            if(AB*BC-AC*BB < 0.0d0) then
                TorF = 0
            else
                TorF = 1
            end if
        end if

        return
    end function CheckPointInternalTriangle


    subroutine CheckEdgeIntersectTriangle(Vertex11,Vertex12, Vertex21, Vertex22, TorF,IntersectCoords)
    implicit none
        double precision, intent(in) :: Vertex11(:), Vertex12(:), Vertex21(:), Vertex22(:)
        integer, intent(out) :: TorF
        double precision, intent(out) :: IntersectCoords(:)

        double precision :: Ksi, Eta, Delta
        double precision :: Lambda, Mu !Coefficient of Line's Gradient

        Ksi = (Vertex22(2)-Vertex21(2))*(Vertex22(1)-Vertex11(1)) - (Vertex22(1)-Vertex21(1))*(Vertex22(2)-Vertex11(2))
        Eta = (Vertex12(1)-Vertex11(1))*(Vertex22(2)-Vertex11(2)) - (Vertex12(2)-Vertex11(2))*(Vertex22(1)-Vertex11(1))
        Delta = (Vertex12(1)-Vertex11(1))*(Vertex22(2)-Vertex21(2)) - (Vertex12(2)-Vertex11(2))*(Vertex22(1)-Vertex21(1))

        Lambda = Ksi / Delta
        Mu = Eta / Delta

        if(Lambda >= 0.0d0 .and. Lambda <= 1.0d0 .and. Mu >=0.0d0 .and. Mu <= 1.0d0) then
            IntersectCoords(1:2) = Vertex11(1:2) + Lambda*(Vertex12(1:2)-Vertex11(1:2))
            IntersectCoords(3) = 0.0d0
            TorF = 1
        else
            TorF = 0
            IntersectCoords = 0.0d0
        end if
    return
    end subroutine CheckEdgeIntersectTriangle


        subroutine CheckOfInterpolateSource
        implicit none
        iCenterX=0
        iCenterY=0
        do iLoop=1,UG%GI%AllCells
            if(MG%NSR%InterpolateFromGlobal(iLoop) == 0) then
                iCenterX = iCenterX+1
            else
                iCenterY = iCenterY+1
            end if
        end do
        write(6,*) UG%GI%AllCells,iCenterX,iCenterY
        stop
        return
        end subroutine CheckOfInterpolateSource


    subroutine LocalGrid2GlobalGrid
    implicit none

        integer :: iGlobalCellNumber

        OSG%FrequencyDistribution = 0

        do  iCell=1, UG%GI%AllCells

            OSG%RelatedU2S(iCell,1:2) = int((MG%NSR%C%NextCoordsInG(iCell,1:2) - Geom%Bound(1,1:2))/(2.0d0*Geom%Width(1,1,1:2)))+1
            OSG%RelatedU2S(iCell,3) = 1

            iGlobalCellNumber = (OSG%RelatedU2S(iCell,2)-1)*Geom%CellNumber(1) + OSG%RelatedU2S(iCell,1)

            !構造格子ijkが非構造格子から参照された回数を記録
            OSG%FrequencyDistribution(iGlobalCellNumber) = OSG%FrequencyDistribution(iGlobalCellNumber) +1

            !最後に構造格子に所属する非構造格子への対応を記録する
            !OSG%RelatedS2U(iGlobalCellNumber,OSG%FrequencyDistribution(iGlobalCellNumber)) = iAdjacentCell !11/21
            OSG%RelatedS2U(iGlobalCellNumber,OSG%FrequencyDistribution(iGlobalCellNumber)) = iCell !11/22
        end do

    return
    end subroutine LocalGrid2GlobalGrid

end subroutine OFindDestinationMK2
