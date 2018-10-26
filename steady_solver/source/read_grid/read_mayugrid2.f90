!***********************************/
!	Name:MayuGrid2形式の格子情報読み込みメソッド
!	Alias:read_mayugrid2
!	Description:バージョン0.9に対応
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.??
!	Update:2018.10.25
!	Other:
!***********************************/
subroutine read_mayugrid2(cFileName, Geom)
    use StructVar_Mod_mk2
    use LoopVar_Mod_mk2
    use FrequentOperation_mk2
    implicit none

    character(len=64), intent(in) :: cFileName
    type(Geometry), intent(inout) :: Geom
    integer :: iXi_max, iEta_max, iZeta_max, iEdgeNum

    real :: rVersion
    character :: cTmp
    integer, allocatable :: iLocalNum(:)

    Geom%iDimension = 2
    iDim = Geom%iDimension !now only-2D
    iZeta_max = 1
    allocate(Geom%CellNumber(3))

    open(unit = 1, file =trim(adjustl(cFileName)), status = 'unknown')
    write(6,*) "mu"
        read(1, *) cTmp, cTmp, cTmp, rVersion ! version
        if(rVersion /= 0.9) then
            write(6, *) "grid file version error"
            stop
        end if

        read(1, *) cTmp, iXi_max
        read(1, *) cTmp, iEta_max

        allocate(Geom%VertexCoords(iDim, iXi_max, iEta_max, iZeta_max))
        !allocate(iPoint2Wall_belong(iXi_max, iEta_max), dPoint2Wall_distance(iXi_max, iEta_max))
        allocate(iLocalNum(2))

        read(1, *) cTmp ! explain file format
        read(1, *) cTmp ! data number
        do iPoint = 1, iXi_max*iEta_max
            call Global2Local(iPoint,iXi_max,iLocalNum)
            read(1, *) Geom%VertexCoords(1, iLocalNum(1), iLocalNum(2), 1), Geom%VertexCoords(2, iLocalNum(1), iLocalNum(2), 1)
        end do

        allocate(Geom%DistanceFromWall(2, iXi_max, iEta_max, 1))   ! 1or2(壁番号/距離), x, y, z
        read(1, *) cTmp
        do iPointY = 1, iEta_max
            do iPointX = 1, iXi_max
                read(1, *) Geom%DistanceFromWall(1, iPointX, iPointY, 1), Geom%DistanceFromWall(2, iPointX, iPointY, 1)

            end do
        end do

        allocate(Geom%BC%PDBW(iXi_max))

        read(1, *) cTmp
        do iEdge = 1, iXi_max
            read(1, *) cTmp, iEdgeNum
            Geom%BC%PDBW(iEdge)%includePoint = iEdgeNum
            allocate(Geom%BC%PDBW(iEdge)%point_id(iEdgeNum, iDim))
            do iPoint = 1, iEdgeNum
                read(1, *) Geom%BC%PDBW(iEdge)%point_id(iPoint, 1), Geom%BC%PDBW(iEdge)%point_id(iPoint, 2)
            end do
        end do
    close(1)

    Geom%CellNumber(1) = iXi_max - 1
    Geom%CellNumber(2) = iEta_max - 1   ! iXi_max, iEta_maxは格子点の数であり，セルはその中心点で定義されるため-1
    Geom%CellNumber(3) = iZeta_max

    return
end subroutine read_mayugrid2

