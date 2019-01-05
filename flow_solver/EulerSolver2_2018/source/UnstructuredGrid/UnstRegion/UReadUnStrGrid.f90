!***********************************/
!	Name:非構造格子の計算領域データを読み込むプログラム
!	Alias:UReadUnStrGrid
!	Description:MAYUGridからデータを読み出すプログラム
!	Type:UG
!	Input:UConf,UG
!	Output:UG
!	Note:必要に応じて凸包の情報を読み取るように変更
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:2018.02.03
!	Other:
!***********************************/
subroutine UReadUnStrGrid(UConf,UCC,UCE,UG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(CellCenter), intent(inout) :: UCC !本プログラム内ではなくUAllocVariables内で大きさを割り当てるためだけに呼び出している
    type(CellEdge), intent(inout) :: UCE !本プログラム内ではなくUAllocVariables内で大きさを割り当てるためだけに呼び出している
    type(UnstructuredGrid), intent(inout) :: UG
    character(len=128) cFileName, cDirectory,cAnnotate

    !cFileName = "UnStrGrid"
    !cFileName = "MiniCircle_Fine.mayu"
    cFileName = "circle_HD.mayu"
    if(UConf%UseJobParallel == 1) cFileName = trim(adjustl(UConf%cGridName))
    if(UConf%CalcEnv == 0) then
        cDirectory = "/mnt/g/Toyota/Data/grid_vtk/NACA4_mayu_HD_course_rev2/"
    else
        cDirectory = ""
    end if
    cFileName = trim(adjustl(cDirectory))//trim(adjustl(cFileName))
    !write(6,*) "Please input file name of *.mayu"
    !read(5,*) cFileName

    !cFileName = "tri_SquareGrid.mayu"
    !write(6,*) trim(adjustl(cFileName))
    UG%InternalRadius = 0.10d0 + epsilon(0.05d0)
!点番号は1始まり
    open(unit=UConf%my_rank+100,file=trim(adjustl(cFileName)),status='unknown')
!格子の基本情報

        read(UConf%my_rank+100,*) cAnnotate, UG%GI%Points
        read(UConf%my_rank+100,*) cAnnotate, UG%GI%Edges
        read(UConf%my_rank+100,*) cAnnotate, UG%GI%RealCells
        read(UConf%my_rank+100,*) cAnnotate, UG%VC%Total
        read(UConf%my_rank+100,*) cAnnotate, UG%GI%OutlineCells

        UG%GI%AllCells = UG%GI%RealCells + UG%VC%Total
        UG%GM%Dimension = 3
        read(UConf%my_rank+100,*) cAnnotate, UG%Tri%Total
        read(UConf%my_rank+100,*) cAnnotate, UG%Quad%Total
        read(UConf%my_rank+100,*) cAnnotate, UG%Line%Total
        read(UConf%my_rank+100,*) cAnnotate, UG%CH%iTotal
        read(UConf%my_rank+100,*) cAnnotate, UG%IO%iTotal

        !非構造格子用の配列確保
        call UAllocVariables(UConf,UG,UCC,UCE)
!三角形要素の幾何的関係
        read(UConf%my_rank+100,*) cAnnotate !"TriPoint"
        do iCell=1, UG%Tri%Total
            read(UConf%my_rank+100,*) (UG%Tri%Point(iCell,iLoop),iLoop=1,3)
        end do

        read(UConf%my_rank+100,*) cAnnotate !"TriEdge"
        do iCell=1, UG%Tri%Total
            read(UConf%my_rank+100,*) (UG%Tri%Edge(iCell,iLoop),iLoop=1,3)
        end do

        read(UConf%my_rank+100,*) cAnnotate !"TriCell"既存のプロジェクトを開く
        do iCell=1, UG%Tri%Total
            read(UConf%my_rank+100,*) (UG%Tri%Cell(iCell,iLoop),iLoop=1,3)
        end do

!四辺形要素の幾何的関係
        if(UG%Quad%Total /= 0) then
            read(UConf%my_rank+100,*) cAnnotate !"coming soon..."
        end if

!線要素の幾何的関係
        read(UConf%my_rank+100,*) cAnnotate !"LinePoint"
        do iEdge=1, UG%Line%Total
            read(UConf%my_rank+100,*) (UG%Line%Point(iEdge,iLoop),iLoop=1,2)
        end do

        read(UConf%my_rank+100,*) cAnnotate !"LineCell"
        do iEdge=1, UG%Line%Total
            read(UConf%my_rank+100,*) ((UG%Line%Cell(iEdge,iSide,iLocalEdge),iLocalEdge=1,2),iSide=1,2)
        end do

!仮想格子の幾何的関係
        read(UConf%my_rank+100,*) cAnnotate !"VCCell"
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            read(UConf%my_rank+100,*) (UG%VC%Cell(iCell,iLoop),iLoop=1,2)

        end do

        read(UConf%my_rank+100,*) cAnnotate !"VCEdge"
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            read(UConf%my_rank+100,*) UG%VC%Edge(iCell)
        end do

!格子点の座標情報
        read(UConf%my_rank+100,*) cAnnotate !"PointC"
        do iPoint=1, UG%GI%Points
            read(UConf%my_rank+100,*) (UG%CD%Point(iPoint,iLoop),iLoop=1,3)
            UG%CD%Point(iPoint,4) = 1.0d0
        end do

!格子境界中心の座標情報
        read(UConf%my_rank+100,*) cAnnotate !"EdgeC"
        do iEdge=1,UG%GI%Edges
            read(UConf%my_rank+100,*) (UG%CD%Edge(iEdge,iLoop),iLoop=1,3)
        end do

!格子要素中心の座標情報
        read(UConf%my_rank+100,*) cAnnotate !"CellC"
        do iCell=1,UG%GI%AllCells
            read(UConf%my_rank+100,*) (UG%CD%Cell(iCell,iLoop),iLoop=1,3)
            UG%CD%Cell(iCell,4) = 1.0d0
        end do

!格子境界の面積
        read(UConf%my_rank+100,*) cAnnotate !"EdgeS"
        do iEdge=1,UG%GI%Edges
            read(UConf%my_rank+100,*) UG%GM%Area(iEdge)
        end do

!格子境界の体積
        read(UConf%my_rank+100,*) cAnnotate !"CellV"
        do iCell=1,UG%GI%RealCells
            read(UConf%my_rank+100,*) UG%GM%Volume(iCell)
        end do

!セル界面における法線ベクトル
        read(UConf%my_rank+100,*) cAnnotate !"EdgeNormal"
        do iEdge=1,UG%GI%Edges
            read(UConf%my_rank+100,*) (UG%GM%Normal(iEdge,iLoop),iLoop=1,3)
        end do

!セル中心から界面中心までの距離ベクトル(中心基準)
        read(UConf%my_rank+100,*) cAnnotate !"Width"
        do iCell=1,UG%GI%RealCells
            read(UConf%my_rank+100,*) ((UG%GM%Width(iCell,iLocalEdge,iLoop),iLoop=1,3),iLocalEdge=1,3)
        end do

        read(UConf%my_rank+100,*) cAnnotate !"BoudnaryCondition of VirtualCell"
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            read(UConf%my_rank+100,*) UG%GM%CellType(iCell,1,1)
        end do

!Radius of Cell's Inscribed-Circle for compute CFL Condition
        read(UConf%my_rank+100,*) cAnnotate
        do iCell=1,UG%GI%RealCells
            read(UConf%my_rank+100,*) UG%InscribedCircle(iCell)
        end do

!内部に含む物体を構成する凸包のデータ
        read(UConf%my_rank+100,*) cAnnotate !"InternalObject"
        if(UG%IO%iTotal /= 0) then
            do iPoint=1, UG%IO%iTotal
                read(UConf%my_rank+100,*) UG%IO%PointNum(iPoint)
            end do
            UG%IO%iMinNumber = minval(UG%IO%PointNum)
        end if

        read(UConf%my_rank+100,*) cAnnotate !"AverageWid"
        !if(UConf%UseMUSCL == 1) then
            do iCell=1,UG%GI%RealCells
                read(UConf%my_rank+100,*) UG%GM%AverageWidth(iCell)
            end do
        !end if
        read(UConf%my_rank+100,*) cAnnotate !"ConvexHull"

        if(UG%CH%iTotal /= 0) then
            do iCell=1,UG%CH%iTotal
                read(UConf%my_rank+100,*) UG%CH%PointNum(iCell)
            end do

        end if

        if(UG%GI%AllCells - UG%GI%RealCells - UG%GI%OutlineCells /= 0) then
            !read(UConf%my_rank+100,*) cAnnotate !"NearestSurfaceBoundaryEdgeNum "
            !do iCell = 1, UG%Tri%Total
                !read(UConf%my_rank+100,*) UG%Tri%Belongs2Wall(iCell)
            !end do

            !read(UConf%my_rank+100,*) cAnnotate !"DistanceFromObjectSurface "
            !do iCell = 1, UG%Tri%Total
                !read(UConf%my_rank+100,*) UG%Tri%Distance(iCell)
            !end do

            read(UConf%my_rank+100,*) cAnnotate !"NearestSurfaceBoundaryEdgeNum4Edge "
            do iEdge = 1, UG%Line%Total
                read(UConf%my_rank+100,*) UG%Line%Belongs2Wall(iEdge)
            end do

            read(UConf%my_rank+100,*) cAnnotate !"DistanceFromObjectSurface4Edge "
            do iEdge = 1, UG%Line%Total
                read(UConf%my_rank+100,*) UG%Line%Distance(iEdge)
            end do

            !read(UConf%my_rank+100,*) cAnnotate, UG%GM%BC%iWallTotal !"Wall2Cell_data "
            if(UConf%TurbulenceModel /= 0) then
                !do iEdge = 1, UG%GM%BC%iWallTotal
                    !read(UConf%my_rank+100,*) UG%GM%BC%VW(iEdge)%iGlobalEdge, UG%GM%BC%VW(iEdge)%iNumberOfMember
                    !allocate(UG%GM%BC%VW(iEdge)%iMemberCell(UG%GM%BC%VW(iEdge)%iNumberOfMember))
                    !do iCell = 1, UG%GM%BC%VW(iEdge)%iNumberOfMember
                        !read(UConf%my_rank+100,*) UG%GM%BC%VW(iEdge)%iMemberCell(iCell)
                    !end do
                !end do

                read(UConf%my_rank+100,*) cAnnotate, UG%GM%BC%iWallTotal !"Wall2Edge_data "
                do iLoop = 1, UG%GM%BC%iWallTotal
                    read(UConf%my_rank+100,*) UG%GM%BC%VW(iLoop)%iGlobalEdge, UG%GM%BC%VW(iLoop)%iNumberOfMemberEdge
                    allocate(UG%GM%BC%VW(iLoop)%iMemberEdge(UG%GM%BC%VW(iLoop)%iNumberOfMemberEdge))
                    do iEdge = 1, UG%GM%BC%VW(iLoop)%iNumberOfMemberEdge
                        read(UConf%my_rank+100,*) UG%GM%BC%VW(iLoop)%iMemberEdge(iEdge)
                    end do
                end do
            end if
        end if

        close(UConf%my_rank+100)

        UG%GM%Bound(1,1) = minval(UG%CD%Point(:,1),1)
        UG%GM%Bound(2,1) = maxval(UG%CD%Point(:,1),1)
        UG%GM%Bound(1,2) = minval(UG%CD%Point(:,2),1)
        UG%GM%Bound(2,2) = maxval(UG%CD%Point(:,2),1)
        UG%GM%Bound(1,3) = minval(UG%CD%Point(:,3),1)
        UG%GM%Bound(2,3) = maxval(UG%CD%Point(:,3),1)

        if(UG%GI%OutlineCells /= UG%VC%Total) then
            UG%InternalBoundary = 1
        else
            UG%InternalBoundary = 0
        end if

    return
end subroutine UReadUnStrGrid
