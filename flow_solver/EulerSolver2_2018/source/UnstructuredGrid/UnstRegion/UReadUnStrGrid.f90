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
    character(len=256) cFileName, cDirectory,cAnnotate
    integer :: iUnit
    !cFileName = "UnStrGrid"
    !cFileName = "MiniCircle_Fine.mayu"
    cFileName = "circle_HD.mayu"
    if(UConf%UseJobParallel == 1 .or. UConf%SwitchProgram > 5) cFileName = trim(adjustl(UConf%cGridName))
    if(UConf%CalcEnv == 0) then
        !cDirectory = "/mnt/g/Toyota/Data/grid_vtk/NACA4_mayu_HD_course_rev2/"
        cDirectory = "/mnt/g/Toyota/Data/grid_vtk/valid/mayu/"
        ! cDirectory = "/mnt/g/Toyota/Data/grid_vtk/NACA4_mayu_course_rev4/"
        cFileName = trim(adjustl(cDirectory))//trim(adjustl(cFileName))
    else
        cDirectory = ""
    end if

    !write(6,*) "Please input file name of *.mayu"
    !read(5,*) cFileName

    !cFileName = "tri_SquareGrid.mayu"
    !write(6,*) trim(adjustl(cFileName))
    UG%InternalRadius = 0.10d0 + epsilon(0.05d0)
!点番号は1始まり
    iUnit = UConf%my_rank + 100
    open(unit=iUnit,file=trim(adjustl(cFileName)),status='unknown')
!格子の基本情報

        read(iUnit,*) cAnnotate, UG%GI%Points
        read(iUnit,*) cAnnotate, UG%GI%Edges
        read(iUnit,*) cAnnotate, UG%GI%RealCells
        read(iUnit,*) cAnnotate, UG%VC%Total
        read(iUnit,*) cAnnotate, UG%GI%OutlineCells

        UG%GI%AllCells = UG%GI%RealCells + UG%VC%Total
        UG%GM%Dimension = 3
        read(iUnit,*) cAnnotate, UG%Tri%Total
        read(iUnit,*) cAnnotate, UG%Quad%Total
        read(iUnit,*) cAnnotate, UG%Line%Total
        read(iUnit,*) cAnnotate, UG%CH%iTotal
        read(iUnit,*) cAnnotate, UG%IO%iTotal

        !非構造格子用の配列確保
        call UAllocVariables(UConf,UG,UCC,UCE)
!三角形要素の幾何的関係
        read(iUnit,*) cAnnotate !"TriPoint"
        do iCell=1, UG%Tri%Total
            read(iUnit,*) (UG%Tri%Point(iCell,iLoop),iLoop=1,3)
        end do

        read(iUnit,*) cAnnotate !"TriEdge"
        do iCell=1, UG%Tri%Total
            read(iUnit,*) (UG%Tri%Edge(iCell,iLoop),iLoop=1,3)
        end do

        read(iUnit,*) cAnnotate !"TriCell"既存のプロジェクトを開く
        do iCell=1, UG%Tri%Total
            read(iUnit,*) (UG%Tri%Cell(iCell,iLoop),iLoop=1,3)
        end do

!四辺形要素の幾何的関係
        if(UG%Quad%Total /= 0) then
            read(iUnit,*) cAnnotate !"coming soon..."
        end if

!線要素の幾何的関係
        read(iUnit,*) cAnnotate !"LinePoint"
        do iEdge=1, UG%Line%Total
            read(iUnit,*) (UG%Line%Point(iEdge,iLoop),iLoop=1,2)
        end do

        read(iUnit,*) cAnnotate !"LineCell"
        do iEdge=1, UG%Line%Total
            read(iUnit,*) ((UG%Line%Cell(iEdge,iSide,iLocalEdge),iLocalEdge=1,2),iSide=1,2)
        end do

!仮想格子の幾何的関係
        read(iUnit,*) cAnnotate !"VCCell"
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            read(iUnit,*) (UG%VC%Cell(iCell,iLoop),iLoop=1,2)

        end do

        read(iUnit,*) cAnnotate !"VCEdge"
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            read(iUnit,*) UG%VC%Edge(iCell)
        end do

!格子点の座標情報
        read(iUnit,*) cAnnotate !"PointC"
        do iPoint=1, UG%GI%Points
            read(iUnit,*) (UG%CD%Point(iPoint,iLoop),iLoop=1,3)
            UG%CD%Point(iPoint,4) = 1.0d0
        end do

!格子境界中心の座標情報
        read(iUnit,*) cAnnotate !"EdgeC"
        do iEdge=1,UG%GI%Edges
            read(iUnit,*) (UG%CD%Edge(iEdge,iLoop),iLoop=1,3)
        end do

!格子要素中心の座標情報
        read(iUnit,*) cAnnotate !"CellC"
        do iCell=1,UG%GI%AllCells
            read(iUnit,*) (UG%CD%Cell(iCell,iLoop),iLoop=1,3)
            UG%CD%Cell(iCell,4) = 1.0d0
        end do

!格子境界の面積
        read(iUnit,*) cAnnotate !"EdgeS"
        do iEdge=1,UG%GI%Edges
            read(iUnit,*) UG%GM%Area(iEdge)
        end do

!格子境界の体積
        read(iUnit,*) cAnnotate !"CellV"
        do iCell=1,UG%GI%RealCells
            read(iUnit,*) UG%GM%Volume(iCell)
        end do

!セル界面における法線ベクトル
        read(iUnit,*) cAnnotate !"EdgeNormal"
        do iEdge=1,UG%GI%Edges
            read(iUnit,*) (UG%GM%Normal(iEdge,iLoop),iLoop=1,3)
        end do

!セル中心から界面中心までの距離ベクトル(中心基準)
        read(iUnit,*) cAnnotate !"Width"
        do iCell=1,UG%GI%RealCells
            read(iUnit,*) ((UG%GM%Width(iCell,iLocalEdge,iLoop),iLoop=1,3),iLocalEdge=1,3)
        end do

        read(iUnit,*) cAnnotate !"BoudnaryCondition of VirtualCell"
        do iCell=UG%GI%RealCells+1, UG%GI%AllCells
            read(iUnit,*) UG%GM%CellType(iCell,1,1)
        end do

!Radius of Cell's Inscribed-Circle for compute CFL Condition
        read(iUnit,*) cAnnotate
        do iCell=1,UG%GI%RealCells
            read(iUnit,*) UG%InscribedCircle(iCell)
        end do

!内部に含む物体を構成する凸包のデータ
        read(iUnit,*) cAnnotate !"InternalObject"
        if(UG%IO%iTotal /= 0) then
            do iPoint=1, UG%IO%iTotal
                read(iUnit,*) UG%IO%PointNum(iPoint)
            end do
            UG%IO%iMinNumber = minval(UG%IO%PointNum)
        end if

        read(iUnit,*) cAnnotate !"AverageWid"
        !if(UConf%UseMUSCL == 1) then
            do iCell=1,UG%GI%RealCells
                read(iUnit,*) UG%GM%AverageWidth(iCell)
            end do
        !end if

        if(UG%CH%iTotal /= 0) then
            read(iUnit,*) cAnnotate !"ConvexHull"
            do iCell=1,UG%CH%iTotal
                read(iUnit,*) UG%CH%PointNum(iCell)
            end do

        end if

        if(UG%GI%AllCells - UG%GI%RealCells - UG%GI%OutlineCells /= 0) then
            !read(iUnit,*) cAnnotate !"NearestSurfaceBoundaryEdgeNum "
            !do iCell = 1, UG%Tri%Total
                !read(iUnit,*) UG%Tri%Belongs2Wall(iCell)
            !end do

            !read(iUnit,*) cAnnotate !"DistanceFromObjectSurface "
            !do iCell = 1, UG%Tri%Total
                !read(iUnit,*) UG%Tri%Distance(iCell)
            !end do

            read(iUnit,*) cAnnotate !"NearestSurfaceBoundaryEdgeNum4Edge "

            do iEdge = 1, UG%Line%Total
                read(iUnit,*) UG%Line%Belongs2Wall(iEdge)
            end do

            read(iUnit,*) cAnnotate !"DistanceFromObjectSurface4Edge "
            do iEdge = 1, UG%Line%Total
                read(iUnit,*) UG%Line%Distance(iEdge)
            end do

            !read(iUnit,*) cAnnotate, UG%GM%BC%iWallTotal !"Wall2Cell_data "
            !if(UConf%TurbulenceModel /= 0) then
                !do iEdge = 1, UG%GM%BC%iWallTotal
                    !read(iUnit,*) UG%GM%BC%VW(iEdge)%iGlobalEdge, UG%GM%BC%VW(iEdge)%iNumberOfMember
                    !allocate(UG%GM%BC%VW(iEdge)%iMemberCell(UG%GM%BC%VW(iEdge)%iNumberOfMember))
                    !do iCell = 1, UG%GM%BC%VW(iEdge)%iNumberOfMember
                        !read(iUnit,*) UG%GM%BC%VW(iEdge)%iMemberCell(iCell)
                    !end do
                !end do

                read(iUnit,*) cAnnotate, UG%GM%BC%iWallTotal !"Wall2Edge_data "
                do iLoop = 1, UG%GM%BC%iWallTotal
                    read(iUnit,*) UG%GM%BC%VW(iLoop)%iGlobalEdge, UG%GM%BC%VW(iLoop)%iNumberOfMemberEdge
                    allocate(UG%GM%BC%VW(iLoop)%iMemberEdge(UG%GM%BC%VW(iLoop)%iNumberOfMemberEdge))
                    do iEdge = 1, UG%GM%BC%VW(iLoop)%iNumberOfMemberEdge
                        read(iUnit,*) UG%GM%BC%VW(iLoop)%iMemberEdge(iEdge)
                    end do
                end do


                read(iUnit,*) cAnnotate, cAnnotate ! "Wall_Curvature ", UG%GM%BC%iWallTotal
                do iLoop = 1, UG%GM%BC%iWallTotal
                    read(iUnit,*) UG%GM%BC%VW(iLoop)%curvature
                end do
                call fixed_curvature(UG)
                call fixed_curvature_naca0012(UG)
            !end if
        end if

        close(iUnit)

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

        !call output_curvature(UG)


    return
contains
    subroutine output_curvature(UG)
        implicit none
        type(UnstructuredGrid), intent(in) :: UG

            iUnit = 1
            open(unit=iUnit,file="curvature.txt",status='unknown')
            write(iUnit,*) UG%GM%BC%iWallTotal
            do iLoop = 1, UG%GM%BC%iWallTotal
                write(iUnit,"(3(2x,f22.14))") UG%CD%Edge(UG%GM%BC%VW(iLoop)%iGlobalEdge, 1), UG%CD%Edge(UG%GM%BC%VW(iLoop)%iGlobalEdge, 2), UG%GM%BC%VW(iLoop)%curvature
            end do
            close(iUnit)

            stop
        return
    end subroutine output_curvature

    subroutine fixed_curvature(UG)
        implicit none
        type(UnstructuredGrid), intent(inout) :: UG
        double precision :: x, a,b,c,d,e
            a = -1.06849984d0
            b = 0.36722201d0
            c = 0.00980996d0
            d = 1.30289925d0
            e = 0.08158373d0

            do iLoop = 1, UG%GM%BC%iWallTotal
                x = UG%GM%BC%VW(iLoop)%curvature
                UG%GM%BC%VW(iLoop)%curvature =  a + b * x + c * x ** 2 + d*exp(e*x)
            end do
        return
    end subroutine fixed_curvature

    subroutine fixed_curvature_naca0012(UG)
        implicit none
        type(UnstructuredGrid), intent(inout) :: UG
        double precision :: x, t = 0.12d0
        double precision :: c_radius, y_t_prime, y_t_doubleprime
            do iLoop = 1, UG%GM%BC%iWallTotal
                x = UG%CD%Edge(UG%GM%BC%VW(iLoop)%iGlobalEdge, 1)
                if (x /= 0.0d0) then
                    y_t_prime = t / 0.2d0 * (0.2969d0 * (0.5d0 / sqrt(x)) - 0.1260d0 - 0.3516d0 * (2.0d0 * x) + 0.2843d0 * (3.0d0 * x ** 2) - 0.1015d0 * (4.0d0 * x ** 3))
                    y_t_doubleprime = t / 0.2d0 * (0.2969d0 * (-0.25d0 / sqrt(x)**3) - 0.3516d0 * (2.0d0) + 0.2843d0 * (6.0d0 * x) - 0.1015d0 * (12.0d0 * x ** 2))
                    if (y_t_doubleprime /= 0.0d0) then
                        c_radius = (1.0d0 + y_t_prime)**(1.5d0) / (y_t_doubleprime)
                        UG%GM%BC%VW(iLoop)%curvature = 1.0d0 / abs(c_radius)
                    else
                        UG%GM%BC%VW(iLoop)%curvature = 0.0d0
                    end if
                else
                    UG%GM%BC%VW(iLoop)%curvature = 0.0d0
                end if
            end do

        return
    end subroutine fixed_curvature_naca0012
end subroutine UReadUnStrGrid
