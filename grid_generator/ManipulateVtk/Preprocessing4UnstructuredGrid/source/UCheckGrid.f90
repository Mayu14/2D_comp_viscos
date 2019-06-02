!***********************************/
!	Name:非構造格子用EulerSolverのための仮想格子並び換えプログラム
!	Alias:UCheckGrid
!	Description:仮想セルの中心座標を基準に仮想セルの番号を振り直す(外部→内部)
!	Type:UnstructuredGrid
!	Input:UG%GM%Dimension,UG%xxx%Point (xxxは使用する種類のセルデータ)
!	Output:UG
!	Note:中心からの距離を算出して，その値の大きい順に番号を振り直す
!	Author:Akitaka Toyota
!	Date:2018.11.07
!	Update:-
!	Other:
!***********************************/
subroutine UCheckGrid(UG)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG
    character(len=256) :: cDirectory,cFileName, cCaseName
    integer :: iCheck
    double precision, allocatable :: curvature(:)

    cDirectory = ""
    cFileName = "chkgrid.vtk"
    cCaseName = "UnstructuredChkGrid" !UConf%CaseName
    open(unit = 1, file =trim(adjustl(cFileName)), status = 'unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,*) cCaseName
        write(1,"('ASCII')")
        write(1,"('DATASET UNSTRUCTURED_GRID')")
        write(1,"('POINTS ',(1x,i7),' double')") UG%GI%Points

        do iPoint=1, UG%GI%Points
            write(1,"(3(1x,f22.17))") UG%CD%Point(iPoint,1),UG%CD%Point(iPoint,2),UG%CD%Point(iPoint,3)
        end do
        write(1,*) ""

        write(1,"('CELLS ',2(1x,i7))") UG%GI%RealCells,UG%GI%RealCells*4
        do iCell=1, UG%GI%RealCells
                write(1,"(4(1x,i7))") 3,(UG%Tri%Point(iCell,iLoop)-1,iLoop=1,3) !Adjusted Point Number
        end do
        write(1,*) ""

        write(1,"('CELL_TYPES ',(1x,i7))") UG%GI%RealCells
        do iCell=1, UG%GI%RealCells
            write(1,"(1x,i1)") 5
        end do

        write(1,"('CELL_DATA ',i7)") UG%GI%RealCells

        write(1,"('SCALARS BoundaryType float')")
        write(1,"('LOOKUP_TABLE default')")

        do iCell=1, UG%GI%RealCells
            iCheck = 0
            do iLoop=1,3
                if(UG%Tri%Cell(iCell,iLoop) > UG%GI%RealCells) then
                    iCheck = iCheck + 1
                    iAdjacentCell = UG%Tri%Cell(iCell,iLoop)
                end if
            end do

            if(iCheck == 0) then
                write(1,"((2x,f22.14))") 0.0d0
            else
                !write(6,*) iCell,dot_product(UCC%ConservedQuantity(2:4,iCell,1,1),-UG%GM%Normal(UG%VC%Edge(iAdjacentCell),1:3))
                write(1,"((2x,f22.14))") float(UG%VC%Type(iAdjacentCell))
            end if
        end do
        !write(1,*) ""

        write(1,"('SCALARS curvature float')")
        write(1,"('LOOKUP_TABLE default')")
        allocate(curvature(UG%GI%RealCells))
        curvature = 0.0d0
        do iEdge=1, UG%GM%BC%iWallTotal
            iCell = UG%Line%Cell(UG%GM%BC%VW(iEdge)%iGlobalEdge, 1, 1)
            !write(6,*) iCell, UG%GI%RealCells
            curvature(iCell) = UG%GM%BC%VW(iEdge)%curvature
        end do
        do iCell = 1, UG%GI%RealCells
            write(1,"((2x,f22.14))") curvature(iCell)
        end do
        !write(1,*) ""
        write(1,"('SCALARS local_edge_id float')")
        write(1,"('LOOKUP_TABLE default')")
        curvature =0.0d0
        do iEdge=1, UG%GM%BC%iWallTotal
            iCell = UG%Line%Cell(UG%GM%BC%VW(iEdge)%iGlobalEdge, 1, 1)
            curvature(iCell) = dble(iEdge)
        end do
        do iCell = 1, UG%GI%RealCells
            write(1,"((2x,f22.14))") curvature(iCell)
        end do

        write(1,"('VECTORS CellDistance float')")
        do iCell=1, UG%GI%RealCells
            write(1, "(3(1x,f22.17))") float(UG%Tri%Belongs2Wall(iCell)), UG%Tri%Distance(iCell), 0.0
        end do
        !write(1,*) ""



    close(1)

return
end subroutine UCheckGrid
