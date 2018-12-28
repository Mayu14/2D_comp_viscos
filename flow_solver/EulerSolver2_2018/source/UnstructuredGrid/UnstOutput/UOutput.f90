subroutine UOutput(UConf,UG,UCC,iStep)
    use StructVar_Mod
    use ConstantVar_Mod
    use LoopVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    integer, intent(in) :: iStep
    character(len=64) :: cDirectory,cFileName, cCaseName
    character(len=32) :: cStep
    integer :: iCheck

    call UConserve2Primitive(UG,UCC)

    write(cStep,*) iStep
    cDirectory = "ResultU/" !UConf%SaveDirectiry
    cFileName = trim(adjustl(cDirectory))//"PrimitiveVariables"//trim(adjustl(cStep))//"th_Result.vtk"
    cCaseName = "UnstructuredShockTube" !UConf%CaseName
    open(unit = 1, file =trim(adjustl(cFileName)), status = 'unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,*) cCaseName
        write(1,"('ASCII')")
        write(1,"('DATASET UNSTRUCTURED_GRID')")
        write(1,"('POINTS ',(1x,i7),' double')") UG%GI%Points

        do iPoint=1, UG%GI%Points
            !write(1,"(3(1x,f22.17))") (UG%CD%Point(iPoint,iLoop),iLoop=1,3) !Adjust Point Number
            write(1,"(3(1x,f22.17))") UG%CD%Point(iPoint,1),UG%CD%Point(iPoint,2),UG%CD%Point(iPoint,3)+1 !for OverSet
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
!密度プロット
        write(1,"('SCALARS Density float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(1,"((2x,f22.14))") UCC%PrimitiveVariable(1,iCell,1,1)
        end do


        write(1,"('VECTORS Velocity float')")
        if(UG%GM%Dimension == 2) then
            do iCell=1, UG%GI%RealCells
                write(1,"((2x,f22.14))") (UCC%PrimitiveVariable(iLoop,iCell,1,1),iLoop=2,3),0.0d0
            end do
        else if(UG%GM%Dimension == 3) then
            do iCell=1, UG%GI%RealCells
                write(1,"((2x,f22.14))") (UCC%PrimitiveVariable(iLoop,iCell,1,1),iLoop=2,4)
            end do
        end if

        write(1,"('SCALARS Pressure float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(1,"((2x,f22.14))") UCC%PrimitiveVariable(UG%GM%Dimension+2,iCell,1,1)
        end do

        if(iStep > 0) then
            write(1,"('SCALARS BoundaryCondition float')")
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
                    write(1,"((2x,f22.14))") dot_product(UCC%ConservedQuantity(2:4,iCell,1,1),-UG%GM%Normal(UG%VC%Edge(iAdjacentCell),1:3))
                end if
            end do

        else
            write(1,"('SCALARS BoundaryCondition float')")
            write(1,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                write(1,"((2x,f22.14))") 0.0d0
            end do
        end if

        write(1,"('SCALARS WallNumber float')")
        write(1,"('LOOKUP_TABLE default')")

        do iCell=1, UG%GI%RealCells
            write(1, "(f22.14)") dble(UG%Line%Belongs2Wall(UG%Tri%Edge(iCell, 1))+UG%Line%Belongs2Wall(UG%Tri%Edge(iCell, 2))+UG%Line%Belongs2Wall(UG%Tri%Edge(iCell, 3))) / 3.0d0
        end do

        write(1,"('SCALARS WallDistance float')")
        write(1,"('LOOKUP_TABLE default')")

        do iCell=1, UG%GI%RealCells
            write(1, "(f22.14)") dble(UG%Line%Distance(UG%Tri%Edge(iCell, 1))+UG%Line%Distance(UG%Tri%Edge(iCell, 2))+UG%Line%Distance(UG%Tri%Edge(iCell, 3))) / 3.0d0
        end do

        if(invicid == .true.) then
            write(1,"('SCALARS LaminarViscosity float')")
            write(1,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                write(1, "(f22.14)") UCC%LaminarViscosity(iCell,1,1)
            end do

            write(1,"('SCALARS EddyViscosity float')")
            write(1,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                write(1, "(f22.14)") UCC%EddyViscosity(iCell,1,1)
            end do
        end if

    close(1)
return
end subroutine UOutput

