subroutine JPUOutput(UConf,UG,UCC,iStep)
    use StructVar_Mod
    use ConstantVar_Mod, gamma => SpecificOfHeatRatio
    use LoopVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    integer, intent(in) :: iStep
    character(len=64) :: cDirectory,cFileName, cCaseName
    character(len=32) :: cStep

    integer :: iCheck, iUnit_num

    call JPUConserve2Primitive(UG,UCC)

    write(cStep,*) iStep
    if(UConf%CalcEnv == 0) then
        cDirectory = "ResultU/" !UConf%SaveDirectiry   ! 研究室PC用
        cFileName = trim(adjustl(cDirectory))//trim(adjustl(UConf%cFileName))//trim(adjustl("_"))//trim(adjustl(cStep))//"th.vtk"
    else if(UConf%CalcEnv == 1) then
        cDirectory = "/work/A/FMa/FMa037/Case1/ResultU/" ! 東北大スパコン用
        cFileName = trim(adjustl(cDirectory))//trim(adjustl(UConf%cFileName))//trim(adjustl("_"))//trim(adjustl(cStep))//"th.vtk"
    end if

    cCaseName = "UnstructuredShockTube" !UConf%CaseName
    iUnit_num = UConf%my_rank + 100

    open(unit = iUnit_num, file =trim(adjustl(cFileName)), status = 'unknown')
        write(iUnit_num,"('# vtk DataFile Version 3.0')")
        write(iUnit_num,*) cCaseName
        write(iUnit_num,"('ASCII')")
        write(iUnit_num,"('DATASET UNSTRUCTURED_GRID')")
        write(iUnit_num,"('POINTS ',(1x,i7),' double')") UG%GI%Points

        do iPoint=1, UG%GI%Points
            !write(iUnit_num,"(3(1x,f22.17))") (UG%CD%Point(iPoint,iLoop),iLoop=1,3) !Adjust Point Number
            write(iUnit_num,"(3(1x,f22.17))") UG%CD%Point(iPoint,1),UG%CD%Point(iPoint,2),UG%CD%Point(iPoint,3)+1 !for OverSet
        end do
        write(iUnit_num,*) ""

        write(iUnit_num,"('CELLS ',2(1x,i7))") UG%GI%RealCells,UG%GI%RealCells*4
        do iCell=1, UG%GI%RealCells
                write(iUnit_num,"(4(1x,i7))") 3,(UG%Tri%Point(iCell,iLoop)-1,iLoop=1,3) !Adjusted Point Number
        end do
        write(iUnit_num,*) ""

        write(iUnit_num,"('CELL_TYPES ',(1x,i7))") UG%GI%RealCells
        do iCell=1, UG%GI%RealCells
            write(iUnit_num,"(1x,i1)") 5
        end do

        write(iUnit_num,"('CELL_DATA ',i7)") UG%GI%RealCells
!密度プロット
        write(iUnit_num,"('SCALARS Density float')")
        write(iUnit_num,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(iUnit_num,"((2x,f22.14))") UCC%PrimitiveVariable(1,iCell,1,1)
        end do

        write(iUnit_num,"('VECTORS Velocity float')")
        do iCell=1, UG%GI%RealCells
            write(iUnit_num,"((2x,f22.14))") (UCC%PrimitiveVariable(iLoop,iCell,1,1),iLoop=2,4)
        end do

        write(iUnit_num,"('SCALARS Pressure float')")
        write(iUnit_num,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(iUnit_num,"((2x,f22.14))") UCC%PrimitiveVariable(5,iCell,1,1)
        end do

        if(iStep > 0) then
            write(iUnit_num,"('SCALARS BoundaryCondition float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                iCheck = 0
                do iLoop=1,3
                    if(UG%Tri%Cell(iCell,iLoop) > UG%GI%RealCells) then
                        iCheck = iCheck + 1
                        iAdjacentCell = UG%Tri%Cell(iCell,iLoop)
                    end if
                end do

                if(iCheck == 0) then
                    write(iUnit_num,"((2x,f22.14))") 0.0d0
                else
                    !write(6,*) iCell,dot_product(UCC%ConservedQuantity(2:4,iCell,1,1),-UG%GM%Normal(UG%VC%Edge(iAdjacentCell),1:3))
                    write(iUnit_num,"((2x,f22.14))") dot_product(UCC%ConservedQuantity(2:4,iCell,1,1),-UG%GM%Normal(UG%VC%Edge(iAdjacentCell),1:3))
                end if
            end do

            write(iUnit_num,"('SCALARS BoundaryType float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                iCheck = 0
                do iLoop=1,3
                    if(UG%Tri%Cell(iCell,iLoop) > UG%GI%RealCells) then
                        iCheck = iCheck + 1
                        iAdjacentCell = UG%Tri%Cell(iCell,iLoop)
                    end if
                end do

                if(iCheck == 0) then
                    write(iUnit_num,"((2x,f22.14))") 0.0d0
                else
                    !write(6,*) iCell,dot_product(UCC%ConservedQuantity(2:4,iCell,1,1),-UG%GM%Normal(UG%VC%Edge(iAdjacentCell),1:3))
                    write(iUnit_num,"((2x,f22.14))") float(UG%GM%CellType(iAdjacentCell, 1, 1))
                end if
            end do
        else
            write(iUnit_num,"('SCALARS BoundaryCondition float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                write(iUnit_num,"((2x,f22.14))") 0.0d0
            end do
        end if

        write(iUnit_num,"('SCALARS Wall float')")
        write(iUnit_num,"('LOOKUP_TABLE default')")

        do iCell=1, UG%GI%RealCells
            write(iUnit_num, "(f22.14)") dble((UG%Line%Belongs2Wall(UG%Tri%Edge(iCell, 1)) + UG%Line%Belongs2Wall(UG%Tri%Edge(iCell, 2)) + UG%Line%Belongs2Wall(UG%Tri%Edge(iCell, 3)))) / 3.0d0
        end do

        write(iUnit_num,"('VECTORS Wall2 float')")
        do iCell=1, UG%GI%RealCells
            iPoint = 0  ! iFlag
            do iLoop=1, UG%GM%BC%iWallTotal
                iAdjacentCell = UG%Line%Cell(UG%GM%BC%VW(iLoop)%iGlobalEdge, 1, 1)
                if(iAdjacentCell == iCell) then
                    iEdge = iAdjacentCell
                    iLocalEdge = iLoop
                    iPoint = 1
                end if
            end do
            if(iPoint == 1) then
                write(iUnit_num, "((2x,f22.14))") dble(iEdge), dble(UCC%debug(iCell,1)), dble(0.0)
            else
                write(iUnit_num, "((2x,f22.14))") dble(0.0), dble(0.0), dble(0.0)
            end if
        end do

        write(iUnit_num,"('SCALARS WallDistance float')")
        write(iUnit_num,"('LOOKUP_TABLE default')")

        do iCell=1, UG%GI%RealCells
            write(iUnit_num, "(f22.14)") dble(UG%Line%Distance(UG%Tri%Edge(iCell, 1))+UG%Line%Distance(UG%Tri%Edge(iCell, 2))+UG%Line%Distance(UG%Tri%Edge(iCell, 3))) / 3.0d0
        end do

        if(invicid == .False.) then
            write(iUnit_num,"('SCALARS AbsVortisity float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                write(iUnit_num, "(f22.14)") UCC%AbsoluteVortisity(iCell,1,1)
            end do

            write(iUnit_num,"('VECTORS StrainRateTensorUxyz float')")
            do iCell=1, UG%GI%RealCells
                write(iUnit_num,"((2x,f22.14))") UCC%StrainRateTensor(1,1,iCell,1,1),UCC%StrainRateTensor(1,2,iCell,1,1),0.0d0
            end do

            write(iUnit_num,"('VECTORS StrainRateTensorVxyz float')")
            do iCell=1, UG%GI%RealCells
                write(iUnit_num,"((2x,f22.14))") UCC%StrainRateTensor(2,1,iCell,1,1),UCC%StrainRateTensor(2,2,iCell,1,1),0.0d0
            end do

            write(iUnit_num,"('SCALARS Temparature float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                write(iUnit_num, "(f22.14)") UCC%Temparature(iCell,1,1)
            end do

            write(iUnit_num,"('SCALARS LaminarViscosity float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                write(iUnit_num, "(f22.14)") UCC%LaminarViscosity(iCell,1,1)
            end do

            write(iUnit_num,"('SCALARS EddyViscosity float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                write(iUnit_num, "(f22.14)") UCC%EddyViscosity(iCell,1,1)
            end do

        end if

        write(iUnit_num,"('SCALARS MachNumber float')")
        write(iUnit_num,"('LOOKUP_TABLE default')")

        do iCell=1, UG%GI%RealCells
            write(iUnit_num, "(f22.14)") sqrt(UCC%PrimitiveVariable(1,iCell,1,1) / (gamma*UCC%PrimitiveVariable(5,iCell,1,1)) &
                                        &   * dot_product(UCC%PrimitiveVariable(2:3,iCell,1,1),UCC%PrimitiveVariable(2:3,iCell,1,1)))
        end do
    close(iUnit_num)

    return
end subroutine JPUOutput
