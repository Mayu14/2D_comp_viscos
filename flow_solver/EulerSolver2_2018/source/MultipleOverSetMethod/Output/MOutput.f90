subroutine MOutput(UConf,UG,NSRP,UCC,iStep,iTargetGrid,MG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(AnotherCoordinateRelation), intent(in) :: NSRP
    type(CellCenter), intent(inout) :: UCC
    integer, intent(in) :: iStep
    integer, intent(in) :: iTargetGrid
    type(MoveGrid), intent(in) :: MG
    character(len=64) :: cDirectory,cFileName, cCaseName
    character(len=32) :: cStep
    character(len=32) :: cGrid

    call UConserve2Primitive(UG,UCC)

    write(cStep,*) iStep
    write(cGrid,*) iTargetGrid

    cDirectory = "Result"//trim(adjustl(cGrid))//"/" !UConf%SaveDirectiry
    cFileName = trim(adjustl(cDirectory))//"PrimitiveVariables"//trim(adjustl(cStep))//"th_Result.vtk"
    cCaseName = "UnstructuredShockTube" !UConf%CaseName
    open(unit = 1, file =trim(adjustl(cFileName)), status = 'unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,*) cCaseName
        write(1,"('ASCII')")
        write(1,"('DATASET UNSTRUCTURED_GRID')")
        write(1,"('POINTS ',(1x,i7),' double')") UG%GI%Points

        do iPoint=1, UG%GI%Points
            write(1,"(3(1x,f22.17))") NSRP%NextCoordsInG(iPoint,1),NSRP%NextCoordsInG(iPoint,2),NSRP%NextCoordsInG(iPoint,3)+1.0/iTargetGrid !for OverSet
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

        write(1,"('VECTORS Limiter float')")
        do iCell=1, UG%GI%RealCells
            write(1,"((2x,f22.14))") UCC%LimiterFunction(1,iCell,1,1),UCC%LimiterFunction(2,iCell,1,1),UCC%LimiterFunction(5,iCell,1,1)
        end do

        write(1,"('SCALARS InfluenceDepth float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(1,"((1x,i5))") MG%IS(iCell)%InfluenceDepth
        end do

        write(1,"('SCALARS Interpolated float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(1,"((1x,i5))") UG%GM%Interpolated(iCell,1,1)
        end do


    close(1)

return
end subroutine MOutput

