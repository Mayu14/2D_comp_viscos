subroutine JPUOutput(UConf,UG,UCC,iStep)
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
    integer :: UseRelativeSpeed = 0
    double precision, allocatable :: RelativeSpeed(:)
    integer :: iCheck

    call JPUConserve2Primitive(UG,UCC)

    write(cStep,*) iStep
    cDirectory = "ResultU/" !UConf%SaveDirectiry   ! 研究室PC用
    !cDirectory = "/work/A/FMa/FMa037/Case1/ResultU/" ! 東北大スパコン用

    cFileName = trim(adjustl(cDirectory))//trim(adjustl(UConf%cFileName))//trim(adjustl("_"))//trim(adjustl(cStep))//"th.vtk"
    cCaseName = "UnstructuredShockTube" !UConf%CaseName
    open(unit = UConf%my_rank + 100, file =trim(adjustl(cFileName)), status = 'unknown')
        write(UConf%my_rank + 100,"('# vtk DataFile Version 3.0')")
        write(UConf%my_rank + 100,*) cCaseName
        write(UConf%my_rank + 100,"('ASCII')")
        write(UConf%my_rank + 100,"('DATASET UNSTRUCTURED_GRID')")
        write(UConf%my_rank + 100,"('POINTS ',(1x,i7),' double')") UG%GI%Points

        do iPoint=1, UG%GI%Points
            !write(UConf%my_rank + 100,"(3(1x,f22.17))") (UG%CD%Point(iPoint,iLoop),iLoop=1,3) !Adjust Point Number
            write(UConf%my_rank + 100,"(3(1x,f22.17))") UG%CD%Point(iPoint,1),UG%CD%Point(iPoint,2),UG%CD%Point(iPoint,3)+1 !for OverSet
        end do
        write(UConf%my_rank + 100,*) ""

        write(UConf%my_rank + 100,"('CELLS ',2(1x,i7))") UG%GI%RealCells,UG%GI%RealCells*4
        do iCell=1, UG%GI%RealCells
                write(UConf%my_rank + 100,"(4(1x,i7))") 3,(UG%Tri%Point(iCell,iLoop)-1,iLoop=1,3) !Adjusted Point Number
        end do
        write(UConf%my_rank + 100,*) ""

        write(UConf%my_rank + 100,"('CELL_TYPES ',(1x,i7))") UG%GI%RealCells
        do iCell=1, UG%GI%RealCells
            write(UConf%my_rank + 100,"(1x,i1)") 5
        end do

        write(UConf%my_rank + 100,"('CELL_DATA ',i7)") UG%GI%RealCells
!密度プロット
        write(UConf%my_rank + 100,"('SCALARS Density float')")
        write(UConf%my_rank + 100,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(UConf%my_rank + 100,"((2x,f22.14))") UCC%PrimitiveVariable(1,iCell,1,1)
        end do

        if(UseRelativeSpeed == 0) then
            write(UConf%my_rank + 100,"('VECTORS Velocity float')")
            if(UG%GM%Dimension == 2) then
                do iCell=1, UG%GI%RealCells
                    write(UConf%my_rank + 100,"((2x,f22.14))") (UCC%PrimitiveVariable(iLoop,iCell,1,1),iLoop=2,3),0.0d0
                end do
            else if(UG%GM%Dimension == 3) then
                do iCell=1, UG%GI%RealCells
                    write(UConf%my_rank + 100,"((2x,f22.14))") (UCC%PrimitiveVariable(iLoop,iCell,1,1),iLoop=2,4)
                end do
            end if

            write(UConf%my_rank + 100,"('SCALARS Pressure float')")
            write(UConf%my_rank + 100,"('LOOKUP_TABLE default')")
            do iCell=1, UG%GI%RealCells
                write(UConf%my_rank + 100,"((2x,f22.14))") UCC%PrimitiveVariable(UG%GM%Dimension+2,iCell,1,1)
            end do

        else
            allocate(RelativeSpeed(3))
            RelativeSpeed(1) = -0.3d0
            RelativeSpeed(2) = 0.0d0
            RelativeSpeed(3) = 0.0d0
            write(UConf%my_rank + 100,"('VECTORS Velocity float')")
            if(UG%GM%Dimension == 2) then
                do iCell=1, UG%GI%RealCells
                    write(UConf%my_rank + 100,"((2x,f22.14))") (UCC%PrimitiveVariable(iLoop,iCell,1,1)-RelativeSpeed(iLoop-1),iLoop=2,3),0.0d0
                end do
            else if(UG%GM%Dimension == 3) then
                do iCell=1, UG%GI%RealCells
                    write(UConf%my_rank + 100,"((2x,f22.14))") (UCC%PrimitiveVariable(iLoop,iCell,1,1)-RelativeSpeed(iLoop-1),iLoop=2,4)
                end do
            end if

            write(UConf%my_rank + 100,"('SCALARS Pressure float')")
            write(UConf%my_rank + 100,"('LOOKUP_TABLE default')")
            do iCell=1, UG%GI%RealCells
                write(UConf%my_rank + 100,"((2x,f22.14))") Gmin1*(UCC%ConservedQuantity(UG%GM%Dimension+2,iCell,1,1)-0.5d0*UCC%ConservedQuantity(1,iCell,1,1)*&
                &(dot_product(UCC%PrimitiveVariable(2:4,iCell,1,1)-RelativeSpeed(1:3),UCC%PrimitiveVariable(2:4,iCell,1,1)-RelativeSpeed(1:3))))
            end do
        end if

        if(iStep > 0) then
            write(UConf%my_rank + 100,"('SCALARS BoundaryCondition float')")
            write(UConf%my_rank + 100,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                iCheck = 0
                do iLoop=1,3
                    if(UG%Tri%Cell(iCell,iLoop) > UG%GI%RealCells) then
                        iCheck = iCheck + 1
                        iAdjacentCell = UG%Tri%Cell(iCell,iLoop)
                    end if
                end do

                if(iCheck == 0) then
                    write(UConf%my_rank + 100,"((2x,f22.14))") 0.0d0
                else
                    !write(6,*) iCell,dot_product(UCC%ConservedQuantity(2:4,iCell,1,1),-UG%GM%Normal(UG%VC%Edge(iAdjacentCell),1:3))
                    write(UConf%my_rank + 100,"((2x,f22.14))") dot_product(UCC%ConservedQuantity(2:4,iCell,1,1),-UG%GM%Normal(UG%VC%Edge(iAdjacentCell),1:3))
                end if
            end do
        else
            write(UConf%my_rank + 100,"('SCALARS BoundaryCondition float')")
            write(UConf%my_rank + 100,"('LOOKUP_TABLE default')")

            do iCell=1, UG%GI%RealCells
                write(UConf%my_rank + 100,"((2x,f22.14))") 0.0d0
            end do
        end if
    close(UConf%my_rank + 100)
return
end subroutine JPUOutput

