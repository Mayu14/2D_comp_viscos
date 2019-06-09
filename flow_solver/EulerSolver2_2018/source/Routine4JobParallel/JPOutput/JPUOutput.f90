subroutine JPUOutput(UConf,UG,UCC,iStep)
    use StructVar_Mod
    use ConstantVar_Mod, gamma => SpecificOfHeatRatio
    use LoopVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    integer, intent(in) :: iStep
    character(len=256) :: cDirectory,cFileName, cCaseName, ctmpDir
    character(len=256) :: cStep
    integer :: access
    integer, allocatable :: iSurfaceEdge(:, :)
    integer :: iCheck, iUnit_num
    integer :: debug = 0, moreplot = 1
    double precision, allocatable :: curvature(:), normal(:,:), vData(:,:), vcoords(:,:)
    call JPUConserve2Primitive(UG,UCC)

    if(UConf%OutputStatus == 1) then
        write(cStep,'("Final")')
        write(ctmpDir,'("Complete/vtk/")')
    else if(UConf%OutputStatus == 2) then
        write(cStep,'("Resume")')
        write(ctmpDir,'("Resume/")')
    else
        write(cStep,*) iStep
        cStep = trim(adjustl(cStep))//"th"
        write(ctmpDir,'("ResultU/")')
    end if
    cDirectory = trim(adjustl(UConf%cDirectory))//trim(adjustl(ctmpDir)) !UConf%SaveDirectiry
    cFileName = trim(adjustl(cDirectory))//trim(adjustl(UConf%cFileName))//"__"//trim(adjustl(cStep))//".vtk"

    cCaseName = UConf%cCaseName
    iUnit_num = UConf%my_rank + 100

    open(unit = iUnit_num, file =trim(adjustl(cFileName)), status = 'unknown')
        write(iUnit_num,"('# vtk DataFile Version 3.0')")
        write(iUnit_num,*) trim(adjustl(cCaseName))
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
            if(isnan(UCC%PrimitiveVariable(1,iCell,1,1))) then
                RetryFlag = 1
            end if
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

        if (moreplot == 1) then
            !本来は定積比熱Cvをかけるべきだが，ここでは省略している
            write(iUnit_num,"('SCALARS Entropy float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")
            do iCell=1, UG%GI%RealCells
                write(iUnit_num, "(f22.14)") log((UCC%PrimitiveVariable(5,iCell,1,1) / UG%GM%BC%InFlowVariable(5)) &
                                                & / ((UCC%PrimitiveVariable(1,iCell,1,1) / UG%GM%BC%InFlowVariable(1)))**gamma )
            end do

            write(iUnit_num,"('SCALARS SoundSpeed float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")
            do iCell=1, UG%GI%RealCells
                write(iUnit_num, "(f22.14)") sqrt(gamma * UCC%PrimitiveVariable(5,iCell,1,1) &
                                                & / (UCC%PrimitiveVariable(1,iCell,1,1)))
            end do

            write(iUnit_num,"('SCALARS Enthalpy float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")
            do iCell=1, UG%GI%RealCells
                write(iUnit_num, "(f22.14)") InverseGmin1 * (gamma * UCC%PrimitiveVariable(5,iCell,1,1) &
                                                & / (UCC%PrimitiveVariable(1,iCell,1,1)))
            end do

            write(iUnit_num,"('SCALARS InternalEnergy float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")
            do iCell=1, UG%GI%RealCells
                write(iUnit_num, "(f22.14)") InverseGmin1 * (UCC%PrimitiveVariable(5,iCell,1,1) &
                                                & / (UCC%PrimitiveVariable(1,iCell,1,1)))
            end do

            write(iUnit_num,"('SCALARS TotalEnergy float')")
            write(iUnit_num,"('LOOKUP_TABLE default')")
            do iCell=1, UG%GI%RealCells
                write(iUnit_num, "(f22.14)") InverseGmin1 * (gamma * UCC%PrimitiveVariable(5,iCell,1,1) &
                                                & / (UCC%PrimitiveVariable(1,iCell,1,1))) &
                                                & + 0.5d0 * dot_product(UCC%PrimitiveVariable(2:4,iCell,1,1), UCC%PrimitiveVariable(2:4,iCell,1,1))
            end do

        end if

    if(debug == 1) then
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

        write(iUnit_num,"('SCALARS TimeStep float')")
        write(iUnit_num,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(iUnit_num, "(f22.14)") UCC%TimeWidth(iCell,1,1)
        end do

        write(iUnit_num,"('SCALARS Sq_SoundSpeed float')")
        write(iUnit_num,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(iUnit_num, "(f22.14)") gamma*Gmin1*(UCC%ConservedQuantity(5,iCell,1,1) / UCC%ConservedQuantity(1,iCell,1,1) &
                            & - 0.5d0 * dot_product(UCC%ConservedQuantity(2:4,iCell,1,1),UCC%ConservedQuantity(2:4,iCell,1,1))/(UCC%ConservedQuantity(1,iCell,1,1)**2))
        end do

        write(iUnit_num,"('SCALARS curvature float')")
        write(iUnit_num,"('LOOKUP_TABLE default')")
        allocate(curvature(UG%GI%RealCells))
        allocate(iSurfaceEdge(UG%GI%RealCells, 3))
        curvature = 0.0d0
        iSurfaceEdge = 0
        do iEdge=1, UG%GM%BC%iWallTotal
            iCell = UG%Line%Cell(UG%GM%BC%VW(iEdge)%iGlobalEdge, 1, 1)
            !write(6,*) iCell, UG%GI%RealCells
            curvature(iCell) = UG%GM%BC%VW(iEdge)%curvature
            iSurfaceEdge(iCell, 1) = iEdge  ! 局所界面番号
            iSurfaceEdge(iCell, 2) = UG%GM%BC%VW(iEdge)%iGlobalEdge ! 大域界面番号
            iSurfaceEdge(iCell, 3) = UG%Line%Cell(UG%GM%BC%VW(iEdge)%iGlobalEdge, 1, 1) ! 隣接実セル大域番号
        end do

        do iCell = 1, UG%GI%RealCells
            write(iUnit_num,"((2x,f22.14))") curvature(iCell)
        end do

        write(iUnit_num,"('VECTORS iSurfE float')")
        do iCell = 1, UG%GI%RealCells
            write(iUnit_num,"(3(2x,f22.14))") float(iSurfaceEdge(iCell,1)),float(iSurfaceEdge(iCell,2)),float(iSurfaceEdge(iCell,3))
        end do

        write(iUnit_num,"('VECTORS Normal float')")
        allocate(normal(UG%GI%RealCells,3))
        normal = 0.0d0
        do iLocalEdge=1, UG%GM%BC%iWallTotal
            iEdge = UG%GM%BC%VW(iLocalEdge)%iGlobalEdge
            iCell = UG%Line%Cell(iEdge,1,1)
            normal(iCell, :) = UG%GM%Normal(iEdge,:)
        end do
        do iCell = 1, UG%GI%RealCells
            write(iUnit_num,"(3(2x,f22.14))") normal(iCell,1),normal(iCell,2),normal(iCell,3)
        end do

        write(iUnit_num,"('VECTORS vData float')")
        allocate(vData(UG%GI%RealCells, 3))
        allocate(vcoords(UG%GI%RealCells,3))
        vData = 0.0d0
        do iLocalEdge=1,UG%GM%BC%iWallTotal
            iEdge = UG%GM%BC%VW(iLocalEdge)%iGlobalEdge
            iCell = UG%Line%Cell(iEdge,1,1)
            iAdjacentCell = UG%Line%Cell(iEdge,2,1)
            vData(iCell, 1) = float(iAdjacentCell)
            vData(iCell, 2) = sqrt(dot_product((UG%GM%Width(iCell, UG%Line%Cell(iEdge,1,2), :)), (UG%GM%Width(iCell, UG%Line%Cell(iEdge,1,2), :))))
            vData(iCell, 3) = sqrt(dot_product(UG%GM%Width(iAdjacentCell, UG%Line%Cell(iEdge,2,2), :), UG%GM%Width(iAdjacentCell, UG%Line%Cell(iEdge,2,2), :)))
            vcoords(iCell,1:3) = UG%CD%Cell(iAdjacentCell,1:3)
        end do
        do iCell = 1, UG%GI%RealCells
            write(iUnit_num,"(3(2x,f22.14))") vData(iCell,1),vData(iCell,2),vData(iCell,3)
        end do
        write(iUnit_num,"('VECTORS vCoords float')")
        do iCell = 1, UG%GI%RealCells
            write(iUnit_num,"(3(2x,f22.14))") vcoords(iCell,1),vcoords(iCell,2),vcoords(iCell,3)
        end do
    end if
    close(iUnit_num)
    !if(RetryFlag == 0) then
        !UCC%PastQuantity = UCC%ConservedQuantity
    !else
        !UCC%ConservedQuantity = UCC%PastQuantity
        !CourantFriedrichsLewyCondition = 0.1d0*CourantFriedrichsLewyCondition
        !RetryFlag = 0
        ! write(6,*) CourantFriedrichsLewyCondition
    !end if
!stop
    if (UConf%OutputStatus /= 0) then  ! Final & Resume出力時にdummyデータが
        write(cStep,'("Resume")')
        write(ctmpDir,'("Running/")')
        cDirectory = trim(adjustl(UConf%cDirectory))//trim(adjustl(ctmpDir)) !UConf%SaveDirectiry
        cFileName = trim(adjustl(cDirectory))//trim(adjustl(UConf%cFileName))//"__"//trim(adjustl(cStep))//".vtk"
        if(access(cFileName, " ") == 0) then    ! 存在していれば
            call system("rm "//trim(adjustl(cFileName)))    ! dummyデータを削除する
        end if
    end if

    return
end subroutine JPUOutput
