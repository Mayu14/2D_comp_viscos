!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:OOverSetS2U
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
subroutine OComputeNextCoords(iStep,UG,MG,UCC)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    integer, intent(in) :: iStep
    type(UnstructuredGrid), intent(inout) :: UG
    type(MoveGrid), intent(inout) :: MG
    type(CellCenter), intent(in) :: UCC

    MG%NSR%C%PresentCoordsInG = MG%NSR%C%NextCoordsInG
    MG%NSR%P%PresentCoordsInG = MG%NSR%P%NextCoordsInG

    do iCell=1, UG%GI%AllCells
        MG%NSR%C%NextCoordsInG(iCell,:) = matmul(MG%TM%Next2Top,MG%RC%Cell(iCell,:))
    end do

    do iCell=1, UG%GI%AllCells
        MG%NSR%C%NextCoordsInN(iCell,:) = matmul(MG%TM%Next2Current,MG%RC%Cell(iCell,:))
    end do

    do iPoint=1, UG%GI%Points
        MG%NSR%P%NextCoordsInG(iPoint,:) = matmul(MG%TM%Next2Top,MG%RC%Point(iPoint,:))
    end do

    do iPoint=1, UG%GI%Points
        MG%NSR%P%NextCoordsInN(iPoint,:) = matmul(MG%TM%Next2Current,MG%RC%Point(iPoint,:))
    end do

!    if(iStep < 200) then
        call TestOutputVTK(iStep)
 !   end if

return
contains

    subroutine TestOutputVTK(iStep)
    implicit none
    integer, intent(in) :: iStep
    character(len=64) :: cDirectory,cFileName, cCaseName
    character(len=32) :: cStep

    write(cStep,*) iStep
    cDirectory = "ResultR/" !UConf%SaveDirectiry
    cFileName = trim(adjustl(cDirectory))//"TestRotate"//trim(adjustl(cStep))//"th_Result.vtk"
    cCaseName = "Rotation" !UConf%CaseName

    open(unit = 1, file =trim(adjustl(cFileName)), status = 'unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,*) cCaseName
        write(1,"('ASCII')")
        write(1,"('DATASET UNSTRUCTURED_GRID')")
        write(1,"('POINTS ',(1x,i7),' double')") UG%GI%Points

        do iPoint=1, UG%GI%Points
            !write(1,"(3(1x,f22.17))") (UG%CD%Point(iPoint,iLoop),iLoop=1,3) !Adjust Point Number
            write(1,"(3(1x,f22.17))") MG%NSR%P%NextCoordsInG(iPoint,1),MG%NSR%P%NextCoordsInG(iPoint,2),MG%NSR%P%NextCoordsInG(iPoint,3)+1 !for OverSet
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

        write(1,"('SCALARS TEST_CELL float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(1,"((2x,f22.14))") UCC%ConservedQuantity(1,iCell,1,1)
        end do

    return
    end subroutine TestOutputVTK

end subroutine OComputeNextCoords
