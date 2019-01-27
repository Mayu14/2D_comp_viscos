!***********************************/
!	Name:数値実験のまとめ
!	Alias:PreNumericalEx
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:
!	Date:
!	Update:
!	Other:
!***********************************/
subroutine UOutputExact(UConf,UG,NE,iStep)
    use StructVar_Mod
    use LoopVar_Mod
    use ExactSol
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(NumericalExperiment), intent(in) :: NE
    integer, intent(in) :: iStep
    character(len=64) :: cDirectory,cFileName, cCaseName
    character(len=32) :: cStep
    double precision :: Width
    integer :: iDir

    Width = 1.0d0/UG%GI%AllCells
    iDir = UConf%UseResume - 1

    write(cStep,*) iStep
    cDirectory = "EXP/" !UConf%SaveDirectiry
    cFileName = trim(adjustl(cDirectory))//"ExactSolution"//trim(adjustl(cStep))//"th_Result.vtk"
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
        write(1,"('SCALARS ExDensity float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            iPoint = nint(UG%CD%Cell(iCell,iDir)/Width)
            write(1,"((2x,f22.14))") NE%STT%PrimitiveVariable(1,iPoint)
        end do

        write(1,"('VECTORS ExVelocity float')")

            do iCell=1, UG%GI%RealCells
                iPoint = nint(UG%CD%Cell(iCell,iDir)/Width)
                write(1,"((2x,f22.14))") (NE%STT%PrimitiveVariable(iLoop,iPoint),iLoop=2,3),0.0d0
            end do


        write(1,"('SCALARS ExPressure float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            iPoint = nint(UG%CD%Cell(iCell,iDir)/Width)
            write(1,"((2x,f22.14))") NE%STT%PrimitiveVariable(4,iPoint)
        end do

    close(1)

    return
end subroutine UOutputExact
