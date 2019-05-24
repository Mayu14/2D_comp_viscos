subroutine UReadInitialCondition(UConf,UG,UCC)
    use StructVar_Mod
    use ConstantVar_Mod
    use LoopVar_Mod
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    character(len=256) :: cFileName, cDirectory, cStep
    character :: cAnnotate,cAnnotate1
    integer :: iUnit

    if(UConf%SwitchProgram == 5) then
        !cFileName = trim(adjustl(UConf%cDirectory))//trim(adjustl("ResultU/"))//UConf%ResumeFileName
        cFileName = UConf%ResumeFileName
    else if(UConf%SwitchProgram /= 7) then
        write(6,*) "Please input VTK-File name of Unstructured Grid for RESUME"
        read(5,*) cFileName
    else
        write(cStep,*) 0
        cDirectory = trim(adjustl(UConf%cDirectory))//trim(adjustl("ResultU/"))
        cFileName = trim(adjustl(cDirectory))//trim(adjustl(UConf%cFileName))//trim(adjustl("_"))//trim(adjustl(cStep))//"th.vtk"
    end if

    iUnit = UConf%my_rank + 100
    open(unit = iUnit, file =trim(adjustl(cFileName)), status = 'old')
        write(6,*) "Find VTK-Data"
        do iLoop=1,5
            read(iUnit,*) cAnnotate !Header
        end do

        do iPoint=1, UG%GI%Points
            read(iUnit,*) cAnnotate,cAnnotate,cAnnotate !Point Coordinate
        end do

        read(iUnit,*) cAnnotate,cAnnotate
        do iCell=1, UG%GI%RealCells
            read(iUnit,*) cAnnotate,cAnnotate,cAnnotate,cAnnotate
        end do

        read(iUnit,*) cAnnotate
        do iCell=1, UG%GI%RealCells
            read(iUnit,*) cAnnotate
        end do

        read(iUnit,*) cAnnotate
!密度プロット
        read(iUnit,*) cAnnotate
        read(iUnit,*) cAnnotate

        do iCell=1, UG%GI%RealCells
            read(iUnit,*) UCC%PrimitiveVariable(1,iCell,1,1)
        end do

        read(iUnit,*) cAnnotate
        do iCell=1, UG%GI%RealCells
            read(iUnit,*) (UCC%PrimitiveVariable(iLoop,iCell,1,1),iLoop=2,4)
        end do

        read(iUnit,*) cAnnotate
        read(iUnit,*) cAnnotate
        do iCell=1, UG%GI%RealCells
            read(iUnit,*) UCC%PrimitiveVariable(UG%GM%Dimension+2,iCell,1,1)
        end do

    close(iUnit)

return
end subroutine UReadInitialCondition

