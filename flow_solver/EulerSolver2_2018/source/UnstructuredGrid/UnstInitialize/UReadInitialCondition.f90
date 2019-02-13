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

    if(UConf%SwitchProgram == 5) then
        cFileName = trim(adjustl(UConf%cDirectory))//trim(adjustl("ResultU/"))//UConf%ResumeFileName
    else if(UConf%SwitchProgram /= 7) then
        write(6,*) "Please input VTK-File name of Unstructured Grid for RESUME"
        read(5,*) cFileName
    else
        write(cStep,*) 0
        cDirectory = trim(adjustl(UConf%cDirectory))//trim(adjustl("ResultU/"))
        cFileName = trim(adjustl(cDirectory))//trim(adjustl(UConf%cFileName))//trim(adjustl("_"))//trim(adjustl(cStep))//"th.vtk"
    end if

    open(unit = 1, file =trim(adjustl(cFileName)), status = 'old')
        write(6,*) "Find VTK-Data"
        do iLoop=1,5
            read(1,*) cAnnotate !Header
        end do

        do iPoint=1, UG%GI%Points
            read(1,*) cAnnotate,cAnnotate,cAnnotate !Point Coordinate
        end do

        read(1,*) cAnnotate,cAnnotate
        do iCell=1, UG%GI%RealCells
            read(1,*) cAnnotate,cAnnotate,cAnnotate,cAnnotate
        end do

        read(1,*) cAnnotate
        do iCell=1, UG%GI%RealCells
            read(1,*) cAnnotate
        end do

        read(1,*) cAnnotate
!密度プロット
        read(1,*) cAnnotate
        read(1,*) cAnnotate

        do iCell=1, UG%GI%RealCells
            read(1,*) UCC%PrimitiveVariable(1,iCell,1,1)
        end do

        read(1,*) cAnnotate
        do iCell=1, UG%GI%RealCells
            read(1,*) (UCC%PrimitiveVariable(iLoop,iCell,1,1),iLoop=2,4)
        end do

        read(1,*) cAnnotate
        read(1,*) cAnnotate
        do iCell=1, UG%GI%RealCells
            read(1,*) UCC%PrimitiveVariable(UG%GM%Dimension+2,iCell,1,1)
        end do

    close(1)

return
end subroutine UReadInitialCondition

