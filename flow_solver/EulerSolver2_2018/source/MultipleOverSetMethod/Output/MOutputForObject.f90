subroutine MOutput4Object(UG,NSRP,iStep,iTargetGrid)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(AnotherCoordinateRelation), intent(in) :: NSRP
    integer, intent(in) :: iStep
    integer, intent(in) :: iTargetGrid
    character(len=64) :: cDirectory,cFileName!, cCaseName
    character(len=32) :: cStep
    character(len=32) :: cGrid

    write(cStep,*) iStep
    write(cGrid,*) iTargetGrid

    cDirectory = "ResultObj"//trim(adjustl(cGrid))//"/" !UConf%SaveDirectiry
    cFileName = trim(adjustl(cDirectory))//"Object"//trim(adjustl(cStep))//"th_Result.vtk"
    !cCaseName = "RigidBody" !UConf%CaseName
    open(unit = 1, file =trim(adjustl(cFileName)), status = 'unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,"('Polygonal Data BODY')")
        write(1,"('ASCII')")
        write(1,"('DATASET POLYDATA')")
        write(1,"('POINTS ',(1x,i7),'double')") UG%IO%iTotal

        do iPoint=1, UG%IO%iTotal
            write(1,"(3(1x,f22.17))") NSRP%NextCoordsInG(UG%IO%PointNum(iPoint),1),NSRP%NextCoordsInG(UG%IO%PointNum(iPoint),2),&
                                    &   NSRP%NextCoordsInG(UG%IO%PointNum(iPoint),3)+2.0d0/iTargetGrid !for OverSet
        end do
        write(1,*) ""

        write(1,"('POLYGONS ',2(1x,i7))") UG%IO%iTotal-2,(UG%IO%iTotal-2)*4
        do iPoint=2, UG%IO%iTotal-1
                write(1,"(4(1x,i7))") 3,UG%IO%PointNum(1)-UG%IO%iMinNumber,&
                    &   UG%IO%PointNum(iPoint)-UG%IO%iMinNumber,UG%IO%PointNum(iPoint+1)-UG%IO%iMinNumber !Adjusted Point Number
        end do
        write(1,*) ""
    close(1)

return
end subroutine MOutput4Object

