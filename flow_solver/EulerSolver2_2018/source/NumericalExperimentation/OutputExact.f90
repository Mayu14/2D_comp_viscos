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
subroutine OutputExact(iStep,NE,Geom) !(iStep,iTargetGrid,CC,UG,MG,Geom)
    use StructVar_Mod
    use LoopVar_Mod
    use ExactSol
    implicit none

    integer,intent(in) :: iStep
    type(NumericalExperiment), intent(inout) :: NE
    type(Geometry),intent(in) :: Geom

    character(len=64) :: cStep,cFileName
    integer, allocatable :: CellNumber(:)

    allocate(CellNumber(2))

    write(cStep,*) iStep

    cFileName = "EXP/ExactSolution"//trim(adjustl(trim(cStep)))//"th_Result.vtk"

    CellNumber(:) = max(1,Geom%CellNumber(2:3))

    open(unit=1,file=cFileName,status='unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,"('Struct_ShockTube')")
        write(1,"('ASCII')")
        write(1,"('DATASET RECTILINEAR_GRID')")
        write(1,"('DIMENSIONS ',3(1x,i5))") Geom%CellNumber(1)+1,CellNumber(1)+1,CellNumber(2)

!座標情報
        write(1,"('X_COORDINATES ',i5,' float')") Geom%CellNumber(1)+1
            do iCenterX=1, Geom%CellNumber(1)+1
                write(1,"(f22.17)") Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * (dble(iCenterX)-1.0d0)
            end do

        write(1,"('Y_COORDINATES ',i5,' float')") CellNumber(1)+1
        do iCenterY=1, CellNumber(1)+1
            if(Geom%Dimension >= 2) then
                write(1,"(f22.17)") Geom%Bound(1,2) + 2.0d0 * Geom%Width(1,1,2) * (dble(iCenterY)-1.0d0)
            else
                write(1,"(f22.17)") 0.0
            end if
        end do

        write(1,"('Z_COORDINATES ',i5,' float')") CellNumber(2)
        do iCenterZ=1, CellNumber(2)
            if(Geom%Dimension == 3) then
                write(1,"(f22.17)") Geom%Bound(1,3) + 2.0d0 * Geom%Width(1,1,3) * (dble(iCenterZ)-0.5d0)
            else
                write(1,"(f22.17)") 0.0
            end if
        end do

        if(Geom%Dimension == 2) then
            write(1,"('CELL_DATA ',i7)") (Geom%CellNumber(1)) * (CellNumber(1))
        else if(Geom%Dimension == 3) then
            write(1,"('CELL_DATA ',i7)") (Geom%CellNumber(1)) * (CellNumber(1)) * (CellNumber(2))
        end if


!密度プロット
        write(1,"('SCALARS Density float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCenterZ=1,CellNumber(2)
            do iCenterY=1,CellNumber(1)
                do iCenterX=1,Geom%CellNumber(1)
                    write(1,"((2x,f22.14))") NE%STT%PrimitiveVariable(1,iCenterX)
                end do
            end do
        end do
!速度プロット
        write(1,"('VECTORS Velocity float')")
        if(Geom%Dimension == 1) then
            do iCenterX=1, Geom%CellNumber(1)
                write(1,"((2x,f22.17))") NE%STT%PrimitiveVariable(2,iCenterX),0.0d0,0.0d0
            end do

        else if(Geom%Dimension == 2) then
            do iCenterY=1, Geom%CellNumber(2)
                do iCenterX=1, Geom%CellNumber(1)
                    write(1,*) NE%STT%PrimitiveVariable(2,iCenterX), &
                        &   NE%STT%PrimitiveVariable(3,iCenterX),0.0d0
                end do
            end do

        else
            write(6,*) "3D cannot use"
            stop
        end if

!圧力プロット
        write(1,"('SCALARS Pressure float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCenterZ=1,CellNumber(2)
            do iCenterY=1,CellNumber(1)
                do iCenterX=1,Geom%CellNumber(1)
                    write(1,"((2x,f22.14))") NE%STT%PrimitiveVariable(4,iCenterX)
                end do
            end do
        end do

    close(1)

    return
end subroutine OutputExact
