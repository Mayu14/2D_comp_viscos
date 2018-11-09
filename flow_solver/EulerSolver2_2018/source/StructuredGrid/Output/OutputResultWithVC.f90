!***********************************/
!	Name:仮想格子まで含んだ計算結果をvtkで吐き出すプログラム
!	Alias:OutputResultWithVC
!	Description:境界条件の代入状況を視覚的に確認できる
!	Type:Geom,iStep
!	Input:Geom,iStep
!	Output:vtk
!	Note:2次元以上にしか対応していない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:VC = VirtualCell
!***********************************/
subroutine OutputResultWithVC(Geom,iStep,CC)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none

    type(Geometry),intent(in) :: Geom
    integer,intent(in) :: iStep
    type(CellCenter),intent(inout) :: CC

    character(len=64) :: cStep,cFileName
    integer,allocatable :: CellNumber(:)

    allocate(CellNumber(2))

    call Conserve2Primitive(Geom,CC)

    write(cStep,*) iStep

    cFileName = "Result/VCPrimitiveVariable"//trim(adjustl(trim(cStep)))//"th_Result.vtk"

    CellNumber(:) = max(1,Geom%CellNumber(2:3))

    open(unit=1,file=cFileName,status='unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,"('Struct_ShockTube')")
        write(1,"('ASCII')")
        write(1,"('DATASET RECTILINEAR_GRID')")
        write(1,"('DIMENSIONS ',3(1x,i5))") Geom%CellNumber(1)+2,CellNumber(1)+2,CellNumber(2)

!座標情報
        write(1,"('X_COORDINATES ',i5,' float')") Geom%CellNumber(1)+2
            do iCenterX=0, Geom%CellNumber(1)+1
                write(1,"(f22.17)") Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * (dble(iCenterX)-0.5d0)
            end do

        write(1,"('Y_COORDINATES ',i5,' float')") CellNumber(1)+2
        do iCenterY=0, CellNumber(1)+1
            if(Geom%Dimension >= 2) then
                write(1,"(f22.17)") Geom%Bound(1,2) + 2.0d0 * Geom%Width(1,1,2) * (dble(iCenterY)-0.5d0)
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

        write(1,"('POINT_DATA ',i7)") (Geom%CellNumber(1)+2) * (CellNumber(1)+2) * CellNumber(2)

!密度プロット
        write(1,"('SCALARS Density float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCenterZ=1,CellNumber(2)
            do iCenterY=0,CellNumber(1)+1
                do iCenterX=0,Geom%CellNumber(1)+1
                    write(1,"((2x,f22.14))") CC%PrimitiveVariable(1,iCenterX,iCenterY,iCenterZ)
                end do
            end do
        end do
!速度プロット
        write(1,"('VECTORS Velocity float')")
        if(Geom%Dimension == 1) then
            do iCenterX=0, Geom%CellNumber(1)+1
                write(1,"((2x,f22.17))") CC%PrimitiveVariable(2,iCenterX,1,1),0.0d0,0.0d0
            end do

        else if(Geom%Dimension == 2) then
            do iCenterY=0, Geom%CellNumber(2)+1
                do iCenterX=0, Geom%CellNumber(1)+1
                    write(1,*) CC%PrimitiveVariable(2,iCenterX,iCenterY,1), &
                        &   CC%PrimitiveVariable(3,iCenterX,iCenterY,1),0.0d0
                end do
            end do

        else
            do iCenterZ=1, Geom%CellNumber(3)
                do iCenterY=1, Geom%CellNumber(2)
                    do iCenterX=1, Geom%CellNumber(1)
                        write(1,"((2x,f22.17))") &
                            &   CC%PrimitiveVariable(2,iCenterX,iCenterY,iCenterZ), &
                            &   CC%PrimitiveVariable(3,iCenterX,iCenterY,iCenterZ), &
                            &   CC%PrimitiveVariable(4,iCenterX,iCenterY,iCenterZ)
                    end do
                end do
            end do
        end if

!圧力プロット
        write(1,"('SCALARS Pressure float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCenterZ=1,CellNumber(2)
            do iCenterY=0,CellNumber(1)+1
                do iCenterX=0,Geom%CellNumber(1)+1
                    write(1,"((2x,f22.14))") CC%PrimitiveVariable(Geom%Dimension+2,iCenterX,iCenterY,iCenterZ)
                end do
            end do
        end do
    close(1)
return
end subroutine OutputResultWithVC
