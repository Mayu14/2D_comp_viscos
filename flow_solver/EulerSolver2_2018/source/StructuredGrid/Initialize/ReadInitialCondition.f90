!***********************************/
!	Name:データから初期条件を読み出して適用するためのプログラム
!	Alias:ReadInitialCondition
!	Description:現状1次元のGnuPlot用出力からの再開にのみ対応
!	Type:Conf,Geom,CC
!	Input:Configulation,Geometory,CellCenter
!	Output:CellCenter
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine ReadInitialCondition(Conf,Geom,CC)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod
    implicit none
    type(Configulation), intent(in) :: Conf
    type(Geometry),intent(in) :: Geom
    type(CellCenter),intent(inout) :: CC
    double precision :: Temporary
    character(len=64) :: cDensity,cVelocity,cPressure,cFileName,cAnnotate
    integer,allocatable :: CellNumber(:)
    allocate(CellNumber(2))

    if(Conf%ResumeFileFormat == 1) then
        write(6,*) "Enter VTK-File Name of Structured Grid's Flow Data"
        read(5,*) cFileName

        open(unit=1,file=cFileName,status='old')
            write(6,*) "Find VTK"
            do iLoop=1,5+(1+Geom%CellNumber(1))+(1+Geom%CellNumber(2))+(1+Geom%CellNumber(3))+3 !読み飛ばす
                read(1,*) cAnnotate
            end do

            do iCenterZ=1,CellNumber(2)
                do iCenterY=1,CellNumber(1)
                    do iCenterX=1,Geom%CellNumber(1)
                        read(1,*) CC%PrimitiveVariable(1,iCenterX,iCenterY,iCenterZ)
                    end do
                end do
            end do

            read(1,*) cAnnotate

            if(Geom%Dimension == 2) then
                do iCenterZ=1, CellNumber(2)
                    do iCenterY=1, CellNumber(1)
                        do iCenterX=1, Geom%CellNumber(1)
                            read(1,*) &
                                &   CC%PrimitiveVariable(2,iCenterX,iCenterY,iCenterZ), &
                                &   CC%PrimitiveVariable(3,iCenterX,iCenterY,iCenterZ), cAnnotate
                        end do
                    end do
                end do
            else
                write(6,*) "1D or 3D can not used"
                stop
            end if

            do iLoop=1,2
                read(1,*) cAnnotate
            end do

            do iCenterZ=1,CellNumber(2)
                do iCenterY=1,CellNumber(1)
                    do iCenterX=1,Geom%CellNumber(1)
                        read(1,*) CC%PrimitiveVariable(Geom%Dimension+2,iCenterX,iCenterY,iCenterZ)
                    end do
                end do
            end do
        close(1)

        call Primitive2Conserve(Geom,CC)


    else if(Conf%ResumeFileFormat == 2) then
        write(cDensity, '("Density/Density",i6.6,".txt")') Conf%ResumeNumber/OutputInterval
        open(unit=1,file=cDensity,status='unknown')
            do iCenterX=1, Geom%CellNumber(1)
                read(1,*) Temporary,CC%PrimitiveVariable(1,iCenterX,1,1)
            end do
        close(1)

    write(cVelocity, '("Velocity/Velocity",i6.6,".txt")') Conf%ResumeNumber/OutputInterval
    open(unit=1,file=cVelocity,status='unknown')
        do iCenterX=1, Geom%CellNumber(1)
            read(1,*) Temporary,CC%PrimitiveVariable(2,iCenterX,1,1)
        end do
    close(1)

    write(cPressure, '("Pressure/Pressure",i6.6,".txt")') Conf%ResumeNumber/OutputInterval
    open(unit=1,file=cPressure,status='unknown')
        do iCenterX=1, Geom%CellNumber(1)
            read(1,*) Temporary,CC%PrimitiveVariable(3,iCenterX,1,1)
        end do
    close(1)

    else if(Conf%ResumeFileFormat == 3) then
        if(Geom%Dimension == 1) then
        open(unit=2, file=Conf%ResumeFileName,status='old')
        !一括で読むためのデータを出力すべき

            do iCenterX=1, Geom%CellNumber(1)
                read(2,*) Temporary,CC%PrimitiveVariable(1,iCenterX,1,1)
            end do

            do iCenterX=1, Geom%CellNumber(1)
                read(2,*) Temporary, CC%PrimitiveVariable(2,iCenterX,1,1)
            end do

            do iCenterX=1, Geom%CellNumber(1)
                read(2,*) Temporary, CC%PrimitiveVariable(3,iCenterX,1,1)
            end do
        close(2)
        end if

    end if

return
end subroutine ReadInitialCondition
