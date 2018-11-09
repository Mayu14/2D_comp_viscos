!***********************************/
!	Name:Gnuplot用のデータを吐き出すプログラム
!	Alias:OutputGnuPlot
!	Description:密度はDensity/以下，速度，圧力もVelocity/,Pressure/以下に吐き出されるためプログラムのあるフォルダが汚れない．
!	Type:Geom,iStep
!	Input:Geom,iStep
!	Output:
!	Note:1次元にしか対応していない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:データを一括で消去するclean.shというバッチファイルと一緒に使おう
!***********************************/
subroutine OutputGnuPlot(Geom,iStep,CC)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none

    type(Geometry),intent(in) :: Geom
    integer,intent(in) :: iStep
    type(CellCenter),intent(inout) :: CC

    character(len=64) :: cStep,cDensity,cVelocity,cPressure
    integer,allocatable :: CellNumber(:)

    allocate(CellNumber(2))

    call Conserve2Primitive(Geom,CC)

    write(cStep,*) iStep
    write(cDensity, '("Density/Density",i6.6,".txt")') iStep

    open(unit=1,file=cDensity,status='unknown')
        do iCenterX=1, Geom%CellNumber(1)
            write(1,*) Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * (dble(iCenterX)-0.5d0),CC%PrimitiveVariable(1,iCenterX,1,1)
        end do
    close(1)

    write(cVelocity, '("Velocity/Velocity",i6.6,".txt")') iStep

    open(unit=1,file=cVelocity,status='unknown')
        do iCenterX=1, Geom%CellNumber(1)
            write(1,*) Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * (dble(iCenterX)-0.5d0),CC%PrimitiveVariable(2,iCenterX,1,1)
        end do
    close(1)

    write(cPressure, '("Pressure/Pressure",i6.6,".txt")') iStep

    open(unit=1,file=cPressure,status='unknown')
        do iCenterX=1, Geom%CellNumber(1)
            write(1,*) Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * (dble(iCenterX)-0.5d0),CC%PrimitiveVariable(3,iCenterX,1,1)
        end do
    close(1)

return
end subroutine OutputGnuPlot
