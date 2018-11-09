!***********************************/
!	Name:Gnuplot用の出力結果を1つにまとめてGif作成用のファイルを作るプログラム
!	Alias:MergeGnuPlot
!	Description:
!	Type:Geom,iStep
!	Input:Density/Density*.txt,Velocity/Velocity*.txt,Pressure/Pressure*.txt
!	Output:
!	Note:これ単体のQMFというプログラムが存在する(Gnuplotスクリプトを吐き出すのでそっちの方が有能)
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine MergeGnuPlot(Geom,iStep)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none

    type(Geometry),intent(in) :: Geom
    integer,intent(in) :: iStep

    integer :: iProgress
    double precision :: Tmp
    double precision, allocatable :: TmpQuantity(:)
    character(len=64), target :: cStep,cDensity,cVelocity,cPressure
!    call Conserve2Primitive(Geom,CC)

    allocate(TmpQuantity(Geom%CellNumber(1)))
    write(cStep,*) iStep
    iProgress = int(dble(iStep)/10.0d0)

    write(6,*) "Density.txt Generating ..."
    open(unit=2,file="Density.txt",status='unknown')
        do iLoop = 0, iStep
            write(cDensity, '("Density/Density",i6.6,".txt")') iLoop

            open(unit=1,file=cDensity,status='unknown')
                do iCenterX=1, Geom%CellNumber(1)
                    read(1,*) Tmp,TmpQuantity(iCenterX)
                end do
            close(1)

            do iCenterX=1, Geom%CellNumber(1)
                write(2,*) Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * (dble(iCenterX)-0.5d0),TmpQuantity(iCenterX)
            end do
            write(2,*) ""
            write(2,*) ""

            if(mod(iLoop,iProgress) == 0) then
                write(6,"((1x,i3),'% Complete')") int(dble(iLoop)/dble(iProgress)*10.0d0)
            end if
        end do
    close(2)

    write(6,*) "Velocity.txt Generating..."
    open(unit=2,file='Velocity.txt',status='unknown')
        do iLoop = 0, iStep
            write(cVelocity, '("Velocity/Velocity",i6.6,".txt")') iLoop

            open(unit=1,file=cVelocity,status='unknown')
                do iCenterX=1, Geom%CellNumber(1)
                    read(1,*) Tmp,TmpQuantity(iCenterX)
                end do
            close(1)

            do iCenterX=1, Geom%CellNumber(1)
                write(2,*) Geom%Bound(1,1) + 2.0d0*Geom%Width(1,1,1)*(dble(iCenterX)-0.5d0), TmpQuantity(iCenterX)
            end do
            write(2,*) ""
            write(2,*) ""
            if(mod(iLoop,iProgress) == 0) then
                write(6,"((1x,i3),'% Complete')") int(dble(iLoop)/dble(iProgress)*10.0d0)
            end if
        end do
    close(2)

    write(6,*) "Pressure.txt Generating...."
    open(unit=2,file='Pressure.txt',status='unknown')
        do iLoop = 0, iStep
            write(cPressure, '("Pressure/Pressure",i6.6,".txt")') iLoop

            open(unit=1,file=cPressure,status='unknown')
                do iCenterX=1, Geom%CellNumber(1)
                    read(1,*) Tmp,TmpQuantity(iCenterX)
                end do
            close(1)

            do iCenterX=1, Geom%CellNumber(1)
                write(2,*) Geom%Bound(1,1) + 2.0d0*Geom%Width(1,1,1)*(dble(iCenterX)-0.5d0), TmpQuantity(iCenterX)
            end do

            write(2,*) ""
            write(2,*) ""
            if(mod(iLoop,iProgress) == 0) then
                write(6,"((1x,i3),'% Complete')") int(dble(iLoop)/dble(iProgress)*10.0d0)
            end if
        end do
    close(2)

return
end subroutine MergeGnuPlot
