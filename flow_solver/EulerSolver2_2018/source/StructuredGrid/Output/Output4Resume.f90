!***********************************/
!	Name:途中から再計算を行うためのデータを吐き出すプログラム
!	Alias:Output4Resume
!	Description:プログラムと同じディレクトリにResume_***th_Resultという名前のファイルが生成される
!	Type:Geom,iStep
!	Input:Geom,iStep
!	Output:
!	Note:1次元にしか対応していない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine Output4Resume(Geom,iStep)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none

    type(Geometry),intent(in) :: Geom
    integer,intent(in) :: iStep

    integer :: iProgress
    double precision :: Tmp
    double precision, allocatable :: TmpQuantity(:)
    character(len=64), target :: cStep,cDensity,cVelocity,cPressure,cFileName
!    call Conserve2Primitive(Geom,CC)

    allocate(TmpQuantity(Geom%CellNumber(1)))
    write(cStep,*) iStep
    iProgress = int(dble(iStep)/10.0d0)

    cFileName = "Resume_"//trim(adjustl(cStep))//"th_Result"

    open(unit=2,file=trim(adjustl(cFileName)),status='unknown')
            write(cDensity, '("Density/Density",i6.6,".txt")') iStep
            open(unit=1,file=cDensity,status='unknown')
                do iCenterX=1, Geom%CellNumber(1)
                    read(1,*) Tmp,TmpQuantity(iCenterX)
                end do
            close(1)

            do iCenterX=1, Geom%CellNumber(1)
                write(2,*) Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * dble(iCenterX),TmpQuantity(iCenterX)
            end do
            write(2,*) ""

            write(cVelocity, '("Velocity/Velocity",i6.6,".txt")') iStep
            open(unit=1,file=cVelocity,status='unknown')
                do iCenterX=1, Geom%CellNumber(1)
                    read(1,*) Tmp,TmpQuantity(iCenterX)
                end do
            close(1)

            do iCenterX=1, Geom%CellNumber(1)
                write(2,*) Geom%Bound(1,1) + 2.0d0*Geom%Width(1,1,1)*dble(iCenterX), TmpQuantity(iCenterX)
            end do
            write(2,*) ""


            write(cPressure, '("Pressure/Pressure",i6.6,".txt")') iStep
            open(unit=1,file=cPressure,status='unknown')
                do iCenterX=1, Geom%CellNumber(1)
                    read(1,*) Tmp,TmpQuantity(iCenterX)
                end do
            close(1)

            do iCenterX=1, Geom%CellNumber(1)
                write(2,*) Geom%Bound(1,1) + 2.0d0*Geom%Width(1,1,1)*dble(iCenterX), TmpQuantity(iCenterX)
            end do
    close(2)
return
end subroutine Output4Resume
