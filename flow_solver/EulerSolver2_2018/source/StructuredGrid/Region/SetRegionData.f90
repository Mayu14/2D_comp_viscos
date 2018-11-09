!***********************************/
!	Name:計算領域の設定をユーザー入力によって行うプログラム&StrGridの出力
!	Alias:SetRegion
!	Description:
!	Type:Geom
!	Input:
!	Output:Geom
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.10
!	Other:
!***********************************/
subroutine SetRegionData(Geom)
use StructVar_Mod
use LoopVar_Mod
implicit none
character :: cAsk
    type(Geometry),intent(inout) :: Geom

        write(6,*) "Please input Dimension"
        read(5,*) Geom%Dimension

    allocate(Geom%CellNumber(3))
    allocate(Geom%Bound(2,Geom%Dimension))
    Geom%CellNumber=0.0d0
    Geom%Bound=0.0d0

    do iLoop = 1, Geom%Dimension
        write(6,*) "Please input Cell Number in Direction",iLoop
        read(5,*) Geom%CellNumber(iLoop)
    end do

    do iLoop = 1, Geom%Dimension
        write(6,*) "Please input Region's Lower/Upper Boundary in Direction",iLoop
        read(5,*) Geom%Bound(1,iLoop),Geom%Bound(2,iLoop)
    end do

    write(6,*) "Do you want to save this setting to ""StrGrid""?[y/n]"
    read(5,*) cAsk

    if(cAsk == "y") then
        open(unit=1,file="StrGrid",status='unknown')
            write(1,*) "次元",Geom%Dimension
            write(1,*) "X分割数",Geom%CellNumber(1)
            write(1,*) "Y分割数",Geom%CellNumber(2)
            write(1,*) "Z分割数",Geom%CellNumber(3)
            write(1,*) "X上下界",Geom%Bound(1,1),Geom%Bound(2,1)
            write(1,*) "Y上下界",Geom%Bound(1,2),Geom%Bound(2,2)
            write(1,*) "Z上下界",Geom%Bound(1,3),Geom%Bound(1,3)
        close(1)
        write(6,*) "Complete to generate StrGrid."
    end if

return
end subroutine SetRegionData
