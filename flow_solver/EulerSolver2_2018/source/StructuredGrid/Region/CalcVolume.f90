!***********************************/
!	Name:計算格子のセル体積を計算するプログラム
!	Alias:CalcEdgeArea
!	Description:1次元の場合格子の幅を用いる
!	Type:Geom
!	Input:Geom%Width
!	Output:Geom＆%Volume
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine CalcVolume(Geom)
use StructVar_Mod
use LoopVar_Mod
implicit none
    type(Geometry),intent(inout) :: Geom

    Geom%Volume(1) = 1.0d0

    do iLoop=1, Geom%Dimension
        Geom%Volume(1) = Geom%Volume(1) * 2.0d0*Geom%Width(1,1,iLoop)
    end do

return
end subroutine CalcVolume

