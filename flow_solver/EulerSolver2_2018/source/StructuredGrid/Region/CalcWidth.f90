!***********************************/
!	Name:計算格子の格子幅を求めるプログラム
!	Alias:CalcWidth
!	Description:分割数と領域の下限から格子を計算する
!	Type:Geom
!	Input:Geom%Bound,Geom%CellNumber
!	Output:Geom%Width
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine CalcWidth(Geom)
use StructVar_Mod
use LoopVar_Mod
implicit none
    type(Geometry),intent(inout) :: Geom

    do iLoop=1, Geom%Dimension
        Geom%Width(1,1,iLoop) = 0.5d0 * dabs(Geom%Bound(2,iLoop) - Geom%Bound(1,iLoop)) / dble(Geom%CellNumber(iLoop)) !セル中心から界面まで(格子幅の半分)
    end do

return
end subroutine CalcWidth
