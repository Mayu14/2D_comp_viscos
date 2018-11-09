!***********************************/
!	Name:ユーザー指定の時間刻みがCFL条件に抵触する場合，自動的に時間刻みを変更するプログラム
!	Alias:CheckCFL4FixedTime
!	Description:変更された時間刻みは，元の時間刻みの整数分の1になるようにして，なるべく統一時間刻みに戻す
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.02.01
!	Update:
!	Other:
!***********************************/
subroutine CheckCFL4FixedTime(Geom,CC,iSplit)
use StructVar_Mod
use ConstantVar_Mod
implicit none

    type(Geometry),intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    integer, intent(inout) :: iSplit

    call CalcTimeWidth(Geom,CC)
    iSplit = ceiling(DefaultTimeStep / minval(CC%TimeWidth(:,1,1),1)) !ceiling:天井= 切り上げ

return
end subroutine CheckCFL4FixedTime
