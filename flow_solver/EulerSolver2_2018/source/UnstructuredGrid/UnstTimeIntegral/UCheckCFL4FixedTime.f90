!***********************************/
!	Name:時間積分制御用プログラム
!	Alias:TimeIntegral
!	Description:Configulationによって時間積分の精度を変更できる
!	Type:Conf,Geom,CellEdge,CellCenter,iStep(RRK2用)
!	Input:Geom,CC,CE,CC,iStep
!	Output:CC%ConservedQuantity
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.06
!	Other:局所時間刻み法はまだ未対応
!***********************************/
subroutine UCheckCFL4FixedTime(UG,UCC,iSplit)
use StructVar_Mod
use ConstantVar_Mod
implicit none

    type(UnstructuredGrid),intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    integer, intent(inout) :: iSplit

    call UCalcTimeWidth(UG,UCC)
    iSplit = ceiling(DefaultTimeStep / minval(UCC%TimeWidth(:,1,1),1))
    !write(6,*) iSplit
return
end subroutine UCheckCFL4FixedTime
