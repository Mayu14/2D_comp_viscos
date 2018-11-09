!***********************************/
!	Name:計算領域の設定をすべて表示するプログラム
!	Alias:TestRegionData
!	Description:おかしいところは目視で確認するに限る
!	Type:Geom
!	Input:
!	Output:Geom
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.10
!	Other:
!***********************************/
subroutine TestRegionData(Geom)
use StructVar_Mod
use LoopVar_Mod
implicit none
    type(Geometry),intent(in) :: Geom

        write(6,*) "Dimension",Geom%Dimension
        write(6,*) "Area",Geom%Area
        write(6,*) "Volume",Geom%Volume
        write(6,*) "Bound",Geom%Bound
        write(6,*) "CellNumber",Geom%CellNumber
        write(6,*) "Width",Geom%Width
!        write(6,*) "Vector",Geom%Vector

return
end subroutine TestRegionData
