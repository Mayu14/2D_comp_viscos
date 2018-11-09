!***********************************/
!	Name:計算格子の界面面積を計算するプログラム
!	Alias:CalcEdgeArea
!	Description:1次元の場合1.0d0を代入する(単位面積だと考える)
!	Type:Geom
!	Input:Geom%Width
!	Output:Geom＆%Area
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine CalcEdgeArea(Geom)
    use StructVar_Mod
    implicit none

    type(Geometry),intent(inout) :: Geom

    if(Geom%Dimension == 1) then
        Geom%Area(1) = 1.0d0 !1次元系の面積(全段面で一様であるから正直なんでもいい)
    end if

    if(Geom%Dimension == 2) then
        Geom%Area(1) = 2.0d0*Geom%Width(1,1,2) !ΔLx = Δy !2次元系の面積=辺の長さ
        Geom%Area(2) = 2.0d0*Geom%Width(1,1,1) !ΔLy = Δx
    end if

    if(Geom%Dimension == 3) then
        Geom%Area(1) = 2.0d0*Geom%Width(1,1,2) * 2.0d0*Geom%Width(1,1,3)  !ΔSx = Δy*Δz
        Geom%Area(2) = 2.0d0*Geom%Width(1,1,3) * 2.0d0*Geom%Width(1,1,1)  !ΔSy = Δz*Δx
        Geom%Area(3) = 2.0d0*Geom%Width(1,1,1) * 2.0d0*Geom%Width(1,1,2)  !ΔSz = Δx*Δy
    end if

return
end subroutine CalcEdgeArea
