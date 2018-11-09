!***********************************/
!	Name:計算領域の設定を制御するためのプログラム
!	Alias:SetRegion
!	Description:
!	Type:Geom
!	Input:
!	Output:Geom
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine SetRegion(Conf,Geom)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none

    type(Configulation),intent(in) :: Conf
    type(Geometry),intent(inout) :: Geom

    if(Conf%UseReadRegion == 1) then
        call ReadRegionData(Geom)
    else if(Conf%UseReadRegion == 0) then
        call SetRegionData(Geom)
    end if

    allocate(Geom%Area(Geom%Dimension))
    allocate(Geom%Width(1,1,Geom%Dimension))
    allocate(Geom%Volume(1))

        iXmax = Geom%CellNumber(1)
        iYmax = max(1,Geom%CellNumber(2)) !セル総数と1の大きい方がiYmax
        iZmax = max(1,Geom%CellNumber(3))
    if(Conf%UseOverSet == 1) then
        !allocate(Geom%CellType(iXmax,iYmax,iZmax))
        allocate(Geom%Interpolated(iXmax,iYmax,iZmax))
    end if

!    allocate(Geom%Vector(Geom%Dimension,Geom%Dimension))

    call CalcWidth(Geom)
    call CalcVolume(Geom)
    call CalcEdgeArea(Geom)
    !call CalcVector(Geom)

return
end subroutine SetRegion
