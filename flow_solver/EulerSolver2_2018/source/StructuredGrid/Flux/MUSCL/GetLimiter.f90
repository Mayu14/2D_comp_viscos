!***********************************/
!	Name:流束制限関数を求めるプログラム
!	Alias:GetLimiter
!	Description:CalcConfigのKindLimiterによって計算する流束制限関数を変更できるようにする予定
!	Type:CellCenter
!	Input:Configulation,Geometory,CellCenter,CellEdge
!	Output:CC%LimiterFunction
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:いまのところVenkatakrishnanLimiterのみに対応
!***********************************/
subroutine GetLimiter(Conf,Geom,CC,CE)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Configulation),intent(in) :: Conf
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CellEdge), intent(inout) :: CE

!隣接セルにおける最大変分と最小変分を検索(1次元の場合は左右の変動量の多い方を見つけるだけ)
!call FindVariation!(Geom,CC)

if(Conf%KindLimiter == 1) then !Venkatakrishnan Limiter
    nullify(iCX,iCY,iCZ,iFX,iFY,iFZ,iVar)
    iVar => iVariable
    iCX => iCenterX
    iFX => iFaceX
    iCY => iCenterY
    iFY => iFaceY
    iCZ => iCenterZ
    iFZ => iFaceZ

    call FindVariation(Geom,CC)

    call GetNormalGradient(Geom,CC,CE)

    call Venkatakrishnan(Geom,CC,CE)

end if

return
end subroutine GetLimiter
