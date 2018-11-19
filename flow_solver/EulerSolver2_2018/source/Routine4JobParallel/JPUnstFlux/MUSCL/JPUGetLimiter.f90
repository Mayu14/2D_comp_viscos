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
subroutine JPUGetLimiter(UConf,UG,UCC,UCE)
    use StructVar_Mod
    use LoopVar_Mod
    use UVenkatakrishnanVar_Mod4OMP
    implicit none
    type(Configulation),intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(CellEdge), intent(inout) :: UCE
    type(PrivateVarOfUVenkatakrishnan) :: PVOUV
!隣接セルにおける最大変分と最小変分を検索(1次元の場合は左右の変動量の多い方を見つけるだけ)
!call FindVariation!(Geom,CC)

    if(UConf%KindLimiter == 1) then !Venkatakrishnan Limiter
        allocate(PVOUV%UFV%Variation(UG%GM%Dimension+2, 3))
        call JPUFindVariation(UG,UCC,PVOUV%UFV)
        call JPUGetNormalGradient(UG,UCC,UCE,PVOUV%UGNG)
        call JPUVenkatakrishnan(UG,UCC,UCE,PVOUV%UV)
    end if

return
end subroutine JPUGetLimiter
