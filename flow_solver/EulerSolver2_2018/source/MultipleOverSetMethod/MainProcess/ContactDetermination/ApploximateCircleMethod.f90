!***********************************/
!	Name:格子を重心を中心とした円として近似して簡易的な接触判定を実施する関数
!	Alias:BoundarySphereMethod
!	Description:入力された2つの格子の粗い半径の和と重心間距離の大小関係を比較し，0または1を返す
!	Type:integer
!	Input:
!	Output:
!	Note:格子が接触していないとき0を，詳細な判定が必要なとき1を返す
!	Author:Akitaka Toyota
!	Date:2017.12.19
!	Update:2017.12.19
!	Other:
!***********************************/
function BoundarySphereMethod(RoughRadius1,RoughRadius2,RC1,RC2,iTimeCombination) result(iToF)
    use StructVar_Mod
!    use LoopVar_Mod
    implicit none
    double precision, intent(in) :: RoughRadius1, RoughRadius2
    type(RelativeCoordinate),intent(in) :: RC1, RC2
    integer, intent(in):: iTimeCombination
    integer :: iToF
    double precision :: GravityCenterDistance

    if(iTimeCombination == 0) then
        GravityCenterDistance = sqrt(sum((RC1%GravityCenter(2,:) - RC2%GravityCenter(2,:))**2))
    else
        GravityCenterDistance = sqrt(sum((RC1%GravityCenter(2,:) - RC2%GravityCenter(1,:))**2))
    end if

    if(GravityCenterDistance > RoughRadius1 + RoughRadius2) then !重心間距離が半径の和より大きい = 明らかに接触していない
        iToF = 0
    else
        iToF = 1
    end if
return
end function BoundarySphereMethod
