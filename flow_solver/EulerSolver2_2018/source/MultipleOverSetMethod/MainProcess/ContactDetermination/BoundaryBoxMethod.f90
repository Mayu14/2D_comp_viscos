!***********************************/
!	Name:非構造セルの簡易接触判定関数
!	Alias:BoundaryBoxMethod
!	Description:名前の通り境界矩形法を用いて判別する
!	Type:integer
!	Input:CountUpVertex
!	Output:iTrueOrFalse
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.26
!	Update:2017.12.26
!	Other:
!***********************************/

function BoundaryBoxMethod(CUV) result(iToF)
use StructVar_Mod
implicit none

    type(CountUpVertex), intent(in) :: CUV

    integer :: iToF
    !double precision :: SquareArea, AreaSum

        !SquareArea = dot_product(CUV%TargetCoords(0,1:3) - CUV%DonorCoords(0,1:3), &
        !                        &   CUV%TargetCoords(0,1:3) - CUV%DonorCoords(0,1:3))

        !AreaSum = CUV%TargetVolume + CUV%DonorVolume

        !if(SquareArea > AreaSum) then !このときセル間距離は明らかに遠く，接触していない
        if(dot_product(CUV%TargetCoords(0,1:3) - CUV%DonorCoords(0,1:3), &
                                &   CUV%TargetCoords(0,1:3) - CUV%DonorCoords(0,1:3)) > &
                                &   CUV%TargetVolume + CUV%DonorVolume) then
            iToF = 0
        else
            iToF = 1
        end if
    return
end function BoundaryBoxMethod

