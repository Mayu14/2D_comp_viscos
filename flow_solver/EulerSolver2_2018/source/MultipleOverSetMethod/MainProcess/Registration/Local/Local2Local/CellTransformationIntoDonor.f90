!***********************************/
!	Name:内挿targetセルをドナー候補の座標系に変換する
!	Alias:CellTransformationIntoDonor
!	Description:iTime=0でドナー格子のN+1座標系へ写す，iTime=1でN座標系へ
!	Type:
!	Input:
!	Output:
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.12.05
!	Update:
!	Other:
!***********************************/
    subroutine CellTransformationIntoDonor(iTime,CUV,Global2Donor)
    use StructVar_Mod
    implicit none
    type(CountUpVertex), intent(inout) :: CUV
    type(TransformationMatrix), intent(in) :: Global2Donor
    integer, intent(in) :: iTime

    if(iTime == 0) then
        CUV%TargetCoords(0,:) = matmul(Global2Donor%Top2Next,CUV%TargetCoords(0,:))
        CUV%TargetCoords(1,:) = matmul(Global2Donor%Top2Next,CUV%TargetCoords(1,:))
        CUV%TargetCoords(2,:) = matmul(Global2Donor%Top2Next,CUV%TargetCoords(2,:))
        CUV%TargetCoords(3,:) = matmul(Global2Donor%Top2Next,CUV%TargetCoords(3,:))
    else
        CUV%TargetCoords(0,:) = matmul(Global2Donor%Top2Current,CUV%TargetCoords(0,:))
        CUV%TargetCoords(1,:) = matmul(Global2Donor%Top2Current,CUV%TargetCoords(1,:))
        CUV%TargetCoords(2,:) = matmul(Global2Donor%Top2Current,CUV%TargetCoords(2,:))
        CUV%TargetCoords(3,:) = matmul(Global2Donor%Top2Current,CUV%TargetCoords(3,:))
    end if

    return
    end subroutine CellTransformationIntoDonor
