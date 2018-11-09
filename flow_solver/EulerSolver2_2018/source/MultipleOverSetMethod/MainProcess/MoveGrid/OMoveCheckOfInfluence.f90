!***********************************/
!	Name:累積変位量による移動可能判定
!	Alias:OMoveCheckOfInfluence
!	Description:累積した重心の並進変位ノルムと，最大の格子幅を比較し，上回った場合モーションウェイトを0にする
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.02.05
!	Update:
!	Other:
!***********************************/
subroutine OMoveCheckOfInfluence(MG1,MG2,UG1,MW1,MW2)
    use StructVar_Mod
    implicit none

    type(MoveGrid), intent(in) :: MG1, MG2
    type(UnstructuredGrid), intent(in) :: UG1
    type(MotionWait), intent(inout) :: MW1, MW2
    double precision :: CenterDistance

    CenterDistance = sqrt(dot_product(MG1%RC%GravityCenter(2,:) - MG2%RC%GravityCenter(2,:),MG1%RC%GravityCenter(2,:) - MG2%RC%GravityCenter(2,:)))

    !if(MG1%MaxInfluenceRange + MG2%MaxInfluenceRange > CenterDistance) then
    if(UG1%RoughRadius + MG2%MaxInfluenceRange*1.2d0 > 0.95d0*CenterDistance) then
        MW1%iMotionWait = 0
        MW2%iMotionWait = 0
    end if

return
end subroutine OMoveCheckOfInfluence
