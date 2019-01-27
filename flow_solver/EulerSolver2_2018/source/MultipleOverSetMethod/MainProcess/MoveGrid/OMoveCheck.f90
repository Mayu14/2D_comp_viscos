!***********************************/
!	Name:累積変位量による移動可能判定
!	Alias:OMoveCheck
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
subroutine OMoveCheck(MW,MinGridWidth,MaxGridWidth)
    use StructVar_Mod
    implicit none

    type(MotionWait), intent(inout) :: MW
    double precision, intent(in) :: MinGridWidth !これと回転角(rad)の積で，重心付近のセルの回転移動量を評価する
    double precision, intent(in) :: MaxGridWidth

    !double precision :: TranslationAmount
    !TranslationAmount = sqrt(dot_product(MW%Accumulate%Translation,MW%Accumulate%Translation))
    !if(TranslationAmount > MaxGridWidth) then

    if(MW%iMotionWait /= 0) then
        if(sqrt(dot_product(MW%Accumulate%Translation,MW%Accumulate%Translation)) + MinGridWidth*MW%Accumulate%Rotation(3) &
                        &   > MaxGridWidth) then
            MW%iMotionWait = 0
        end if

    else !Motion Wait = 0 but still accumulate distance is not enough
        if(sqrt(dot_product(MW%Accumulate%Translation,MW%Accumulate%Translation)) + MinGridWidth*MW%Accumulate%Rotation(3) &
                        &   < MaxGridWidth) then
            MW%iMotionWait = 1
        end if
    end if


return
end subroutine OMoveCheck
