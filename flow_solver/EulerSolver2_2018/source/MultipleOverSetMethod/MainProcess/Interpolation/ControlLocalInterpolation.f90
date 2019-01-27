!***********************************/
!	Name:ローカル格子の内挿制御ルーチン
!	Alias:ControlLocalInterpolation
!	Description:格子が移動したとき全セルの初期値を内挿，格子が移動していないとき境界および上位格子セルとかぶったセルのみ内挿を実施
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.01.10
!	Update:
!	Update:
!	Other:
!***********************************/
    subroutine ControlLocalInterpolation(iGrid,LG,CC,Geom,iStep,MConf)
    use StructVar_Mod
    implicit none
    integer, intent(in) :: iGrid
    type(LocalGrids), intent(inout) :: LG
    type(CellCenter), intent(in) :: CC
    type(Geometry), intent(in) :: Geom
    integer, intent(in) :: iStep
    type(Configulation), intent(in) :: MConf

    type(LocalInterpolate) :: LI
    integer :: iTargetCell

    allocate(LI%Moment(5),LI%ArithmeticAverage(5))
    allocate(LI%UpwindMoment(5),LI%DownwindMoment(5),LI%TmpMoment(5))

    LG%CC(iGrid)%InterpolatedQuantity = LG%CC(iGrid)%ConservedQuantity !物理量をコピーする(自格子の前段階の値からも内挿する可能性があるため)

    !if(iStep == 1) then !初回のみ流束制限関数が未計算であるため用意しておく !ShockConverge~を使う時だけ必要
        !call UConserve2Primitive(LG%UG(iGrid),LG%CC(iGrid))
        !call UGetGradient(LG%UG(iGrid),LG%CC(iGrid))
        !call UGetLimiter(MConf,LG%UG(iGrid),LG%CC(iGrid),LG%CE(iGrid))
    !end if


    if(LG%MW(iGrid)%iMotionWait == 0) then !格子移動を表す何か
        !全セルについて補間を実行
        do iTargetCell=1, LG%UG(iGrid)%GI%AllCells
            call InterpolationOfLocalGrid(iTargetCell,iGrid,LG,CC,Geom,MConf,LI)
            !call ShockConvergeLocalInterpolateMain(iTargetCell,iGrid,LG,CC,Geom,MConf,LI)
        end do

        LG%MW(iGrid)%iEstimateWait = LG%MG(iGrid)%iStepFromLastMove
        LG%MG(iGrid)%iStepFromLastMove = 0

    else !格子が移動していないとき
        do iTargetCell=1, LG%UG(iGrid)%GI%AllCells
            if(LG%UG(iGrid)%GM%Interpolated(iTargetCell,1,1) == 1 .or. iTargetCell > LG%UG(iGrid)%GI%RealCells) then !Only Overlapped with SUPERIOR Cell
                !call ShockConvergeLocalInterpolateMain(iTargetCell,iGrid,LG,CC,Geom,MConf,LI)
                call InterpolationOfLocalGrid(iTargetCell,iGrid,LG,CC,Geom,MConf,LI)
            end if
        end do

    end if

    LG%CC(iGrid)%ConservedQuantity = LG%CC(iGrid)%InterpolatedQuantity !内挿が終わったデータを元の変数に格納する


    if(LG%MW(iGrid)%iMotionWait == 0) then
        do iTargetCell = 1, LG%UG(iGrid)%GI%RealCells
           if(LG%MG(iGrid)%IS(iTargetCell)%BoundaryAdjacent(1) /= 0) then
                if(LG%MG(iGrid)%IS(iTargetCell)%BoundaryAdjacent(3) == 0) then
                    LG%CC(iGrid)%ConservedQuantity(:,iTargetCell,1,1) = &
                        &   0.5d0*(LG%CC(iGrid)%ConservedQuantity(:,iTargetCell,1,1) + LG%CC4MB(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1))
                end if
           end if
        end do

    end if

    return
    end subroutine ControlLocalInterpolation
