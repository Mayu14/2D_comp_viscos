!***********************************/
!	Name:移動先の要素の初期値を決定するための補間処理(流束制限関数を用いて衝撃波面の維持を試みたケース)
!	Alias:ShockConvergeLocalInterpolateMain
!	Description:
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2018.01.10
!	Update:2018.01.10
!	Other:
!***********************************/
    subroutine ShockConvergeLocalInterpolateMain(iTargetCell,iGrid,LG,CC,Geom,MConf,LI)
    use StructVar_Mod
    use FrequentOperation
    implicit none
    integer, intent(in) :: iTargetCell
    integer, intent(in) :: iGrid
    type(LocalGrids), intent(inout) :: LG
    type(CellCenter), intent(in) :: CC
    type(Geometry), intent(in) :: Geom
    type(Configulation), intent(in) :: MConf
    type(DonorElement), pointer :: DE

    type(LocalInterpolate), intent(inout) :: LI
        LI%Threshold = 0.000d0
        LI%AreaSum = 0.0d0
        LI%Moment = 0.0d0
        LI%MinLimiterValue = 1.0d0
        LI%ArithmeticAverage = 0.0d0
        LI%SummationCellNumber = 0
        DE => LG%MG(iGrid)%ADE(iTargetCell)%DE !rootへ移動

        do while(associated(DE)) !次のレコードがなくなるまで続行

            !if(LG%MW(iGrid)%iMotionWait /= 0 .and. DE%Grid == iGrid) then !のとき内挿対象とみなさない
                !何もせず次のレコードへ
            !else
                LI%AreaSum = LI%AreaSum + DE%SharedArea !面積の総和計算

                if(DE%Grid > 0) then !0番はGlobalGridであり，ConservedQuantityの配列構造が異なる，詳細はRegisterGlobal2Local参照のこと

                    LI%MinLimiterValue = min(LI%MinLimiterValue,minval(LG%CC(iGrid)%LimiterFunction(:,DE%Cell,1,1)))

                    LI%Moment = LI%Moment + LG%CC(DE%Grid)%ConservedQuantity(:,DE%Cell,1,1) * DE%SharedArea !面積×保存量の総和計算
                    LI%ArithmeticAverage = LI%ArithmeticAverage + LG%CC(DE%Grid)%ConservedQuantity(:,DE%Cell,1,1)


                else
                    call Global2Local(DE%Cell,Geom%CellNumber(1),LI%iCenter(1:2))
                    LI%Moment(1:3) = LI%Moment(1:3) + CC%ConservedQuantity(1:3,LI%iCenter(1),LI%iCenter(2),1) * DE%SharedArea !面積×保存量の総和計算
                    LI%Moment(4) = 0.0d0
                    LI%Moment(5) = LI%Moment(5) + CC%ConservedQuantity(4,LI%iCenter(1),LI%iCenter(2),1) * DE%SharedArea !面積×保存量の総和計算

                    LI%ArithmeticAverage(1:3) = CC%ConservedQuantity(1:3,LI%iCenter(1),LI%iCenter(2),1)
                    LI%ArithmeticAverage(4) = 0.0d0
                    LI%ArithmeticAverage(5) = CC%ConservedQuantity(4,LI%iCenter(1),LI%iCenter(2),1)
                end if

                LI%SummationCellNumber = LI%SummationCellNumber + 1
            !end if
            DE => DE%Next !次のレコードへ

        end do

        if(LI%MinLimiterValue >= LI%Threshold) then !interpolated cell not including shock
            LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = LI%Moment/LI%AreaSum

        else !interpolated cell including shock
            LI%ArithmeticAverage = LI%ArithmeticAverage/dble(LI%SummationCellNumber)
            DE => LG%MG(iGrid)%ADE(iTargetCell)%DE !rootへ移動

            LI%UpwindMoment = 0.0d0 !Initialize
            LI%DownwindMoment = 0.0d0
            LI%UpwindArea= 0.0d0
            LI%DownwindArea = 0.0d0


            do while(associated(DE)) !次のレコードがなくなるまで続行
                if(DE%Grid /= 0) then
                    LI%TmpMoment = LG%CC(DE%Grid)%ConservedQuantity(:,DE%Cell,1,1)

                else
                    call Global2Local(DE%Cell,Geom%CellNumber(1),LI%iCenter(1:2))
                    LI%TmpMoment(1:3) = CC%ConservedQuantity(1:3,LI%iCenter(1),LI%iCenter(2),1)
                    LI%TmpMoment(4) = 0.0d0
                    LI%TmpMoment(5) = CC%ConservedQuantity(4,LI%iCenter(1),LI%iCenter(2),1)

                end if

                if(LI%TmpMoment(1) >= LI%ArithmeticAverage(1)) then !Upwind side
                    LI%UpwindMoment = LI%UpwindMoment + LI%TmpMoment * DE%SharedArea
                    LI%UpwindArea = LI%UpwindArea + DE%SharedArea

                else !Downwind Side
                    LI%DownwindMoment = LI%DownwindMoment + LI%TmpMoment * DE%SharedArea
                    LI%DownwindArea = LI%DownwindArea + DE%SharedArea
                end if

                DE => DE%Next !次のレコードへ
            end do

            if(LI%UpwindArea > LI%DownwindArea) then !interpolated cell near to shocks's UPWIND
                LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = LI%UpwindMoment/LI%UpwindArea

            else
                LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = LI%DownwindMoment/LI%DownwindArea
            end if

        end if

    !end do

    !LG%CC(iGrid)%ConservedQuantity = LG%CC(iGrid)%InterpolatedQuantity


    return
    end subroutine ShockConvergeLocalInterpolateMain
