!***********************************/
!	Name:移動先の格子を補間するための補間データ提供元格子検索
!	Alias:InterpolationOfLocalGrid
!	Description:
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.11.16
!	Update:2017.11.16
!	Other:
!***********************************/
subroutine InterpolationOfLocalGrid(iTargetCell,iGrid,LG,CC,Geom,MConf,LI)
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
    integer, allocatable :: Local(:)
    integer :: iHalfDepth, iExponent
    integer :: iLoop, iCheck
    double precision, allocatable :: SekoQuantity(:)
    allocate(Local(2), SekoQuantity(5))
        iExponent = 1


        LI%AreaSum = 0.0d0
        LI%Moment = 0.0d0
        DE => LG%MG(iGrid)%ADE(iTargetCell)%DE !rootへ移動

        do while(associated(DE)) !次のレコードがなくなるまで続行
            if(LG%MW(iGrid)%iMotionWait /= 0 .and. DE%Grid == iGrid) then !のとき内挿対象とみなさない
                !何もせず次のレコードへ
                write(6,*) "nnaaaaaaahhh"
            else

                LI%AreaSum = LI%AreaSum + DE%SharedArea !面積の総和計算

                if(DE%Grid > 0) then !0番はGlobalGridであり，ConservedQuantityの配列構造が異なる，詳細はRegisterGlobal2Local参照のこと

                    LI%Moment = LI%Moment + LG%CC(DE%Grid)%ConservedQuantity(:,DE%Cell,1,1) * DE%SharedArea !面積×保存量の総和計算

                else
                    call Global2Local(DE%Cell,Geom%CellNumber(1),LI%iCenter(1:2))

                    LI%Moment(1:3) = LI%Moment(1:3) + CC%ConservedQuantity(1:3,LI%iCenter(1),LI%iCenter(2),1) * DE%SharedArea !面積×保存量の総和計算
                    LI%Moment(4) = 0.0d0
                    LI%Moment(5) = LI%Moment(5) + CC%ConservedQuantity(4,LI%iCenter(1),LI%iCenter(2),1) * DE%SharedArea !面積×保存量の総和計算
                end if
            end if
            DE => DE%Next !次のレコードへ
        end do

        if(LG%UG(iGrid)%InternalBoundary == 0) then !内部に物体を持たないとき
            call NormalInterpolate

        else !物体があるとき
            if(iTargetCell > LG%UG(iGrid)%GI%RealCells) then !仮想セルは標準の処理
                call NormalInterpolate

            else !実セル
                if(LG%MG(iGrid)%IS(iTargetCell)%InfluenceDepth == 0) then !影響域でないとき標準の処理
                    iCheck = 0 !if an Adjacent-Cell is in the influence area then Use Sigmoid Merge
                    SekoQuantity = 0.0d0
                    do iLoop = 1, 3
                        if(LG%UG(iGrid)%Tri%Cell(iTargetCell,iLoop) < LG%UG(iGrid)%GI%RealCells+1) exit !skip virtual cell

                        if(LG%MG(iGrid)%IS(LG%UG(iGrid)%Tri%Cell(iTargetCell,iLoop))%InfluenceDepth /= 0) then !adjacent cell is in the influence area
                            call SpreadInterpolate(1,LG%UG(iGrid)%Tri%Cell(iTargetCell,iLoop),0)
                        else
                            call SpreadInterpolate(2,LG%UG(iGrid)%Tri%Cell(iTargetCell,iLoop),0)
                            iCheck = iCheck + 1 !Normal Cell Number
                            LG%MG(iGrid)%IS(iTargetCell)%BoundaryAdjacent(iCheck) = LG%UG(iGrid)%Tri%Cell(iTargetCell,iLoop)
                        end if

                    end do

                    if(iCheck == 3) then
                        call NormalInterpolate
                    else
                        call SpreadInterpolate(3,iTargetCell,3)
                        call NormalInterpolate
                        !call InterpolateWithInfluenceMK2
                        !call SigmoidMerge
                    end if

                else
                    !call InterpolateWithInfluenceRegion
                    call InterpolateWithInfluenceMK2
                    !call SigmoidMerge
                end if

            end if

        end if

    return
contains
    subroutine SpreadInterpolate(iSwitch,iAdjacent,iCount) !this is SEKO program
    implicit none
    integer, intent(in) :: iSwitch, iAdjacent, iCount

        if(iSwitch == 1) then
            SekoQuantity = SekoQuantity + LG%CC4MB(iGrid)%ConservedQuantity(:,iAdjacent,1,1)

        else if(iSwitch == 2) then
            SekoQuantity = SekoQuantity + LG%CC(iGrid)%InterpolatedQuantity(:,iAdjacent,1,1)

        else if(iSwitch == 3) then
            LG%CC4MB(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = SekoQuantity/dble(iCount)

        end if

    return
    end subroutine SpreadInterpolate


    subroutine NormalInterpolate
    implicit none
        LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = LI%Moment/LI%AreaSum !標準の処理
    return
    end subroutine NormalInterpolate


    subroutine InterpolateWithInfluenceRegion
    implicit none

        !※重み付き平均
        LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = &
                    &   (LG%CC4MB(LG%UG(iGrid)%InternalBoundary)%ConservedQuantity(:,iTargetCell,1,1)*dble(LG%MG(iGrid)%IS(iTargetCell)%InfluenceDepth**iExponent) &
                    & + (LI%Moment/LI%AreaSum)*dble(LG%MG(iGrid)%iStepFromLastMove - LG%MG(iGrid)%IS(iTargetCell)%InfluenceDepth)**iExponent) &
                    & / (dble(LG%MG(iGrid)%iStepFromLastMove**iExponent))

    return
    end subroutine InterpolateWithInfluenceRegion

    subroutine InterpolateWithInfluenceMK2
    implicit none

        LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = &
                    &   LG%CC4MB(LG%UG(iGrid)%InternalBoundary)%ConservedQuantity(:,iTargetCell,1,1)

    return
    end subroutine InterpolateWithInfluenceMK2


    subroutine SigmoidMerge
    use ConstantVar_Mod, Gain => SigmoidGain
    implicit none
    !マージの重みもシグモイド関数にしてみる
    double precision :: Weight, Distance, SigmoidX !SigmoidRange is defined in the Constant Variables module
    double precision :: MergeGain

        !write(6,*) Gain, LG%MG(iGrid)%IS(iTargetCell)%InfluenceDepth, LG%MG(iGrid)%iStepFromLastMove
        Distance = sqrt(dot_product(LG%MG(iGrid)%RC%Cell(iTargetCell,1:3),LG%MG(iGrid)%RC%Cell(iTargetCell,1:3)))
        MergeGain = -2.0d0/(2.0d0*LG%UG(iGrid)%GM%AverageWidth(iTargetCell))*log((3.0d0-sqrt(3.0d0))/(3.0d0+sqrt(3.0d0)))

        if(iTargetCell < LG%UG(iGrid)%GI%RealCells + 1) then
            SigmoidX = MergeGain*(Distance - (LG%MG(iGrid)%MaxInfluenceRange + (LG%UG(iGrid)%GM%AverageWidth(iTargetCell))))
        else
            SigmoidX = MergeGain*(Distance - (LG%MG(iGrid)%MaxInfluenceRange + (LG%UG(iGrid)%GM%AverageWidth(LG%UG(iGrid)%VC%Cell(iTargetCell,1)))))
        end if


            if(SigmoidX < - SigmoidRange) then !Error procedure 1
                Weight = 1.0d0 - 10.0d0**(-15)

            else if(SigmoidX > SigmoidRange) then !Error procedure 2
                Weight = 10.0d0**(-15)

            else !Default procedure
                Weight = -1.0d0 / (1.0d0 + exp(-SigmoidX)) + 1.0d0
            end if

        !Weight = - 1.0d0/(1.0d0 + exp(-Gain*dble(LG%MG(iGrid)%IS(iTargetCell)%InfluenceDepth - LG%MG(iGrid)%iStepFromLastMove)))
        !write(6,*) Weight
        LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = &
                    &   (LG%CC4MB(LG%UG(iGrid)%InternalBoundary)%ConservedQuantity(:,iTargetCell,1,1)*Weight) &
                    & + (LI%Moment/LI%AreaSum)*(1.0d0 - Weight)
!    write(6,*) LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1)
!stop
    return
    end subroutine SigmoidMerge

end subroutine InterpolationOfLocalGrid
