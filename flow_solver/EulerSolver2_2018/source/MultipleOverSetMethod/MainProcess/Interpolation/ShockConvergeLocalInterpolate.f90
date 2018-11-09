!***********************************/
!	Name:移動先の格子を補間するための補間データ提供元格子検索
!	Alias:ShockConvergeLocalInterpolate
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
    subroutine ShockConvergeLocalInterpolate(iGrid,LG,CC,Geom,iStep,MConf)
    use StructVar_Mod
    implicit none
    integer, intent(in) :: iGrid
    type(LocalGrids), intent(inout) :: LG
    type(CellCenter), intent(in) :: CC
    type(Geometry), intent(in) :: Geom
    integer, intent(in) :: iStep
    type(Configulation), intent(in) :: MConf
    type(DonorElement), pointer :: DE
    integer :: iTargetCell
    integer :: iCenterX, iCenterY
    double precision :: AreaSum
    double precision, allocatable :: Moment(:)

    double precision :: MinLimiterValue
    double precision :: Threshold = 0.0005
    double precision, allocatable :: ArithmeticAverage(:)
    integer :: SummationCellNumber
    double precision :: UpwindArea, DownwindArea
    double precision, allocatable :: UpwindQuantity(:),DownwindQuantity(:),TmpQuantity(:)

!    write(6,*) "Cell value interpolation of Local Grid"
    allocate(Moment(5),ArithmeticAverage(5))
    allocate(UpwindQuantity(5),DownwindQuantity(5),TmpQuantity(5))

    LG%CC(iGrid)%InterpolatedQuantity = LG%CC(iGrid)%ConservedQuantity
    if(iStep == 1) then
        call UConserve2Primitive(LG%UG(iGrid),LG%CC(iGrid))
        call UGetGradient(LG%UG(iGrid),LG%CC(iGrid))
        call UGetLimiter(MConf,LG%UG(iGrid),LG%CC(iGrid),LG%CE(iGrid))
    end if

    !全セルについて補間を実行
    do iTargetCell=1, LG%UG(iGrid)%GI%RealCells

        AreaSum = 0.0d0
        Moment = 0.0d0
        MinLimiterValue = 1.0d0
        ArithmeticAverage = 0.0d0
        SummationCellNumber = 0
        DE => LG%MG(iGrid)%ADE(iTargetCell)%DE !rootへ移動

        do while(associated(DE)) !次のレコードがなくなるまで続行

            AreaSum = AreaSum + DE%SharedArea !面積の総和計算

            if(DE%Grid > 0) then !0番はGlobalGridであり，ConservedQuantityの配列構造が異なる，詳細はRegisterGlobal2Local参照のこと

                MinLimiterValue = min(MinLimiterValue,minval(LG%CC(iGrid)%LimiterFunction(:,DE%Cell,1,1)))

                Moment = Moment + LG%CC(DE%Grid)%ConservedQuantity(:,DE%Cell,1,1) * DE%SharedArea !面積×保存量の総和計算
                ArithmeticAverage = ArithmeticAverage + LG%CC(DE%Grid)%ConservedQuantity(:,DE%Cell,1,1)


            else

                iCenterX = DE%Cell - int(DE%Cell/Geom%CellNumber(1))*Geom%CellNumber(1)
                iCenterY = int(DE%Cell/Geom%CellNumber(1)) + 1

                Moment(1:3) = Moment(1:3) + CC%ConservedQuantity(1:3,iCenterX,iCenterY,1) * DE%SharedArea !面積×保存量の総和計算
                Moment(4) = 0.0d0
                Moment(5) = Moment(5) + CC%ConservedQuantity(4,iCenterX,iCenterY,1) * DE%SharedArea !面積×保存量の総和計算

                ArithmeticAverage(1:3) = CC%ConservedQuantity(1:3,iCenterX,iCenterY,1)
                ArithmeticAverage(4) = 0.0d0
                ArithmeticAverage(5) = CC%ConservedQuantity(4,iCenterX,iCenterY,1)
            end if

            SummationCellNumber = SummationCellNumber + 1
            DE => DE%Next !次のレコードへ

        end do

        if(MinLimiterValue > Threshold) then !interpolated cell not including shock
            LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = Moment/AreaSum

        else !interpolated cell including shock
            ArithmeticAverage = ArithmeticAverage/dble(SummationCellNumber)
            DE => LG%MG(iGrid)%ADE(iTargetCell)%DE !rootへ移動

            UpwindQuantity = 0.0d0 !Initialize
            DownwindQuantity = 0.0d0
            UpwindArea= 0.0d0
            DownwindArea = 0.0d0


            do while(associated(DE)) !次のレコードがなくなるまで続行
                if(DE%Grid /= 0) then
                    TmpQuantity = LG%CC(DE%Grid)%ConservedQuantity(:,DE%Cell,1,1)

                else
                    iCenterX = DE%Cell - int(DE%Cell/Geom%CellNumber(1))*Geom%CellNumber(1)
                    iCenterY = int(DE%Cell/Geom%CellNumber(1)) + 1
                    TmpQuantity(1:3) = CC%ConservedQuantity(1:3,iCenterX,iCenterY,1)
                    TmpQuantity(4) = 0.0d0
                    TmpQuantity(5) = CC%ConservedQuantity(4,iCenterX,iCenterY,1)

                end if

                if(TmpQuantity(1) >= ArithmeticAverage(1)) then !Upwind side
                    UpwindQuantity = UpwindQuantity + TmpQuantity * DE%SharedArea
                    UpwindArea = UpwindArea + DE%SharedArea

                else !Downwind Side
                    DownwindQuantity = DownwindQuantity + TmpQuantity * DE%SharedArea
                    DownwindArea = DownwindArea + DE%SharedArea
                end if

                DE => DE%Next !次のレコードへ
            end do

            if(UpwindArea > DownwindArea) then !interpolated cell near to shocks's UPWIND
                LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = UpwindQuantity/UpwindArea

            else
                LG%CC(iGrid)%InterpolatedQuantity(:,iTargetCell,1,1) = DownwindQuantity/DownwindArea
            end if

        end if

    end do

    LG%CC(iGrid)%ConservedQuantity = LG%CC(iGrid)%InterpolatedQuantity


    return
    end subroutine ShockConvergeLocalInterpolate
