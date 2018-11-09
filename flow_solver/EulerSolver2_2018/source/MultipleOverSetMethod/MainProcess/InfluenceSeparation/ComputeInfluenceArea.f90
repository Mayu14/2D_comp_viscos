!***********************************/
!	Name:影響域分離法で用いる影響域の計算
!	Alias:ComputeInfluenceArea
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.31
!	Update:
!	Other:
!***********************************/
    subroutine ComputeInfluenceArea(UG,CC,MW,MG)
    use StructVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none

    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(in) :: CC
    type(MotionWait), intent(in) :: MW
    type(MoveGrid), intent(inout) :: MG
    type(ComputeInfluenceAreaWithOMP) :: CIA
    integer :: iCell
    integer :: iLoop
    double precision :: MaxInfluenceRange, InfluenceRange
    !MG%IS(iCell)%InfluenceDepth = 0 !Initialize routine

    allocate(CIA%FlowVelcity(3),CIA%BodyVelocity(3),CIA%UnitDirectionVector(3))

        !物体の移動速度ベクトルを格納する
        CIA%BodyVelocity = MW%StepDistance(1:3)/FixedTimeStep
        MaxInfluenceRange = 0.0d0

        if(MG%iStepFromLastMove == 1) then !1st step
            do iCell = UG%GI%RealCells+1, UG%GI%AllCells
                if(UG%GM%CellType(iCell,1,1) == 2) then !Wall Boundary = Internal Object Surface's Virtual Cell
                    MG%IS(UG%VC%Cell(iCell,1))%InfluenceDepth = 1
                    do iLoop = 1, 3
                        if(UG%Tri%Cell(UG%VC%Cell(iCell,1),iLoop) < UG%GI%RealCells + 1) then !Real Adjacent-Cell
                            MG%IS(UG%Tri%Cell(UG%VC%Cell(iCell,1),iLoop))%InfluenceDepth = 1
                        end if
                    end do
                end if
            end do

            MG%MaxInfluenceRange = UG%InternalRadius + maxval(UG%InscribedCircle,1) + 2.0d0*maxval(UG%GM%AverageWidth,1)

        else !not first step
!$omp parallel num_threads(CoreNumberOfCPU),shared(UG,CC,MW,MG),firstprivate(iCell,CIA,InfluenceRange)
!$omp do reduction(max:MaxInfluenceRange)
        do iCell = 1, UG%GI%RealCells
            if(MG%IS(iCell)%InfluenceDepth /= 0) then
                MG%IS(iCell)%InfluenceDepth = MG%IS(iCell)%InfluenceDepth + 1 !影響域の深度を上げる(内挿時の重み付けに影響)

            else
                !物体中心からの距離を求める．
                CIA%Distance = sqrt(dot_product(MG%RC%Cell(iCell,1:3),MG%RC%Cell(iCell,1:3))) ! - UG%InscribedCircle(iCell)
                !単位方向ベクトルを求める
                CIA%UnitDirectionVector = MG%RC%Cell(iCell,1:3) / CIA%Distance
                !当該セル方向のサポート写像を求める
                !CIA%ContourOfBody = SupportFunction(UG%CH,CIA%UnitDirectionVector)
                CIA%ContourOfBody = UG%InternalRadius + UG%InscribedCircle(iCell) !円の場合のみ直打ち込みで対応可能

                !当該セルの主流速度ベクトルを格納する
                CIA%FlowVelcity = CC%ConservedQuantity(2:4,iCell,1,1) / CC%ConservedQuantity(1,iCell,1,1)

                !当該セル内の音速を求める
                CIA%LocalSoundSpeed = sqrt(Gamma*Gmin1 * (CC%ConservedQuantity(5,iCell,1,1) / CC%ConservedQuantity(1,iCell,1,1) &
                                    & - 0.5d0 * dot_product(CIA%FlowVelcity,CIA%FlowVelcity))) !局所音速
                !当該セル内の(主流+物体)速度と単位方向ベクトルとの内積を取る
                CIA%NormalVelocity = dot_product(-CIA%BodyVelocity + CIA%FlowVelcity,CIA%UnitDirectionVector)
                !CIA%NormalVelocity = dot_product(-CIA%BodyVelocity + CIA%FlowVelcity,-CIA%BodyVelocity + CIA%FlowVelcity)
                !CIA%NormalVelocity = dot_product(abs(CIA%BodyVelocity) + abs(CIA%FlowVelcity),abs(CIA%BodyVelocity) + abs(CIA%FlowVelcity))
                !CIA%NormalVelocity = max(dot_product(-CIA%BodyVelocity + CIA%FlowVelcity,CIA%UnitDirectionVector),dot_product(-CIA%BodyVelocity + CIA%FlowVelcity,-CIA%BodyVelocity + CIA%FlowVelcity))
                !判定を行う
                InfluenceRange = CIA%ContourOfBody + (CIA%NormalVelocity + CIA%LocalSoundSpeed)*FixedTimeStep*dble(MG%iStepFromLastMove + 1) !&
                            !&     + min(2.0d0*UG%GM%AverageWidth(iCell)*dble(MG%iStepFromLastMove)/dble(MW%iEstimateWait),2.0d0*UG%GM%AverageWidth(iCell))

                MaxInfluenceRange = max(MaxInfluenceRange,InfluenceRange)

                if(InfluenceRange - (CIA%Distance  - 6.0d0*UG%GM%AverageWidth(iCell)) > 0.0d0) then !in the influence region
                    !write(6,*) "Influence"
                    MG%IS(iCell)%InfluenceDepth = 1
                end if

            end if
        end do
!$omp end do
!$omp end parallel

    MG%MaxInfluenceRange = MaxInfluenceRange
    end if

    return
    end subroutine ComputeInfluenceArea
