!***********************************/
!	Name:次のモーションウェイト量の計算
!	Alias:OCalculateMotionWait
!	Description:衝撃波が1セル分通過するまで格子移動にウェイトを与える
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.10
!	Update:
!	Other:
!***********************************/
subroutine OCalculateMotionWait(CC,UG,TimeStepWidth,MW,Geom)
    use StructVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, CFL => CourantFriedrichsLewyCondition
    implicit none

    type(CellCenter), intent(in) :: CC
    type(UnstructuredGrid), intent(in) :: UG
    double precision, intent(in) :: TimeStepWidth
    type(MotionWait), intent(inout) :: MW
    type(Geometry), intent(in) :: Geom
    integer  :: iCell
    type(CalculateMotionWaitWithOMP) :: CMW
    integer :: iRequiredWait !衝撃波が通過するために必要なウェイト

    iRequiredWait = 0

!$omp parallel num_threads(CoreNumberOfCPU),shared(UG,CC,TimeStepWidth,Geom),firstprivate(iCell,CMW,iRequiredWait)
!$omp do reduction(max:iRequiredWait)
    do iCell = 1, UG%GI%RealCells
        CMW%AbsoluteVelocity2 = dot_product(CC%ConservedQuantity(2:4,iCell,1,1),CC%ConservedQuantity(2:4,iCell,1,1)) !速度L2ノルム

        CMW%LocalSoundSpeed = sqrt(Gamma*(Gmin1 * CC%ConservedQuantity(5,iCell,1,1) / CC%ConservedQuantity(1,iCell,1,1)) &
                            & + 0.5d0 * CMW%AbsoluteVelocity2) !局所音速

        CMW%AbsoluteVelocity = sqrt(CMW%AbsoluteVelocity2) !速度絶対値

        !CMW%ShockPassingTime = CFL * UG%DeltaXofCFLCondition / (CMW%LocalSoundSpeed + CMW%AbsoluteVelocity) !衝撃波が通過する速度
        !CMW%ShockPassingTime = CFL * sqrt(3.0d0)*UG%GM%AverageWidth(iCell) / (CMW%LocalSoundSpeed + CMW%AbsoluteVelocity) !衝撃波が通過する速度
        !CMW%ShockPassingTime = CFL * sqrt(3.0d0)*sqrt(Geom%Area(1)**2+Geom%Area(2)**2 / (CMW%LocalSoundSpeed + CMW%AbsoluteVelocity) !衝撃波が通過する速度
        CMW%ShockPassingTime = CFL * 2.0d0*max(Geom%Area(1),Geom%Area(2),UG%InscribedCircle(iCell)) / (CMW%LocalSoundSpeed + CMW%AbsoluteVelocity) !衝撃波が通過する速度

        iRequiredWait = max(iRequiredWait, int(CMW%ShockPassingTime/TimeStepWidth))

    end do
!$omp end do
!$omp end parallel

    MW%iMotionWait = iRequiredWait

return
end subroutine OCalculateMotionWait
