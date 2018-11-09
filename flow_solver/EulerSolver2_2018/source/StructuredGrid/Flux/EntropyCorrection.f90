subroutine EntropyCorrection(Conf,Geom,RA,FM) !没(？)未完成
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Epsilon => EntropyCorrection
    implicit none

    type(Configulation), intent(in) :: Conf
    type(Geometry), intent(in) :: Geom
    type(RoeAverage), intent(in) :: RA
    type(FDSMatrix), intent(inout) :: FM

    !if(Conf%KindEntropyCorrect == 1)
        do iLoop=1,Geom%Dimension+2
            if(FM%AbsEigenValues(iLoop,iLoop) < Epsilon) then
                FM%AbsEigenValues(iLoop,iLoop) = 0.5d0*(FM%AbsEigenValues(iLoop,iLoop)**2 + Epsilon**2)/Epsilon
            end if
        end do
    !else if(Conf%KindEntropyCorrect == 2)
        do iLoop=1,Geom%Dimension+2
            !音速求める
            !u-c,u,u+c求める
            if(RA%SideQuantity)
        end do

return
end subroutine EntropyCorrection
