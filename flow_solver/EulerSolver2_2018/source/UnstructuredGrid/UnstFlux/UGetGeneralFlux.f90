!***********************************/
!	Name:Roeスキームによって中心差分の流束を計算するプログラム
!	Alias:GetGeneralFlux
!	Description:計算結果はRoe流束の第1項に相当する
!	Type:RoeAverage
!	Input:Geometory,RoeAverage
!	Output:RoeAverage
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.22
!	Other:1次元と2次元に対応
!***********************************/
subroutine UGetGeneralFlux(iEdge,UG,RA,P4O)
    use StructVar_Mod
    !use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, Epsilon => EntropyCorrection
    implicit none
    integer, intent(in) :: iEdge
    type(UnstructuredGrid), intent(in) :: UG
    type(RoeAverage), intent(inout) :: RA
    type(PrivateVar4OMP), intent(inout) :: P4O

    P4O%GGF%InverseDensity(:) = 1.0d0/RA%SideQuantity(1,:)
    P4O%GGF%Velocity(1,:) = RA%SideQuantity(2,:)*P4O%GGF%InverseDensity
    P4O%GGF%Velocity(2,:) = RA%SideQuantity(3,:)*P4O%GGF%InverseDensity
    P4O%GGF%Velocity(3,:) = RA%SideQuantity(4,:)*P4O%GGF%InverseDensity

    P4O%GGF%KinemeticEnergy(:) = 0.5d0*(P4O%GGF%Velocity(1,:)**2+P4O%GGF%Velocity(2,:)**2+P4O%GGF%Velocity(3,:)**2)

    P4O%GGF%NormalVelocity(:) = (P4O%GGF%Velocity(1,:)*UG%GM%Normal(iEdge,1) + P4O%GGF%Velocity(2,:)*UG%GM%Normal(iEdge,2) &
        & + P4O%GGF%Velocity(3,:)*UG%GM%Normal(iEdge,3))

    P4O%GGF%TotalEnergy(:) = Gamma*RA%SideQuantity(5,:) + P4O%GGF%KinemeticEnergy(:)

    P4O%GGF%FluxJacobianMatrix(1,1,:) = 0.0d0
    P4O%GGF%FluxJacobianMatrix(1,2,:) = UG%GM%Normal(iEdge,1)
    P4O%GGF%FluxJacobianMatrix(1,3,:) = UG%GM%Normal(iEdge,2)
    P4O%GGF%FluxJacobianMatrix(1,4,:) = UG%GM%Normal(iEdge,3)
    P4O%GGF%FluxJacobianMatrix(1,5,:) = 0.0d0

    P4O%GGF%FluxJacobianMatrix(2,1,:) = Gmin1*P4O%GGF%KinemeticEnergy(:)*UG%GM%Normal(iEdge,1) - P4O%GGF%Velocity(1,:)*P4O%GGF%NormalVelocity(:)
    P4O%GGF%FluxJacobianMatrix(2,2,:) = P4O%GGF%NormalVelocity(:) - (Gamma-2.0d0)*P4O%GGF%Velocity(1,:)*UG%GM%Normal(iEdge,1)
    P4O%GGF%FluxJacobianMatrix(2,3,:) = P4O%GGF%Velocity(1,:)*UG%GM%Normal(iEdge,2) - Gmin1*P4O%GGF%Velocity(2,:)*UG%GM%Normal(iEdge,1)
    P4O%GGF%FluxJacobianMatrix(2,4,:) = P4O%GGF%Velocity(1,:)*UG%GM%Normal(iEdge,3) - Gmin1*P4O%GGF%Velocity(3,:)*UG%GM%Normal(iEdge,1)
    P4O%GGF%FluxJacobianMatrix(2,5,:) = Gmin1*UG%GM%Normal(iEdge,1)

    P4O%GGF%FluxJacobianMatrix(3,1,:) = Gmin1*P4O%GGF%KinemeticEnergy*UG%GM%Normal(iEdge,2) - P4O%GGF%Velocity(2,:)*P4O%GGF%NormalVelocity(:)
    P4O%GGF%FluxJacobianMatrix(3,2,:) = P4O%GGF%Velocity(2,:)*UG%GM%Normal(iEdge,1) - Gmin1*P4O%GGF%Velocity(1,:)*UG%GM%Normal(iEdge,2)
    P4O%GGF%FluxJacobianMatrix(3,3,:) = P4O%GGF%NormalVelocity(:)-(Gamma-2.0d0)*P4O%GGF%Velocity(2,:)*UG%GM%Normal(iEdge,2)
    P4O%GGF%FluxJacobianMatrix(3,4,:) = P4O%GGF%Velocity(2,:)*UG%GM%Normal(iEdge,3) - Gmin1*P4O%GGF%Velocity(3,:)*UG%GM%Normal(iEdge,2)
    P4O%GGF%FluxJacobianMatrix(3,5,:) = Gmin1*UG%GM%Normal(iEdge,2)

    P4O%GGF%FluxJacobianMatrix(4,1,:) = GMin1*P4O%GGF%KinemeticEnergy*UG%GM%Normal(iEdge,3) - P4O%GGF%Velocity(3,:)*P4O%GGF%NormalVelocity(:)
    P4O%GGF%FluxJacobianMatrix(4,2,:) = P4O%GGF%Velocity(3,:)*UG%GM%Normal(iEdge,1) - Gmin1*P4O%GGF%Velocity(1,:)*UG%GM%Normal(iEdge,3)
    P4O%GGF%FluxJacobianMatrix(4,3,:) = P4O%GGF%Velocity(3,:)*UG%GM%Normal(iEdge,2) - Gmin1*P4O%GGF%Velocity(2,:)*UG%GM%Normal(iEdge,3)
    P4O%GGF%FluxJacobianMatrix(4,4,:) = P4O%GGF%NormalVelocity(:) - (Gamma-2.0d0)*P4O%GGF%Velocity(3,:)*UG%GM%Normal(iEdge,3)
    P4O%GGF%FluxJacobianMatrix(4,5,:) = Gmin1*UG%GM%Normal(iEdge,3)

    P4O%GGF%FluxJacobianMatrix(5,1,:) = (Gmin1*P4O%GGF%KinemeticEnergy(:) - P4O%GGF%TotalEnergy(:))*P4O%GGF%NormalVelocity(:)
    P4O%GGF%FluxJacobianMatrix(5,2,:) = P4O%GGF%TotalEnergy(:)*UG%GM%Normal(iEdge,1) - Gmin1*P4O%GGF%Velocity(1,:)*P4O%GGF%NormalVelocity(:)
    P4O%GGF%FluxJacobianMatrix(5,3,:) = P4O%GGF%TotalEnergy(:)*UG%GM%Normal(iEdge,2) - Gmin1*P4O%GGF%Velocity(2,:)*P4O%GGF%NormalVelocity(:)
    P4O%GGF%FluxJacobianMatrix(5,4,:) = P4O%GGF%TotalEnergy(:)*UG%GM%Normal(iEdge,3) - Gmin1*P4O%GGF%Velocity(3,:)*P4O%GGF%NormalVelocity(:)
    P4O%GGF%FluxJacobianMatrix(5,5,:) = Gamma*P4O%GGF%NormalVelocity(:)

    RA%GeneralFlux(:,1,1) = matmul(P4O%GGF%FluxJacobianMatrix(:,:,1),RA%SideQuantity(:,1))
    RA%GeneralFlux(:,1,2) = matmul(P4O%GGF%FluxJacobianMatrix(:,:,2),RA%SideQuantity(:,2))

    return
end subroutine UGetGeneralFlux
