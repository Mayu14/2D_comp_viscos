!***********************************/
!	Name:FluxDifferenceSplitting法で用いる左右の固有行列を計算するためのプログラム
!	Alias:GetEigenMatrix
!	Description:特性速度の絶対値を成分に持つ対角行列に左右の固有行列を作用させるといわゆる|A|が得られ，特性速度がどちらに向いていても風上法で扱えるようになる
!	Type:FDSMatrix
!	Input:Geometory,RoeAverage,FDSMatrix
!	Output:FDSMatrix
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.22
!	Other:1次元と2次元に対応
!***********************************/
    subroutine UGetEigenMatrix(iEdge,UG,RA,FM,P4O)
    use StructVar_Mod
    !use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, Epsilon => EntropyCorrection
    implicit none
    integer, intent(in) :: iEdge
    type(UnstructuredGrid), intent(in) :: UG
    type(RoeAverage), intent(inout) :: RA
    type(FDSMatrix), intent(inout) :: FM
    type(PrivateVar4OMP), intent(inout) :: P4O

            P4O%GEM%KinemeticEnergy = 0.5d0*(sum(RA%RoeAverage(2:4)**2))
            P4O%GEM%NormalVelocity = sum(UG%GM%Normal(iEdge,1:3)*RA%RoeAverage(2:4))
            P4O%GEM%InverseSSS = 1.0d0/RA%RoeAverage(1)**2


    if(abs(UG%GM%Normal(iEdge,1)) > abs(UG%GM%Normal(iEdge,2))) then
        call AvoidSingularityOfX
    else
        call AvoidSingularityOfY
    !else
        !call AvoidSingularitiesOfZ !法線ベクトルの絶対値が最大となる成分が分母にくるように呼び出す固有行列を変える
    end if
        FM%AbsEigenValues = 0.0d0
        FM%AbsEigenValues(1,1) = dabs(P4O%GEM%NormalVelocity - RA%RoeAverage(1))
        FM%AbsEigenValues(2,2) = dabs(P4O%GEM%NormalVelocity)
        FM%AbsEigenValues(3,3) = dabs(P4O%GEM%NormalVelocity + RA%RoeAverage(1))
        FM%AbsEigenValues(4,4) = dabs(P4O%GEM%NormalVelocity)
        FM%AbsEigenValues(5,5) = dabs(P4O%GEM%NormalVelocity)

    return
contains
    subroutine AvoidSingularityOfX
    implicit none

        FM%RightEigenMatrix(1,1) = 1.0d0
        FM%RightEigenMatrix(1,2) = 1.0d0
        FM%RightEigenMatrix(1,3) = 1.0d0
        FM%RightEigenMatrix(1,4) = 0.0d0
        FM%RightEigenMatrix(1,5) = 0.0d0

        FM%RightEigenMatrix(2,1) = RA%RoeAverage(2) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,1)
        FM%RightEigenMatrix(2,2) = RA%RoeAverage(2)
        FM%RightEigenMatrix(2,3) = RA%RoeAverage(2) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,1)
        FM%RightEigenMatrix(2,4) = UG%GM%Normal(iEdge,2)
        FM%RightEigenMatrix(2,5) = - UG%GM%Normal(iEdge,3)

        FM%RightEigenMatrix(3,1) = RA%RoeAverage(3) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,2)
        FM%RightEigenMatrix(3,2) = RA%RoeAverage(3)
        FM%RightEigenMatrix(3,3) = RA%RoeAverage(3) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,2)
        FM%RightEigenMatrix(3,4) = - UG%GM%Normal(iEdge,1)
        FM%RightEigenMatrix(3,5) = 0.0d0

        FM%RightEigenMatrix(4,1) = RA%RoeAverage(4) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,3)
        FM%RightEigenMatrix(4,2) = RA%RoeAverage(4)
        FM%RightEigenMatrix(4,3) = RA%RoeAverage(4) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,3)
        FM%RightEigenMatrix(4,4) = 0.0d0
        FM%RightEigenMatrix(4,5) = UG%GM%Normal(iEdge,1)

        FM%RightEigenMatrix(5,1) = (RA%RoeAverage(5)+P4O%GEM%KinemeticEnergy) - RA%RoeAverage(1)*P4O%GEM%NormalVelocity
        FM%RightEigenMatrix(5,2) = P4O%GEM%KinemeticEnergy
        FM%RightEigenMatrix(5,3) = (RA%RoeAverage(5)+P4O%GEM%KinemeticEnergy) + RA%RoeAverage(1)*P4O%GEM%NormalVelocity
        FM%RightEigenMatrix(5,4) = RA%RoeAverage(2)*UG%GM%Normal(iEdge,2) - RA%RoeAverage(3)*UG%GM%Normal(iEdge,1)
        FM%RightEigenMatrix(5,5) = RA%RoeAverage(4)*UG%GM%Normal(iEdge,1) - RA%RoeAverage(2)*UG%GM%Normal(iEdge,3)

            P4O%GEM%CalcAssistB2 = Gmin1*P4O%GEM%InverseSSS
            P4O%GEM%CalcAssistB1 = P4O%GEM%KinemeticEnergy * P4O%GEM%CalcAssistB2

        FM%LeftEigenMatrix(1,1) = 0.5d0*P4O%GEM%InverseSSS*(Gmin1*P4O%GEM%KinemeticEnergy + RA%RoeAverage(1)*P4O%GEM%NormalVelocity)
        FM%LeftEigenMatrix(1,2) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(2) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,1))
        FM%LeftEigenMatrix(1,3) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(3) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,2))
        FM%LeftEigenMatrix(1,4) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(4) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,3))
        FM%LeftEigenMatrix(1,5) = 0.5d0*P4O%GEM%CalcAssistB2

        FM%LeftEigenMatrix(2,1) = 1.0d0-P4O%GEM%CalcAssistB1
        FM%LeftEigenMatrix(2,2) = P4O%GEM%CalcAssistB2 * RA%RoeAverage(2)
        FM%LeftEigenMatrix(2,3) = P4O%GEM%CalcAssistB2 * RA%RoeAverage(3)
        FM%LeftEigenMatrix(2,4) = P4O%GEM%CalcAssistB2 * RA%RoeAverage(4)
        FM%LeftEigenMatrix(2,5) = - P4O%GEM%CalcAssistB2

        FM%LeftEigenMatrix(3,1) = 0.5d0*P4O%GEM%InverseSSS*(Gmin1*P4O%GEM%KinemeticEnergy - RA%RoeAverage(1)*P4O%GEM%NormalVelocity)
        FM%LeftEigenMatrix(3,2) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(2) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,1))
        FM%LeftEigenMatrix(3,3) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(3) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,2))
        FM%LeftEigenMatrix(3,4) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(4) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,3))
        FM%LeftEigenMatrix(3,5) = 0.5d0*P4O%GEM%CalcAssistB2

        FM%LeftEigenMatrix(4,1) = (RA%RoeAverage(3) - P4O%GEM%NormalVelocity*UG%GM%Normal(iEdge,2))/UG%GM%Normal(iEdge,1)
        FM%LeftEigenMatrix(4,2) = UG%GM%Normal(iEdge,2)
        FM%LeftEigenMatrix(4,3) = (UG%GM%Normal(iEdge,2)**2-1.0d0) / UG%GM%Normal(iEdge,1)
        FM%LeftEigenMatrix(4,4) = UG%GM%Normal(iEdge,2)*UG%GM%Normal(iEdge,3)/UG%GM%Normal(iEdge,1)
        FM%LeftEigenMatrix(4,5) = 0.0d0

        FM%LeftEigenMatrix(5,1) = (P4O%GEM%NormalVelocity*UG%GM%Normal(iEdge,3) - RA%RoeAverage(4))/UG%GM%Normal(iEdge,1)
        FM%LeftEigenMatrix(5,2) = - UG%GM%Normal(iEdge,3)
        FM%LeftEigenMatrix(5,3) = - UG%GM%Normal(iEdge,2)*UG%GM%Normal(iEdge,3)/UG%GM%Normal(iEdge,1)
        FM%LeftEigenMatrix(5,4) = (1.0d0 - UG%GM%Normal(iEdge,3)**2) / UG%GM%Normal(iEdge,1)
        FM%LeftEigenMatrix(5,5) = 0.0d0

        return
    end subroutine AvoidSingularityOfX

    subroutine AvoidSingularityOfY
    implicit none

        FM%RightEigenMatrix(1,1) = 1.0d0
        FM%RightEigenMatrix(1,2) = 1.0d0
        FM%RightEigenMatrix(1,3) = 1.0d0
        FM%RightEigenMatrix(1,4) = 0.0d0
        FM%RightEigenMatrix(1,5) = 0.0d0

        FM%RightEigenMatrix(2,1) = RA%RoeAverage(2) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,1)
        FM%RightEigenMatrix(2,2) = RA%RoeAverage(2)
        FM%RightEigenMatrix(2,3) = RA%RoeAverage(2) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,1)
        FM%RightEigenMatrix(2,4) = UG%GM%Normal(iEdge,2)
        FM%RightEigenMatrix(2,5) = 0.0d0

        FM%RightEigenMatrix(3,1) = RA%RoeAverage(3) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,2)
        FM%RightEigenMatrix(3,2) = RA%RoeAverage(3)
        FM%RightEigenMatrix(3,3) = RA%RoeAverage(3) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,2)
        FM%RightEigenMatrix(3,4) = - UG%GM%Normal(iEdge,1)
        FM%RightEigenMatrix(3,5) = UG%GM%Normal(iEdge,3)

        FM%RightEigenMatrix(4,1) = RA%RoeAverage(4) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,3)
        FM%RightEigenMatrix(4,2) = RA%RoeAverage(4)
        FM%RightEigenMatrix(4,3) = RA%RoeAverage(4) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,3)
        FM%RightEigenMatrix(4,4) = 0.0d0
        FM%RightEigenMatrix(4,5) = - UG%GM%Normal(iEdge,2)


        FM%RightEigenMatrix(5,1) = (RA%RoeAverage(5)+P4O%GEM%KinemeticEnergy) - RA%RoeAverage(1)*P4O%GEM%NormalVelocity
        FM%RightEigenMatrix(5,2) = P4O%GEM%KinemeticEnergy
        FM%RightEigenMatrix(5,3) = (RA%RoeAverage(5)+P4O%GEM%KinemeticEnergy) + RA%RoeAverage(1)*P4O%GEM%NormalVelocity
        FM%RightEigenMatrix(5,4) = RA%RoeAverage(2)*UG%GM%Normal(iEdge,2) - RA%RoeAverage(3)*UG%GM%Normal(iEdge,1)
        FM%RightEigenMatrix(5,5) = RA%RoeAverage(3)*UG%GM%Normal(iEdge,3) - RA%RoeAverage(4)*UG%GM%Normal(iEdge,2)

            P4O%GEM%CalcAssistB2 = Gmin1*P4O%GEM%InverseSSS
            P4O%GEM%CalcAssistB1 = P4O%GEM%KinemeticEnergy * P4O%GEM%CalcAssistB2

        FM%LeftEigenMatrix(1,1) = 0.5d0*P4O%GEM%InverseSSS*(Gmin1*P4O%GEM%KinemeticEnergy + RA%RoeAverage(1)*P4O%GEM%NormalVelocity)
        FM%LeftEigenMatrix(1,2) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(2) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,1))
        FM%LeftEigenMatrix(1,3) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(3) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,2))
        FM%LeftEigenMatrix(1,4) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(4) - RA%RoeAverage(1)*UG%GM%Normal(iEdge,3))
        FM%LeftEigenMatrix(1,5) = 0.5d0*P4O%GEM%CalcAssistB2

        FM%LeftEigenMatrix(2,1) = 1.0d0-P4O%GEM%CalcAssistB1
        FM%LeftEigenMatrix(2,2) = P4O%GEM%CalcAssistB2 * RA%RoeAverage(2)
        FM%LeftEigenMatrix(2,3) = P4O%GEM%CalcAssistB2 * RA%RoeAverage(3)
        FM%LeftEigenMatrix(2,4) = P4O%GEM%CalcAssistB2 * RA%RoeAverage(4)
        FM%LeftEigenMatrix(2,5) = - P4O%GEM%CalcAssistB2

        FM%LeftEigenMatrix(3,1) = 0.5d0*P4O%GEM%InverseSSS*(Gmin1*P4O%GEM%KinemeticEnergy - RA%RoeAverage(1)*P4O%GEM%NormalVelocity)
        FM%LeftEigenMatrix(3,2) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(2) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,1))
        FM%LeftEigenMatrix(3,3) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(3) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,2))
        FM%LeftEigenMatrix(3,4) = 0.5d0*P4O%GEM%InverseSSS*(-Gmin1*RA%RoeAverage(4) + RA%RoeAverage(1)*UG%GM%Normal(iEdge,3))
        FM%LeftEigenMatrix(3,5) = 0.5d0*P4O%GEM%CalcAssistB2

        FM%LeftEigenMatrix(4,1) = (P4O%GEM%NormalVelocity*UG%GM%Normal(iEdge,1) - RA%RoeAverage(2))/UG%GM%Normal(iEdge,2)
        FM%LeftEigenMatrix(4,2) = (1.0d0 - UG%GM%Normal(iEdge,1)**2)/UG%GM%Normal(iEdge,2)
        FM%LeftEigenMatrix(4,3) = - UG%GM%Normal(iEdge,1)
        FM%LeftEigenMatrix(4,4) = - UG%GM%Normal(iEdge,1)*UG%GM%Normal(iEdge,3)/UG%GM%Normal(iEdge,2)
        FM%LeftEigenMatrix(4,5) = 0.0d0

        FM%LeftEigenMatrix(5,1) = (RA%RoeAverage(4) - P4O%GEM%NormalVelocity*UG%GM%Normal(iEdge,3))/UG%GM%Normal(iEdge,2)
        FM%LeftEigenMatrix(5,2) = UG%GM%Normal(iEdge,1)*UG%GM%Normal(iEdge,3)/UG%GM%Normal(iEdge,2)
        FM%LeftEigenMatrix(5,3) = UG%GM%Normal(iEdge,3)
        FM%LeftEigenMatrix(5,4) = (UG%GM%Normal(iEdge,3)**2 - 1.0d0) / UG%GM%Normal(iEdge,2)
        FM%LeftEigenMatrix(5,5) = 0.0d0

        return
    end subroutine AvoidSingularityOfY


end subroutine UGetEigenMatrix
