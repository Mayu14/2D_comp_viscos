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
!	Update:2017.11.09
!	Other:1次元と2次元に対応
!***********************************/
    subroutine GetEigenMatrix(Geom,RA,FM,SP4O)
    use StructVar_Mod
    !use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, Epsilon => EntropyCorrection
    implicit none
    type(Geometry), intent(in) :: Geom
    type(RoeAverage), intent(in) :: RA
    type(FDSMatrix), intent(inout) :: FM
    type(PrivateVar4OMP), intent(inout) :: SP4O

    if(Geom%Dimension == 1) then
        FM%RightEigenMatrix(1,1) = 1.0d0
        FM%RightEigenMatrix(1,2) = 1.0d0
        FM%RightEigenMatrix(1,3) = 1.0d0
        FM%RightEigenMatrix(2,1) = RA%RoeAverage(2) - RA%RoeAverage(1)
        FM%RightEigenMatrix(2,2) = RA%RoeAverage(2)
        FM%RightEigenMatrix(2,3) = RA%RoeAverage(2) + RA%RoeAverage(1)
        FM%RightEigenMatrix(3,1) = RA%RoeAverage(3) - RA%RoeAverage(2)*RA%RoeAverage(1)
        FM%RightEigenMatrix(3,2) = 0.5d0 * RA%RoeAverage(2)**2
        FM%RightEigenMatrix(3,3) = RA%RoeAverage(3) + RA%RoeAverage(2)*RA%RoeAverage(1)

        SP4O%GEM%CalcAssistB2 = Gmin1/RA%RoeAverage(1)**2
        SP4O%GEM%CalcAssistB1 = 0.5d0*RA%RoeAverage(2)**2 * SP4O%GEM%CalcAssistB2

        FM%LeftEigenMatrix(1,1) = 0.5d0*(SP4O%GEM%CalcAssistB1 + RA%RoeAverage(2)/RA%RoeAverage(1))
        FM%LeftEigenMatrix(1,2) = -0.5d0*(1.0d0/RA%RoeAverage(1) + SP4O%GEM%CalcAssistB2*RA%RoeAverage(2))
        FM%LeftEigenMatrix(1,3) = 0.5d0*(SP4O%GEM%CalcAssistB2)
        FM%LeftEigenMatrix(2,1) = 1.0d0 - SP4O%GEM%CalcAssistB1
        FM%LeftEigenMatrix(2,2) = SP4O%GEM%CalcAssistB2 * RA%RoeAverage(2)
        FM%LeftEigenMatrix(2,3) = -SP4O%GEM%CalcAssistB2
        FM%LeftEigenMatrix(3,1) = 0.5d0*(SP4O%GEM%CalcAssistB1 - RA%RoeAverage(2)/RA%RoeAverage(1))
        FM%LeftEigenMatrix(3,2) = 0.5d0*(1.0d0/RA%RoeAverage(1) - SP4O%GEM%CalcAssistB2*RA%RoeAverage(2))
        FM%LeftEigenMatrix(3,3) = 0.5d0*(SP4O%GEM%CalcAssistB2)

        FM%AbsEigenValues = 0.0d0
        FM%AbsEigenValues(1,1) = dabs(RA%RoeAverage(2) - RA%RoeAverage(1))
        FM%AbsEigenValues(2,2) = dabs(RA%RoeAverage(2))
        FM%AbsEigenValues(3,3) = dabs(RA%RoeAverage(2) + RA%RoeAverage(1))

    else if(Geom%Dimension == 2) then
        if(RA%Direction == 1) then
            FM%RightEigenMatrix(1,1) = 1.0d0
            FM%RightEigenMatrix(1,2) = 1.0d0
            FM%RightEigenMatrix(1,3) = 1.0d0
            FM%RightEigenMatrix(1,4) = 0.0d0
            FM%RightEigenMatrix(2,1) = RA%RoeAverage(2) - RA%RoeAverage(1)
            FM%RightEigenMatrix(2,2) = RA%RoeAverage(2)
            FM%RightEigenMatrix(2,3) = RA%RoeAverage(2) + RA%RoeAverage(1)
            FM%RightEigenMatrix(2,4) = 0.0d0
            FM%RightEigenMatrix(3,1) = RA%RoeAverage(3)
            FM%RightEigenMatrix(3,2) = RA%RoeAverage(3)
            FM%RightEigenMatrix(3,3) = RA%RoeAverage(3)
            FM%RightEigenMatrix(3,4) = 1.0d0
            FM%RightEigenMatrix(4,1) = RA%RoeAverage(4) - RA%RoeAverage(2)*RA%RoeAverage(1)
            FM%RightEigenMatrix(4,2) = 0.5d0 * (RA%RoeAverage(2)**2 + RA%RoeAverage(3)**2)
            FM%RightEigenMatrix(4,3) = RA%RoeAverage(4) + RA%RoeAverage(2)*RA%RoeAverage(1)
            FM%RightEigenMatrix(4,4) = RA%RoeAverage(3)

            SP4O%GEM%CalcAssistB2 = Gmin1/RA%RoeAverage(1)**2
            SP4O%GEM%CalcAssistB1 = 0.5d0*(RA%RoeAverage(2)**2 + RA%RoeAverage(3)**2) * SP4O%GEM%CalcAssistB2

            FM%LeftEigenMatrix(1,1) = 0.5d0*(SP4O%GEM%CalcAssistB1 + RA%RoeAverage(2)/RA%RoeAverage(1))
            FM%LeftEigenMatrix(1,2) = -0.5d0*(1.0d0/RA%RoeAverage(1) + SP4O%GEM%CalcAssistB2*RA%RoeAverage(2))
            FM%LeftEigenMatrix(1,3) = -0.5d0*(SP4O%GEM%CalcAssistB2*RA%RoeAverage(3))
            FM%LeftEigenMatrix(1,4) = 0.5d0*(SP4O%GEM%CalcAssistB2)
            FM%LeftEigenMatrix(2,1) = 1.0d0 - SP4O%GEM%CalcAssistB1
            FM%LeftEigenMatrix(2,2) = SP4O%GEM%CalcAssistB2 * RA%RoeAverage(2)
            FM%LeftEigenMatrix(2,3) = SP4O%GEM%CalcAssistB2 * RA%RoeAverage(3)
            FM%LeftEigenMatrix(2,4) = - SP4O%GEM%CalcAssistB2
            FM%LeftEigenMatrix(3,1) = 0.5d0*(SP4O%GEM%CalcAssistB1 - RA%RoeAverage(2)/RA%RoeAverage(1))
            FM%LeftEigenMatrix(3,2) = 0.5d0*(1.0d0/RA%RoeAverage(1) - SP4O%GEM%CalcAssistB2*RA%RoeAverage(2))
            FM%LeftEigenMatrix(3,3) = -0.5d0*(SP4O%GEM%CalcAssistB2*RA%RoeAverage(3))
            FM%LeftEigenMatrix(3,4) = 0.5d0*(SP4O%GEM%CalcAssistB2)
            FM%LeftEigenMatrix(4,1) = -RA%RoeAverage(3)
            FM%LeftEigenMatrix(4,2) = 0.0d0
            FM%LeftEigenMatrix(4,3) = 1.0d0
            FM%LeftEigenMatrix(4,4) = 0.0d0

            FM%AbsEigenValues = 0.0d0
            FM%AbsEigenValues(1,1) = dabs(RA%RoeAverage(2) - RA%RoeAverage(1))
            FM%AbsEigenValues(2,2) = dabs(RA%RoeAverage(2))
            FM%AbsEigenValues(3,3) = dabs(RA%RoeAverage(2) + RA%RoeAverage(1))
            FM%AbsEigenValues(4,4) = dabs(RA%RoeAverage(2))

        else if(RA%Direction == 2) then
            FM%RightEigenMatrix(1,1) = 1.0d0
            FM%RightEigenMatrix(1,2) = 1.0d0
            FM%RightEigenMatrix(1,3) = 1.0d0
            FM%RightEigenMatrix(1,4) = 0.0d0
            FM%RightEigenMatrix(2,1) = RA%RoeAverage(2)
            FM%RightEigenMatrix(2,2) = RA%RoeAverage(2)
            FM%RightEigenMatrix(2,3) = RA%RoeAverage(2)
            FM%RightEigenMatrix(2,4) = 1.0d0
            FM%RightEigenMatrix(3,1) = RA%RoeAverage(3) - RA%RoeAverage(1)
            FM%RightEigenMatrix(3,2) = RA%RoeAverage(3)
            FM%RightEigenMatrix(3,3) = RA%RoeAverage(3) + RA%RoeAverage(1)
            FM%RightEigenMatrix(3,4) = 0.0d0
            FM%RightEigenMatrix(4,1) = RA%RoeAverage(4) - RA%RoeAverage(3)*RA%RoeAverage(1)
            FM%RightEigenMatrix(4,2) = 0.5d0 * (RA%RoeAverage(2)**2 + RA%RoeAverage(3)**2)
            FM%RightEigenMatrix(4,3) = RA%RoeAverage(4) + RA%RoeAverage(3)*RA%RoeAverage(1)
            FM%RightEigenMatrix(4,4) = RA%RoeAverage(2)

            SP4O%GEM%CalcAssistB2 = Gmin1/RA%RoeAverage(1)**2
            SP4O%GEM%CalcAssistB1 = 0.5d0*(RA%RoeAverage(2)**2 + RA%RoeAverage(3)**2) * SP4O%GEM%CalcAssistB2

            FM%LeftEigenMatrix(1,1) = 0.5d0*(SP4O%GEM%CalcAssistB1 + RA%RoeAverage(3)/RA%RoeAverage(1))
            FM%LeftEigenMatrix(1,2) = -0.5d0*(SP4O%GEM%CalcAssistB2*RA%RoeAverage(2))
            FM%LeftEigenMatrix(1,3) = -0.5d0*(1.0d0/RA%RoeAverage(1) + SP4O%GEM%CalcAssistB2*RA%RoeAverage(3))
            FM%LeftEigenMatrix(1,4) = 0.5d0*(SP4O%GEM%CalcAssistB2)
            FM%LeftEigenMatrix(2,1) = 1.0d0 - SP4O%GEM%CalcAssistB1
            FM%LeftEigenMatrix(2,2) = SP4O%GEM%CalcAssistB2 * RA%RoeAverage(2)
            FM%LeftEigenMatrix(2,3) = SP4O%GEM%CalcAssistB2 * RA%RoeAverage(3)
            FM%LeftEigenMatrix(2,4) = - SP4O%GEM%CalcAssistB2
            FM%LeftEigenMatrix(3,1) = 0.5d0*(SP4O%GEM%CalcAssistB1 - RA%RoeAverage(3)/RA%RoeAverage(1))
            FM%LeftEigenMatrix(3,2) = -0.5d0*(SP4O%GEM%CalcAssistB2*RA%RoeAverage(2))
            FM%LeftEigenMatrix(3,3) = 0.5d0*(1.0d0/RA%RoeAverage(1) - SP4O%GEM%CalcAssistB2*RA%RoeAverage(3))
            FM%LeftEigenMatrix(3,4) = 0.5d0*(SP4O%GEM%CalcAssistB2)
            FM%LeftEigenMatrix(4,1) = -RA%RoeAverage(2)
            FM%LeftEigenMatrix(4,2) = 1.0d0
            FM%LeftEigenMatrix(4,3) = 0.0d0
            FM%LeftEigenMatrix(4,4) = 0.0d0

            FM%AbsEigenValues = 0.0d0
            FM%AbsEigenValues(1,1) = dabs(RA%RoeAverage(3) - RA%RoeAverage(1))
            FM%AbsEigenValues(2,2) = dabs(RA%RoeAverage(3))
            FM%AbsEigenValues(3,3) = dabs(RA%RoeAverage(3) + RA%RoeAverage(1))
            FM%AbsEigenValues(4,4) = dabs(RA%RoeAverage(3))
        end if
    end if
    return
end subroutine GetEigenMatrix
