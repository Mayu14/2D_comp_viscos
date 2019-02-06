!***********************************/
!	Name:ひずみ速度テンソルを計算するためのプログラム
!	Alias:UGetStrainRateTensor
!	Description:
!	Type:CellEdge
!	Input:Configulation,Geometory,CellCenter,CellEdge
!	Output:CE%NormalFluxDiff
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:1次元と2次元,1次精度と2次精度のみに対応
!***********************************/
subroutine UGetStrainRateTensor_mk2(UConf, UG, UCC, UCE)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    use ConstantVar_Mod, ci => ImaginaryNumber
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(CellEdge), intent(in) :: UCE
    integer :: iAdjC1, iAdjC2, iAdjC3
    double precision :: MyLength, AdjLength
    double precision :: totalLength
    integer :: iVelDir, iDirectionW !面基準でのの向き(1 or 2)
    double precision :: DirectionN  !法線ベクトルの向き(Normalにこれを掛けると外向き法線になる)
    double precision, allocatable :: weight_ave_value(:)    !1~3:velocty, 4:Temparature; 距離の加重平均による速度と温度
    double precision :: MyTemparature, AdjTemparature
    double precision, allocatable :: dUdi(:), dVdi(:), dTdi(:)    ! 速度・温度の微分(1:x微分,2:y微分)


    iDim = UG%GM%Dimension
    allocate(weight_ave_value(4))
    allocate(dUdi(2), dVdi(2), dTdi(2))

    do iCell = 1, UG%GI%RealCells
        dUdi = 0.0d0
        dVdi = 0.0d0
        dTdi = 0.0d0
        do iLocalEdge = 1, 3
            iAdjacentCell = UG%Tri%Cell(iCell, iLocalEdge)
            iEdge = UG%Tri%Edge(iCell, iLocalEdge)

            call CheckBackOrFront(iCell, iAdjacentCell, iDirectionW, DirectionN)

            MyLength = AbsVector(UG%GM%Width(iCell, iLocalEdge, :))
            AdjLength = AbsVector(UG%GM%Width(iAdjacentCell, UG%Line%Cell(iEdge, iDirectionW, 2), :))

            !if(UConf%UseMUSCL == 0) then
                weight_ave_value(1:3) = (MyLength * UCC%PrimitiveVariable(2:4,iCell,1,1) &
                                        & + AdjLength * UCC%PrimitiveVariable(2:4,iAdjacentCell,1,1)) &
                                        & / (MyLength + AdjLength)

                MyTemparature = gamma * Gmin1 * UCC%PrimitiveVariable(5,iCell,1,1)/UCC%PrimitiveVariable(1,iCell,1,1)
                AdjTemparature = gamma * Gmin1 * UCC%PrimitiveVariable(5,iAdjacentCell,1,1)/UCC%PrimitiveVariable(1,iAdjacentCell,1,1)
            !else
            !end if
            weight_ave_value(4) = (MyLength * MyTemparature + AdjLength * AdjTemparature) / (MyLength + AdjLength)



            dUdi(1:2) = dUdi(1:2) + weight_ave_value(1) * (DirectionN * UG%GM%Normal(iEdge,1:2)) * UG%GM%Area(iEdge)
            dVdi(1:2) = dVdi(1:2) + weight_ave_value(2) * (DirectionN * UG%GM%Normal(iEdge,1:2)) * UG%GM%Area(iEdge)
            dTdi(1:2) = dTdi(1:2) + weight_ave_value(4) * (DirectionN * UG%GM%Normal(iEdge,1:2)) * UG%GM%Area(iEdge)
        end do
        UCC%StrainRateTensor(1,1:2,iCell,1,1) = dUdi(1:2)
        UCC%StrainRateTensor(2,1:2,iCell,1,1) = dUdi(1:2)
        UCC%TemparatureGrad(1:2,iCell,1,1) = dTdi(1:2)

    end do

    UCC%AbsoluteVortisity(iCell, 1, 1) = abs(UCC%StrainRateTensor(1, 2, iCell, 1, 1) &
                                         & - UCC%StrainRateTensor(2, 1, iCell, 1, 1))

    end do

    return
contains

    function GetDeterminantOfA(len1, len2, len3) result(detA)
        implicit none
        double precision, intent(in) :: len1, len2, len3
        double precision :: detA

        detA = (((len1) * (len2 ** 2) * (len3 ** 3)) + ((len2 ** 3) * (len1 ** 2) * (len3)) + ((len1 ** 3) * (len2) * (len3 ** 2))) &
            & -(((len3) * (len2 ** 2) * (len1 ** 3)) + ((len3 ** 3) * (len2 ** 2) * (len1)) + ((len3 ** 3) * (len2) * (len1 ** 2)))

        return
    end function GetDeterminantOfA

end subroutine UGetStrainRateTensor_mk2
