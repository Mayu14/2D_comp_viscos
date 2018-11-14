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
subroutine UGetStrainRateTensor(UG,UCC)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    use ConstantVar_Mod, ci => ImaginaryNumber
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    integer :: iAdjC1, iAdjC2, iAdjC3
    complex(kind(0d0)) :: alpha1, alpha2, alpha3, alpha0, tmpDerivetive
    double precision :: totalLength
    integer :: iVelDir


    iDim = UG%GM%Dimension
    do iCell = 1, UG%GI%RealCells
        call GetWeightOfCenterDifferenceScheme(iCell, iAdjC1, iAdjC2, iAdjC3, alpha1, alpha2, alpha3, totalLength)
        alpha0 = alpha1 + alpha2 + alpha3
        do iVelDir = 1, iDim
            tmpDerivetive = (alpha1 * UCC%PrimitiveVariable(1 + iVelDir, iAdjC1, 1, 1) &
                         & + alpha2 * UCC%PrimitiveVariable(1 + iVelDir, iAdjC2, 1, 1) &
                         & + alpha3 * UCC%PrimitiveVariable(1 + iVelDir, iAdjC3, 1, 1) &
                         & - alpha0 * UCC%PrimitiveVariable(1 + iVelDir, iCell, 1, 1)) / totalLength

            UCC%StrainRateTensor(iVelDir, 1, iCell, 1, 1) = real(tmpDerivetive)
            UCC%StrainRateTensor(iVelDir, 2, iCell, 1, 1) = aimag(tmpDerivetive)
        end do
    end do


return
contains
    subroutine GetWeightOfCenterDifferenceScheme(iCell, iAdjC1, iAdjC2, iAdjC3, alpha1, alpha2, alpha3, totalLength)
        implicit none
        integer, intent(in) :: iCell
        integer, intent(out) :: iAdjC1, iAdjC2, iAdjC3  ! iAdjCN: N番の局所面番号でiCellと隣接するセルの番号
        complex(kind(0d0)), intent(out) :: alpha1, alpha2, alpha3   ! alphaN: N番の局所面番号でiCellと隣接するセルの重み
        double precision, intent(out) :: totalLength
        integer :: iAdjE1, iAdjE2, iAdjE3   ! iAdjEN: iCellのN番目の局所面番号を大域面番号にしたもの
        integer :: iDirectionW1, iDirectionW2, iDirectionW3 ! 面基準での裏表(1 or 2)
        double precision :: DirectionN1, DirectionN2, DirectionN3   ! 法線ベクトルの向き(-1 or 1)
        double precision :: length1, length2, length3, theta1, theta2, theta3, determinantA ! lengthN: iCellと局所面N番の隣接セルの中心間距離，thetaN: 法線ベクトルとx軸のなす角度，重みを計算する式の行列式


        iAdjC1 = UG%Tri%Cell(iCell, 1)
        iAdjC2 = UG%Tri%Cell(iCell, 2)
        iAdjC3 = UG%Tri%Cell(iCell, 3)

        iAdjE1 = UG%Tri%Edge(iCell, 1)
        iAdjE2 = UG%Tri%Edge(iCell, 2)
        iAdjE3 = UG%Tri%Edge(iCell, 3)

        call CheckBackOrFront(iCell, iAdjC1, iDirectionW1, DirectionN1)
        call CheckBackOrFront(iCell, iAdjC2, iDirectionW2, DirectionN2)
        call CheckBackOrFront(iCell, iAdjC3, iDirectionW3, DirectionN3)

        length1 = AbsVector(UG%GM%Width(iCell, 1, :)) + AbsVector(UG%GM%Width(iAdjC1, UG%Line%Cell(iAdjE1, iDirectionW1, 2), :))
        length2 = AbsVector(UG%GM%Width(iCell, 2, :)) + AbsVector(UG%GM%Width(iAdjC2, UG%Line%Cell(iAdjE1, iDirectionW1, 2), :))
        length3 = AbsVector(UG%GM%Width(iCell, 3, :)) + AbsVector(UG%GM%Width(iAdjC3, UG%Line%Cell(iAdjE1, iDirectionW1, 2), :))

        theta1 = atan2(DirectionN1*UG%GM%Normal(iAdjE1, 2), DirectionN1*UG%GM%Normal(iAdjE1, 1))
        theta2 = atan2(DirectionN2*UG%GM%Normal(iAdjE2, 2), DirectionN2*UG%GM%Normal(iAdjE2, 1))
        theta3 = atan2(DirectionN3*UG%GM%Normal(iAdjE3, 2), DirectionN3*UG%GM%Normal(iAdjE3, 1))

        determinantA = GetDeterminantOfA(length1, length2, length3)

        alpha1 = (length2 * length3) ** 2 * (length3 - length2) / determinantA * (exp(ci * (theta2 + theta3)) + exp(ci * (theta2 + theta3 + 0.5d0 * dPi)))
        alpha2 = (length3 * length1) ** 2 * (length1 - length3) / determinantA * (exp(ci * (theta3 + theta1)) + exp(ci * (theta3 + theta1 + 0.5d0 * dPi)))
        alpha3 = (length1 * length2) ** 2 * (length2 - length1) / determinantA * (exp(ci * (theta1 + theta2)) + exp(ci * (theta1 + theta2 + 0.5d0 * dPi)))

        totalLength = length1 + length2 + length3
        return
    end subroutine

    subroutine CheckBackOrFront(iMainCell, iAdjCell, iBackOrFront12, iBackOrFrontm11)
        implicit none
        integer, intent(in) :: iMainCell, iAdjCell
        integer, intent(out) :: iBackOrFront12
        double precision, intent(out) :: iBackOrFrontm11

        if(iMainCell > iAdjCell) then ! Main側の数字が大きい = Main側が表
            iBackOrFront12 = 1  ! GM%Width用
            iBackOrFrontm11 = -1.0d0    ! GM%Normalの向き調整用
        else    ! Main側の数字が小さい = Main側が裏
            iBackOrFront12 = 2
            iBackOrFrontm11 = 1.0d0
        end if
        return
    end subroutine CheckBackOrFront

    function GetDeterminantOfA(len1, len2, len3) result(detA)
        implicit none
        double precision, intent(in) :: len1, len2, len3
        double precision :: detA

        detA = (((len1) * (len2 ** 2) * (len3 ** 3)) + ((len2 ** 3) * (len1 ** 2) * (len3)) + ((len1 ** 3) * (len2) * (len3 ** 2))) &
            & -(((len3) * (len2 ** 2) * (len1 ** 3)) + ((len3 ** 3) * (len2 ** 2) * (len1)) + ((len3 ** 3) * (len2) * (len1 ** 2)))

        return
    end function GetDeterminantOfA

end subroutine UGetStrainRateTensor
