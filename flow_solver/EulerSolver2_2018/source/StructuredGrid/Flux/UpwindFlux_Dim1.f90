!***********************************/
!	Name:FDS法による1次元風上流束を求めるプログラム
!	Alias:UpwindFlux_Dim1
!	Description:FDS法を用いて1次元系における風上流束を求める
!	Type:CE
!	Input:Geometory,CellEdge,CellCenter
!	Output:CellEdge
!	Note:次元ごとにサブルーチンを分けているのは1つだけでも割と長いプログラムになってしまったため
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:1次元と2次元に対応
!***********************************/
    subroutine UpwindFlux_Dim1(Geom,CE,CC) !MUSCL経由の場合CEのみ，1次精度の場合CCを渡す
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, Epsilon => EntropyCorrection
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellEdge), intent(inout) :: CE
    type(CellCenter), intent(in), optional :: CC

    type(FDSMatrix) :: FM
    type(RoeAverage) :: RA
    double precision, allocatable :: DeltaQuantity(:)
    double precision, allocatable :: RoeFlux(:)
    integer iTestMode
    iTestMode = 0
    allocate(RA%GeneralFlux(3,1,2), DeltaQuantity(3),RoeFlux(3))
    allocate(FM%RightEigenMatrix(3,3),FM%LeftEigenMatrix(3,3),FM%AbsEigenValues(3,3))
    allocate(RA%RoeAverage(3),RA%SideQuantity(3,2))

    !Roeスキームの流束を求める
    do iFaceX=0,Geom%CellNumber(1) !すべての界面についての総和
        if( present(CC)) then !保存変数からの流束計算
            RA%SideQuantity(:,1) = CC%ConservedQuantity(:,iCX,1,1) !負側の保存量
            RA%SideQuantity(:,2) = CC%ConservedQuantity(:,iCX+1,1,1) !正側の保存量
            RA%TypeOfQuantity = 1
        else  !基礎変数からの流束計算
            RA%SideQuantity(1,1) = CE%RebuildQunatity(1,iCX,1,1,2)
            RA%SideQuantity(2,1) = CE%RebuildQunatity(2,iCX,1,1,2)*CE%RebuildQunatity(1,iCX,1,1,2)
            RA%SideQuantity(3,1) = InverseGmin1*CE%RebuildQunatity(3,iCX,1,1,2) &
                &   + 0.5d0*CE%RebuildQunatity(1,iCX,1,1,2)*CE%RebuildQunatity(2,iCX,1,1,2)**2

            RA%SideQuantity(1,2) = CE%RebuildQunatity(1,iCX+1,1,1,1)
            RA%SideQuantity(2,2) = CE%RebuildQunatity(2,iCX+1,1,1,1)*CE%RebuildQunatity(1,iCX+1,1,1,1)
            RA%SideQuantity(3,2) = InverseGmin1*CE%RebuildQunatity(3,iCX+1,1,1,1)  &
                &   + 0.5d0*CE%RebuildQunatity(1,iCX+1,1,1,1)*CE%RebuildQunatity(2,iCX+1,1,1,1)**2
            RA%TypeOfQuantity = 1
        end if

        DeltaQuantity(:) = RA%SideQuantity(:,2) - RA%SideQuantity(:,1)
        call GetGeneralFlux(Geom,RA)

        !call GetRoeAverage(Geom,RA)

!左右の固有ベクトルから風上化流束を求める
        !call GetEigenMatrix(Geom,RA,FM)

!エントロピー補正
        do iLoop =1,3
        if(FM%AbsEigenValues(iLoop,iLoop) < Epsilon) then
            FM%AbsEigenValues(iLoop,iLoop) = 0.5d0*(FM%AbsEigenValues(iLoop,iLoop)**2 + Epsilon**2)/Epsilon
        end if
        end do
        !TestCode
        iTestMode = 0
        if(iTestMode == 0) then
        RoeFlux = matmul(matmul(matmul(FM%RightEigenMatrix,FM%AbsEigenValues),FM%LeftEigenMatrix),DeltaQuantity)
        else
            call RoeMultiply
        end if
        CE%NormalFluxDiff(:,iCX,1,1,1) = &
            &   0.5d0 * ((RA%GeneralFlux(:,1,1)+RA%GeneralFlux(:,1,2)) - RoeFlux(:))

    end do
    return
contains
    subroutine RoeMultiply
    implicit none
    double precision :: CalcAssistB2,CalcAssistB1
        CalcAssistB2 = Gmin1/RA%RoeAverage(1)**2
        CalcAssistB1 = 0.5d0*RA%RoeAverage(2)**2 * CalcAssistB2

        RoeFlux(1) = &
        & (0.5d0*(CalcAssistB1 + RA%RoeAverage(2)/RA%RoeAverage(1)) &
        &+(1.0d0 - CalcAssistB1) &
        &+0.5d0*(CalcAssistB1 - RA%RoeAverage(2)/RA%RoeAverage(1)))*DeltaQuantity(1) &
        &*dabs(RA%RoeAverage(2) - RA%RoeAverage(1))

        RoeFlux(2) = &
        & ((RA%RoeAverage(2) - RA%RoeAverage(1))&
        &*(-0.5d0*(1.0d0/RA%RoeAverage(1) + CalcAssistB2*RA%RoeAverage(2))) &
        &+RA%RoeAverage(2)*CalcAssistB2 * RA%RoeAverage(2) &
        &+(RA%RoeAverage(2) + RA%RoeAverage(1))&
        &*0.5d0*(1.0d0/RA%RoeAverage(1) - CalcAssistB2*RA%RoeAverage(2))) &
        &*DeltaQuantity(2)*dabs(RA%RoeAverage(2))

        RoeFlux(3) = &
        & ((RA%RoeAverage(3) - RA%RoeAverage(2)*RA%RoeAverage(1))*0.5d0*CalcAssistB2 &
        &+0.5d0 * RA%RoeAverage(2)**2 *(-CalcAssistB2) &
        &+(RA%RoeAverage(3) + RA%RoeAverage(2)*RA%RoeAverage(1))*0.5d0*CalcAssistB2) &
        &*DeltaQuantity(3)*dabs(RA%RoeAverage(2) + RA%RoeAverage(1))
    return
    end subroutine RoeMultiply

    end subroutine UpwindFlux_Dim1
