!***********************************/
!	Name:FDS法によって2次元風上流束を求めるプログラム
!	Alias:UpwindFlux_Dim2
!	Description:FDS法を用いて2次元系における風上流束を求める
!	Type:CE
!	Input:Geometory,CellEdge,CellCenter
!	Output:CellEdge
!	Note:次元ごとにサブルーチンを分けているのは1つだけでも割と長いプログラムになってしまったため
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.12.27
!	Other:下の方でx方向とy方向を分ているのは法線ベクトルを計算していないため !12/27 UGetRoeAverageの引数からUGを削除
!***********************************/
    subroutine UUpwindFlux_Dim2(UG,UCE,UCC) !MUSCL経由の場合CEのみ，1次精度の場合CCを渡す
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, Epsilon => EntropyCorrection
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellEdge), intent(inout) :: UCE
    type(CellCenter), intent(in), optional :: UCC

    type(FDSMatrix) :: FM
    type(RoeAverage) :: RA
    type(PrivateVar4OMP) :: P4O
    iDim = UG%GM%Dimension
    allocate(RA%GeneralFlux(iDim+2,iDim,2), P4O%UF%DeltaQuantity(iDim+2),P4O%UF%RoeFlux(iDim+2))
    allocate(FM%RightEigenMatrix(iDim+2,iDim+2),FM%LeftEigenMatrix(iDim+2,iDim+2),FM%AbsEigenValues(iDim+2,iDim+2))
    allocate(RA%RoeAverage(iDim+2),RA%SideQuantity(iDim+2,2))

        allocate(P4O%GRA%SqrtDensity(2),P4O%GRA%Enthalpy(2))
        allocate(P4O%GRA%Velocity(UG%GM%Dimension,2))
        allocate(P4O%GGF%KinemeticEnergy(2),P4O%GGF%NormalVelocity(2),P4O%GGF%InverseDensity(2), &
            &   P4O%GGF%TotalEnergy(2),P4O%GGF%FluxJacobianMatrix(5,5,2))
        allocate(P4O%GGF%Velocity(3,2))

    !Roeスキームの流束を求める
!$omp parallel num_threads(CoreNumberOfCPU),shared(UG,UCE,UCC),firstprivate(iEdge,FM,RA,P4O,iLoop)
!$omp do
    do iEdge=1, UG%GI%Edges !すべての界面について
        if(present(UCC)) then
            iCell = UG%Line%Cell(iEdge,1,1) !表側セル
            iAdjacentCell = UG%Line%Cell(iEdge,2,1) !裏側セル
            RA%SideQuantity(:,2) = UCC%ConservedQuantity(:,iCell,1,1) !保存変数の場合 表側セルの保存量
            RA%SideQuantity(:,1) = UCC%ConservedQuantity(:,iAdjacentCell,1,1) !裏側セルの保存量
        else !基礎変数からの変換
            RA%SideQuantity(1,1) = UCE%RebuildQunatity(1,1,1,2,iEdge)
            RA%SideQuantity(2,1) = UCE%RebuildQunatity(1,1,1,2,iEdge)*UCE%RebuildQunatity(2,1,1,2,iEdge)
            RA%SideQuantity(3,1) = UCE%RebuildQunatity(1,1,1,2,iEdge)*UCE%RebuildQunatity(3,1,1,2,iEdge)
            RA%SideQuantity(4,1) = UCE%RebuildQunatity(1,1,1,2,iEdge)*UCE%RebuildQunatity(4,1,1,2,iEdge)
            RA%SideQuantity(5,1) = InverseGmin1*UCE%RebuildQunatity(5,1,1,2,iEdge) &
                & + 0.5d0*UCE%RebuildQunatity(1,1,1,2,iEdge) * sum(UCE%RebuildQunatity(2:4,1,1,2,iEdge)**2)

            RA%SideQuantity(1,2) = UCE%RebuildQunatity(1,1,1,1,iEdge)
            RA%SideQuantity(2,2) = UCE%RebuildQunatity(1,1,1,1,iEdge)*UCE%RebuildQunatity(2,1,1,1,iEdge)
            RA%SideQuantity(3,2) = UCE%RebuildQunatity(1,1,1,1,iEdge)*UCE%RebuildQunatity(3,1,1,1,iEdge)
            RA%SideQuantity(4,2) = UCE%RebuildQunatity(1,1,1,1,iEdge)*UCE%RebuildQunatity(4,1,1,1,iEdge)
            RA%SideQuantity(5,2) = InverseGmin1*UCE%RebuildQunatity(5,1,1,1,iEdge) &
                & + 0.5d0*UCE%RebuildQunatity(1,1,1,1,iEdge) * sum(UCE%RebuildQunatity(2:4,1,1,1,iEdge)**2)

        end if

        P4O%UF%DeltaQuantity(:) = RA%SideQuantity(:,2) - RA%SideQuantity(:,1) !面を挟んだ保存量の差
        P4O%Cell1 = UG%Line%Cell(iEdge,1,1) !表側セル
        P4O%Cell2 = UG%Line%Cell(iEdge,2,1) !裏側セル

        call UGetRoeAverage(RA,P4O) !Roe平均を求める

        call UGetGeneralFlux(iEdge,UG,RA,P4O)

        call UGetEigenMatrix(iEdge,UG,RA,FM,P4O)

        do iLoop =1,iDim+2
            if(FM%AbsEigenValues(iLoop,iLoop) < Epsilon) then
                FM%AbsEigenValues(iLoop,iLoop) = 0.5d0*(FM%AbsEigenValues(iLoop,iLoop)**2 + Epsilon**2)/Epsilon
            end if
        end do

        P4O%UF%RoeFlux = matmul(matmul(matmul(FM%RightEigenMatrix,FM%AbsEigenValues),FM%LeftEigenMatrix),P4O%UF%DeltaQuantity)

        UCE%NormalFluxDiff(:,1,1,1,iEdge) = &
        &   + 0.5d0*(RA%GeneralFlux(:,1,1)+RA%GeneralFlux(:,1,2) - P4O%UF%RoeFlux(:))

    end do
!$omp end do
!$omp end parallel

    return
end subroutine UUpwindFlux_Dim2
