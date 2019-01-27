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
!	Update:2017.11.09
!	Other:下の方でx方向とy方向を分ているのは法線ベクトルを計算していないため
!***********************************/
    subroutine UpwindFlux_Dim2(Geom,CE,CC) !MUSCL経由の場合CEのみ，1次精度の場合CCを渡す
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, Epsilon => EntropyCorrection
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellEdge), intent(inout) :: CE
    type(CellCenter), intent(in), optional :: CC

    type(FDSMatrix) :: FM
    type(RoeAverage) :: RA
    type(PrivateVar4OMP) :: SP4O
    integer :: iSkip

    iDim = Geom%Dimension
    allocate(RA%GeneralFlux(iDim+2,1,2), SP4O%UF%DeltaQuantity(iDim+2),SP4O%UF%RoeFlux(iDim+2))
    allocate(FM%RightEigenMatrix(iDim+2,iDim+2),FM%LeftEigenMatrix(iDim+2,iDim+2),FM%AbsEigenValues(iDim+2,iDim+2))
    allocate(RA%RoeAverage(iDim+2),RA%SideQuantity(iDim+2,2))

        allocate(SP4O%GRA%SqrtDensity(2),SP4O%GRA%Enthalpy(2))
        allocate(SP4O%GRA%Velocity(Geom%Dimension,2))


    !Roeスキームの流束を求める

    RA%Direction = 1
!$omp parallel num_threads(CoreNumberOfCPU),shared(Geom,CE,CC),firstprivate(FM,RA,SP4O,iLoop,iSkip,iCenterX,iCenterY,iFaceX,iFaceY)
!$omp do
    do iCenterY=1,Geom%CellNumber(2) !y方向はセル中心位置として
        do iFaceX=0,Geom%CellNumber(1) !すべてのX界面についての総和を取る
!___________________!call directionX!________________________________
            if(present(CC)) then !保存変数からの流束計算
                RA%SideQuantity(:,1) = CC%ConservedQuantity(:,iFaceX,iCenterY,1) !負側の保存量
                RA%SideQuantity(:,2) = CC%ConservedQuantity(:,iFaceX+1,iCenterY,1) !正側の保存量
                RA%TypeOfQuantity = 1
            else  !基礎変数からの流束計算

                RA%SideQuantity(1,1) = CE%RebuildQunatity(1,iFaceX,iCenterY,1,2) !負側の保存量
                RA%SideQuantity(2,1) = CE%RebuildQunatity(2,iFaceX,iCenterY,1,2) * CE%RebuildQunatity(1,iFaceX,iCenterY,1,2)
                RA%SideQuantity(3,1) = CE%RebuildQunatity(3,iFaceX,iCenterY,1,2) * CE%RebuildQunatity(1,iFaceX,iCenterY,1,2)
                RA%SideQuantity(4,1) = InverseGmin1*CE%RebuildQunatity(4,iFaceX,iCenterY,1,2) &
                &   + 0.5d0*CE%RebuildQunatity(1,iFaceX,iCenterY,1,2) &
                &   * (CE%RebuildQunatity(2,iFaceX,iCenterY,1,2)**2 + CE%RebuildQunatity(3,iFaceX,iCenterY,1,2)**2)

                RA%SideQuantity(1,2) = CE%RebuildQunatity(1,iFaceX+1,iCenterY,1,1)
                RA%SideQuantity(2,2) = CE%RebuildQunatity(2,iFaceX+1,iCenterY,1,1)*CE%RebuildQunatity(1,iFaceX+1,iCenterY,1,1)
                RA%SideQuantity(3,2) = CE%RebuildQunatity(3,iFaceX+1,iCenterY,1,1)*CE%RebuildQunatity(1,iFaceX+1,iCenterY,1,1)

                RA%SideQuantity(4,2) = InverseGmin1*CE%RebuildQunatity(4,iFaceX+1,iCenterY,1,1) &
                &   + 0.5d0*CE%RebuildQunatity(1,iFaceX+1,iCenterY,1,1) &
                &   * (CE%RebuildQunatity(2,iFaceX+1,iCenterY,1,1)**2 + CE%RebuildQunatity(3,iFaceX+1,iCenterY,1,1)**2)
                RA%TypeOfQuantity = 1
            end if
!___________________!call common2DFDS!________________________________

        call GetGeneralFlux(Geom,RA)
        SP4O%UF%DeltaQuantity(:) = RA%SideQuantity(:,2) - RA%SideQuantity(:,1)
        iSkip = 0
        do iVariable=1,4
            if(SP4O%UF%DeltaQuantity(iVariable) == 0.0d0) then
                iSkip = iSkip + 1
            end if
        end do
        if(iSkip /= 4) then
            call GetRoeAverage(Geom,RA,SP4O)
            !左右の固有ベクトルから風上化流束を求める
            call GetEigenMatrix(Geom,RA,FM,SP4O)
            !エントロピー補正
            do iLoop =1,iDim+2
                if(FM%AbsEigenValues(iLoop,iLoop) < Epsilon) then
                    FM%AbsEigenValues(iLoop,iLoop) = 0.5d0*(FM%AbsEigenValues(iLoop,iLoop)**2 + Epsilon**2)/Epsilon
                end if
            end do
            SP4O%UF%RoeFlux = matmul(matmul(matmul(FM%RightEigenMatrix,FM%AbsEigenValues),FM%LeftEigenMatrix),SP4O%UF%DeltaQuantity)
        else
            SP4O%UF%RoeFlux = 0.0d0
        end if
!___________________end common2DFDS________________________________

        CE%NormalFluxDiff(:,iFaceX,iCenterY,1,1) = &
        &   0.5d0 * ((RA%GeneralFlux(:,1,1)+RA%GeneralFlux(:,1,2)) - SP4O%UF%RoeFlux(:))

!___________________end directionX________________________________
        end do
    end do
!$omp end do
!$omp barrier
    RA%Direction = 2
!$omp do
    do iCenterX=1,Geom%CellNumber(1) !x方向はセル中心位置として
        do iFaceY=0,Geom%CellNumber(2) !すべてのy界面についての総和を取る
!___________________call directionY________________________________
        if(present(CC)) then !保存変数からの流束計算
            RA%SideQuantity(:,1) = CC%ConservedQuantity(:,iCenterX,iFaceY,1) !負側の保存量
            RA%SideQuantity(:,2) = CC%ConservedQuantity(:,iCenterX,iFaceY+1,1) !正側の保存量
            RA%TypeOfQuantity = 1
        else  !基礎変数からの流束計算
            RA%SideQuantity(1,1) = CE%RebuildQunatity(1,iCenterX,iFaceY,1,4) !負側の保存量
            RA%SideQuantity(2,1) = CE%RebuildQunatity(2,iCenterX,iFaceY,1,4) * CE%RebuildQunatity(1,iCenterX,iFaceY,1,4)
            RA%SideQuantity(3,1) = CE%RebuildQunatity(3,iCenterX,iFaceY,1,4) * CE%RebuildQunatity(1,iCenterX,iFaceY,1,4)
            RA%SideQuantity(4,1) = InverseGmin1*CE%RebuildQunatity(4,iCenterX,iFaceY,1,4) &
            &   + 0.5d0*CE%RebuildQunatity(1,iCenterX,iFaceY,1,4) &
            &   * (CE%RebuildQunatity(2,iCenterX,iFaceY,1,4)**2 + CE%RebuildQunatity(3,iCenterX,iFaceY,1,4)**2)
                RA%SideQuantity(1,2) = CE%RebuildQunatity(1,iCenterX,iFaceY+1,1,3)

            RA%SideQuantity(2,2) = CE%RebuildQunatity(2,iCenterX,iFaceY+1,1,3)*CE%RebuildQunatity(1,iCenterX,iFaceY+1,1,3)
            RA%SideQuantity(3,2) = CE%RebuildQunatity(3,iCenterX,iFaceY+1,1,3)*CE%RebuildQunatity(1,iCenterX,iFaceY+1,1,3)
            RA%SideQuantity(4,2) = InverseGmin1*CE%RebuildQunatity(4,iCenterX,iFaceY+1,1,3) &
            &   + 0.5d0*CE%RebuildQunatity(1,iCenterX,iFaceY+1,1,3) &
            &   * (CE%RebuildQunatity(2,iCenterX,iFaceY+1,1,3)**2 + CE%RebuildQunatity(3,iCenterX,iFaceY+1,1,3)**2)
            RA%TypeOfQuantity = 1
        end if

!___________________call common2DFDS________________________________
        call GetGeneralFlux(Geom,RA)
        SP4O%UF%DeltaQuantity(:) = RA%SideQuantity(:,2) - RA%SideQuantity(:,1)
        iSkip = 0
        do iVariable=1,4
            if(SP4O%UF%DeltaQuantity(iVariable) == 0.0d0) then
                iSkip = iSkip + 1
            end if
        end do
        if(iSkip /= 4) then
            call GetRoeAverage(Geom,RA,SP4O)
            !左右の固有ベクトルから風上化流束を求める
            call GetEigenMatrix(Geom,RA,FM,SP4O)
            !エントロピー補正
            do iLoop =1,iDim+2
                if(FM%AbsEigenValues(iLoop,iLoop) < Epsilon) then
                    FM%AbsEigenValues(iLoop,iLoop) = 0.5d0*(FM%AbsEigenValues(iLoop,iLoop)**2 + Epsilon**2)/Epsilon
                end if
            end do
            SP4O%UF%RoeFlux = matmul(matmul(matmul(FM%RightEigenMatrix,FM%AbsEigenValues),FM%LeftEigenMatrix),SP4O%UF%DeltaQuantity)
        else
            SP4O%UF%RoeFlux = 0.0d0
        end if
!___________________end common2DFDS________________________________

        CE%NormalFluxDiff(:,iCenterX,iFaceY,1,2) = &
        &   0.5d0 * ((RA%GeneralFlux(:,1,1)+RA%GeneralFlux(:,1,2)) - SP4O%UF%RoeFlux(:))

!___________________end directionY________________________________
        end do
    end do
!$omp end do
!$omp end parallel
    return

contains

    subroutine directionX !X方向界面流束
    implicit none
        if( present(CC)) then !保存変数からの流束計算
            RA%SideQuantity(:,1) = CC%ConservedQuantity(:,iFaceX,iCenterY,1) !負側の保存量
            RA%SideQuantity(:,2) = CC%ConservedQuantity(:,iFaceX+1,iCenterY,1) !正側の保存量
            RA%TypeOfQuantity = 1
        else  !基礎変数からの流束計算

            RA%SideQuantity(1,1) = CE%RebuildQunatity(1,iFaceX,iCenterY,1,2) !負側の保存量
            RA%SideQuantity(2,1) = CE%RebuildQunatity(2,iFaceX,iCenterY,1,2) * CE%RebuildQunatity(1,iFaceX,iCenterY,1,2)
            RA%SideQuantity(3,1) = CE%RebuildQunatity(3,iFaceX,iCenterY,1,2) * CE%RebuildQunatity(1,iFaceX,iCenterY,1,2)
            RA%SideQuantity(4,1) = InverseGmin1*CE%RebuildQunatity(4,iFaceX,iCenterY,1,2) &
            &   + 0.5d0*CE%RebuildQunatity(1,iFaceX,iCenterY,1,2) &
            &   * (CE%RebuildQunatity(2,iFaceX,iCenterY,1,2)**2 + CE%RebuildQunatity(3,iFaceX,iCenterY,1,2)**2)

            RA%SideQuantity(1,2) = CE%RebuildQunatity(1,iFaceX+1,iCenterY,1,1)

            RA%SideQuantity(2,2) = CE%RebuildQunatity(2,iFaceX+1,iCenterY,1,1)*CE%RebuildQunatity(1,iFaceX+1,iCenterY,1,1)

            RA%SideQuantity(3,2) = CE%RebuildQunatity(3,iFaceX+1,iCenterY,1,1)*CE%RebuildQunatity(1,iFaceX+1,iCenterY,1,1)

            RA%SideQuantity(4,2) = InverseGmin1*CE%RebuildQunatity(4,iFaceX+1,iCenterY,1,1) &
            &   + 0.5d0*CE%RebuildQunatity(1,iFaceX+1,iCenterY,1,1) &
            &   * (CE%RebuildQunatity(2,iFaceX+1,iCenterY,1,1)**2 + CE%RebuildQunatity(3,iFaceX+1,iCenterY,1,1)**2)

            RA%TypeOfQuantity = 1
        end if

        call common2DFDS

        CE%NormalFluxDiff(:,iFaceX,iCenterY,1,1) = &
        &   0.5d0 * ((RA%GeneralFlux(:,1,1)+RA%GeneralFlux(:,1,2)) - SP4O%UF%RoeFlux(:))

    return
    end subroutine directionX
!________________________________________________________

    subroutine directionY !Y方向界面流束
    implicit none

        if(present(CC)) then !保存変数からの流束計算
            RA%SideQuantity(:,1) = CC%ConservedQuantity(:,iCenterX,iFaceY,1) !負側の保存量
            RA%SideQuantity(:,2) = CC%ConservedQuantity(:,iCenterX,iFaceY+1,1) !正側の保存量
            RA%TypeOfQuantity = 1
        else  !基礎変数からの流束計算
            RA%SideQuantity(1,1) = CE%RebuildQunatity(1,iCenterX,iFaceY,1,4) !負側の保存量
            RA%SideQuantity(2,1) = CE%RebuildQunatity(2,iCenterX,iFaceY,1,4) * CE%RebuildQunatity(1,iCenterX,iFaceY,1,4)
            RA%SideQuantity(3,1) = CE%RebuildQunatity(3,iCenterX,iFaceY,1,4) * CE%RebuildQunatity(1,iCenterX,iFaceY,1,4)
            RA%SideQuantity(4,1) = InverseGmin1*CE%RebuildQunatity(4,iCenterX,iFaceY,1,4) &
            &   + 0.5d0*CE%RebuildQunatity(1,iCenterX,iFaceY,1,4) &
            &   * (CE%RebuildQunatity(2,iCenterX,iFaceY,1,4)**2 + CE%RebuildQunatity(3,iCenterX,iFaceY,1,4)**2)
                RA%SideQuantity(1,2) = CE%RebuildQunatity(1,iCenterX,iFaceY+1,1,3)

            RA%SideQuantity(2,2) = CE%RebuildQunatity(2,iCenterX,iFaceY+1,1,3)*CE%RebuildQunatity(1,iCenterX,iFaceY+1,1,3)
            RA%SideQuantity(3,2) = CE%RebuildQunatity(3,iCenterX,iFaceY+1,1,3)*CE%RebuildQunatity(1,iCenterX,iFaceY+1,1,3)
            RA%SideQuantity(4,2) = InverseGmin1*CE%RebuildQunatity(4,iCenterX,iFaceY+1,1,3) &
            &   + 0.5d0*CE%RebuildQunatity(1,iCenterX,iFaceY+1,1,3) &
            &   * (CE%RebuildQunatity(2,iCenterX,iFaceY+1,1,3)**2 + CE%RebuildQunatity(3,iCenterX,iFaceY+1,1,3)**2)
            RA%TypeOfQuantity = 1
        end if

        call common2DFDS
        CE%NormalFluxDiff(:,iCenterX,iFaceY,1,2) = &
        &   0.5d0 * ((RA%GeneralFlux(:,1,1)+RA%GeneralFlux(:,1,2)) - SP4O%UF%RoeFlux(:))

    return
    end subroutine directionY
!________________________________________________________

    subroutine common2DFDS
    implicit none

        call GetGeneralFlux(Geom,RA)
        SP4O%UF%DeltaQuantity(:) = RA%SideQuantity(:,2) - RA%SideQuantity(:,1)
        iSkip = 0
        do iVariable=1,4
            if(SP4O%UF%DeltaQuantity(iVariable) == 0.0d0) then
                iSkip = iSkip + 1
            end if
        end do
        if(iSkip /= 4) then
            call GetRoeAverage(Geom,RA,SP4O)
            !左右の固有ベクトルから風上化流束を求める
            call GetEigenMatrix(Geom,RA,FM,SP4O)
            !エントロピー補正
            do iLoop =1,iDim+2
                if(FM%AbsEigenValues(iLoop,iLoop) < Epsilon) then
                    FM%AbsEigenValues(iLoop,iLoop) = 0.5d0*(FM%AbsEigenValues(iLoop,iLoop)**2 + Epsilon**2)/Epsilon
                end if
            end do
            SP4O%UF%RoeFlux = matmul(matmul(matmul(FM%RightEigenMatrix,FM%AbsEigenValues),FM%LeftEigenMatrix),SP4O%UF%DeltaQuantity)
        else
            SP4O%UF%RoeFlux = 0.0d0
        end if

        return
    end subroutine common2DFDS

end subroutine UpwindFlux_Dim2
