subroutine PrimaryAccuracyFlux(Geom,CC,CE)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, Epsilon => EntropyCorrection
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter),intent(inout) :: CC
    type(CellEdge),intent(inout) :: CE

    call Conserve2Primitive(Geom,CC)

    nullify(iVar,iFX,iFY,iFZ,iCX,iCY,iCZ)
    iVar => iVariable
    iFX => iFaceX
    iFY => iFaceY
    iFZ => iFaceZ
    iCX => iFaceX
    iCY => iFaceY
    iCZ => iFaceZ

    if(Geom%Dimension == 1) call FluxDim1_UpWind
    if(Geom%Dimension == 2) call FluxDim2
    !if(Geom%Dimension == 3) call FluxDim3
    nullify(iVar,iFX,iFY,iFZ,iCX,iCY,iCZ)
    return

contains
    subroutine FluxDim1_Upwind
    implicit none
    double precision, allocatable :: GeneralFlux(:,:), DeltaQuantity(:)
    double precision, allocatable :: RoeFlux(:)
    type(FDSMatrix) :: FM

    type(RoeAverage) :: RA
    double precision :: CalcAssistB1,CalcAssistB2
    allocate(GeneralFlux(3,2), DeltaQuantity(3),RoeFlux(3))
    allocate(FM%RightEigenMatrix(3,3),FM%LeftEigenMatrix(3,3),FM%AbsEigenValues(3,3))
    allocate(RA%RoeAverage(3),RA%SideQuantity(3,2))
    !Roeスキームの流束を求める
    do iFaceX=0,Geom%CellNumber(1) !すべての界面についての総和
        RA%SideQuantity(1,1) = CC%ConservedQuantity(1,iCX,1,1) !負側の保存量
        RA%SideQuantity(2,1) = CC%ConservedQuantity(2,iCX,1,1)
        RA%SideQuantity(3,1) = CC%ConservedQuantity(3,iCX,1,1)

        RA%SideQuantity(1,2) = CC%ConservedQuantity(1,iCX+1,1,1) !正側の保存量
        RA%SideQuantity(2,2) = CC%ConservedQuantity(2,iCX+1,1,1)
        RA%SideQuantity(3,2) = CC%ConservedQuantity(3,iCX+1,1,1)

        DeltaQuantity(:) = RA%SideQuantity(:,2) - RA%SideQuantity(:,1)

        GeneralFlux(1,:) = RA%SideQuantity(2,:)

        GeneralFlux(2,:) = RA%SideQuantity(2,:)**2/RA%SideQuantity(1,:) &
            &   + Gmin1*(RA%SideQuantity(3,:)-0.5d0*RA%SideQuantity(2,:)**2/RA%SideQuantity(1,:))

        GeneralFlux(3,:) = RA%SideQuantity(2,:)/RA%SideQuantity(1,:) &
            &   * (Gamma*RA%SideQuantity(3,:)-0.5d0*Gmin1*RA%SideQuantity(2,:)**2/RA%SideQuantity(1,:))

        RA%TypeOfQuantity = 1
        call GetRoeAverage(Geom,RA)

!左右の固有ベクトルから風上化流束を求める
!        call GetRightEigenMatrix(Geom,RoeAverage,RightEigenMatrix)


        FM%RightEigenMatrix(1,1) = 1.0d0
        FM%RightEigenMatrix(1,2) = 1.0d0
        FM%RightEigenMatrix(1,3) = 1.0d0
        FM%RightEigenMatrix(2,1) = RA%RoeAverage(2) - RA%RoeAverage(1)
        FM%RightEigenMatrix(2,2) = RA%RoeAverage(2)
        FM%RightEigenMatrix(2,3) = RA%RoeAverage(2) + RA%RoeAverage(1)
        FM%RightEigenMatrix(3,1) = RA%RoeAverage(3) - RA%RoeAverage(2)*RA%RoeAverage(1)
        FM%RightEigenMatrix(3,2) = 0.5d0 * RA%RoeAverage(2)**2
        FM%RightEigenMatrix(3,3) = RA%RoeAverage(3) + RA%RoeAverage(2)*RA%RoeAverage(1)


!       call GetLeftEigenMatrix(Geom,RoeAverage,FM%LeftEigenMatrix)
        CalcAssistB2 = Gmin1/RA%RoeAverage(1)**2
        CalcAssistB1 = 0.5d0*RA%RoeAverage(2)**2 * CalcAssistB2

        FM%LeftEigenMatrix(1,1) = 0.5d0*(CalcAssistB1 + RA%RoeAverage(2)/RA%RoeAverage(1))
        FM%LeftEigenMatrix(1,2) = -0.5d0*(1.0d0/RA%RoeAverage(1) + CalcAssistB2*RA%RoeAverage(2))
        FM%LeftEigenMatrix(1,3) = 0.5d0*(CalcAssistB2)
        FM%LeftEigenMatrix(2,1) = 1.0d0 - CalcAssistB1
        FM%LeftEigenMatrix(2,2) = CalcAssistB2 * RA%RoeAverage(2)
        FM%LeftEigenMatrix(2,3) = -CalcAssistB2
        FM%LeftEigenMatrix(3,1) = 0.5d0*(CalcAssistB1 - RA%RoeAverage(2)/RA%RoeAverage(1))
        FM%LeftEigenMatrix(3,2) = 0.5d0*(1.0d0/RA%RoeAverage(1) - CalcAssistB2*RA%RoeAverage(2))
        FM%LeftEigenMatrix(3,3) = 0.5d0*(CalcAssistB2)

        FM%AbsEigenValues = 0.0d0
        FM%AbsEigenValues(1,1) = dabs(RA%RoeAverage(2) - RA%RoeAverage(1))
        FM%AbsEigenValues(2,2) = dabs(RA%RoeAverage(2))
        FM%AbsEigenValues(3,3) = dabs(RA%RoeAverage(2) + RA%RoeAverage(1))
        do iLoop =1,3
        if(FM%AbsEigenValues(iLoop,iLoop) < Epsilon) then
            FM%AbsEigenValues(iLoop,iLoop) = 0.5d0*(FM%AbsEigenValues(iLoop,iLoop)**2 + Epsilon**2)/Epsilon
        end if
        end do

        RoeFlux = matmul(matmul(matmul(FM%RightEigenMatrix,FM%AbsEigenValues),FM%LeftEigenMatrix),DeltaQuantity)
        !write(6,*) RoeFlux
        !write(6,*) FM%RightEigenMatrix
        !write(6,*) FM%LeftEigenMatrix
 !       write(6,*) DeltaQuantity
        CE%NormalFluxDiff(:,iCX,1,1,1) = &
            &   0.5d0 * ((GeneralFlux(:,1)+GeneralFlux(:,2)) - RoeFlux(:))


    end do
    return
    end subroutine FluxDim1_Upwind





    subroutine FluxDim1_UpWind2!没
    implicit none
    double precision, allocatable :: ConservedQuantity(:,:), GeneralFlux(:,:), DeltaFlux(:),DeltaQuantity(:)
    allocate(ConservedQuantity(3,2), GeneralFlux(3,2), DeltaFlux(3),DeltaQuantity(3))
    do iFaceX=0,Geom%CellNumber(1)
        ConservedQuantity(1,1)=CC%ConservedQuantity(1,iCX,1,1)
        ConservedQuantity(2,1)=CC%ConservedQuantity(2,iCX,1,1)
        ConservedQuantity(3,1)=CC%ConservedQuantity(3,iCX,1,1)

        ConservedQuantity(1,2)=CC%ConservedQuantity(1,iCX+1,1,1)
        ConservedQuantity(2,2)=CC%ConservedQuantity(2,iCX+1,1,1)
        ConservedQuantity(3,2)=CC%ConservedQuantity(3,iCX+1,1,1)

        GeneralFlux(1,:) = ConservedQuantity(2,:)

        GeneralFlux(2,:) = ConservedQuantity(2,:)**2/ConservedQuantity(1,:) &
            &   + Gmin1*(ConservedQuantity(3,:)-0.5d0*ConservedQuantity(2,:)**2/ConservedQuantity(1,:))

        GeneralFlux(3,:) = ConservedQuantity(2,:)/ConservedQuantity(1,:) &
            &   * (Gamma*ConservedQuantity(3,:)-0.5d0*Gmin1*ConservedQuantity(2,:)**2/ConservedQuantity(1,:))
!        write(6,*) GeneralFlux(1,1),GeneralFlux(1,2)
!        write(6,*) GeneralFlux(2,1),GeneralFlux(2,2)
        !write(6,*) GeneralFlux(3,1),GeneralFlux(3,2)

        DeltaFlux(:) = GeneralFlux(:,2) - GeneralFlux(:,1)
        DeltaQuantity(:) = ConservedQuantity(:,2) - ConservedQuantity(:,1)

        do iVariable=1,3
            DeltaFlux(iVar) = dsign(dabs(DeltaFlux(iVar)),DeltaQuantity(iVar))
                if(DeltaQuantity(iVar) /= 0.0d0) then
                    write(6,*) "flag 1"
                    if(DeltaFlux(iVar) < Epsilon) then
                        write(6,*) "flag 2"
                        DeltaFlux(iVar) = 0.5d0 * (DeltaFlux(iVar)**2 / Epsilon + Epsilon)
                    end if
                end if

                CE%NormalFluxDiff(iVar,iCX,1,1,1) = &
                    &   0.5d0 * ((GeneralFlux(iVar,1)+GeneralFlux(iVar,2)) - DeltaFlux(iVar))
        end do
        write(6,*) DeltaFlux,iCX
!        write(6,*) DeltaFlux
!        write(6,*) CE%NormalFluxDiff(:,iCX,1,1,1)
    end do
    return
    end subroutine FluxDim1_UpWind2


    subroutine FluxDim1_Quick!没
    implicit none
    double precision :: Q1,Q2,Q3
    do iFaceX=0,Geom%CellNumber(1)
        if(iFaceX == 0) then
            Q1=CC%ConservedQuantity(1,iCX,1,1)
            Q2=CC%ConservedQuantity(2,iCX,1,1)
            Q3=CC%ConservedQuantity(3,iCX,1,1)
        else
            Q1=0.375d0*(CC%ConservedQuantity(1,iCX+1,1,1)+0.75d0*CC%ConservedQuantity(1,iCX,1,1))&
                &   -0.125d0*CC%ConservedQuantity(1,iCX-1,1,1)
            Q2=0.375d0*(CC%ConservedQuantity(2,iCX+1,1,1)+0.75d0*CC%ConservedQuantity(2,iCX,1,1))&
                &   -0.125d0*CC%ConservedQuantity(2,iCX-1,1,1)
            Q3=0.375d0*(CC%ConservedQuantity(3,iCX+1,1,1)+0.75d0*CC%ConservedQuantity(3,iCX,1,1))&
                &   -0.125d0*CC%ConservedQuantity(3,iCX-1,1,1)
        end if

        CE%NormalFluxDiff(1,iFaceX,1,1,1) = Q2

        CE%NormalFluxDiff(2,iFaceX,1,1,1) = Q2**2/Q1+Gmin1*(Q3-0.5d0*Q2**2/Q1)

        CE%NormalFluxDiff(3,iFaceX,1,1,1) = Q2/Q1*(Gamma*Q3-0.5d0*Gmin1*Q2**2/Q1)
    end do
    return
    end subroutine FluxDim1_Quick



    subroutine FluxDim1_CenterSpace!没
    implicit none
    double precision :: Q1,Q2,Q3
    do iFaceX=0,Geom%CellNumber(1)
        Q1=0.5d0*(CC%ConservedQuantity(1,iCX+1,1,1)+CC%ConservedQuantity(1,iCX,1,1))
        Q2=0.5d0*(CC%ConservedQuantity(2,iCX+1,1,1)+CC%ConservedQuantity(2,iCX,1,1))
        Q3=0.5d0*(CC%ConservedQuantity(3,iCX+1,1,1)+CC%ConservedQuantity(3,iCX,1,1))

        CE%NormalFluxDiff(1,iFaceX,1,1,1) = Q2

        CE%NormalFluxDiff(2,iFaceX,1,1,1) = Q2**2/Q1+Gmin1*(Q3-0.5d0*Q2**2/Q1)

        CE%NormalFluxDiff(3,iFaceX,1,1,1) = Q2/Q1*(Gamma*Q3-0.5d0*Gmin1*Q2**2/Q1)
    end do
    return
    end subroutine FluxDim1_CenterSpace

    subroutine FluxDim14!没
    implicit none

    do iFaceX=0, Geom%CellNumber(1)
        CE%NormalFluxDiff(1,iFaceX,1,1,1) = CC%ConservedQuantity(2,iCX+1,1,1) - CC%ConservedQuantity(2,iCX,1,1)

        CE%NormalFluxDiff(2,iFaceX,1,1,1) = (CC%ConservedQuantity(2,iCX+1,1,1)*CC%PrimitiveVariable(2,iCX+1,1,1) &
            &  + Gmin1*(CC%ConservedQuantity(3,iCX+1,1,1) &
            &           - 0.5d0*CC%ConservedQuantity(2,iCX+1,1,1)*CC%PrimitiveVariable(2,iCX+1,1,1)))&
            &  - (CC%ConservedQuantity(2,iCX,1,1)*CC%PrimitiveVariable(2,iCX,1,1) &
            &  + Gmin1*(CC%ConservedQuantity(3,iCX,1,1)-0.5d0*CC%ConservedQuantity(2,iCX,1,1)*CC%PrimitiveVariable(2,iCX,1,1)))

        CE%NormalFluxDiff(3,iFaceX,1,1,1) = &
            &    (CC%PrimitiveVariable(2,iCX+1,1,1)*(Gamma*CC%ConservedQuantity(3,iCX+1,1,1) &
            &  - Gmin1*0.5d0*CC%ConservedQuantity(2,iCX+1,1,1)*CC%PrimitiveVariable(2,iCX+1,1,1))) &
            &  - (CC%PrimitiveVariable(2,iCX,1,1)*(Gamma*CC%ConservedQuantity(3,iCX,1,1) &
            &  - Gmin1*0.5d0*CC%ConservedQuantity(2,iCX,1,1)*CC%PrimitiveVariable(2,iCX,1,1)))

        do iVariable=1,3
            if(CE%NormalFluxDiff(iVar,iFaceX,1,1,1) /= CE%NormalFluxDiff(iVar,iFaceX,1,1,1)) call ErrorMessage(3)
        end do
    end do
    return
    end subroutine FluxDim14


    subroutine FluxDim2!没
    implicit none
    do iFaceY = 0, Geom%CellNumber(2)
        do iFaceX = 0, Geom%CellNumber(1)
            CE%NormalFluxDiff(1,iFX,iFY,1,1) = CC%ConservedQuantity(2,iCX+1,iCY,1) - CC%ConservedQuantity(2,iCX,iCY,1)

            CE%NormalFluxDiff(2,iFX,iFY,1,1) = &
                &     (CC%ConservedQuantity(2,iCX+1,iCY,1)*CC%PrimitiveVariable(2,iCX+1,iCY,1) &
                &   + CC%PrimitiveVariable(4,iCX+1,iCY,1)) &
                &   - (CC%ConservedQuantity(2,iCX,iCY,1)*CC%PrimitiveVariable(2,iCX,iCY,1) &
                &   + CC%PrimitiveVariable(4,iCX,iCY,1))

            CE%NormalFluxDiff(3,iFX,iFY,1,1) = &
                &    (CC%ConservedQuantity(1,iCX+1,iCY,1) &
                &   * CC%PrimitiveVariable(2,iCX+1,iCY,1) * CC%PrimitiveVariable(3,iCX+1,iCY,1)) &
                &   -(CC%ConservedQuantity(1,iCX,iCY,1) &
                &   * CC%PrimitiveVariable(2,iCX,iCY,1) * CC%PrimitiveVariable(3,iCX,iCY,1))

            CE%NormalFluxDiff(4,iFX,iFY,1,1) = &
                &    ((CC%ConservedQuantity(4,iCX+1,iCY,1) + CC%PrimitiveVariable(4,iCX+1,iCY,1)) &
                &   * CC%PrimitiveVariable(2,iCX+1,iCY,1)) &
                &   -((CC%ConservedQuantity(4,iCX,iCY,1) + CC%PrimitiveVariable(4,iCX,iCY,1)) &
                &   * CC%PrimitiveVariable(2,iCX,iCY,1))


            CE%NormalFluxDiff(1,iFX,iFY,1,2) = CC%ConservedQuantity(3,iCX,iCY+1,1) - CC%ConservedQuantity(3,iCX,iCY,1)

            CE%NormalFluxDiff(2,iFX,iFY,1,2) =  &
                &    (CC%ConservedQuantity(1,iCX,iCY+1,1) &
                &   * CC%PrimitiveVariable(2,iCX,iCY+1,1) * CC%PrimitiveVariable(3,iCX,iCY+1,1)) &
                &   -(CC%ConservedQuantity(1,iCX,iCY,1) &
                &   * CC%PrimitiveVariable(2,iCX,iCY,1) * CC%PrimitiveVariable(3,iCX,iCY,1))

            CE%NormalFluxDiff(3,iFX,iFY,1,2) = &
                &    (CC%ConservedQuantity(3,iCX,iCY+1,1)*CC%PrimitiveVariable(3,iCX,iCY+1,1) &
                &   + CC%PrimitiveVariable(4,iCX,iCY+1,1)) &
                &   -(CC%ConservedQuantity(3,iCX,iCY,1)*CC%PrimitiveVariable(3,iCX,iCY,1) &
                &   + CC%PrimitiveVariable(4,iCX,iCY,1))

            CE%NormalFluxDiff(4,iFX,iFY,1,2) = &
                &    ((CC%ConservedQuantity(4,iCX,iCY+1,1) + CC%PrimitiveVariable(4,iCX,iCY+1,1)) &
                &   * CC%PrimitiveVariable(3,iCX,iCY+1,1)) &
                &   -((CC%ConservedQuantity(4,iCX,iCY,1) + CC%PrimitiveVariable(4,iCX,iCY,1)) &
                &   * CC%PrimitiveVariable(3,iCX,iCY,1))
        do iSide=1,2
        do iVariable=1,4
            if(CE%NormalFluxDiff(iVar,iFX,iFY,1,iSide) /= CE%NormalFluxDiff(iVar,iFX,iFY,1,iSide)) call ErrorMessage(4)
        end do
        end do
        end do
    end do
    return
    end subroutine FluxDim2





end subroutine PrimaryAccuracyFlux
