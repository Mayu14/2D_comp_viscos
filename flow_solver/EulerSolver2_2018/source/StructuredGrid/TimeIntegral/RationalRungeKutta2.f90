!***********************************/
!	Name:2段階有理ルンゲクッタ法による時間積分プログラム
!	Alias:RationalRungeKutta2
!	Description:時間2次精度
!	Type:Geom,CellEdge,CellCenter,iSwitch(1段目と2段目の切り替え用)
!	Input:Geom,CC,CE
!	Output:CC%ConservedQuantity
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.06
!	Other:現在，ExplicitEulerMethodと比較して1/3の時間進展速度になるという問題が確認されている
!***********************************/
subroutine RationalRungeKutta2(Geom,CE,CC,iSwitch)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod,RKPB1 => RungeKuttaParameterB1, RKPB2 => RungeKuttaParameterB2, RKPC2 => RungeKuttaParameterC2
implicit none
    type(Geometry),intent(in) :: Geom
    type(CellEdge), intent(in) :: CE
    type(CellCenter), intent(inout) :: CC
    integer, intent(in) :: iSwitch

    double precision :: InverseVolume
    double precision,allocatable :: FluxSum(:)
    type(RungeKuttaParameter) :: RKP
    allocate(FluxSum(Geom%Dimension+2))
    InverseVolume = 1.0d0/Geom%Volume(1)

    if(iSwitch == 1) then !1段階目
        call FirstStepOfRRK2

    else if(iSwitch == 0) then !2段階目"
        call SecondStepOfRRK2
    end if
return

contains

    subroutine FirstStepOfRRK2
    implicit none
        CC%PreviousQuantity = CC%ConservedQuantity
        if(Geom%Dimension == 1) then
            do iCenterX=1, Geom%CellNumber(1)
                FluxSum(:) = CE%NormalFluxDiff(:,iFX,1,1,1) - CE%NormalFluxDiff(:,iFX-1,1,1,1)
                CC%RungeKuttaG1(:,iCX,1,1) =  - InverseVolume*CC%TimeWidth(iTX,iTY,iTZ)*FluxSum
            end do

        else if(Geom%Dimension == 2)then
            do iCenterY=1,Geom%CellNumber(2)
                do iCenterX=1, Geom%CellNumber(1)
                    FluxSum(:) = (CE%NormalFluxDiff(:,iFX,iCY,1,1) - CE%NormalFluxDiff(:,iFX-1,iCY,1,1)) * Geom%Area(1) &
                            &   + (CE%NormalFluxDiff(:,iCX,iFY,1,2) - CE%NormalFluxDiff(:,iCX,iFY-1,1,2)) * Geom%Area(2)
                    CC%RungeKuttaG1(:,iCX,iCY,1) =  - InverseVolume*CC%TimeWidth(iTX,iTY,iTZ)*FluxSum
                end do
            end do

        else if(Geom%Dimension == 3) then
            do iCenterZ=1,Geom%CellNumber(3)
                do iCenterY=1,Geom%CellNumber(2)
                    do iCenterX=1, Geom%CellNumber(1)
                        FluxSum(:) = (CE%NormalFluxDiff(:,iFX,iCY,iCZ,1) - CE%NormalFluxDiff(:,iFX-1,iCY,iCZ,1))*Geom%Area(1) &
                                &   + (CE%NormalFluxDiff(:,iCX,iFY,iCZ,2) - CE%NormalFluxDiff(:,iCX,iFY-1,iCZ,2))*Geom%Area(2) &
                                &   + (CE%NormalFluxDiff(:,iCX,iCY,iFZ,3) - CE%NormalFluxDiff(:,iCX,iCY,iFZ-1,3))*Geom%Area(3)
                        CC%RungeKuttaG1(:,iCX,iCY,iCZ) =  - InverseVolume*CC%TimeWidth(iTX,iTY,iTZ)*FluxSum
                    end do
                end do
            end do
        end if
        CC%ConservedQuantity = CC%ConservedQuantity + RKPC2 * CC%RungeKuttaG1

    return
    end subroutine FirstStepOfRRK2

!_________________________________________________

    subroutine SecondStepOfRRK2
    implicit none
        if(Geom%Dimension == 1) then
            do iCenterX=1, Geom%CellNumber(1)
                FluxSum(:) = CE%NormalFluxDiff(:,iFX,1,1,1) - CE%NormalFluxDiff(:,iFX-1,1,1,1)
                CC%RungeKuttaG3(:,iCX,1,1) = RKPB1 * CC%RungeKuttaG1(:,iCX,1,1) &
                    &   - RKPB2*(-InverseVolume*CC%TimeWidth(iTX,iTY,iTZ)*FluxSum) !RungeKuttaG2
            end do

            RKP%G11 = sum(sum(sum(sum(CC%RungeKuttaG1*CC%RungeKuttaG1,4),3),2),1)
            RKP%G13 = sum(sum(sum(sum(CC%RungeKuttaG1*CC%RungeKuttaG3,4),3),2),1)
            RKP%InvG33 = 1.0d0 / sum(sum(sum(sum(CC%RungeKuttaG3*CC%RungeKuttaG3,4),3),2),1)

            do iCenterX=1, Geom%CellNumber(1)
                CC%ConservedQuantity(:,iCX,1,1) = CC%PreviousQuantity(:,iCX,1,1) &
                    &   + (2.0d0*CC%RungeKuttaG1(:,iCX,1,1) * RKP%G13 &
                    & - CC%RungeKuttaG3(:,iCX,1,1)*RKP%G11) * RKP%InvG33
            end do

    else if(Geom%Dimension == 2) then
            do iCenterY=1, Geom%CellNumber(2)
                do iCenterX=1, Geom%CellNumber(1)
                    FluxSum(:) = (CE%NormalFluxDiff(:,iFX,iCY,1,1) - CE%NormalFluxDiff(:,iFX-1,iCY,1,1)) * Geom%Area(1) &
                            &   + (CE%NormalFluxDiff(:,iCX,iFY,1,2) - CE%NormalFluxDiff(:,iCX,iFY-1,1,2)) * Geom%Area(2)
                    CC%RungeKuttaG3(:,iCX,iCY,1) = RKPB1 * CC%RungeKuttaG1(:,iCX,iCY,1) &
                    &   - RKPB2*(-InverseVolume*CC%TimeWidth(iTX,iTY,iTZ)*FluxSum) !RungeKuttaG2
                end do
            end do

            RKP%G11 = sum(sum(sum(sum(CC%RungeKuttaG1*CC%RungeKuttaG1,4),3),2),1)
            RKP%G13 = sum(sum(sum(sum(CC%RungeKuttaG1*CC%RungeKuttaG3,4),3),2),1)
            RKP%InvG33 = 1.0d0 / sum(sum(sum(sum(CC%RungeKuttaG3*CC%RungeKuttaG3,4),3),2),1)

            do iCenterY=1, Geom%CellNumber(2)
                do iCenterX=1, Geom%CellNumber(1)
                    CC%ConservedQuantity(:,iCX,iCY,1) = CC%PreviousQuantity(:,iCX,iCY,1) &
                        &   + (2.0d0*CC%RungeKuttaG1(:,iCX,iCY,1) * RKP%G13 &
                        & - CC%RungeKuttaG3(:,iCX,iCY,1)*RKP%G11) * RKP%InvG33
                end do
            end do

    else if(Geom%Dimension == 3) then
            do iCenterZ=1, Geom%CellNumber(3)
                do iCenterY=1, Geom%CellNumber(2)
                    do iCenterX=1, Geom%CellNumber(1)
                        FluxSum(:) = (CE%NormalFluxDiff(:,iFX,iCY,iCZ,1) - CE%NormalFluxDiff(:,iFX-1,iCY,iCZ,1))*Geom%Area(1) &
                                &   + (CE%NormalFluxDiff(:,iCX,iFY,iCZ,2) - CE%NormalFluxDiff(:,iCX,iFY-1,iCZ,2))*Geom%Area(2) &
                               &   + (CE%NormalFluxDiff(:,iCX,iCY,iFZ,3) - CE%NormalFluxDiff(:,iCX,iCY,iFZ-1,3))*Geom%Area(3)
                        CC%RungeKuttaG3(:,iCX,iCY,iCZ) = RKPB1 * CC%RungeKuttaG1(:,iCX,iCY,iCZ) &
                        &   - RKPB2*(-InverseVolume*CC%TimeWidth(iTX,iTY,iTZ)*FluxSum) !RungeKuttaG2
                    end do
                end do
            end do

            RKP%G11 = sum(sum(sum(sum(CC%RungeKuttaG1*CC%RungeKuttaG1,4),3),2),1)
            RKP%G13 = sum(sum(sum(sum(CC%RungeKuttaG1*CC%RungeKuttaG3,4),3),2),1)
            RKP%InvG33 = 1.0d0 / sum(sum(sum(sum(CC%RungeKuttaG3*CC%RungeKuttaG3,4),3),2),1)

            do iCenterZ=1, Geom%CellNumber(3)
                do iCenterY=1, Geom%CellNumber(2)
                    do iCenterX=1, Geom%CellNumber(1)
                        CC%ConservedQuantity(:,iCX,iCY,iCZ) = CC%PreviousQuantity(:,iCX,iCY,iCZ) &
                            &   + (2.0d0*CC%RungeKuttaG1(:,iCX,iCY,iCZ) * RKP%G13 &
                            & - CC%RungeKuttaG3(:,iCX,iCY,iCZ)*RKP%G11) * RKP%InvG33
                    end do
                end do
            end do

        end if
        return
    end subroutine SecondStepOfRRK2

end subroutine RationalRungeKutta2
