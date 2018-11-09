!***********************************/
!	Name:
!	Alias:VelocityCorrection
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.31
!	Update:
!	Other:
!***********************************/
    subroutine VelocityCorrection(iSwitch,UG,MG,CC4MB,MW,Timestep,MaxInfluenceRange)
    use StructVar_Mod
    use ConstantVar_Mod
    implicit none

    integer, intent(in) :: iSwitch
    double precision, intent(in) :: Timestep
    type(UnstructuredGrid), intent(in) :: UG
    type(MoveGrid), intent(in) :: MG
    type(MotionWait), intent(in) :: MW
    type(CellCenter), intent(inout) :: CC4MB
    double precision, intent(in) :: MaxInfluenceRange

    integer :: iCell
    double precision, allocatable :: ObjectVelocity(:),Velocity(:), AngularVelocity(:)
    double precision :: InverseTimeStep
    double precision :: Distance, InternalWindWeight, Gain!, InverseRadius
    double precision, allocatable :: RotationalMatrix(:,:)
    allocate(ObjectVelocity(3),Velocity(3))

    !InverseTimeStep = 1.0d0 / Timestep
    !ObjectVelocity = MW%StepDistance(1:3) * InverseTimeStep
    ObjectVelocity = MW%StepDistance(1:3)/Timestep
    !InverseRadius = -1.0d0 / (UG%RoughRadius)

    call UConserve2Primitive(UG,CC4MB)

    if(iSwitch == 2) ObjectVelocity = - ObjectVelocity

    !call SigmoidWeight
    !call LinearWeight
    call NoWeight

    return
contains
    subroutine NoWeight
    implicit none

        do iCell=1, UG%GI%AllCells

            CC4MB%ConservedQuantity(2:4,iCell,1,1) = CC4MB%ConservedQuantity(1,iCell,1,1)&
                &  * (CC4MB%PrimitiveVariable(2:4,iCell,1,1) - ObjectVelocity(1:3))

            CC4MB%ConservedQuantity(5,iCell,1,1) = &
                & InverseGmin1*CC4MB%PrimitiveVariable(5,iCell,1,1) + 0.5d0*CC4MB%ConservedQuantity(1,iCell,1,1) &
                & * (dot_product(CC4MB%PrimitiveVariable(2:4,iCell,1,1) - ObjectVelocity(1:3),CC4MB%PrimitiveVariable(2:4,iCell,1,1) - ObjectVelocity(1:3)))

        end do

    return
    end subroutine NoWeight


    subroutine LinearWeight
    implicit none
    double precision :: InverseRoughRadius

    InverseRoughRadius = 1.0d0/UG%RoughRadius

    do iCell=1, UG%GI%AllCells
            Distance = sqrt(dot_product(MG%RC%Cell(iCell,1:3),MG%RC%Cell(iCell,1:3)))

            InternalWindWeight = (-InverseRoughRadius*Distance + 1.0d0)

            Velocity = InternalWindWeight*ObjectVelocity

            CC4MB%ConservedQuantity(2:4,iCell,1,1) = CC4MB%ConservedQuantity(1,iCell,1,1)&
                &  * (CC4MB%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3))

            CC4MB%ConservedQuantity(5,iCell,1,1) = &
                & InverseGmin1*CC4MB%PrimitiveVariable(5,iCell,1,1) + 0.5d0*CC4MB%ConservedQuantity(1,iCell,1,1) &
                & * (dot_product(CC4MB%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3),CC4MB%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3)))


    end do

        return
    end subroutine LinearWeight


    subroutine SigmoidWeight
    implicit none
    double precision :: SigmoidX !SigmoidRange is defined in the Constant Variables module
    double precision :: GainMargin

        !Gain = SigmoidGain
        GainMargin = (6.0d0*min(1.0d0, dble(MG%iStepFromLastMove)/dble(MW%iEstimateWait))-2.0d0)

        do iCell=1, UG%GI%AllCells
            Distance = sqrt(dot_product(MG%RC%Cell(iCell,1:3),MG%RC%Cell(iCell,1:3)))

            if(iCell < UG%GI%RealCells+1) then
                Gain = -2.0d0/(UG%GM%AverageWidth(iCell))*log((3.0d0-sqrt(3.0d0))/(3.0d0+sqrt(3.0d0)))
                !SigmoidX = Gain*(Distance - (MaxInfluenceRange + 2.0d0*UG%GM%AverageWidth(iCell)) )
                SigmoidX = Gain*(Distance - (MaxInfluenceRange + GainMargin*UG%GM%AverageWidth(iCell)) )
                !SigmoidX = Gain*(Distance - (MaxInfluenceRange + GainMargin*UG%GM%AverageWidth(iCell) &
                !        &  * (6.0d0*min(1.0d0, dble(MG%iStepFromLastMove)/dble(MW%iEstimateWait))-2.0d0)))
            else
                Gain = -2.0d0/(UG%GM%AverageWidth(UG%VC%Cell(iCell,1)))*log((3.0d0-sqrt(3.0d0))/(3.0d0+sqrt(3.0d0)))
                SigmoidX = Gain*(Distance - (MaxInfluenceRange + GainMargin*UG%GM%AverageWidth(UG%VC%Cell(iCell,1))))
                !SigmoidX = Gain*(Distance - (MaxInfluenceRange + GainMargin*UG%GM%AverageWidth(UG%VC%Cell(iCell,1)) &
                !        &  * (6.0d0*min(1.0d0, dble(MG%iStepFromLastMove)/dble(MW%iEstimateWait))-2.0d0)))
            end if

            !Error procedure
                if(SigmoidX < - SigmoidRange) then !Error1
                    InternalWindWeight = 1.0d0 - 10.0d0**(-15)

                else if(SigmoidX > SigmoidRange) then !Error2
                    InternalWindWeight = 0.0d0

                else !Sigmoid
                    InternalWindWeight = -1.0d0 / (1.0d0 + exp(-SigmoidX)) + 1.0d0

                end if

                !InternalWindWeight = -InverseRadius*Distance+1.0d0

            Velocity = InternalWindWeight*ObjectVelocity

            CC4MB%ConservedQuantity(2:4,iCell,1,1) = CC4MB%ConservedQuantity(1,iCell,1,1)&
                &  * (CC4MB%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3))

            !CC4MB%ConservedQuantity(5,iCell,1,1) = &
            !    & InverseGmin1*CC4MB%PrimitiveVariable(5,iCell,1,1) + 0.5d0*CC4MB%ConservedQuantity(1,iCell,1,1) &
            !    & * (dot_product(CC4MB%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3),CC4MB%PrimitiveVariable(2:4,iCell,1,1) - Velocity(1:3)))

        end do

    return
    end subroutine SigmoidWeight


    subroutine MakeRotationMatrixOfZ(RotationAngle,RotationalMatrix)
    implicit none
        double precision, intent(in) :: RotationAngle
        double precision, intent(inout) :: RotationalMatrix(:,:)

        RotationalMatrix(1,1) = cos(RotationAngle)
        RotationalMatrix(1,2) =-sin(RotationAngle)
        RotationalMatrix(1,3) = 0.0d0
        RotationalMatrix(2,1) = sin(RotationAngle)
        RotationalMatrix(2,2) = cos(RotationAngle)
        RotationalMatrix(2,3) = 0.0d0
        RotationalMatrix(3,1) = 0.0d0
        RotationalMatrix(3,2) = 0.0d0
        RotationalMatrix(3,3) = 1.0d0

    return
    end subroutine MakeRotationMatrixOfZ

    end subroutine VelocityCorrection
