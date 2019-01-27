!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:OStoreDisplacement
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.01.10
!	Update:
!	Other:
!***********************************/
subroutine OStoreDisplacement(MW,iStep,iGrid)
    use StructVar_Mod
    use ConstantVar_Mod
    implicit none

    integer, intent(in) :: iStep,iGrid
    type(MotionWait), intent(inout) :: MW

    type(GravityCenterDisplacement) :: Input
    double precision :: deltaX,deltaY,deltaR,Zero
    integer :: iModulo, iRound
    allocate(Input%Translation(3),Input%Rotation(3))

    !call Sliding(1)
    !call DemoCylinder
    call DemoCylinder2
    !call DemoComparison
    !call DemoSQ
     !call MoveDemonstration
    !call FourCircleRotation
!Store input displacement (written in global coordinates)
    call Scaling(0.3d0)
    !call Scaling(0.0d0)
    !call StopGrid
    !Input%Rotation = 0.000001d0 * Input%Rotation
    Input%Rotation = 0.0d0
    MW%Accumulate%Translation = MW%Accumulate%Translation + Input%Translation
    MW%Accumulate%Rotation = MW%Accumulate%Rotation + Input%Rotation
    MW%StepDistance(1:3) = Input%Translation(1:3)
    MW%StepDistance(0) = Input%Rotation(3)

return
contains
    subroutine DemoCylinder
    implicit none

        call Sliding(1)
        if(iGrid == 2) then
            call Scaling(-1.0d0)
        end if

        return
    end subroutine DemoCylinder

    subroutine DemoCylinder2
    implicit none

        call Sliding(1)
        if(iGrid == 1) then
            call Scaling(-1.0d0)
        end if

        return
    end subroutine DemoCylinder2

    subroutine DemoSQ
    implicit none

        if(iGrid == 1) then
            call Sliding(1)
            if(dble(iStep)*FixedTimeStep > 1.0d0) then
                call Scaling(-1.0d0)
            end if

        else if(iGrid == 2) then
            call Sliding(2)
            if(dble(iStep)*FixedTimeStep < 1.0d0) then
                call Scaling(-1.0d0)
            end if

        else if(iGrid == 3) then
            call Sliding(1)
            if(dble(iStep)*FixedTimeStep < 1.0d0) then
                call Scaling(-1.0d0)
            end if

        else if(iGrid == 4) then
            call Sliding(2)
            if(dble(iStep)*FixedTimeStep > 1.0d0) then
                call Scaling(-1.0d0)
            end if
        end if

        call Scaling(0.1d0)
        return
    end subroutine DemoSQ


    subroutine StopGrid
    implicit none
        Input%Translation = 0.0d0
        Input%Rotation = 0.0d0
        return
    end subroutine


    subroutine Scaling(Scale)
    implicit none
        double precision, intent(in) :: Scale
        Input%Translation = Scale*Input%Translation
        Input%Rotation = Scale*Input%Rotation
    end subroutine Scaling


    subroutine Moving2D(Input,X,Y,R)
    use StructVar_Mod
    implicit none
    type(GravityCenterDisplacement), intent(inout) :: Input
    double precision, intent(in) :: X,Y,R
        Input%Translation(1) = X
        Input%Translation(2) = Y
        Input%Translation(3) = 0.0d0

        Input%Rotation(1) =  0.0d0
        Input%Rotation(2) =  0.0d0
        Input%Rotation(3) =  R
    return
    end subroutine


    subroutine Sliding(dir)
    implicit none
    integer, intent(in) :: dir

        deltaX = FixedTimeStep*dble(iStep) - FixedTimeStep*dble(iStep-1)
        deltaR = -2.0d0*dPi*(FixedTimeStep*dble(iStep) - FixedTimeStep*dble(iStep-1))
        !if(iGrid == 2) deltaX = 1.5d0*deltaX
    if(dir == 1) then
        call Moving2D(Input,deltaX,0.0d0,deltaR)
    else
        call Moving2D(Input,0.0d0,deltaX,deltaR)
    end if
    return
    end subroutine


    subroutine InputRotation(Input)
    use StructVar_Mod
    implicit none
    type(GravityCenterDisplacement), intent(inout) :: Input

        Input%Translation(1) = 0.0d0
        Input%Translation(2) = 0.0d0
        Input%Translation(3) = 0.0d0

        Input%Rotation(1) =  0.0d0
        Input%Rotation(2) =  0.0d0
        Input%Rotation(3) =  -20.0d0 !dPi/6.0d0*cos(dble(iStep)*dPi/4.0d0)
    return
    end subroutine InputRotation

    subroutine InputTranslation(Input)
    use StructVar_Mod
    implicit none
    type(GravityCenterDisplacement), intent(inout) :: Input

        Input%Translation(1) = 0.005d0
        Input%Translation(2) = 0.000d0
        Input%Translation(3) = 0.0d0

        Input%Rotation(1) =  0.0d0
        Input%Rotation(2) =  0.0d0
        Input%Rotation(3) =  0.0d0
    return
    end subroutine InputTranslation

    subroutine YokodaoshiHachinoji
    implicit none

        deltaX = 0.25d0*(cos(dPi*(1.0d0 + 2.0d0*FixedTimeStep*dble(iStep))) - cos(dPi*(1.0d0 + 2.0d0*FixedTimeStep*dble(iStep-1))))
        deltaY = 0.25d0*(sin(2.0d0*dPi*FixedTimeStep*dble(iStep)) - sin(2.0d0*dPi*FixedTimeStep*dble(iStep-1)))

        if(int(1.0d0 / FixedTimeStep) > iStep) then ! 1.0secまで
        else if(int(2.0d0 / FixedTimeStep) > iStep) then ! 2.0secまで
            deltaX = -deltaX
        else if(int(3.0d0 / FixedTimeStep) > iStep) then ! 3.0secまで
            deltaY = -deltaY
        else if(int(4.0d0 / FixedTimeStep) > iStep) then ! 4.0secまで
            deltaX = -deltaX
            deltaY = -deltaY
        end if

        call Moving2D(Input,deltaX,deltaY,0.0d0)
    return
    end subroutine

    subroutine FourCircleRotation
    implicit none
    double precision :: Radius
    Radius = 1.0d0
    Zero=0.0d0
        deltaX = Radius*(cos(dPi*(1.0d0 + 2.0d0*FixedTimeStep*dble(iStep))) - cos(dPi*(1.0d0 + 2.0d0*FixedTimeStep*dble(iStep-1))))
        deltaY = Radius*(sin(2.0d0*dPi*FixedTimeStep*dble(iStep)) - sin(2.0d0*dPi*FixedTimeStep*dble(iStep-1)))
    if(iGrid == 1) then
        call Moving2D(Input,deltaX,deltaY,0.0d0)
    else if(iGrid == 2) then
        call Moving2D(Input,deltaY,-deltaX,0.0d0)
    else if(iGrid == 3) then
        call Moving2D(Input,-deltaX,-deltaY,0.0d0)
    else if(iGrid == 4) then
        call Moving2D(Input,-deltaY,deltaX,0.0d0)
    end if
    return
    end subroutine FourCircleRotation

    subroutine MoveDemonstration
    implicit none
    Zero = 0.0d0
    deltaX = FixedTimeStep !0.002d0
    deltaY = FixedTimeStep !0.002d0
    deltaR = FixedTimeStep*dPi !0.0003d0 !0.05d0

    if(iGrid == 1) then
        iRound = int(1.0d0/FixedTimeStep)*8
    else if(iGrid == 2) then
        deltaX = -deltaX
        deltaY = -deltaY
        deltaR = FixedTimeStep !0.0006d0 !0.006d0
        iRound = int(1.0d0/FixedTimeStep)*8
    else if(iGrid == 3) then
        deltaR = deltaX
        deltaX = deltaY
        deltaY = deltaR
        deltaR = 0.0d0
        iRound = int(1.0d0/FixedTimeStep)*8
    else if(iGrid == 4) then
        deltaR = deltaX
        deltaX = - deltaY
        deltaY = - deltaR
        deltaR = 0.0d0
        iRound = int(1.0d0/FixedTimeStep)*8
    end if
        iModulo = mod(iStep,iRound)
        if(iModulo < iRound/8) then
            if(iGrid <= 2) call Moving2D(Input,Zero,deltaY,-deltaR)
            if(iGrid > 2) call Moving2D(Input,deltaX,Zero,-deltaR)
        else if(iModulo < iRound/4) then
            if(iGrid <= 2) call Moving2D(Input,deltaX,Zero,-deltaR)
            if(iGrid > 2) call Moving2D(Input,Zero,deltaY,-deltaR)
        else if (iModulo < iRound/8*3) then
            if(iGrid <= 2) call Moving2D(Input,Zero,-deltaY,-deltaR)
            if(iGrid > 2) call Moving2D(Input,-deltaX,Zero,-deltaR)
        else if(iModulo < iRound/8*5) then
            if(iGrid <= 2) call Moving2D(Input,-deltaX,Zero,-deltaR)
            if(iGrid > 2) call Moving2D(Input,Zero,-deltaY,-deltaR)
        else if (iModulo < iRound/4*3) then
            if(iGrid <= 2) call Moving2D(Input,Zero,-deltaY,-deltaR)
            if(iGrid > 2) call Moving2D(Input,-deltaX,Zero,-deltaR)
        else
            call Moving2D(Input,deltaX/2.0d0,deltaY/2.0d0,-deltaR)
        end if
        !call Moving2D(Input,deltaR*cos(iModulo*dPi/200),deltaR*sin(iModulo*dPi/200),-deltaR)

    !call InputRotation(Input)
    !call InputTranslation(Input)
    return
    end subroutine MoveDemonstration

end subroutine OStoreDisplacement
