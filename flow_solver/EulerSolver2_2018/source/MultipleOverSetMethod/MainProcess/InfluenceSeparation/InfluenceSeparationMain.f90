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
    subroutine InfluenceSeparationMain(iStep,MConf,MW,UG,MG,CC,CC4MB,CE4MB,SW)
    use StructVar_Mod
    use ConstantVar_Mod
!$ use omp_lib
    implicit none

    integer, intent(in) :: iStep
    type(Configulation), intent(in) :: MConf
    type(MotionWait), intent(inout) :: MW
    type(UnstructuredGrid), intent(in) :: UG
    type(MoveGrid), intent(inout) :: MG
    type(CellCenter), intent(in) :: CC
    type(CellCenter), intent(inout) :: CC4MB
    type(CellEdge), intent(inout) :: CE4MB
    type(StopWatch), intent(inout) :: SW
    integer :: iSplit, iLoop
    double precision :: InfluenceRange
    double precision, allocatable :: CopyDistance(:)

!$ SW%LapTime(1,9) = omp_get_wtime()
        allocate(CopyDistance(3))
        CopyDistance = MW%StepDistance
        !Initialize
        if(iStep == 1 .or. MW%iMotionWait == 0) then
            if(iStep /= 0) call MOutput4InfluenceArea(MConf,UG,MG%NSR%P,CC4MB,iStep,1,MG)
            call InitializeInfluenceSeparation(UG,MG,CC,CC4MB)
        end if

        !if(iStep /= 0) call ModifyQuantity(CC4MB,CC,MG,UG)
        !if(mod(iStep,2) == 0) call MOutput4InfluenceArea(MConf,UG,MG%NSR%P,CC4MB,iStep,1,MG) !For Test
        !call MOutput(MConf,UG,MG%NSR%P,CC,1,1,MG)
        !call MOutput4InfluenceArea(MConf,UG,MG%NSR%P,CC4MB,2,1,MG) !For Test
        !!速度順補正
        !call VelocityCorrection(1,UG,MG,CC4MB,MW,FixedTimeStep,InfluenceRange) !for influence separation

        !運動量の座標変換(ローカル系での値に修正)
        call OTransferInflowMomentum(UG,CC4MB,MG)

        !ローカル格子でのEulerEQ
        if(MConf%UseLocalTimeStep == 0) then
            call UCheckCFL4FixedTime(UG,CC4MB,iSplit)
            FixedTimeStep = FixedTimeStep / dble(iSplit)
            MW%StepDistance = MW%StepDistance / dble(iSplit)
        else
            iSplit = 1
        end if
!$ SW%LapTime(2,9) = omp_get_wtime()
!$ SW%LapTime(3,9) = SW%LapTime(3,9) + SW%LapTime(2,9) - SW%LapTime(1,9)
        do iLoop=1, iSplit

            MG%iStepFromLastMove = MG%iStepFromLastMove + 1
            !影響域の計算
!$ SW%LapTime(1,9) = omp_get_wtime()
            call ComputeInfluenceArea(UG,CC4MB,MW,MG)
            InfluenceRange = MG%MaxInfluenceRange

!$ SW%LapTime(2,9) = omp_get_wtime()
            !速度順補正
            if(iSplit == 1) call TransferDistance(1)
            call VelocityCorrection(1,UG,MG,CC4MB,MW,FixedTimeStep,InfluenceRange) !for influence separation
            !call MOutput4InfluenceArea(MConf,UG,MG%NSR%P,CC4MB,iStep,1,MG) !For Test
!$ SW%LapTime(3,9) = SW%LapTime(3,9) + SW%LapTime(2,9) - SW%LapTime(1,9)
            !Euler Solver
            call OUnstEuler(iStep,MConf,UG,CC4MB,CE4MB,SW)

            !速度逆補正
            call VelocityCorrection(2,UG,MG,CC4MB,MW,FixedTimeStep,InfluenceRange) !inverse correction    return
            !call TransferDistance(2)
            !VelocityCorrectionMK2
            !call MOutput4InfluenceArea(MConf,UG,MG%NSR%P,CC4MB,4,1,MG) !For Test
        end do

!$ SW%LapTime(1,9) = omp_get_wtime()
        FixedTimeStep = DefaultTimeStep

        !運動量の座標変換(グローバル系での値に修正)
        call OTransferOutflowMomentum(UG,CC4MB,MG)

        !!速度逆補正
        !!call VelocityCorrection(2,UG,MG,CC4MB,MW,FixedTimeStep,InfluenceRange) !inverse correction    return

!$ SW%LapTime(2,9) = omp_get_wtime()
!$ SW%LapTime(3,9) = SW%LapTime(3,9) + SW%LapTime(2,9) - SW%LapTime(1,9)

    return
contains
    subroutine TransferDistance(iSwitch)
    implicit none
    integer, intent(in) :: iSwitch

        if(iSwitch == 1) then !inflow
            MW%StepDistance(1:3) = matmul(MG%TM%Top2Next(1:3,1:3),MW%StepDistance(1:3))
        else
            MW%StepDistance = CopyDistance
        end if

        return
    end subroutine TransferDistance

    end subroutine InfluenceSeparationMain
