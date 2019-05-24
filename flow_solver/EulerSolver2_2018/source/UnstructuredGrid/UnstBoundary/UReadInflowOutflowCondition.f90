!***********************************/
!	Name:流入条件・流出条件の読み込み
!	Alias:UReadInflowOutflowCondition
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.20
!	Update:
!	Other:
!***********************************/
subroutine UReadInflowOutflowCondition(UG, Uconf)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
use FrequentOperation
implicit none
    type(UnstructuredGrid), intent(inout) :: UG
    type(Configulation), intent(in) :: Uconf
    character(len=64) :: cFileName
    character :: cAnnotate

    !write(6,*) "Please input the name of In/Outflow condition's datafile"
    !read(5,*) cFileName
    cFileName = "InOutFlowCondition"
    if(.not. allocated(UConf%ResumeInFlowVars)) then
        open(unit=Uconf%my_rank+100,file=trim(adjustl(cFileName)),status='unknown')
            read(Uconf%my_rank+100,*) cAnnotate
            do iLoop=1,5
                read(Uconf%my_rank+100,*) cAnnotate, UG%GM%BC%InFlowVariable(iLoop)
            end do

            read(Uconf%my_rank+100,*) cAnnotate
            do iLoop=1,5
                read(Uconf%my_rank+100,*) cAnnotate, UG%GM%BC%OutFlowVariable(iLoop)
            end do
        close(Uconf%my_rank+100)
    else
        UG%GM%BC%InFlowVariable = Uconf%ResumeInFlowVars
        UG%GM%BC%OutFlowVariable = Uconf%ResumeOutFlowVars
    end if
    MachNumber = AbsVector(UG%GM%BC%InFlowVariable(2:4))
    UG%GM%BC%InFlowVariable(5) = UG%GM%BC%InFlowVariable(1)/Gamma !無次元化
    !UG%GM%BC%InFlowVariable(5) = &
    !&   UG%GM%BC%InFlowVariable(5)/Gmin1 + 0.5d0 * UG%GM%BC%InFlowVariable(1) &
    !&   * dot_product(UG%GM%BC%InFlowVariable(2:4),UG%GM%BC%InFlowVariable(2:4))

    !UG%GM%BC%InFlowVariable(2:4) = -UG%GM%BC%InFlowVariable(2:4)*UG%GM%BC%InFlowVariable(1) !速度→運動量 & 座標変換! 従来は物体の進行速度を入力していたため正負が反転

    UG%GM%BC%OutFlowVariable(5) = UG%GM%BC%OutFlowVariable(1)/Gamma !圧力の無次元化
    !UG%GM%BC%OutFlowVariable(5) = &
    !&   UG%GM%BC%OutFlowVariable(5)/Gmin1 + 0.5d0 * UG%GM%BC%OutFlowVariable(1) &
    !&   * dot_product(UG%GM%BC%OutFlowVariable(2:4),UG%GM%BC%OutFlowVariable(2:4)) !圧力→エネルギー

    !UG%GM%BC%OutFlowVariable(2:4) = -UG%GM%BC%OutFlowVariable(2:4)*UG%GM%BC%OutFlowVariable(1) !速度→運動量 & 座標変換(符号反転はそのため)

    call ChangeAngle(Uconf%dAttackAngle, UG%GM%BC%InFlowVariable(2:4))
    call ChangeAngle(Uconf%dAttackAngle, UG%GM%BC%OutFlowVariable(2:4))

    return

contains

    subroutine ChangeAngle(AttackAngle, Velocity)
        implicit none
        double precision, intent(in) :: AttackAngle
        double precision, intent(inout) :: Velocity(:)

        Velocity(1) = Velocity(1) * cos(AttackAngle) - Velocity(2) * sin(AttackAngle)
        Velocity(2) = Velocity(1) * sin(AttackAngle) + Velocity(2) * cos(AttackAngle)
        return
    end subroutine ChangeAngle

end subroutine UReadInflowOutflowCondition
