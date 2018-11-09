!***********************************/
!	Name:数値実験のまとめ
!	Alias:NumericalExperimentation
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:
!	Date:
!	Update:
!	Other:
!***********************************/
subroutine NumericalExperimentation(iStep,MConf,NE,CC,Geom,UG)!(iStep,iTargetGrid,CC,UG,MG,Geom)
    use StructVar_Mod
    use LoopVar_Mod
    use ExactSol
    use ConstantVar_Mod
    implicit none

    integer, intent(in) :: iStep
    type(Configulation), intent(in) :: MConf
    type(NumericalExperiment), intent(inout) :: NE
    type(CellCenter), intent(inout) :: CC
    type(Geometry), intent(in) :: Geom
    type(UnstructuredGrid), intent(in), optional :: UG
    integer :: iStep4Plot

    if(MConf%SwitchProgram == 0) then !Struct Euler
        if(MConf%UseResume < 5) then
            !厳密解呼出
            call ExactSolver(dble(iStep)*FixedTimeStep,NE%STT%Position,NE%STT%PrimitiveVariable(2,:),NE%STT%PrimitiveVariable(4,:),NE%STT%PrimitiveVariable(1,:))
            NE%STT%PrimitiveVariable(3,:) = 0.0d0
            call CalcResidual(CC,NE,Geom)

            if(mod(iStep,OutputInterval) == 0) then
                iStep4Plot = iStep/ OutputInterval
                call OutputExact(iStep4Plot,NE,Geom)
            end if

        end if

    else !unst Euler
        if(present(UG)) then
            if(MConf%UseResume < 5) then
                call ExactSolver(dble(iStep)*FixedTimeStep,NE%STT%Position,NE%STT%PrimitiveVariable(2,:),NE%STT%PrimitiveVariable(4,:),NE%STT%PrimitiveVariable(1,:))
                NE%STT%PrimitiveVariable(3,:) = 0.0d0
                call UCalcResidual(MConf,CC,NE,UG)

                if(mod(iStep,OutputInterval) == 0) then
                    iStep4Plot = iStep/ OutputInterval
                    call UOutputExact(MConf,UG,NE,iStep)
                end if
            end if
        else
            write(6,*) "UG is not loaded, at NumericalExperimentation Line53"
            stop
        end if
    end if



return
end subroutine NumericalExperimentation
