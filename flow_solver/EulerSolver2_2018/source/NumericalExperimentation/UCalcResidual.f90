!***********************************/
!	Name:数値実験のまとめ
!	Alias:PreNumericalEx
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
subroutine UCalcResidual(Conf,CC,NE,UG) !(iStep,iTargetGrid,CC,UG,MG,Geom)
    use StructVar_Mod
    use LoopVar_Mod
    use ExactSol
    implicit none

    type(Configulation), intent(in) :: Conf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: CC
    type(NumericalExperiment), intent(inout) :: NE
    integer :: iDir
    double precision :: Width

    call UConserve2Primitive(UG,CC)
    NE%STT%ErrorNorm(:) = 0.0d0

    if(Conf%UseResume == 2) then
        iDir = 1
    else if(Conf%UseResume == 3) then
        iDir = 2
    end if

    Width = 1.0d0/dble(UG%GI%AllCells)
    do iCell = 1, UG%GI%RealCells
        iPoint = nint(UG%CD%Cell(iCell,iDir)/Width)
        NE%STT%ErrorNorm(:) = NE%STT%ErrorNorm(:)+sqrt((CC%PrimitiveVariable(:,iCell,1,1) - NE%STT%PrimitiveVariable(:,iPoint))**2)
    end do

    NE%STT%ErrorNorm = NE%STT%ErrorNorm/dble(UG%GI%RealCells)

    do iLoop = 1,4
        write(100+iLoop,'(f22.14)') NE%STT%ErrorNorm(iLoop)
    end do

    return
end subroutine UCalcResidual
