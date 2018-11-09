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
subroutine PreNumericalEx(maxWidth,NE) !(iStep,iTargetGrid,CC,UG,MG,Geom)
    use StructVar_Mod
    use LoopVar_Mod
    use ExactSol
    implicit none

    integer, intent(in) :: maxWidth
    type(NumericalExperiment), intent(inout) :: NE

    allocate(NE%STT%Position(maxWidth))
    allocate(NE%STT%PrimitiveVariable(4,maxWidth))
    allocate(NE%STT%ErrorNorm(4))

    return
end subroutine PreNumericalEx
