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
subroutine CalcResidual(CC,NE,Geom) !(iStep,iTargetGrid,CC,UG,MG,Geom)
    use StructVar_Mod
    use LoopVar_Mod
    use ExactSol
    implicit none

    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(NumericalExperiment), intent(inout) :: NE

    call Conserve2Primitive(Geom,CC)
    NE%STT%ErrorNorm(:) = 0.0d0

    do iCenterY = 1, Geom%CellNumber(2)
        do iCenterX = 1,  Geom%CellNumber(1)
            NE%STT%ErrorNorm(:) = NE%STT%ErrorNorm(:)+sqrt((CC%PrimitiveVariable(:,iCenterX,iCenterY,1) - NE%STT%PrimitiveVariable(:,iCenterX))**2)
        end do
    end do

    NE%STT%ErrorNorm = NE%STT%ErrorNorm/(Geom%CellNumber(1)*Geom%CellNumber(2))

    do iLoop = 1,4
        write(100+iLoop,'(f22.14)') NE%STT%ErrorNorm(iLoop)
    end do

    return
end subroutine CalcResidual
