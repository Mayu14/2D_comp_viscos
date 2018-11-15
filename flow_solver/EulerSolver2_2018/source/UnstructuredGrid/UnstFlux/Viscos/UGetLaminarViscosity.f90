!***********************************/
!	Name:Sutherland's lawに基づき粘度を求めるプログラム
!	Alias:UGetLaminarViscosity
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.11.13
!	Update:-
!	Other:
!***********************************/
subroutine UGetLaminarViscosity(UConf, UG,UCC)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    use ConstantVar_Mod, gamma => SpecificOfHeatRatio
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    double precision :: Mu0, STS, SC1, Mach2!, Tref

    if(UConf%UseSutherlandLaw == 0) then
        UCC%LaminarViscosity = 1.0d0

    else
        Mu0 = ReferenceViscosity_Mu0
        !Tref = ReferenceTemperature_Tref
        STS = SutherlandTemperature_S
        SC1 = SutherlandCoefficient1
        Mach2 = MachNumber ** 2

        do iCell = 1, UG%GI%RealCells
            UCC%Temparature(iCell, 1, 1) = gamma * UCC%PrimitiveVariable(5, iCell, 1, 1) / UCC%PrimitiveVariable(1, iCell, 1, 1) * Mach2
            UCC%LaminarViscosity(iCell, 1, 1) = SC1 * UCC%Temparature(iCell, 1, 1) ** 1.5d0 / (UCC%Temparature(iCell, 1, 1) + STS)
        end do

    end if

    return
end subroutine UGetLaminarViscosity
