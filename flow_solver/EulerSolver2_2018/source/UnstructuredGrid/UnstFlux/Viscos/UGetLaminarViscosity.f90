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
subroutine UGetLaminarViscosity(UG,UCC)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    use ConstantVar_Mod, Mu0 => ReferenceViscosity_Mu0, Tref => ReferenceTemperature_Tref, SS => SutherlandTemperature_S, SC1 = SutherlandCoefficient1
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC



return
end subroutine UGetLaminarViscosity
