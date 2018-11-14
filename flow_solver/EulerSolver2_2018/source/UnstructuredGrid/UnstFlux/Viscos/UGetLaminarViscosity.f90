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
    use ConstantVar_Mod, Mu0 => ReferenceViscosity_Mu0  !&
    !&   , Tref => ReferenceTemperature_Tref, SS => SutherlandTemperature_S, SC1 = SutherlandCoefficient1
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC

    if(UConf%UseSutherlandLaw == 0) then
        UCC%LaminarViscosity = 1.0d0    !
    else
        !call sutherland_law
    end if

return
contains
    subroutine Sutherland_Law
        implicit none

        return
    end subroutine Sutherland_Law

end subroutine UGetLaminarViscosity
