!***********************************/
!	Name:Sutherland's lawに基づき粘度を求めるプログラム
!	Alias:UGetLaminarViscosity_mk2
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
subroutine UGetLaminarViscosity_mk2(UConf, UG, UCC, UCE)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    use ConstantVar_Mod, gamma => SpecificOfHeatRatio
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(CellEdge), intent(inout) :: UCE
    double precision :: Mu0, STS, SC1, Mach2, Tinf
    double precision :: EdgeTemparature
    logical :: debug = .true.

    Mach2 = MachNumber ** 2
    if(UConf%UseSutherlandLaw == 0) then
        UCE%LaminarViscosity = 1.0d0
        do iCell = 1, UG%GI%RealCells
            UCC%Temparature(iCell, 1, 1) = gamma * UCC%PrimitiveVariable(5, iCell, 1, 1) / UCC%PrimitiveVariable(1, iCell, 1, 1) * Mach2 ! Calculate "NON"-Dimensional Value
        end do

    else
        Mu0 = ReferenceViscosity_Mu0
        STS = SutherlandTemperature_S
        SC1 = SutherlandCoefficient1
        Tinf = InfinityTemperature

        do iCell = 1, UG%GI%RealCells
            UCC%Temparature(iCell, 1, 1) = gamma * UCC%PrimitiveVariable(5, iCell, 1, 1) / UCC%PrimitiveVariable(1, iCell, 1, 1) * Mach2 * Tinf ! Calculate Dimensional Value
        end do

        do iEdge = 1, UG%GI%Edges
            EdgeTemparature = 0.5d0 * (UCC%Temparature(UG%Line%Cell(iEdge, 1, 1), 1, 1) + UCC%Temparature(UG%Line%Cell(iEdge, 2, 1), 1, 1))
            UCE%LaminarViscosity(iEdge,1,1) = SC1 * EdgeTemparature ** 1.5d0 / (EdgeTemparature + STS) / InfinityStaticViscosity ! Calculate Non-Dimensional Value
        end do
    end if

    if(debug == .true.) then
        do iCell = 1, UG%GI%RealCells
            UCC%LaminarViscosity(iCell,1,1) = (UCE%LaminarViscosity(UG%Tri%Edge(iCell,1),1,1)+UCE%LaminarViscosity(UG%Tri%Edge(iCell,2),1,1)+UCE%LaminarViscosity(UG%Tri%Edge(iCell,3),1,1)) / 3.0
            UCC%EddyViscosity(iCell,1,1) = (UCE%EddyViscosity(UG%Tri%Edge(iCell,1),1,1)+UCE%EddyViscosity(UG%Tri%Edge(iCell,2),1,1)+UCE%EddyViscosity(UG%Tri%Edge(iCell,3),1,1)) / 3.0
        end do
    end if

    return
end subroutine UGetLaminarViscosity_mk2
