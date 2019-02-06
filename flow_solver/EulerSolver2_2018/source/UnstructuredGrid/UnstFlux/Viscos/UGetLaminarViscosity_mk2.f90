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
    double precision :: Mu0, SC1, Mach2, Tinflow
    double precision :: EdgeTemparature
    double precision :: FrontLength, BackLength
    logical :: debug = .true.

    Mach2 = MachNumber ** 2
    if(UConf%UseSutherlandLaw == 0) then
        UCE%LaminarViscosity = 1.0d0
        if(debug == .true.) then
            do iCell = 1, UG%GI%RealCells
                UCC%LaminarViscosity(iCell,1,1) = (UCE%LaminarViscosity(UG%Tri%Edge(iCell,1),1,1)+UCE%LaminarViscosity(UG%Tri%Edge(iCell,2),1,1)+UCE%LaminarViscosity(UG%Tri%Edge(iCell,3),1,1)) / 3.0
                UCC%EddyViscosity(iCell,1,1) = (UCE%EddyViscosity(UG%Tri%Edge(iCell,1),1,1)+UCE%EddyViscosity(UG%Tri%Edge(iCell,2),1,1)+UCE%EddyViscosity(UG%Tri%Edge(iCell,3),1,1)) / 3.0
            end do
        end if

    else
        Mu0 = 1.0d0
        Tinflow = gamma * UG%GM%BC%InFlowVariable(5) / UG%GM%BC%InFlowVariable(1)
        SC1 = SutherlandTemperature_S / ReferenceTemperature_Tref

        !　セル中心で温度・粘性・速度ノルム求める
        do iCell = 1, UG%GI%AllCells
            UCC%Temparature(iCell, 1, 1) = gamma * UCC%PrimitiveVariable(5, iCell, 1, 1) / UCC%PrimitiveVariable(1, iCell, 1, 1) * Gmin1
            UCC%LaminarViscosity(iCell,1,1) = SutherlandViscosity(UCC%Temparature(iCell,1,1), Tinflow)
        end do

        do iEdge = 1, UG%GI%Edges
            iBackCell = UG%Line%Cell(iEdge,2,1)
            if(iBackCell > UG%GI%RealCells) then
                UCE%LaminarViscosity(iEdge,1,1) = 0.5d0 * (UCC%LaminarViscosity(iFrontCell,1,1) + UCC%LaminarViscosity(iBackCell,1,1))
            else
                call GetLengthBetweenEdge(UG,iEdge,iFrontCell,iBackCell,FrontLength,BackLength)
                UCE%LaminarViscosity(iEdge,1,1) = (FrontLength * UCC%LaminarViscosity(iFrontCell,1,1) &
                                             &  +  BackLength * UCC%LaminarViscosity(iBackCell,1,1))  &
                                             &  / (FrontLength + BackLength)

            end if
        end do
    end if

    return
end subroutine UGetLaminarViscosity_mk2
