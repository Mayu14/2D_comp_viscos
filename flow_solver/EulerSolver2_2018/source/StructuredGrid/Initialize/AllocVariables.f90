!***********************************/
!	Name:格子情報に基づいて動的割当を確保するプログラム
!	Alias:AllocTest
!	Description:
!	Type:Conf,Geom,CC,CE
!	Input:Configulation,Geometory,CellEdge,CellCenter
!	Output:
!	Note:1~3次元まですべて同じプログラムで実行するためにiXmaxなどを導入している
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2018.02.01
!	Other:
!***********************************/
subroutine AllocVariables(Conf,Geom,CC,CE)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none

    type(Configulation), intent(in) :: Conf
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CellEdge), intent(inout) :: CE

    integer, allocatable :: iMin(:)
    integer :: iYmax1,iZmax1

    iDim = Geom%Dimension

    allocate(iMin(3))
    iMin = 0
        if(Geom%CellNumber(2) == 0) iMin(2) = 1
        if(Geom%CellNumber(3) == 0) iMin(3) = 1

        iXmax = Geom%CellNumber(1)
        iYmax = max(1,Geom%CellNumber(2)) !セル総数と1の大きい方がiYmax
        iZmax = max(1,Geom%CellNumber(3))
        iYmax1 = Geom%CellNumber(2) + 1 !1次元のとき1,多次元のときiYmax+1
        iZmax1 = Geom%CellNumber(3) + 1 !2次元以下のとき1,3次元のときiZmax+1

    allocate(CC%ConservedQuantity(iDim+2,0:Geom%CellNumber(1)+1,iMin(2):Geom%CellNumber(2)+1,iMin(3):Geom%CellNumber(3)+1))
    allocate(CC%PrimitiveVariable(iDim+2,0:Geom%CellNumber(1)+1,iMin(2):Geom%CellNumber(2)+1,iMin(3):Geom%CellNumber(3)+1))
    allocate(CE%NormalFluxDiff(iDim+2,0:Geom%CellNumber(1),iMin(2):iYmax,iMin(3):iZmax,iDim))
    CC%ConservedQuantity = 0.0d0
    CC%PrimitiveVariable = 0.0d0
    CE%NormalFluxDiff = 0.0d0

    if(Conf%UseLocalTimeStep == 1 .or. Conf%UseOverSet == 1) then
        allocate(CC%TimeWidth(iXmax,iYmax,iZmax))
        allocate(CC%TmpTimeWidth(iMin(1):iXmax+1,iMin(2):iYmax1,iMin(3):iZmax1,iDim))
        CC%TmpTimeWidth = 0.0d0
    else
        allocate(CC%TimeWidth(1,1,1))
    end if
        CC%TimeWidth = 0.0d0

    if(Conf%UseRRK2 == 1) then
        allocate(CC%RungeKuttaG1(iDim+2,0:Geom%CellNumber(1)+1,iMin(2):Geom%CellNumber(2)+1,iMin(3):Geom%CellNumber(3)+1))
        allocate(CC%RungeKuttaG3(iDim+2,0:Geom%CellNumber(1)+1,iMin(2):Geom%CellNumber(2)+1,iMin(3):Geom%CellNumber(3)+1))
        allocate(CC%PreviousQuantity(iDim+2,0:Geom%CellNumber(1)+1,iMin(2):Geom%CellNumber(2)+1,iMin(3):Geom%CellNumber(3)+1))
        CC%RungeKuttaG1 = 0.0d0
        CC%RungeKuttaG3 = 0.0d0
        CC%PreviousQuantity = 0.0d0
    end if

    if(Conf%UseMUSCL == 1) then
        allocate(CC%GradientOfVariable(iDim+2,Geom%CellNumber(1),iYmax,iZmax,iDim))
        allocate(CC%VariationOfVariable(2,iDim+2,Geom%CellNumber(1),iYmax,iZmax))
        allocate(CC%LimiterFunction(iDim+2,Geom%CellNumber(1),iMin(2):iYmax,iMin(3):iZmax))
        CC%GradientOfVariable = 0.0d0
        CC%VariationOfVariable = 0.0d0
        CC%LimiterFunction = 0.0d0


        allocate(CE%NormalGradient(iDim+2,0:Geom%CellNumber(1)+1,iMin(2):iYmax1,iMin(3):iZmax1,2*iDim))
        allocate(CE%RebuildQunatity(iDim+2,0:Geom%CellNumber(1)+1,iMin(2):iYmax1,iMin(3):iZmax1,2*iDim)) !修正1101(セル中心からみた面番号に変更)
        allocate(CE%TmpLimiterFunction(iDim+2,0:Geom%CellNumber(1)+1,iMin(2):iYmax1,iMin(3):iZmax1,2*iDim))
        CE%NormalGradient = 0.0d0
        CE%RebuildQunatity = 0.0d0
        CE%TmpLimiterFunction = 0.0d0
    end if

    if(Conf%UseOverSet == 1) then
        allocate(CC%InterpolatedQuantity(iDim+2,0:Geom%CellNumber(1)+1,iMin(2):Geom%CellNumber(2)+1,iMin(3):Geom%CellNumber(3)+1))
    end if

return
end subroutine AllocVariables
