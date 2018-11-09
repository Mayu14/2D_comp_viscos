!***********************************/
!	Name:VenkatakrishnanリミッターのΔmax，Δminを求めるプログラム
!	Alias:FindVariation
!	Description:(隣接セルの基本変数)-(自セルの基本変数)を計算して，変数ごとに最小値と最大値を評価する
!	Type:CellCenter
!	Input:CC%PrimitiveVariable
!	Output:CC%VariationOfVariable
!	Note:法線ベクトルについての考慮は式自体に内包してあるのでいらない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine UFindVariation(UG,UCC,UFV)
    use StructVar_Mod
    !use LoopVar_Mod
    use ConstantVar_Mod, only:CoreNumberOfCPU
    use UVenkatakrishnanVar_Mod4OMP
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(UFindVariationWithOMP), intent(inout) :: UFV
    integer :: iCell,iLocalEdge

!    double precision, allocatable :: UFV%Variation(:,:)

!隣接セルにおける最大変分と最小変分を検索(1次元の場合は左右の変動量の多い方を見つけるだけ)
!call FindUFV%Variation!(Geom,CC)
    !allocate(UFV%Variation(UG%GM%Dimension+2, 3)) !各変数,界面の数
!$omp parallel num_threads(CoreNumberOfCPU),shared(UG,UCC),firstprivate(iCell,iLocalEdge,UFV)
!$omp do
    do iCell=1,UG%GI%RealCells
        do iLocalEdge=1,3
            UFV%iAdjacentCell = UG%Tri%Cell(iCell,iLocalEdge)
            UFV%Variation(:,iLocalEdge) = UCC%PrimitiveVariable(:,UFV%iAdjacentCell,1,1) - UCC%PrimitiveVariable(:,iCell,1,1)
        end do
        UCC%VariationOfVariable(1,:,iCell,1,1) = maxval(UFV%Variation,2)
        UCC%VariationOfVariable(2,:,iCell,1,1) = minval(UFV%Variation,2)

    end do
!$omp end do
!$omp end parallel

    return
end subroutine UFindVariation
