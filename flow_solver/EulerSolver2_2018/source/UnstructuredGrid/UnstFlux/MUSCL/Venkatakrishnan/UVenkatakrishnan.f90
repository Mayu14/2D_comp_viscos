!***********************************/
!	Name:Venkatakrishnanリミッターを求めるプログラム
!	Alias:Venkatakrishnan
!	Description:既に求めてあるΔmax,Δmin,Δ_を元に各格子における流束制限関数を求める
!	Type:CellCenter
!	Input:Geometory,CC%VariationOfVariable,CC%GradientOfVariable,CE%NormalGradient
!	Output:CC%LimiterFunction
!	Note:法線ベクトルについての考慮は式自体に内包してあるのでいらない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2018.01.11
!	Other:
!***********************************/
subroutine UVenkatakrishnan(UG,UCC,UCE,UV)
    use StructVar_Mod
    !use LoopVar_Mod
    use ConstantVar_Mod, VPK => VenkatakrishnanParameterK
    use UVenkatakrishnanVar_Mod4OMP
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(CellEdge), intent(inout) :: UCE
    integer :: iEdge, iSide, iVariable, iCell
    type(UVenkatakrishnanWithOMP), intent(inout) :: UV

!初期値として1.0(制限関数の最大値)を与える
        UCE%TmpLimiterFunction = 1.0d0

!$omp parallel num_threads(CoreNumberOfCPU),shared(UG,UCE,UCC),firstprivate(iCell,iEdge,iSide,iVariable,UV)
!$omp do
        do iEdge = 1, UG%GI%Edges
            do iSide=1,2
                do iVariable=1,5
                if(UCE%NormalGradient(iVariable,1,1,iSide,iEdge) == 0.0d0) then
                    UCE%TmpLimiterFunction(iVariable,1,1,iSide,iEdge) = 1.0d0

                else
                    iCell = UG%Line%Cell(iEdge,iSide,1) !界面の隣接要素を与えるが
                    if(iCell > UG%GI%RealCells) then
                        UCE%TmpLimiterFunction(iVariable,1,1,iSide,iEdge) = 1.0d0
                    else
                        if(UCE%NormalGradient(iVariable,1,1,iSide,iEdge) > 0.0d0) then
                            UV%iMaxMin => iOne !Δmax
                        else !if(UCE%NormalGradient(iVariable,1,1,iSide,iEdge) <= 0.0d0) then
                            UV%iMaxMin => iTwo !Δmin
                        end if

                        UV%VE2 = (VPK * 2.0d0*UG%GM%AverageWidth(iCell))**3

                        UCE%TmpLimiterFunction(iVariable,1,1,iSide,iEdge) = &
                        &   (UCC%VariationOfVariable(UV%iMaxMin,iVariable,iCell,1,1)**2 + UV%VE2 &
                            &   + 2.0d0 * UCE%NormalGradient(iVariable,1,1,iSide,iEdge) &
                            &   * UCC%VariationOfVariable(UV%iMaxMin,iVariable,iCell,1,1)) &
                        & / (UCC%VariationOfVariable(UV%iMaxMin,iVariable,iCell,1,1)**2  &
                            &   + 2.0d0*UCE%NormalGradient(iVariable,1,1,iSide,iEdge)**2 &
                            &   + UCC%VariationOfVariable(UV%iMaxMin,iVariable,iCell,1,1) &
                            &   * UCE%NormalGradient(iVariable,1,1,iSide,iEdge) + UV%VE2)

                    end if
                end if
                end do
            end do
        end do
!$omp end do
!$omp barrier
!$omp do
        do iCell = 1, UG%GI%RealCells
            UCC%LimiterFunction(:,iCell,1,1) = &
            & min(UCE%TmpLimiterFunction(:,1,1,1,UG%Tri%Edge(iCell,1)), UCE%TmpLimiterFunction(:,1,1,2,UG%Tri%Edge(iCell,1)), &
            &     UCE%TmpLimiterFunction(:,1,1,1,UG%Tri%Edge(iCell,2)), UCE%TmpLimiterFunction(:,1,1,2,UG%Tri%Edge(iCell,2)), &
            &     UCE%TmpLimiterFunction(:,1,1,1,UG%Tri%Edge(iCell,3)), UCE%TmpLimiterFunction(:,1,1,2,UG%Tri%Edge(iCell,3)))

        end do
!$omp end do
!$omp end parallel
    return
end subroutine UVenkatakrishnan
