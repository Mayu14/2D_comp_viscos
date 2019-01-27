!***********************************/
!	Name:VenkatakrishnanリミッターのΔ_を求めるプログラム
!	Alias:GetNormalGradient
!	Description:セル内の基礎変数勾配にセル中心，界面中心間の位置ベクトルとの要素積を取る
!	Type:CellEdge
!	Input:Geometory,CC%GradientOfVariable,
!	Output:CE%NormalGradient
!	Note:法線ベクトルについての考慮は式自体に内包してあるのでいらない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2018.01.10
!	Other:
!***********************************/
subroutine UGetNormalGradient(UG,UCC,UCE,UGNG)
    use StructVar_Mod
    !use LoopVar_Mod
    use ConstantVar_Mod, only:CoreNumberOfCPU
    use UVenkatakrishnanVar_Mod4OMP
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(in) :: UCC
    type(CellEdge), intent(inout) :: UCE
    integer :: iEdge, iVariable
    type(UGetNormalGradientWithOMP) :: UGNG

!$omp parallel num_threads(CoreNumberOfCPU),shared(UG,UCE,UCC),firstprivate(iEdge,iVariable,UGNG)
!$omp do
    do iEdge=1, UG%GI%Edges !すべての界面について
        UGNG%iFrontCell = UG%Line%Cell(iEdge,1,1)
        UGNG%iBackCell = UG%Line%Cell(iEdge,2,1)
        do iVariable=1, 5
!            if(UGNG%iFrontCell > UG%GI%RealCells) then
                !UCE%NormalGradient(iVariable,1,1,1,iEdge) = (-1.0d0)*&
                !    &   dot_product(UCC%GradientOfVariable(iVariable,UGNG%iBackCell,1,1,:),UG%GM%Normal(iEdge,:)) !表側面
                !UCE%NormalGradient(iVariable,1,1,2,iEdge) = &
                !    &   dot_product(-UCC%GradientOfVariable(iVariable,UGNG%iBackCell,1,1,:),UG%GM%Normal(iEdge,:)) !表側面


            if(UGNG%iBackCell > UG%GI%RealCells) then
                UCE%NormalGradient(iVariable,1,1,1,iEdge) = &
                    &   dot_product(UCC%GradientOfVariable(iVariable,UGNG%iFrontCell,1,1,:),UG%GM%Width(UGNG%iFrontCell,UG%Line%Cell(iEdge,1,2),:)) !表側面

                UCE%NormalGradient(iVariable,1,1,2,iEdge) = 0.0d0 !仮想セルは再構築とかいらなさそうだし... !というかこの値は使わないはず
                !UCE%NormalGradient(iVariable,1,1,2,iEdge) = (-1.0d0)*&
                !    &   dot_product(UCC%GradientOfVariable(iVariable,UGNG%iFrontCell,1,1,:),UG%GM%Width(UGNG%iFrontCell,UG%Line%Cell(iEdge,2,2),:)) !表側面

            else
                UCE%NormalGradient(iVariable,1,1,1,iEdge) = &
                    &   dot_product(UCC%GradientOfVariable(iVariable,UGNG%iFrontCell,1,1,:),UG%GM%Width(UGNG%iFrontCell,UG%Line%Cell(iEdge,1,2),:)) !表側面

                UCE%NormalGradient(iVariable,1,1,2,iEdge) = &
                    &   dot_product(UCC%GradientOfVariable(iVariable,UGNG%iBackCell,1,1,:),UG%GM%Width(UGNG%iBackCell,UG%Line%Cell(iEdge,2,2),:)) !表側面
            end if

        end do
    end do
!$omp end do
!$omp end parallel

    return
end subroutine UGetNormalGradient
