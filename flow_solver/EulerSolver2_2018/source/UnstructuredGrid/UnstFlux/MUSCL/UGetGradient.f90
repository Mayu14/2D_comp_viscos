!***********************************/
!	Name:セル中心における基礎変数勾配を求めるプログラム
!	Alias:GetGradient
!	Description:※勾配計算の式はスカラーに対するガウスの発散定理から導出される
!	Type:CellCenter
!	Input:Geometory,CCPrimitibeVariable
!	Output:CC%GradientOfVariable
!	Note:法線ベクトルについての考慮は式自体に内包してあるのでいらない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2018.01.11
!	Other:体積で除してなかったと思われる
!***********************************/
subroutine UGetGradient(UG,UCC)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, only:CoreNumberOfCPU
    use FrequentOperation
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    double precision, allocatable :: EdgeNormal(:)
    double precision :: MyLength, AdjLength
    allocate(EdgeNormal(3))

!面ごとに保存量と法線×面積...の定義式を計算
!総和を計算したあとに体積で除す
!本来はセル中心と界面の距離を元に荷重平均を取るべきではあるが，有意な差が出るとも思えないため単純に平均値を取る
!Cell内における勾配を求める
    UCC%GradientOfVariable = 0.0d0 !初期化

!たぶんここOMP入れたら値が不正になる
!!$omp parallel num_threads(CoreNumberOfCPU),shared(UG,UCC),firstprivate(iCell,iLocalEdge,iAdjacentCell,iVariable,EdgeNormal)
!!$omp do
    do iCell=1, UG%GI%RealCells !すべての実セルについて
        do iLocalEdge=1,3        !すべての局所界面についての総和を取る
            iEdge = UG%Tri%Edge(iCell,iLocalEdge) !セル界面の面番号
            iAdjacentCell = UG%Tri%Cell(iCell,iLocalEdge)   !隣接セルの番号

            if(iAdjacentCell > iCell) then !自セル位置の表裏を判定して法線ベクトルの向きを変え
                EdgeNormal(:) = UG%GM%Normal(iEdge,:)
                call GetLengthBetweenEdge(UG, iEdge, iCell, iAdjacentCell, MyLength, AdjLength)
            else
                EdgeNormal(:) = - UG%GM%Normal(iEdge,:)
                call GetLengthBetweenEdge(UG, iEdge,iAdjacentCell, iCell, AdjLength, MyLength)
            end if

            do iVariable = 1, 5 !全変数について勾配ベクトルを計算する
                !UCC%GradientOfVariable(iVariable,iCell,1,1,:) = UCC%GradientOfVariable(iVariable,iCell,1,1,:) &
                !& + 0.5d0 * (UCC%PrimitiveVariable(iVariable,iCell,1,1) + UCC%PrimitiveVariable(iVariable,iAdjacentCell,1,1)) &
                !& * UG%GM%Area(iEdge) * EdgeNormal(:)
                ! 界面上では距離に応じた荷重平均値を与える
                UCC%GradientOfVariable(iVariable,iCell,1,1,:) = UCC%GradientOfVariable(iVariable,iCell,1,1,:) &
                & + (MyLength * UCC%PrimitiveVariable(iVariable,iCell,1,1) + AdjLength * UCC%PrimitiveVariable(iVariable,iAdjacentCell,1,1)) &
                & / (MyLength + AdjLength) * UG%GM%Area(iEdge) * EdgeNormal(:)

            end do

        end do

        UCC%GradientOfVariable(:,iCell,1,1,:) = UCC%GradientOfVariable(:,iCell,1,1,:)/UG%GM%Volume(iCell)

    end do
!!$omp end do
!!$omp end parallel

return
end subroutine UGetGradient
