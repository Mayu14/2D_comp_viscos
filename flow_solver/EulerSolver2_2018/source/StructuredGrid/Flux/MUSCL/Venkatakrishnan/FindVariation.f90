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
!	Update:2018.01.11
!	Other:
!***********************************/
subroutine FindVariation(Geom,CC)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, only:CoreNumberOfCPU
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    double precision, allocatable :: Variation(:,:)
!隣接セルにおける最大変分と最小変分を検索(1次元の場合は左右の変動量の多い方を見つけるだけ)
!call FindVariation!(Geom,CC)
    allocate(Variation(Geom%Dimension+2, 2*Geom%Dimension)) !各変数,xyzの両側
    nullify(iCX,iCY,iCZ,iFX,iFY,iFZ,iVar)
    iVar => iVariable
    iCX => iCenterX
    iCY => iCenterY
    iCZ => iCenterZ


        if(Geom%Dimension == 1) then
            do iCenterX = 1, Geom%CellNumber(1)
                Variation(:,1) = CC%PrimitiveVariable(:,iCX-1,1,1) - CC%PrimitiveVariable(:,iCX,1,1)
                Variation(:,2) = CC%PrimitiveVariable(:,iCX+1,1,1) - CC%PrimitiveVariable(:,iCX,1,1)
                CC%VariationOfVariable(1,:,iCX,1,1) = maxval(Variation,2)
                CC%VariationOfVariable(2,:,iCX,1,1) = minval(Variation,2)
            end do

        else if(Geom%Dimension == 2) then
            !全実セルについてループ
!!$omp parallel num_threads(CoreNumberOfCPU),shared(Geom,CC),firstprivate(iCenterX,iCenterY,Variation)
!!$omp do
            do iCenterY = 1, Geom%CellNumber(2)
                do iCenterX = 1, Geom%CellNumber(1)
                    Variation(:,1) = CC%PrimitiveVariable(:,iCenterX-1,iCenterY,1) - CC%PrimitiveVariable(:,iCenterX,iCenterY,1)
                    Variation(:,2) = CC%PrimitiveVariable(:,iCenterX+1,iCenterY,1) - CC%PrimitiveVariable(:,iCenterX,iCenterY,1)
                    Variation(:,3) = CC%PrimitiveVariable(:,iCenterX,iCenterY-1,1) - CC%PrimitiveVariable(:,iCenterX,iCenterY,1)
                    Variation(:,4) = CC%PrimitiveVariable(:,iCenterX,iCenterY+1,1) - CC%PrimitiveVariable(:,iCenterX,iCenterY,1)
                    CC%VariationOfVariable(1,:,iCenterX,iCenterY,1) = maxval(Variation,2)
                    CC%VariationOfVariable(2,:,iCenterX,iCenterY,1) = minval(Variation,2)
                end do
            end do
!!$omp end do
!!$omp end parallel

        else if(Geom%Dimension == 3) then
            do iCenterZ = 1, Geom%CellNumber(3)
                do iCenterY = 1, Geom%CellNumber(2)
                    do iCenterX = 1, Geom%CellNumber(1)
                        Variation(:,1) = CC%PrimitiveVariable(:,iCX-1,iCY,iCZ) - CC%PrimitiveVariable(:,iCX,iCY,iCZ)
                        Variation(:,2) = CC%PrimitiveVariable(:,iCX+1,iCY,iCZ) - CC%PrimitiveVariable(:,iCX,iCY,iCZ)
                        Variation(:,3) = CC%PrimitiveVariable(:,iCX,iCY-1,iCZ) - CC%PrimitiveVariable(:,iCX,iCY,iCZ)
                        Variation(:,4) = CC%PrimitiveVariable(:,iCX,iCY+1,iCZ) - CC%PrimitiveVariable(:,iCX,iCY,iCZ)
                        Variation(:,5) = CC%PrimitiveVariable(:,iCX,iCY,iCZ-1) - CC%PrimitiveVariable(:,iCX,iCY,iCZ)
                        Variation(:,6) = CC%PrimitiveVariable(:,iCX,iCY,iCZ+1) - CC%PrimitiveVariable(:,iCX,iCY,iCZ)
                        CC%VariationOfVariable(1,:,iCX,iCY,1) = maxval(Variation,2)
                        CC%VariationOfVariable(2,:,iCX,iCY,1) = minval(Variation,2)
                    end do
                end do
            end do

        end if
        deallocate(Variation)
    return
end subroutine FindVariation

