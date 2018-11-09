!***********************************/
!	Name:保存変数から基礎変数へ変換するプログラム
!	Alias:Conserve2Primitive
!	Description:定義通りに計算するだけ
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.06
!	Other:
!***********************************/
subroutine Conserve2Primitive(Geom,CC)
    use StructVar_Mod
!    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter),intent(inout) :: CC
    integer :: iDim
    integer :: iCenterX,iCenterY,iCenterZ


    iDim = Geom%Dimension

    if(iDim == 1) then
            do iCenterX=0, Geom%CellNumber(1)+1

                CC%PrimitiveVariable(1,iCenterX,1,1) = CC%ConservedQuantity(1,iCenterX,1,1)


                CC%PrimitiveVariable(2:iDim+1,iCenterX,1,1) = &
                        &   CC%ConservedQuantity(2:iDim+1,iCenterX,1,1)/CC%ConservedQuantity(1,iCenterX,1,1)

                CC%PrimitiveVariable(iDim+2,iCenterX,1,1) = &
                    &   Gmin1 * (CC%ConservedQuantity(iDim+2,iCenterX,1,1) - 0.5d0 * CC%ConservedQuantity(1,iCenterX,1,1) &
                    &   * dot_product(CC%PrimitiveVariable(2:iDim+1,iCenterX,1,1),CC%PrimitiveVariable(2:iDim+1,iCenterX,1,1)))

            end do


    else if(iDim == 2) then
!$omp parallel num_threads(CoreNumberOfCPU),shared(Geom,CC,iDim),firstprivate(iCenterX,iCenterY)
!$omp do
            do iCenterY=0, Geom%CellNumber(2)+1
                do iCenterX=0, Geom%CellNumber(1)+1

                CC%PrimitiveVariable(1,iCenterX,iCenterY,1) = CC%ConservedQuantity(1,iCenterX,iCenterY,1)

                CC%PrimitiveVariable(2:iDim+1,iCenterX,iCenterY,1) = &
                    &   CC%ConservedQuantity(2:iDim+1,iCenterX,iCenterY,1)/CC%ConservedQuantity(1,iCenterX,iCenterY,1)

                CC%PrimitiveVariable(iDim+2,iCenterX,iCenterY,1) = &
                    &   Gmin1 * (CC%ConservedQuantity(iDim+2,iCenterX,iCenterY,1) - 0.5d0 * CC%ConservedQuantity(1,iCenterX,iCenterY,1) &
                    &   * dot_product(CC%PrimitiveVariable(2:iDim+1,iCenterX,iCenterY,1),CC%PrimitiveVariable(2:iDim+1,iCenterX,iCenterY,1)))

                end do
            end do
!$omp end do
!$omp end parallel

    else
        do iCenterZ=0, Geom%CellNumber(3)+1
            do iCenterY=0, Geom%CellNumber(2)+1
                do iCenterX=0, Geom%CellNumber(1)+1

                CC%PrimitiveVariable(1,iCenterX,iCenterY,iCenterZ) = CC%ConservedQuantity(1,iCenterX,iCenterY,iCenterZ)


                CC%PrimitiveVariable(2:iDim+1,iCenterX,iCenterY,iCenterZ) = &
                    &   CC%ConservedQuantity(2:iDim+1,iCenterX,iCenterY,iCenterZ)/CC%ConservedQuantity(1,iCenterX,iCenterY,iCenterZ)


                CC%PrimitiveVariable(iDim+2,iCenterX,iCenterY,iCenterZ) = &
                    &   Gmin1 * (CC%ConservedQuantity(iDim+2,iCenterX,iCenterY,iCenterZ) - 0.5d0 * CC%ConservedQuantity(1,iCenterX,iCenterY,iCenterZ) &
                    &   * dot_product(CC%PrimitiveVariable(2:iDim+1,iCenterX,iCenterY,iCenterZ),CC%PrimitiveVariable(2:iDim+1,iCenterX,iCenterY,iCenterZ)))

                end do
            end do
        end do

    end if



return
end subroutine Conserve2Primitive
