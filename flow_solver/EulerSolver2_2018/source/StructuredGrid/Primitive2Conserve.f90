!***********************************/
!	Name:基礎変数から保存変数へ変換するプログラム
!	Alias:Primitive2Conserve
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
subroutine Primitive2Conserve(Geom,CC)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter),intent(inout) :: CC

    nullify(iVar,iCX,iCY,iCZ)
    iVar => iVariable
    iCX => iCenterX
    iCY => iCenterY
    iCZ => iCenterZ
    iDim = Geom%Dimension

        if(iDim <= 2) iCZ => iOne
        if(iDim == 1) iCY => iOne

    do iCenterZ=0, Geom%CellNumber(3)+1
        do iCenterY=0, Geom%CellNumber(2)+1
            do iCenterX=0, Geom%CellNumber(1)+1
                CC%ConservedQuantity(1,iCX,iCY,iCZ) = CC%PrimitiveVariable(1,iCX,iCY,iCZ)

                do iVariable = 2, iDim+1
                    CC%ConservedQuantity(iVar,iCX,iCY,iCZ) = &
                        &   CC%PrimitiveVariable(iVar,iCX,iCY,iCZ)*CC%PrimitiveVariable(1,iCX,iCY,iCZ)
                end do

                CC%ConservedQuantity(iDim+2,iCX,iCY,iCZ) = &
                    &   InverseGmin1 * (CC%PrimitiveVariable(iDim+2,iCX,iCY,iCZ) + 0.5d0 * CC%PrimitiveVariable(1,iCX,iCY,iCZ) &
                    &   * dot_product(CC%PrimitiveVariable(2:iDim+1,iCX,iCY,iCZ),CC%PrimitiveVariable(2:iDim+1,iCX,iCY,iCZ)))

                !if(CC%PrimitiveVariable(1,iCX,iCY,iCZ) <= 0.0d0) call ErrorMessage(1)
                do iVariable=1,iDim+2
                if(CC%PrimitiveVariable(iVariable,iCX,iCY,iCZ) <= 0.0d0)then
                    write(6,*) CC%ConservedQuantity(1,iCX,iCY,iCZ),iCX,iCY,iCZ
                    write(6,*) CC%PrimitiveVariable(1,iCX,iCY,iCZ)
                    write(6,*) lbound(CC%ConservedQuantity)
                    write(6,*) ubound(CC%ConservedQuantity)
                   stop
               end if
               end do
               ! if(CC%PrimitiveVariable(iDim+2,iCX,iCY,iCZ) <= 0.0d0) call ErrorMessage(2)
            end do
            if(iDim == 1) exit
        end do
        if(iDim <= 2) exit
    end do

return
end subroutine Primitive2Conserve
