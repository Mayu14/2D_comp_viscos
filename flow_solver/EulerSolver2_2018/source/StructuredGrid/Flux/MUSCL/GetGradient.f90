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
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine GetGradient(Geom,CC)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, only:CoreNumberOfCPU
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    double precision :: InverseVolume

    nullify(iVar,iCX,iFX,iCY,iFY,iCZ,iFZ)
    iVar => iVariable
    iCX => iCenterX
    iCY => iCenterY
    iCZ => iCenterZ

    InverseVolume = 1.0d0 / Geom%Volume(1)

    if(Geom%Dimension == 1) then
        do iCenterX = 1, Geom%CellNumber(1)
            CC%GradientOfVariable(:,iCX,1,1,1) = 0.5d0*InverseVolume *  &
            &   (CC%PrimitiveVariable(:,iCX+1,1,1) - CC%PrimitiveVariable(:,iCX-1,1,1))
        end do

    else if(Geom%Dimension == 2) then
!$omp parallel num_threads(CoreNumberOfCPU),shared(Geom,CC,InverseVolume),private(iCenterX,iCenterY)
!$omp do
        do iCenterY = 1 , Geom%CellNumber(2)
            do iCenterX = 1, Geom%CellNumber(1)
                CC%GradientOfVariable(:,iCenterX,iCenterY,1,1) = 0.5d0*InverseVolume * Geom%Area(1) * &
                &   (CC%PrimitiveVariable(:,iCenterX+1,iCenterY,1) - CC%PrimitiveVariable(:,iCenterX-1,iCenterY,1))

                CC%GradientOfVariable(:,iCenterX,iCenterY,1,2) = 0.5d0*InverseVolume * Geom%Area(2) * &
                &   (CC%PrimitiveVariable(:,iCenterX,iCenterY+1,1) - CC%PrimitiveVariable(:,iCenterX,iCenterY-1,1))
            end do
        end do
!$omp end do
!$omp end parallel

    else if(Geom%Dimension == 3) then
        do iCenterZ = 1, Geom%CellNumber(3)
            do iCenterY = 1 , Geom%CellNumber(2)
                do iCenterX = 1, Geom%CellNumber(1)
                    CC%GradientOfVariable(:,iCX,iCY,iCZ,1) = 0.5d0*InverseVolume * Geom%Area(1) * &
                    &   (CC%PrimitiveVariable(:,iCX+1,iCY,iCZ) - CC%PrimitiveVariable(:,iCX-1,iCY,iCZ))

                    CC%GradientOfVariable(:,iCX,iCY,iCZ,2) = 0.5d0*InverseVolume * Geom%Area(2) * &
                    &   (CC%PrimitiveVariable(:,iCX,iCY+1,iCZ) - CC%PrimitiveVariable(:,iCX,iCY-1,iCZ))

                    CC%GradientOfVariable(:,iCX,iCY,iCZ,3) = 0.5d0*InverseVolume * Geom%Area(3) * &
                    &   (CC%PrimitiveVariable(:,iCX,iCY,iCZ+1) - CC%PrimitiveVariable(:,iCX,iCY,iCZ-1))
                end do
            end do
        end do
    end if

return
end subroutine GetGradient
