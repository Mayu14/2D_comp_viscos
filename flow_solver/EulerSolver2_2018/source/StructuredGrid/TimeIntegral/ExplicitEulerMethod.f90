!***********************************/
!	Name:オイラー陽解法による時間積分プログラム
!	Alias:ExplicitEulerMethod
!	Description:時間1次精度(基本これで十分)
!	Type:Geom,CellEdge,CellCenter
!	Input:Geom,CC,CE
!	Output:CC%ConservedQuantity
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.06
!	Other:
!***********************************/
subroutine ExplicitEulerMethod(Geom,CE,CC)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod
implicit none
    type(Geometry),intent(in) :: Geom
    type(CellEdge), intent(in) :: CE
    type(CellCenter), intent(inout) :: CC
    double precision :: InverseVolume
    double precision,allocatable :: FluxSum(:)

        allocate(FluxSum(Geom%Dimension+2))

        InverseVolume = 1.0d0/Geom%Volume(1)

    if(Geom%Dimension == 1) then

        do iCenterX = 1, Geom%CellNumber(1)
            FluxSum(:) = CE%NormalFluxDiff(:,iFX,1,1,1)-CE%NormalFluxDiff(:,iFX-1,1,1,1)
            CC%ConservedQuantity(:,iCX,1,1) = CC%ConservedQuantity(:,iCX,1,1) - InverseVolume &
                &   * CC%TimeWidth(iTX,iTY,iTZ)*FluxSum
        end do


    else if(Geom%Dimension == 2) then
        !if(allocated(Geom%Interpolated)) then
        if(1 == 0) then
            do iCenterY = 1, Geom%CellNumber(2)
                do iCenterX = 1, Geom%CellNumber(1)
                    if(Geom%Interpolated(iCenterX,iCenterY,1) == 0) then
                        FluxSum(:) = (CE%NormalFluxDiff(:,iFX,iFY,1,1) - CE%NormalFluxDiff(:,iFX-1,iFY,1,1)) * Geom%Area(1) &
                            &   + (CE%NormalFluxDiff(:,iFX,iFY,1,2) - CE%NormalFluxDiff(:,iFX,iFY-1,1,2)) * Geom%Area(2)
                        CC%ConservedQuantity(:,iCX,iCY,1) = CC%ConservedQuantity(:,iCX,iCY,1) - InverseVolume &
                            &   * CC%TimeWidth(iTX,iTY,iTZ)*FluxSum

                    end if
                end do
            end do

        else

            do iCenterY = 1, Geom%CellNumber(2)
                do iCenterX = 1, Geom%CellNumber(1)
                        FluxSum(:) = (CE%NormalFluxDiff(:,iFX,iFY,1,1) - CE%NormalFluxDiff(:,iFX-1,iFY,1,1)) * Geom%Area(1) &
                            &   + (CE%NormalFluxDiff(:,iFX,iFY,1,2) - CE%NormalFluxDiff(:,iFX,iFY-1,1,2)) * Geom%Area(2)
                        CC%ConservedQuantity(:,iCX,iCY,1) = CC%ConservedQuantity(:,iCX,iCY,1) - InverseVolume &
                            &   * CC%TimeWidth(iTX,iTY,iTZ)*FluxSum
                end do
            end do
        end if

    else
        do iCenterZ = 1, Geom%CellNumber(3)
            do iCenterY = 1, Geom%CellNumber(2)
                do iCenterX = 1, Geom%CellNumber(1)
                    FluxSum(:) = (CE%NormalFluxDiff(:,iFX,iFY,iFZ,1) - CE%NormalFluxDiff(:,iFX-1,iFY,iFZ,1)) * Geom%Area(1) &
                            &   + (CE%NormalFluxDiff(:,iFX,iFY,iFZ,2) - CE%NormalFluxDiff(:,iFX,iFY-1,iFZ,2)) * Geom%Area(2) &
                            &   + (CE%NormalFluxDiff(:,iFX,iFY,iFZ,3) - CE%NormalFluxDiff(:,iFX,iFY,iFZ-1,3)) * Geom%Area(3)
                    CC%ConservedQuantity(:,iCX,iCY,iCZ) = CC%ConservedQuantity(:,iCX,iCY,iCZ) - InverseVolume &
                        &   * CC%TimeWidth(iTX,iTY,iTZ)*FluxSum
                end do
            end do
        end do
    end if

return
end subroutine ExplicitEulerMethod
