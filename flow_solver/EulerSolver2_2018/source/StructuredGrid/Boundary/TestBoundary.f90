!***********************************/
!	Name:境界条件が正しく代入されているか確認するためのプログラム
!	Alias:TestBoundary
!	Description:現状，勾配なしとすべり壁(反射境界)のみ，NaNの判定を行う
!	Type:CellCenter
!	Input:Geometry%Dimension,Geom%CellNumber, CellCenter%ConservedQuantity
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.??
!	Update:2017.11.01
!	Other:
!***********************************/
subroutine TestBoundary(Geom,CC)
use StructVar_Mod
use LoopVar_Mod
implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC

    iXmax = Geom%CellNumber(1)
    iYmax = max(1,Geom%CellNumber(2))
    iZmax = max(1,Geom%CellNumber(3))
    iDim = Geom%Dimension
!Density/Energy : 境界において密度勾配なし

do iLoop = 1, iDim+2,iDim+1
    do iCenterZ=1, iZmax
        do iCenterY=1,iYmax
            do iCenterX=1,iXmax
                if(CC%ConservedQuantity(iLoop,0,iCenterY,iCenterZ) /= CC%ConservedQuantity(iLoop,1,iCenterY,iCenterZ)) then
                    call ErrorMessage(5)
                end if
                if(CC%ConservedQuantity(iLoop,iXmax+1,iCenterY,iCenterZ)/= CC%ConservedQuantity(iLoop,iXmax,iCenterY,iCenterZ)) then
                    call ErrorMessage(5)
                end if
            if(iDim >= 2) then
                if(CC%ConservedQuantity(iLoop,iCenterX,0,iCenterZ) /= CC%ConservedQuantity(iLoop,iCenterX,1,iCenterZ)) then
                    call ErrorMessage(5)
                end if
                if(CC%ConservedQuantity(iLoop,iCenterX,iYmax+1,iCenterZ)/= CC%ConservedQuantity(iLoop,iCenterX,iYmax,iCenterZ)) then
                    call ErrorMessage(5)
                end if
            end if
            if(iDim == 3) then
                if(CC%ConservedQuantity(iLoop,iCenterX,iCenterY,0) /= CC%ConservedQuantity(iLoop,iCenterX,iCenterZ,1)) then
                    call ErrorMessage(5)
                end if
                if(CC%ConservedQuantity(iLoop,iCenterX,iCenterY,iZmax+1)/= CC%ConservedQuantity(iLoop,iCenterX,iCenterY,iZmax)) then
                    call ErrorMessage(5)
                end if
            end if
            end do
        end do
    end do
end do


!Momentum：境界において速度(運動量)反射
    do iLoop = 2, iDim+1
        do iCenterZ=1, iZmax
            do iCenterY=1,iYmax
                do iCenterX=1,iXmax
        if(CC%ConservedQuantity(iLoop,0,iCenterY,iCenterZ)/= - CC%ConservedQuantity(iLoop,1,iCenterY,iCenterZ)) then
            call ErrorMessage(5)
        end if
        if(CC%ConservedQuantity(iLoop,iXmax+1,iCenterY,iCenterZ)/= - CC%ConservedQuantity(iLoop,iXmax,iCenterY,iCenterZ)) then
            call ErrorMessage(5)
        end if

        if(iDim >= 2) then
            if(CC%ConservedQuantity(iLoop,iCenterX,0,iCenterZ)/= - CC%ConservedQuantity(iLoop,iCenterX,1,iCenterZ)) then
                call ErrorMessage(5)
            end if
            if(CC%ConservedQuantity(iLoop,iCenterX,iYmax+1,iCenterZ)/= - CC%ConservedQuantity(iLoop,iCenterX,iYmax,iCenterZ)) then
                call ErrorMessage(5)
            end if
        end if
            if(iDim == 3) then
                if(CC%ConservedQuantity(iLoop,iCenterX,iCenterY,0)/= - CC%ConservedQuantity(iLoop,iCenterX,iCenterY,1)) then
                    call ErrorMessage(5)
                end if
                if(CC%ConservedQuantity(iLoop,iCenterX,iCenterY,iZmax+1)/=-CC%ConservedQuantity(iLoop,iCenterX,iCenterY,iZmax)) then
                    call ErrorMessage(5)
                end if
            end if
                end do
            end do
        end do
    end do
return
end subroutine TestBoundary
