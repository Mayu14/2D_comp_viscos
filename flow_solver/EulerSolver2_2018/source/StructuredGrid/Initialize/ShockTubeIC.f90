!***********************************/
!	Name:1~2次元衝撃波管問題の初期条件を与えるプログラム
!	Alias:ShockTubeIC
!	Description:X軸に並行，Y軸に並行，XYの対角線に沿った初期条件の3つから選べる
!	Type:Geom,CC
!	Input:Geometory,CellCenter
!	Output:
!	Note:1次元の場合選択肢はない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine ShockTubeIC(iIC,Geom,CC)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none
    integer, intent(in) :: iIC
    type(Geometry),intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC

    iDim = Geom%Dimension
    nullify(iVar,iCX,iCY,iCZ)
    iVar => iVariable
    iCX => iCenterX
    iCY => iCenterY
    iCZ => iCenterZ

    if(Geom%Dimension == 1) then
        call Initial_YZ_Symmetry
    else
        !write(6,*) "Please choose initial condition."
        !write(6,*) "1:parallel2XAxis,2:parallel2YAxis,3:DiagonalXY"
        !read(5,*) iIC
        !iIC=1
            if(iIC == 3) then
                call Initial_XZ_Symmetry
            else if(iIC == 4) then
                call Initial_Diag
            else if(iIC == 99) then
                call TestInput
            else
                call Initial_YZ_Symmetry
            end if
    end if
return

contains
!_________________________________________________________
    subroutine TestInput
    implicit none
        CC%ConservedQuantity(1,:,:,:) = 1.0d0
        CC%ConservedQuantity(4,:,:,:) = 1.0d0
        CC%ConservedQuantity(2:3,:,:,:) = 0.0d0
        CC%ConservedQuantity(2,Geom%CellNumber(1)/2,Geom%CellNumber(2)/2,1) = 0.5d0
    return
    end subroutine TestInput


!_________________________________________________________
    subroutine Initial_YZ_Symmetry
    implicit none
        iXmax = Geom%CellNumber(1)
        iYmax = max(1,Geom%CellNumber(2)) !セル総数と1の大きい方がiYmax
        iZmax = max(1,Geom%CellNumber(3))
        CC%ConservedQuantity = 1.0d0
        do iCenterZ = 1, iZmax
            do iCenterY = 1, iYmax
                do iCenterX = 1, iXmax/2
                        CC%ConservedQuantity(1,iCX,iCY,iCZ) = 1.0d0
                    do iVariable = 2, iDim+1
                        CC%ConservedQuantity(iVar,iCX,iCY,iCZ) = 0.0d0
                    end do
                        CC%ConservedQuantity(iDim+2,iCX,iCY,iCZ) = 1.0d0/(Gamma - 1.0d0)/Gamma
                end do

                do iCenterX =iXmax/2+1, iXmax
                        CC%ConservedQuantity(1,iCX,iCY,iCZ) = 0.1d0
                    do iVariable = 2, iDim+1
                        CC%ConservedQuantity(iVar,iCX,iCY,iCZ) = 0.0d0
                    end do
                        CC%ConservedQuantity(iDim+2,iCX,iCY,iCZ) = 0.1d0/(Gamma - 1.0d0)/Gamma
                end do
            end do
        end do
    return
    end subroutine Initial_YZ_Symmetry

!_________________________________________________________

    subroutine Initial_XZ_Symmetry
    implicit none
        iXmax = Geom%CellNumber(1)
        iYmax = max(1,Geom%CellNumber(2)) !セル総数と1の大きい方がiYmax
        iZmax = max(1,Geom%CellNumber(3))
        CC%ConservedQuantity = 1.0d0
        do iCenterZ = 1, iZmax
            do iCenterX = 1, iXmax
                do iCenterY = 1, iYmax/2
                        CC%ConservedQuantity(1,iCX,iCY,iCZ) = 1.0d0
                    do iVariable = 2, iDim+1
                        CC%ConservedQuantity(iVar,iCX,iCY,iCZ) = 0.0d0
                    end do
                        CC%ConservedQuantity(iDim+2,iCX,iCY,iCZ) = 1.0d0/(Gamma - 1.0d0)/Gamma
                end do

                do iCenterY = iYmax/2+1, iYmax
                        CC%ConservedQuantity(1,iCX,iCY,iCZ) = 0.1d0
                    do iVariable = 2, iDim+1
                        CC%ConservedQuantity(iVar,iCX,iCY,iCZ) = 0.0d0
                    end do
                        CC%ConservedQuantity(iDim+2,iCX,iCY,iCZ) = 0.1d0/(Gamma - 1.0d0)/Gamma
                end do
            end do
        end do
    return
    end subroutine Initial_XZ_Symmetry

!_________________________________________________________

    subroutine Initial_Diag
    implicit none
    double precision :: CoordX,CoordY
    iXmax = Geom%CellNumber(1)
    iYmax = Geom%CellNumber(2)
    iZmax = max(1,Geom%CellNumber(3))
    CC%ConservedQuantity = 1.0d0

    do iCenterZ = 1, iZmax
        do iCenterY = 1, iYmax
            do iCenterX = 1, iXmax
                CoordX = Geom%Bound(1,1) + 2.0d0 * Geom%Width(1,1,1) * dble(iCX)
                CoordY = Geom%Bound(1,2) + 2.0d0 * Geom%Width(1,1,2) * dble(iCY)
                if(CoordX >= CoordY) then
                        CC%ConservedQuantity(1,iCX,iCY,iCZ) = 1.0d0
                    do iVariable = 2, iDim+1
                        CC%ConservedQuantity(iVar,iCX,iCY,iCZ) = 0.0d0
                    end do
                        CC%ConservedQuantity(iDim+2,iCX,iCY,iCZ) = 1.0d0/(Gamma - 1.0d0)/Gamma
                else
                        CC%ConservedQuantity(1,iCX,iCY,iCZ) = 0.1d0
                    do iVariable = 2, iDim+1
                        CC%ConservedQuantity(iVar,iCX,iCY,iCZ) = 0.0d0
                    end do
                        CC%ConservedQuantity(iDim+2,iCX,iCY,iCZ) = 0.1d0/(Gamma - 1.0d0)/Gamma
                end if
            end do
        end do
    end do
    return
    end subroutine Initial_Diag

end subroutine ShockTubeIC
