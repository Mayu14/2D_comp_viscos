!***********************************/
!	Name:仮想格子に境界条件を代入するプログラム
!	Alias:SetBoundary
!	Description:現状，勾配なしとすべり壁(反射境界)のみ
!	Type:CellCenter
!	Input:Geometry%Dimension,Geom%CellNumber, CellCenter%ConservedQuantity
!	Output:CellCenter%ConservedQuantity
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.??
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine SetBoundary(Geom,CC)
use StructVar_Mod
use LoopVar_Mod
implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    integer :: iWall
    iXmax = Geom%CellNumber(1)
    iYmax = max(1,Geom%CellNumber(2))
    iZmax = max(1,Geom%CellNumber(3))
    iDim = Geom%Dimension

!以下のプログラムでは計算領域の頂点セルがx,y,zすべてで上書きされるため不具合が起きそうだが，
!計算領域の頂点セル(x,y,z)=(0,0,0)のセルは計算式に現れないためセーフ
    !Density/Energy : 境界において密度勾配なし
    !速度も反転しない部分があるため一括で=にしとく
    do iWall=1, 4
        if(Geom%WallType(iWall) == 1) then
            call InOutFlowBoundary(iWall)

        else if(Geom%WallType(iWall) == 2) then
            call WallBoundary(iWall)

        else if(Geom%WallType(iWall) == 3) then
            call ReflectBoundary(iWall)

        end if
    end do



return
contains

    subroutine InOutFlowBoundary(iWallNumber)
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none
    integer, intent(in) :: iWallNumber
    integer :: iInOrOut, iSSOrSub !1:In or 2:Out , 1:SuperSonic or 2:Subsonic
    double precision, allocatable :: LocalVelocity(:), OuterNormal(:)
    double precision :: EdgeNormalVelocity, LocalSoundSpeed
    integer, pointer :: iEnd, iCX1, iCY1
    integer, target :: iXmax1, iYmax1, iBCLoop, iXmaxP, iYmaxP
    !0th step: Prepare pointer
        iXmaxP = iXmax
        iYmaxP = iYmax
        iXmax1 = iXmax + 1
        iYmax1 = iYmax + 1
        allocate(LocalVelocity(2),OuterNormal(2))

        if(iWallNumber == 1) then
            OuterNormal(1) = -1.0d0
            OuterNormal(2) = 0.0d0
            iCX => iZero !iCXは外側，iCX1は内側
            iCX1 => iOne
            iCY => iBCLoop
            iCY1 => iBCLoop
            iEnd => iYmaxP
        else if(iWallNumber == 2) then
            OuterNormal(1) = 1.0d0
            OuterNormal(2) = 0.0d0
            iCX => iXmax1
            iCX1 => iXmaxP
            iCY => iBCLoop
            iCY1 => iBCLoop
            iEnd => iYmaxP
        else if(iWallNumber == 3) then
            OuterNormal(1) = 0.0d0
            OuterNormal(2) = -1.0d0
            iCX => iBCLoop
            iCX1 => iBCLoop
            iCY => iZero
            iCY1 => iOne
            iEnd => iXmaxP
        else if(iWallNumber == 4) then
            OuterNormal(1) = 0.0d0
            OuterNormal(2) = 1.0d0
            iCX => iBCLoop
            iCX1 => iBCLoop
            iCY => iYmax1
            iCY1 => iYmaxP
            iEnd => iXmaxP
        end if

        do iBCLoop = 1, iEnd
        !1st step : judge of flow
            LocalVelocity = CC%ConservedQuantity(2:3,iCX1,iCY1,1) / CC%ConservedQuantity(1,iCX1,iCY1,1)

            EdgeNormalVelocity = dot_product(LocalVelocity,OuterNormal)

            if(EdgeNormalVelocity < 0.0d0) then
                iInOrOut = 1 !隣接実セルの持っている速度ベクトルが格子の外向き法線ベクトルに対して逆方向であるとき流入
            else
                iInOrOut = 2 !隣接実セルの保有する速度ベクトルと外向き法線ベクトルの向きが同方向であるとき流出
            end if

            LocalSoundSpeed = sqrt(Gamma*Gmin1 * (CC%ConservedQuantity(4,iCX1,iCY1,1) / CC%ConservedQuantity(1,iCX1,iCY1,1) &
                                & - 0.5d0 * dot_product(LocalVelocity,LocalVelocity))) !局所音速

            if(abs(EdgeNormalVelocity) > LocalSoundSpeed) then !法線方向速度の絶対値
                iSSOrSub = 1 !supersonic flow
            else
                iSSOrSub = 2 !subsonic flow
            end if

        !2nd step : Input Data
            if(iInOrOut == 1) then !Inflow
                if(iSSOrSub == 1) then !supersonic inflow
                    CC%ConservedQuantity(1,iCX,iCY,1) = Geom%BC%InFlowVariable(1)
                    CC%ConservedQuantity(2:3,iCX,iCY,1) = Geom%BC%InFlowVariable(1) * Geom%BC%InFlowVariable(2:3)
                    CC%ConservedQuantity(4,iCX,iCY,1) = &
                    &   Gmin1 * (Geom%BC%InFlowVariable(4) - 0.5d0 * Geom%BC%InFlowVariable(1) &
                    & * dot_product(Geom%BC%InFlowVariable(2:3),Geom%BC%InFlowVariable(2:3)))

                else !subsonic inflow
                    CC%PrimitiveVariable(4,iCX,iCY,1) = &
                    &   Gmin1 * (CC%ConservedQuantity(4,iCX1,iCY1,1) - 0.5d0 * Geom%BC%InFlowVariable(1) &
                    &   * dot_product(Geom%BC%InFlowVariable(2:3),Geom%BC%InFlowVariable(2:3)))

                    CC%ConservedQuantity(1,iCX,iCY,1) = Geom%BC%InFlowVariable(1)
                    CC%ConservedQuantity(2:3,iCX,iCY,1) = Geom%BC%InFlowVariable(1) * Geom%BC%InFlowVariable(2:3)
                    CC%ConservedQuantity(4,iCX,iCY,1) = &
                    &   InverseGmin1 * CC%PrimitiveVariable(4,iCX,iCY,1)  &
                    & + 0.5d0 * Geom%BC%InFlowVariable(1) * dot_product(Geom%BC%InFlowVariable(2:3),Geom%BC%InFlowVariable(2:3))

                end if

            else !OutFlow
                if(iSSOrSub == 1) then !supersonic outflow
                    CC%ConservedQuantity(:,iCX,iCY,1) = CC%ConservedQuantity(:,iCX1,iCY1,1)

                else !subsonic outflow
                    CC%ConservedQuantity(1:3,iCX,iCY,1) = CC%ConservedQuantity(1:3,iCX1,iCY1,1)
                    CC%ConservedQuantity(4,iCX,iCY,1) = &
                    &   Gmin1 * (Geom%BC%InFlowVariable(4) - 0.5d0 * CC%ConservedQuantity(1,iCX1,iCY1,1) &
                    & * dot_product(CC%ConservedQuantity(2:3,iCX1,iCY1,1),CC%ConservedQuantity(2:3,iCX1,iCY1,1)))

                end if
            end if
        end do

        return
    end subroutine InOutFlowBoundary

!__________________________

    subroutine WallBoundary(iWallNumber)
    implicit none
    integer, intent(in) :: iWallNumber !1:west, 2:east, 3:south, 4:north, 5:bottom, 6:top

    if(iWallNumber == 1) then
        CC%ConservedQuantity(:,0,:,:) = CC%ConservedQuantity(:,1,:,:)
        CC%ConservedQuantity(2,0,:,:) = 0.0d0 !Momentum：X境界においてX速度(運動量)のみ反射，他の速度成分は勾配なし

    else if(iWallNumber == 2) then
        CC%ConservedQuantity(:,iXmax+1,:,:) = CC%ConservedQuantity(:,iXmax,:,:)
        CC%ConservedQuantity(2,iXmax+1,:,:) = 0.0d0

    else if(iWallNumber == 3) then
        CC%ConservedQuantity(:,:,0,:) = CC%ConservedQuantity(:,:,1,:)
        CC%ConservedQuantity(3,:,0,:) = 0.0d0

    else if(iWallNumber == 4) then
        CC%ConservedQuantity(:,:,iYmax+1,:) = CC%ConservedQuantity(:,:,iYmax,:)
        CC%ConservedQuantity(3,:,iYmax+1,:) = 0.0d0

    else if(iWallNumber == 5) then
        CC%ConservedQuantity(:,:,:,0) = CC%ConservedQuantity(:,:,:,1)
        CC%ConservedQuantity(4,:,:,0) = 0.0d0

    else if(iWallNumber == 6) then
        CC%ConservedQuantity(:,:,:,iZmax+1) = CC%ConservedQuantity(:,:,:,iZmax)
        CC%ConservedQuantity(4,:,:,iZmax+1) = 0.0d0
    end if

    return
    end subroutine WallBoundary



    subroutine ReflectBoundary(iWallNumber)
    implicit none
    integer, intent(in) :: iWallNumber !1:west, 2:east, 3:south, 4:north, 5:bottom, 6:top
    !Density/Energy : 境界において密度勾配なし
    !速度も反転しない部分があるため一括で=にしとく
    if(iWallNumber == 1) then
        CC%ConservedQuantity(:,0,:,:) = CC%ConservedQuantity(:,1,:,:)
        CC%ConservedQuantity(2,0,:,:) = - CC%ConservedQuantity(2,1,:,:) !Momentum：X境界においてX速度(運動量)のみ反射，他の速度成分は勾配なし

    else if(iWallNumber == 2) then
        CC%ConservedQuantity(:,iXmax+1,:,:) = CC%ConservedQuantity(:,iXmax,:,:)
        CC%ConservedQuantity(2,iXmax+1,:,:) = - CC%ConservedQuantity(2,iXmax,:,:)

    else if(iWallNumber == 3) then
        CC%ConservedQuantity(:,:,0,:) = CC%ConservedQuantity(:,:,1,:)
        CC%ConservedQuantity(3,:,0,:) = - CC%ConservedQuantity(3,:,1,:)

    else if(iWallNumber == 4) then
        CC%ConservedQuantity(:,:,iYmax+1,:) = CC%ConservedQuantity(:,:,iYmax,:)
        CC%ConservedQuantity(3,:,iYmax+1,:) = - CC%ConservedQuantity(3,:,iYmax,:)

    else if(iWallNumber == 5) then
        CC%ConservedQuantity(:,:,:,0) = CC%ConservedQuantity(:,:,:,1)
        CC%ConservedQuantity(4,:,:,0) = - CC%ConservedQuantity(4,:,:,1)

    else if(iWallNumber == 6) then
        CC%ConservedQuantity(:,:,:,iZmax+1) = CC%ConservedQuantity(:,:,:,iZmax)
        CC%ConservedQuantity(4,:,:,iZmax+1) = - CC%ConservedQuantity(4,:,:,iZmax)
    end if

    return
    end subroutine ReflectBoundary


!________________________not use!
    subroutine OldSetBoundary(Geom,CC)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC

    iXmax = Geom%CellNumber(1)
    iYmax = max(1,Geom%CellNumber(2))
    iZmax = max(1,Geom%CellNumber(3))
    iDim = Geom%Dimension

!以下のプログラムでは計算領域の頂点セルがx,y,zすべてで上書きされるため不具合が起きそうだが，
!計算領域の頂点セル(x,y,z)=(0,0,0)のセルは計算式に現れないためセーフ
    !Density/Energy : 境界において密度勾配なし
    !速度も反転しない部分があるため一括で=にしとく

        CC%ConservedQuantity(:,0,:,:) = CC%ConservedQuantity(:,1,:,:)
        CC%ConservedQuantity(:,iXmax+1,:,:) = CC%ConservedQuantity(:,iXmax,:,:)

        if(iDim >= 2) then
            CC%ConservedQuantity(:,:,0,:) = CC%ConservedQuantity(:,:,1,:)
            CC%ConservedQuantity(:,:,iYmax+1,:) = CC%ConservedQuantity(:,:,iYmax,:)

            if(iDim == 3) then
                CC%ConservedQuantity(:,:,:,0) = CC%ConservedQuantity(:,:,:,1)
                CC%ConservedQuantity(:,:,:,iZmax+1) = CC%ConservedQuantity(:,:,:,iZmax)
            end if
        end if

    !Momentum：X境界においてX速度(運動量)のみ反射，他の速度成分は勾配なし
        CC%ConservedQuantity(2,0,:,:) = - CC%ConservedQuantity(2,1,:,:)
        CC%ConservedQuantity(2,iXmax+1,:,:) = - CC%ConservedQuantity(2,iXmax,:,:)

        if(iDim >= 2) then
            CC%ConservedQuantity(3,:,0,:) = - CC%ConservedQuantity(3,:,1,:)
            CC%ConservedQuantity(3,:,iYmax+1,:) = - CC%ConservedQuantity(3,:,iYmax,:)

            if(iDim == 3) then
                CC%ConservedQuantity(4,:,:,0) = - CC%ConservedQuantity(4,:,:,1)
                CC%ConservedQuantity(4,:,:,iZmax+1) = - CC%ConservedQuantity(4,:,:,iZmax)
            end if
        end if


    return
    end subroutine OldSetBoundary

end subroutine SetBoundary
