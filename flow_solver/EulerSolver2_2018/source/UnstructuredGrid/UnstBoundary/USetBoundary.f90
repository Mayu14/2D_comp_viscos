!***********************************/
!	Name:仮想格子に境界条件を代入するプログラム
!	Alias:USetBoundary
!	Description:現状，勾配なしとすべり壁(反射境界)のみ
!	Type:CellCenter
!	Input:Geometry%Dimension,Geom%CellNumber, CellCenter%ConservedQuantity
!	Output:CellCenter%ConservedQuantity
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.11.11
!	Update:2017.11.11
!	Other:
!***********************************/
subroutine USetBoundary(UG,UCC)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC

    iDim = UG%GM%Dimension

    do iCell = UG%GI%RealCells+1, UG%GI%AllCells

        if(UG%GM%CellType(iCell,1,1) == 1) then !Inflow and Outflow
            call InOutFlowBoundary

        else if(UG%GM%CellType(iCell,1,1) == 2) then !Wall
            if(invicid == .true.) then
                call WallBoundary
            else
                call NonSlipBoundary
            end if

        else if(UG%GM%CellType(iCell,1,1) == 3) then !Reflect Wall
            call ReflectBoundary

        else if(UG%GM%CellType(iCell,1,1) == 4) then !Overset Boundary
            !Not Use(Interpolate from other grid)

        else if(UG%GM%CellType(iCell,1,1) == 5) then !Overset inflow Boundary (only moving grid)
            !OversetPsuedoFlow

        end if

    end do

!write(6,*)
return
contains

    subroutine InOutFlowBoundary
    implicit none
    integer :: iInOrOut, iSSOrSub !1:In or 2:Out , 1:SuperSonic or 2:Subsonic
    double precision, allocatable :: LocalVelocity(:)
    double precision :: EdgeNormalVelocity, LocalSoundSpeed, VelocityEpsilon !VE is SEKO variable
    !1st step : judge of flow
        allocate(LocalVelocity(3))
        iAdjacentCell = UG%VC%Cell(iCell,1) !隣接する実セル
        iAdjacentEdge = UG%VC%Edge(iCell) !共有界面
        LocalVelocity = UCC%ConservedQuantity(2:iDim+1,iAdjacentCell,1,1) / UCC%ConservedQuantity(1,iAdjacentCell,1,1)
        !EdgeNormalVelocity = dot_product(LocalVelocity,-UG%GM%Normal(iAdjacentEdge,1:iDim))
        !SEKO procedure
        EdgeNormalVelocity = dot_product(UG%GM%BC%InFlowVariable(2:4),-UG%GM%Normal(iAdjacentEdge,1:iDim))/(sqrt(dot_product(UG%GM%BC%InFlowVariable(2:4),UG%GM%BC%InFlowVariable(2:4))))

        VelocityEpsilon = 0.02d0
        !VelocityEpsilon = dPi/6.0d0 !SEKO??
        !write(6,*) iAdjacentEdge,EdgeNormalVelocity
        if(EdgeNormalVelocity < - VelocityEpsilon) then

            iInOrOut = 1 !隣接実セルの持っている速度ベクトルが格子の外向き法線ベクトルに対して逆方向であるとき流入
        else if(EdgeNormalVelocity > VelocityEpsilon) then

            iInOrOut = 2 !隣接実セルの保有する速度ベクトルと外向き法線ベクトルの向きが同方向であるとき流出
        else
            iInOrOut = 2
            !iInOrOut = 3 !SEKO
            !call WallBoundary !SEKO
            !call ReflectBoundary !SEKO
        end if

        if(iInOrOut /= 3) then
        LocalSoundSpeed = sqrt(Gamma*Gmin1 * (UCC%ConservedQuantity(iDim+2,iAdjacentCell,1,1) / UCC%ConservedQuantity(1,iAdjacentCell,1,1) &
                            & - 0.5d0 * dot_product(LocalVelocity,LocalVelocity))) !局所音速

        if(abs(EdgeNormalVelocity) > LocalSoundSpeed) then !法線方向速度の絶対値
            iSSOrSub = 1 !supersonic flow
        else
            iSSOrSub = 2 !subsonic flow
        end if

        !2nd step : Input Data
        if(iInOrOut == 1) then !Inflow
            if(iSSOrSub == 1) then !supersonic inflow
                UCC%ConservedQuantity(1,iCell,1,1) = UG%GM%BC%InFlowVariable(1)
                UCC%ConservedQuantity(2:4,iCell,1,1) = UG%GM%BC%InFlowVariable(1) * UG%GM%BC%InFlowVariable(2:4)
                UCC%ConservedQuantity(5,iCell,1,1) = InverseGmin1*UG%GM%BC%InFlowVariable(5) + 0.5d0 *  UG%GM%BC%InFlowVariable(1) &
                    &  * dot_product(UG%GM%BC%InFlowVariable(2:4),UG%GM%BC%InFlowVariable(2:4))

            else !subsonic inflow
                UCC%PrimitiveVariable(5,iAdjacentCell,1,1) = &
                &   Gmin1 * (UCC%ConservedQuantity(5,iAdjacentCell,1,1) - 0.5d0 * UCC%ConservedQuantity(1,iAdjacentCell,1,1) &
                &   * dot_product(UCC%PrimitiveVariable(2:4,iAdjacentCell,1,1),UCC%PrimitiveVariable(2:4,iAdjacentCell,1,1)))

                UCC%ConservedQuantity(1,iCell,1,1) = UG%GM%BC%InFlowVariable(1)
                UCC%ConservedQuantity(2:4,iCell,1,1) = UG%GM%BC%InFlowVariable(1) * UG%GM%BC%InFlowVariable(2:4)
                UCC%ConservedQuantity(5,iCell,1,1) = InverseGmin1*UCC%PrimitiveVariable(5,iAdjacentCell,1,1) &
                    & + 0.5d0 *  UG%GM%BC%InFlowVariable(1) * dot_product(UG%GM%BC%InFlowVariable(2:4),UG%GM%BC%InFlowVariable(2:4))

            end if

        else !OutFlow
            if(iSSOrSub == 1) then !supersonic outflow
                UCC%ConservedQuantity(:,iCell,1,1) = UCC%ConservedQuantity(:,iAdjacentCell,1,1)

            else !subsonic outflow
                UCC%ConservedQuantity(1:4,iCell,1,1) = UCC%ConservedQuantity(1:4,iAdjacentCell,1,1)
                UCC%ConservedQuantity(iDim+2,iCell,1,1) = InverseGmin1*UG%GM%BC%OutFlowVariable(5) &
                    & + 0.5d0 * UCC%ConservedQuantity(1,iAdjacentCell,1,1) * dot_product(UCC%ConservedQuantity(2:4,iAdjacentCell,1,1),UCC%ConservedQuantity(2:4,iAdjacentCell,1,1))

            end if
        end if

        end if

    return
    end subroutine InOutFlowBoundary


    subroutine WallBoundary !速度の法線方向成分を0にする
    implicit none

        iAdjacentCell = UG%VC%Cell(iCell,1)
        iAdjacentEdge = UG%VC%Edge(iCell)
        UCC%ConservedQuantity(1,iCell,1,1) = UCC%ConservedQuantity(1,iAdjacentCell,1,1) !Grad 0
        UCC%ConservedQuantity(iDim+2,iCell,1,1) = UCC%ConservedQuantity(iDim+2,iAdjacentCell,1,1)
        !write(6,*) UCC%ConservedQuantity(1,iCell,1,1),iCell,UG%GI%RealCells
        UCC%ConservedQuantity(2:iDim+1,iCell,1,1) = UCC%ConservedQuantity(2:iDim+1,iAdjacentCell,1,1) &
        & - 1.0d0 * UG%GM%Normal(iAdjacentEdge,1:iDim) &
        & * sum(UCC%ConservedQuantity(2:iDim+1,iAdjacentCell,1,1)*UG%GM%Normal(iAdjacentEdge,1:iDim)) !Normal Direction Velocity is 0

    return
    end subroutine WallBoundary

    subroutine NonSlipBoundary ! すべての速度を0にする
    implicit none

        iAdjacentCell = UG%VC%Cell(iCell,1)
        iAdjacentEdge = UG%VC%Edge(iCell)
        UCC%ConservedQuantity(1,iCell,1,1) = UCC%ConservedQuantity(1,iAdjacentCell,1,1) !Grad 0
        UCC%ConservedQuantity(iDim+2,iCell,1,1) = UCC%ConservedQuantity(iDim+2,iAdjacentCell,1,1)   !Energy Gradient 0 !これ正しいの？
        UCC%ConservedQuantity(2:iDim+1,iCell,1,1) = 0.0d0

    return
    end subroutine NonSlipBoundary

    subroutine ReflectBoundary
    implicit none

        iAdjacentCell = UG%VC%Cell(iCell,1)
        iAdjacentEdge = UG%VC%Edge(iCell)
        UCC%ConservedQuantity(:,iCell,1,1) = UCC%ConservedQuantity(:,iAdjacentCell,1,1) !Grad 0
        UCC%ConservedQuantity(2:iDim+1,iCell,1,1) = UCC%ConservedQuantity(2:iDim+1,iAdjacentCell,1,1) &
        & - 2.0d0 * UG%GM%Normal(iAdjacentEdge,1:iDim) &
        & * sum(UCC%ConservedQuantity(2:iDim+1,iAdjacentCell,1,1)*UG%GM%Normal(iAdjacentEdge,1:iDim)) !Reflect

    return
    end subroutine ReflectBoundary


end subroutine USetBoundary
