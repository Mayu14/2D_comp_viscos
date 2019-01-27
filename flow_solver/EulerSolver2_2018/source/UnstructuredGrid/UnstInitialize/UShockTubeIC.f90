!***********************************/
!	Name:1~2次元衝撃波管問題の初期条件を与えるプログラム
!	Alias:ShockTubeIC
!	Description:X軸に並行，Y軸に並行，XYの対角線に沿った初期条件の3つから選べる
!	Type:Geom,CC
!	Input:Geometory,CellCenter
!	Output:
!	Note:1次元の場合選択肢はない
!	Author:Akitaka Toyota
!	Date:2017.11.11
!	Update:2017.12.27
!	Other:mid周辺書き換え
!***********************************/
subroutine UShockTubeIC(UG,UCC)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none
    integer :: iIC
    type(UnstructuredGrid),intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    double precision :: mid
    iDim = UG%GM%Dimension

       ! write(6,*) "Please choose initial condition."
        !write(6,*) "1:parallel2XAxis, 2:parallel2YAxis, 3:DiagonalXY"
        !read(5,*) iIC
        iIC=1
            UCC%ConservedQuantity(1,:,:,:) = 1.0d0 !初期条件プロット時のConserve2Primitive対策
            UCC%ConservedQuantity(2:iDim+1,:,:,:) = 0.0d0 !Velocity can be neglect since all velocities are 0 in all region.
            if(iIC == 2) then
                !H:0<=X<0.5 L:0.5<x<=1.0
                mid = 0.5d0*(UG%GM%Bound(1,2)+UG%GM%Bound(2,2))
                call UInitial_XZ_Symmetry(mid)
            else if(iIC == 3) then
                call UInitial_Diag
            else
                !H:0<=X<0.5 L:0.5<x<=1.0
                mid = 0.5d0*(UG%GM%Bound(1,1)+UG%GM%Bound(2,1))
                call UInitial_YZ_Symmetry(mid)
            end if
    return

contains
!_________________________________________________________
    subroutine UHighPressureArea
    implicit none
        UCC%ConservedQuantity(1,iCell,1,1) = 1.0d0
        UCC%ConservedQuantity(iDim+2,iCell,1,1) = 1.0d0/(Gamma)/(Gamma - 1.0d0) !無次元化のため
    return
    end subroutine UHighPressureArea
!_________________________________________________________
    subroutine ULowPressureArea
    implicit none
        UCC%ConservedQuantity(1,iCell,1,1) = 0.1d0
        UCC%ConservedQuantity(iDim+2,iCell,1,1) = 0.1d0/(Gamma)/(Gamma - 1.0d0)
    return
    end subroutine ULowPressureArea
!_________________________________________________________
    subroutine UInitial_YZ_Symmetry(midX)
    implicit none
    double precision, intent(in) :: midX

        do iCell=1, UG%GI%RealCells
            if(UG%CD%Cell(iCell,1) < midX) then !H side
                call UHighPressureArea
            else
                call ULowPressureArea
            end if
        end do
    return
    end subroutine UInitial_YZ_Symmetry
!_________________________________________________________

    subroutine UInitial_XZ_Symmetry(midY)
    implicit none
    double precision, intent(in) :: midY

        do iCell=1, UG%GI%RealCells
            if(UG%CD%Cell(iCell,2) < midY) then !H side
                call UHighPressureArea
            else
                call ULowPressureArea
            end if
        end do
    return
    end subroutine UInitial_XZ_Symmetry
!_________________________________________________________

    subroutine UInitial_Diag
    implicit none
        !H:X<Y L:Y<=X
        do iCell=1, UG%GI%RealCells
            if(UG%CD%Cell(iCell,1) < UG%CD%Cell(iCell,2)) then !H side
                call UHighPressureArea
            else
                call ULowPressureArea
            end if
        end do
    return
    end subroutine UInitial_Diag

end subroutine UShockTubeIC
