!***********************************/
!	Name:初期条件の代入を制御するためのプログラム
!	Alias:UInitialize
!	Description:データからの読み出し(途中からの計算の再開)や複数の初期条件の切り替え等に対応
!	Type:Conf,Geom,CC,CE
!	Input:Configulation,Geometory,CellEdge,CellCenter
!	Output:
!	Note:ShockTube以外の初期条件はまだない
!	Author:Akitaka Toyota
!	Date:2017.11.11
!	Update:2017.11.11
!	Other:
!***********************************/
subroutine UInitialize(UConf,UG,UCC)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Configulation),intent(in) :: UConf
    type(UnstructuredGrid),intent(in) :: UG
    type(CellCenter),intent(inout) :: UCC

    if(UConf%UseResume == 1) then
        write(6,*) "Resume..."
        call UReadInitialCondition(UConf,UG,UCC)
        call UPrimitive2Conserve(UG,UCC)

    else if (UConf%UseResume == 5) then
        call UFlatIC

    else if (UConf%UseResume == 6) then
        call MFlatIC2

    else if (UConf%UseResume == 7) then
        call MFlatIC3

    else if (UConf%UseResume == 8) then
        call MFlatICIN
    else
        call UShockTubeIC(UG,UCC)

    end if
    !call AllocTest(CC,CE,Geom)
    call ZeroStepBoundary

return
contains
    subroutine UFlatIC
    implicit none

        UCC%ConservedQuantity(1,:,1,1) = 1.0d0
        UCC%ConservedQuantity(2,:,1,1) = 0.0d0
        UCC%ConservedQuantity(3,:,1,1) = 0.0d0
        UCC%ConservedQuantity(4,:,1,1) = 0.0d0
        UCC%ConservedQuantity(5,:,1,1) = (1.0d0/1.4d0/0.4d0)+0.5d0*(dot_product(UCC%ConservedQuantity(2:4,1,1,1),UCC%ConservedQuantity(2:4,1,1,1)))

    return
    end subroutine

    subroutine ZeroStepBoundary
    implicit none

        do iCell = UG%GI%RealCells+1, UG%GI%AllCells
            UCC%ConservedQuantity(:,iCell,1,1) = UCC%ConservedQuantity(:,UG%VC%Cell(iCell,1),1,1)
        end do

        return
    end subroutine ZeroStepBoundary


    subroutine MFlatIC2
    implicit none

        UCC%ConservedQuantity(1,:,1,1) = 1.0d0
        UCC%ConservedQuantity(2,:,1,1) = -0.3d0
        UCC%ConservedQuantity(3,:,1,1) = 0.0d0
        UCC%ConservedQuantity(4,:,1,1) = 0.0d0
        UCC%ConservedQuantity(5,:,1,1) = (1.0d0/1.4d0/0.4d0)+0.5d0*(dot_product(UCC%ConservedQuantity(2:4,1,1,1),UCC%ConservedQuantity(2:4,1,1,1)))

    return
    end subroutine


    subroutine MFlatIC3
    implicit none

        UCC%ConservedQuantity(1,:,1,1) = 1.0d0
        UCC%ConservedQuantity(2,:,1,1) = -1.5d0
        UCC%ConservedQuantity(3,:,1,1) = 0.0d0
        UCC%ConservedQuantity(4,:,1,1) = 0.0d0
        UCC%ConservedQuantity(5,:,1,1) = (1.0d0/1.4d0/0.4d0)+0.5d0*(dot_product(UCC%ConservedQuantity(2:4,1,1,1),UCC%ConservedQuantity(2:4,1,1,1)))

    return
    end subroutine

subroutine MFlatICIN
    implicit none

        UCC%ConservedQuantity(1,:,1,1) = UG%GM%BC%InFlowVariable(1)
        UCC%ConservedQuantity(2,:,1,1) = UG%GM%BC%InFlowVariable(1)*UG%GM%BC%InFlowVariable(2)
        UCC%ConservedQuantity(3,:,1,1) = UG%GM%BC%InFlowVariable(1)*UG%GM%BC%InFlowVariable(3)
        UCC%ConservedQuantity(4,:,1,1) = UG%GM%BC%InFlowVariable(1)*UG%GM%BC%InFlowVariable(4)
        UCC%ConservedQuantity(5,:,1,1) = (1.0d0/1.4d0/0.4d0)+0.5d0*(dot_product(UG%GM%BC%InFlowVariable(2:4),UG%GM%BC%InFlowVariable(2:4)))

    return
    end subroutine



end subroutine UInitialize
