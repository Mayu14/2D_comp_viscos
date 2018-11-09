!***********************************/
!	Name:初期条件の代入を制御するためのプログラム
!	Alias:UInitialize
!	Description:データからの読み出し(途中からの計算の再開)や複数の初期条件の切り替え等に対応
!	Type:Conf,Geom,CC,CE
!	Input:Configulation,Geometory,CellEdge,CellCenter
!	Output:
!	Note:ShockTube以外の初期条件はまだない
!	Author:Akitaka Toyota
!	Date:2017.12.27
!	Update:
!	Other:
!***********************************/
subroutine MInitialize(UConf,UG,MG,UCC,Geom)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Geometry), intent(in) :: Geom
    type(MoveGrid), intent(in) :: MG
    type(Configulation),intent(in) :: UConf
    type(UnstructuredGrid),intent(in) :: UG
    type(CellCenter),intent(inout) :: UCC

    if(UConf%UseResume == 1) then
        write(6,*) "Resume..."
        call UReadInitialCondition(UConf,UG,UCC)
        call UPrimitive2Conserve(UG,UCC)

    else if (UConf%UseResume == 5) then
        call MFlatIC

    else if (UConf%UseResume == 6) then
        call MFlatIC2

    else
        call MShockTubeIC(UConf%UseResume,UG,MG,UCC,Geom)

    end if
    !call AllocTest(CC,CE,Geom)
return
contains

    subroutine MFlatIC
    implicit none

        UCC%ConservedQuantity(1,:,1,1) = 1.0d0
        UCC%ConservedQuantity(2,:,1,1) = 0.0d0
        UCC%ConservedQuantity(3,:,1,1) = 0.0d0
        UCC%ConservedQuantity(4,:,1,1) = 0.0d0
        UCC%ConservedQuantity(5,:,1,1) = (1.0d0/1.4d0/0.4d0)+0.5d0*(dot_product(UCC%ConservedQuantity(2:4,1,1,1),UCC%ConservedQuantity(2:4,1,1,1)))

    return
    end subroutine


    subroutine MFlatIC2
    implicit none

        UCC%ConservedQuantity(1,:,1,1) = 1.0d0
        UCC%ConservedQuantity(2,:,1,1) = -0.3d0
        UCC%ConservedQuantity(3,:,1,1) = 0.0d0
        UCC%ConservedQuantity(4,:,1,1) = 0.0d0
        UCC%ConservedQuantity(5,:,1,1) = (1.0d0/1.4d0/0.4d0)+0.5d0*(dot_product(UCC%ConservedQuantity(2:4,1,1,1),UCC%ConservedQuantity(2:4,1,1,1)))

    return
    end subroutine

end subroutine MInitialize
