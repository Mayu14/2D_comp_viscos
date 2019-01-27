!***********************************/
!	Name:初期条件の代入を制御するためのプログラム
!	Alias:Initialize
!	Description:データからの読み出し(途中からの計算の再開)や複数の初期条件の切り替え等に対応
!	Type:Conf,Geom,CC,CE
!	Input:Configulation,Geometory,CellEdge,CellCenter
!	Output:
!	Note:ShockTube以外の初期条件はまだない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine Initialize(Conf,Geom,CC,CE)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Configulation),intent(in) :: Conf
    type(Geometry),intent(in) :: Geom
    type(CellCenter),intent(inout) :: CC
    type(CellEdge),intent(inout) :: CE

    call AllocVariables(Conf,Geom,CC,CE)
    if(Conf%UseResume == 0) then
        write(6,*) "Please put VTK-file of last calculate on current directory"
        call ReadInitialCondition(Conf,Geom,CC)
        call Primitive2Conserve(Geom,CC)

    else if (Conf%UseResume == 5) then
        call SFlatIC

    else if (Conf%UseResume == 8) then
        call SFlatIC2

    else
        call ShockTubeIC(Conf%UseResume,Geom,CC)
    end if

    !call AllocTest(CC,CE,Geom)

return
contains
    subroutine SFlatIC
    implicit none

        CC%ConservedQuantity(1,:,:,1) = 1.0d0
        CC%ConservedQuantity(2,:,:,1) = 0.0d0
        CC%ConservedQuantity(3,:,:,1) = 0.0d0
        CC%ConservedQuantity(4,:,:,1) = (1.0d0/1.4d0/0.4d0)+0.5d0*(CC%ConservedQuantity(2,1,1,1)**2+CC%ConservedQuantity(3,1,1,1)**2)

    return
    end subroutine SFlatIC

    subroutine SFlatIC2
    implicit none

        CC%ConservedQuantity(1,:,:,1) = Geom%BC%InFlowVariable(1)
        CC%ConservedQuantity(2,:,:,1) = Geom%BC%InFlowVariable(1)*Geom%BC%InFlowVariable(2)
        CC%ConservedQuantity(3,:,:,1) = Geom%BC%InFlowVariable(1)*Geom%BC%InFlowVariable(3)
        CC%ConservedQuantity(4,:,:,1) = (1.0d0/1.4d0/0.4d0)+0.5d0*(Geom%BC%InFlowVariable(2)**2+Geom%BC%InFlowVariable(3)**2)

    return
    end subroutine SFlatIC2


end subroutine Initialize
