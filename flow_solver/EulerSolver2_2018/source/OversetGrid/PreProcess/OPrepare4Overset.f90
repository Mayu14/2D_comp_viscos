!***********************************/
!	Name:重合格子法で用いる配列のallocate
!	Alias:OPrepare4Overset
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.11.16
!	Update:2017.12.19
!	Other:
!***********************************/
subroutine OPrepare4Overset(UG,Geom,MG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(Geometry), intent(inout) :: Geom
    type(MoveGrid), intent(inout) :: MG

    integer :: iTotalStructCells
    double precision :: AverageVolume

    AverageVolume = sum(UG%GM%Volume)/dble(UG%GI%RealCells)

    iXmax = max(1,Geom%CellNumber(1))
    iYmax = max(1,Geom%CellNumber(2))
    iZmax = max(1,Geom%CellNumber(3))

    iTotalStructCells = iXmax * iYmax * iZmax



!以下廃止予定
    allocate(MG%GCD%Translation(3),MG%GCD%Rotation(3))
!以下実装予定
    allocate(MG%TM%Current2Next(4,4), MG%TM%Current2Top(4,4), MG%TM%Next2Current(4,4))
    allocate(MG%TM%Next2Top(4,4),MG%TM%Top2Current(4,4), MG%TM%Top2Next(4,4))
        MG%TM%Current2Next = 0.0d0
        do iLoop=1,4
            MG%TM%Current2Next(iLoop,iLoop) = 1.0d0
        end do
        !単位行列で初期化
        MG%TM%Current2Top = MG%TM%Current2Next
        MG%TM%Next2Current= MG%TM%Current2Next
        MG%TM%Next2Top = MG%TM%Current2Next
        MG%TM%Top2Current = MG%TM%Current2Next
        MG%TM%Top2Next = MG%TM%Current2Next

    allocate(MG%RC%Cell(UG%GI%AllCells,4))
    allocate(MG%RC%Point(UG%GI%Points,4))
    allocate(MG%RC%GravityCenter(2,4))

    allocate(MG%ASG%Width(3),MG%ASG%CellNumber(3))
    allocate(MG%ASG%Bound(2,3))

    allocate(MG%NSR%Previous2Next(UG%GI%AllCells,9))
    allocate(MG%NSR%InterpolateFromGlobal(UG%GI%AllCells))
    allocate(MG%NSR%C%PresentCoordsInG(UG%GI%AllCells,4))
    allocate(MG%NSR%C%NextCoordsInG(UG%GI%AllCells,4), MG%NSR%C%NextCoordsInN(UG%GI%AllCells,4))
    allocate(MG%NSR%P%PresentCoordsInG(UG%GI%Points,4))
    allocate(MG%NSR%P%NextCoordsInG(UG%GI%Points,4), MG%NSR%P%NextCoordsInN(UG%GI%Points,4))

    allocate(CopyOfUCC%ConservedQuantity(5,UG%GI%AllCells,1,1))

    allocate(Geom%Interpolated(iXmax,iYmax,iZmax)) !この処理は1回だけでいいのでこのプログラムの外に出す


return
end subroutine OPrepare4Overset
