!***********************************/
!	Name:MultipleOversetのLocalGrid用の初期条件代入ルーチン
!	Alias:MInitialize
!	Description:GlobalGridの初期条件との整合性を取る
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.26
!	Update:2017.12.27
!	Other:
!***********************************/
    subroutine MInitialize
    use StructVar_Mod
    implicit none
    !type(GlobalGrid), intent(inout) :: GG
    !type(Geometry), intent(in) :: Geom

    call MShockTubeIC

    return

    contains
    subroutine MShockTubeIC
    use StructVar_Mod
    implicit none
    type(CellCenter), intent(in) :: CC
    type(UnstructuredGrid), intent(in) :: UG
    type(MoveGrid), intent(in) :: MG

    integer :: iCell

        CC%ConservedQuantity(2:4,:,1,1) = 0.0d0     !common

    call YZ_Symmetry

    return
    end subroutine MShockTubeIC


    end subroutine MInitialize

