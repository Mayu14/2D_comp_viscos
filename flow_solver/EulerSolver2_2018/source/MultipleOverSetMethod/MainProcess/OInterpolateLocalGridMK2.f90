!***********************************/
!	Name:移動先の格子を補間するための補間データ提供元格子検索
!	Alias:OInterpolatedLocalGrid
!	Description:
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.11.16
!	Update:2017.11.16
!	Other:
!***********************************/
subroutine OInterpolateLocalGridMK2!(UG,OSG,MG,CC,UCC,CopyOfUCC)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(OverSetGrid), intent(in) :: OSG
    type(MoveGrid), intent(inout) :: MG
    type(CellCenter), intent(in) :: CC
    type(CellCenter), intent(inout) :: UCC, CopyOfUCC !ConservedQuantity only

    double precision, allocatable :: Moment(:)
    double precision :: InverseDistanceSum, InverseDistance, Distance
    double precision :: RadiasOfCircumCircle
    integer :: iAuxiliaryGridNumber, iContent

    allocate(Moment(5))


end subroutine OInterpolateLocalGrid
