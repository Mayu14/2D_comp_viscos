!***********************************/
!	Name:移動先の格子を補間するための補間データ提供元格子検索
!	Alias:OFindDestination
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
subroutine InterpolationOfGlobalGrid(CC,Geom,LG)
    use StructVar_Mod
    use FrequentOperation
    implicit none
    type(CellCenter), intent(inout) :: CC
    type(Geometry), intent(in) :: Geom
    !type(UnstructuredGrid), intent(in) :: UG
    type(LocalGrids), intent(in) :: LG
    type(DonorElement), pointer :: DE
!    type(MoveGrid), intent(inout) :: MG

    integer :: iGlobalCellNumber
    integer :: iCenterX, iCenterY
    double precision :: AreaSum
    double precision, allocatable :: Moment(:)

!    write(6,*) "Cell value interpolation of Global Grid"
    allocate(Moment(4))
    CC%InterpolatedQuantity = CC%ConservedQuantity !補間処理が入らない格子はすべて前段階の値で上書き

    do iCenterY=1, Geom%CellNumber(2)
        do iCenterX=1, Geom%CellNumber(1)
            iGlobalCellNumber = Local2Global(iCenterX,iCenterY,Geom%CellNumber(1))

            DE => LG%GG%GDE(iGlobalCellNumber)%DE
            if(associated(DE)) then !そもそも補間の必要がない格子に関してはこのあとの処理はスルーで

                if(DE%Grid > 0) then
                    AreaSum = 0.0d0 !初期化
                    Moment = 0.0d0
                    do while(associated(DE)) !次のレコードがなくなるまで続行
                        AreaSum = AreaSum + DE%SharedArea !面積の総和計算
                        Moment(1:3) = Moment(1:3) + LG%CC(DE%Grid)%ConservedQuantity(1:3,DE%Cell,1,1) * DE%SharedArea !面積×保存量の総和計算
                        Moment(4) = Moment(4) + LG%CC(DE%Grid)%ConservedQuantity(5,DE%Cell,1,1) * DE%SharedArea !面積×保存量の総和計算 !StructEulerは4番にエネルギー入れてるため
                        DE => DE%Next !次のレコードへ

                    end do
                    CC%InterpolatedQuantity(:,iCenterX,iCenterY,1) = Moment/AreaSum

                else
                    !CC%InterpolatedQuantity(:,iCenterX,iCenterY,1) = 0.0d0 !物体の内部
                end if
            end if
        end do
    end do


    CC%ConservedQuantity = CC%InterpolatedQuantity

    return
end subroutine InterpolationOfGlobalGrid
