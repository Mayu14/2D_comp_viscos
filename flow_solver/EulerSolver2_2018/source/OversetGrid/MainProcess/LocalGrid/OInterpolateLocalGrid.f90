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
subroutine OInterpolateLocalGrid(UG,OSG,MG,CC,UCC,CopyOfUCC)
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

    CopyOfUCC%ConservedQuantity = UCC%ConservedQuantity

    do iCell=1, UG%GI%AllCells !ローカル格子のすべてのセルについて
        if(iCell <= UG%GI%RealCells) then
            if(MG%NSR%InterpolateFromGlobal(iCell) == 0) then !interpolate from local cell

                InverseDistanceSum = 0.0d0
                Moment = 0.0d0

                !本来は補間対象が三角形格子の面積中に存在しているという条件の下で加重平均を取るべきだが，計算が複雑になりすぎるので外接円内に存在するセルなら大丈夫とする
                !加重平均の基準に使うための外接円半径の計算
                RadiasOfCircumCircle = 0.25d0 / UG%GM%Volume(iCell) &
                & * UG%GM%Area(UG%Tri%Edge(iCell,1))*UG%GM%Area(UG%Tri%Edge(iCell,2))*UG%GM%Area(UG%Tri%Edge(iCell,3))

                do iLoop=1,9 !被る可能性のある補助格子すべてについて
                    iAuxiliaryGridNumber = MG%NSR%Previous2Next(iCell,iLoop) !iCellと被る可能性のある補助格子のうちiLoop番目に登録された格子のグローバル番号

                    do iContent=1, MG%ASG%ContentsAS2U(iAuxiliaryGridNumber) !その補助格子中に含まれるすべてのセルの数iContentについて
!                        iAdjacentCell = MG%ASG%RelatedAS2U(iAuxiliaryGridNumber,iContent) !補助格子iAuxiliaryGridNumber中のiContent番目のセル番号を取り出し
                        Distance = sqrt(sum((MG%NSR%C%PresentCoordsInG(iAdjacentCell,1:3) - MG%NSR%C%NextCoordsInG(iCell,1:3))**2)) !

                        if(Distance <= RadiasOfCircumCircle) then !セル中心間距離が外接円半径以下であるとき加重平均の計算対象にする（それ以外はスルー）
                            InverseDistance = 1.0d0/Distance
                            InverseDistanceSum = InverseDistanceSum + InverseDistance
                            !(距離の逆数)*(物理量)!!(距離) 近いほど影響が大きい
                            Moment = Moment + InverseDistance * CopyOfUCC%ConservedQuantity(:,iAdjacentCell,1,1)
                        end if
                    end do
                end do

                UCC%ConservedQuantity(:,iCell,1,1) = Moment/InverseDistanceSum

            else !interpolate from global cell

                call InterpolateFromGlobal

            end if

        else

            call InterpolateFromGlobal

        end if

    end do

return
contains

    subroutine InterpolateFromGlobal
    implicit none

        UCC%ConservedQuantity(1:3,iCell,1,1) = &
        &   CC%ConservedQuantity(1:3,OSG%RelatedU2S(iCell,1),OSG%RelatedU2S(iCell,2),1)

        UCC%ConservedQuantity(5,iCell,1,1) = &
        &   CC%ConservedQuantity(4,OSG%RelatedU2S(iCell,1),OSG%RelatedU2S(iCell,2),1)

    return
    end subroutine InterpolateFromGlobal
end subroutine OInterpolateLocalGrid
