!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:OExclusionObviousMistake
!	Description:明らかに遠いセルを補間対象の検討から除外する
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.12.05
!	Update:
!	Other:not parallelize
!***********************************/
    subroutine OExclusionObviousMistake(UG,MG,RCN,EOM)
    use StructVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(MoveGrid), intent(in) :: MG
    type(RelatedCellNumber), intent(in) :: RCN
    type(ExclusionObviousMistake), intent(inout) :: EOM

    EOM%SquareArea = dot_product(MG%RC%Cell(RCN%iPartnerCell,1:3)-MG%NSR%C%NextCoordsInN(RCN%iMyCell,1:3), &
                            &   MG%RC%Cell(RCN%iPartnerCell,1:3)-MG%NSR%C%NextCoordsInN(RCN%iMyCell,1:3))

    EOM%AreaSum = UG%GM%Area(RCN%iMyCell) + UG%GM%Area(RCN%iPartnerCell)

    if(EOM%SquareArea > EOM%AreaSum) then
        EOM%iTrueOrFalse = 0
    else
        EOM%iTrueOrFalse = 1
    end if

    return
    end subroutine OExclusionObviousMistake
