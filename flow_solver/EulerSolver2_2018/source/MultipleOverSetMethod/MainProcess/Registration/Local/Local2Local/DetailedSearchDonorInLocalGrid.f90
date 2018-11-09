!***********************************/
!	Name:内挿targetセルをドナー候補の座標系に変換する
!	Alias:
!	Description:iTime=0でドナー格子のN+1座標系へ写す，iTime=1でN座標系へ
!	Type:
!	Input:
!	Output:
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.12.05
!	Update:
!	Other:
!***********************************/
    subroutine DetailedSearchDonorInLocalGrid(D_NSR,D_UG,iTime,iUCell,CUV,T_MG,iFindDonor)
    use StructVar_Mod
    implicit none
    type(NextStepRelation), intent(in) :: D_NSR
    type(UnstructuredGrid), intent(in) :: D_UG
    integer, intent(in) :: iTime, iUCell !=TempS2U%Cell

    type(CountUpVertex), intent(inout) :: CUV
    type(MoveGrid), intent(inout) :: T_MG
    integer, intent(inout) :: iFindDonor
    integer :: BoundaryBoxMethod

    !第1：距離の２乗和と面積和
    !接触判定を行う相手セルをCUVに登録
        !call InputCUV(2,iTime,iUCell,D_UG,D_NSR,CUV)

        if(BoundaryBoxMethod(CUV) == 1) then !補助格子に登録された要素がtarget要素と接触してそうな距離にある(=極端に離れていない)とき1
        !第2：詳細な接触判定と，共有部分の頂点探索
        !write(6,*) CUV%iTargetNumber(2)
        !if(CUV%TargetCoords(3,CUV%iTargetNumber(2)) > 10.0d0**6) then
        !    write(6,*) CUV%TargetCoords(3,CUV%iTargetNumber(2))
        !    stop
        !end if
            call OCountUpVertexFromPoint(CUV) !内点の判定
            call OCountUpVertexFromEdge(CUV) !辺の交差判定

            if(CUV%iTotalVertexNum > 2) then !共有部分の頂点が3つ以上ある =点列が求積可能な多角形をなしている時
                call OSortVertexByAngle(CUV) !共有面積の算出

                if(CUV%SharedArea > 0.0d0) then !共有面積を計算できた格子の番号と共有面積を単連結リストに登録
                    call AddSharedDonorCellNumber(D_UG,T_MG,CUV)
                    iFindDonor = 1
                end if

            end if

        end if

    return
    end subroutine DetailedSearchDonorInLocalGrid
