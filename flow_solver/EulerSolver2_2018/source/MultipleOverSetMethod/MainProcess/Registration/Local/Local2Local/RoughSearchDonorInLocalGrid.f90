!***********************************/
!	Name:内挿targetセルをドナー候補の座標系に変換する
!	Alias:RoughSearchDonorInLocalGrid
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
    subroutine RoughSearchDonorInLocalGrid(iTimeCombination,D_MG,D_UG,T_MG,CUV,iFindDonor,iAuxiliary,RSDILG)
    use StructVar_Mod
    use RegistVar_Mod4OMP
    implicit none
    integer, intent(in) :: iTimeCombination
    type(MoveGrid), intent(in) :: D_MG
    type(UnstructuredGrid), intent(in) :: D_UG
    type(MoveGrid), intent(inout) :: T_MG
    type(CountUpVertex), intent(inout) :: CUV
    integer, intent(inout) :: iFindDonor
    integer, intent(inout) :: iAuxiliary
    type(RoughSearchDonorInLocalGridWithOMP), intent(inout) :: RSDILG

    !type(Structured2Unstructured), pointer :: RSDILG%TempS2U

        !要素の頂点座標・中心座標をグローバルに戻したうえでUG2(Donor)系に変換する
        call CellTransformationIntoDonor(1,CUV,D_MG%TM)
        !座標から所属する補助格子の番号を返す(補助格子の値が不正な場合は例外0を返す)
        !要素の面積と補助格子の面積から隣接9つのセル番号を取得する(補助格子の値が不正な場合は例外0を取得する)
        call MatchingWithDonorInAuxiliaryGrid(D_MG,CUV,RSDILG%MWDIAG)

        if(CUV%AuxiliaryCellNumber(5) /= 0) then !適切な補助格子が見つかった場合次のループが実行される，見つからなかった場合は次の格子へ
            !補助格子内の要素とターゲット要素との接触判定を行う

            do iAuxiliary = 1, 9
                if(CUV%AuxiliaryCellNumber(iAuxiliary) /= 0) then !適切な番号が割り振られた補助格子が存在するとき

                    RSDILG%TempS2U => D_MG%ASG%AS2U(CUV%AuxiliaryCellNumber(iAuxiliary))%S2U !保有要素すべてについてのループ
                    do while(associated(RSDILG%TempS2U))
                        !call DetailedSearchDonorInLocalGrid(D_MG%NSR,D_UG,iTimeCombination,RSDILG%TempS2U%Cell,CUV,T_MG,iFindDonor)

                        call SimpleSearchDonorInLocalGrid(D_MG%NSR,D_UG,iTimeCombination,RSDILG%TempS2U%Cell,CUV,T_MG,iFindDonor,D_MG,RSDILG%SSDIG)

                        RSDILG%TempS2U => RSDILG%TempS2U%Next

                    end do

                end if
            end do

        end if

    return
    end subroutine RoughSearchDonorInLocalGrid
