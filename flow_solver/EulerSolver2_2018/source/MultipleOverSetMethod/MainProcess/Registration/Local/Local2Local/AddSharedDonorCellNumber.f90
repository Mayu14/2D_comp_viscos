!***********************************/
!	Name:共有面積が計算できた要素について，Donor格子番号，要素番号，共有面積情報を単連結リストに登録
!	Alias:AddSharedDonorCellNumber
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.23
!	Update:2017.01.06
!	Other:revise Initialize (nullify Pointer)
!***********************************/
subroutine AddSharedDonorCellNumber(D_UG,T_MG,CUV,ASDCN)
    use StructVar_Mod
    use RegistVar_Mod4OMP
    implicit none
    type(UnstructuredGrid), intent(in) :: D_UG !Donor
    type(MoveGrid), intent(inout) :: T_MG !Target
    type(CountUpVertex), intent(in) :: CUV

    type(AddSharedDonorCellNumberWithOMP), intent(inout) :: ASDCN

    !type(DonorElement), pointer :: ASDCN%TempDE

        allocate(ASDCN%TempDE) !新規ノードの確保
        ASDCN%TempDE%Cell = CUV%iDonorNumber(0) !DonorGrid's CellNumber !ノードへの記録
        ASDCN%TempDE%Grid = D_UG%Number !DonorGridNumber
        ASDCN%TempDE%SharedArea = CUV%SharedArea

        if(CUV%iTotalDonor == 1) then !hajimete
            nullify(ASDCN%TempDE%Next)
        else
            ASDCN%TempDE%Next => T_MG%ADE(CUV%iTargetNumber(0))%DE !前のノードへのリンク作成
        end if

        T_MG%ADE(CUV%iTargetNumber(0))%DE => ASDCN%TempDE !データの保存


return
end subroutine AddSharedDonorCellNumber
