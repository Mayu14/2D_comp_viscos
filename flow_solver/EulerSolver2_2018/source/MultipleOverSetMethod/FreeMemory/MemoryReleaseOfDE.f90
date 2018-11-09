!***********************************/
!	Name:毎回計算しなおす必要がある部分で使用しているポインタの解放
!	Alias:FreeMemory
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:本来は単連結List専用のメソッドとして用意すべきだが，object的実装法に慣れていないためsubroutineとして書いている
!	Author:Akitaka Toyota
!	Date:2017.12.26
!	Update:2017.12.26
!	Other:
!***********************************/

subroutine MemoryReleaseOfGDE(GG,iTotalCellNumber)
    use StructVar_Mod
    implicit none

    type(GlobalGrid), intent(inout) :: GG
    integer, intent(in) :: iTotalCellNumber
    type(DonorElement), pointer :: DE
    type(DonorElement), pointer :: TempDE
    integer :: iCell

        do iCell = 1, iTotalCellNumber
            DE => GG%GDE(iCell)%DE
            do while(associated(DE))
                TempDE => DE%Next
                    nullify(DE%Next) !?
                deallocate(DE)
                DE => TempDE
            end do
                nullify(GG%GDE(iCell)%DE) !?
        end do
    return
end subroutine MemoryReleaseOfGDE


