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
subroutine FreeMemory(iTotalCell,MG)
    use StructVar_Mod
    implicit none
    integer :: iTotalCell
    integer :: iCell
    type(MoveGrid), intent(inout) :: MG
    type(DonorElement), pointer :: DE
    type(DonorElement), pointer :: TempDE

    !write(6,*) "Free Pointer"
    !do iGrid = 1, iTotalGridNum外側で回す,グローバル格子は別途回す
        do iCell=1, iTotalCell
            DE => MG%ADE(iCell)%DE
            do while(associated(DE))
                TempDE => DE%Next
                    nullify(DE%Next) !?
                deallocate(DE)
                DE => TempDE
            end do
                nullify(MG%ADE(iCell)%DE) !?
        end do
    !end do

    return
end subroutine FreeMemory


