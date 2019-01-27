!***********************************/
!	Name:非構造格子の計算領域データを確認するプログラム
!	Alias:UTestVTK
!	Description:MAYUGridから読み出したデータの格子点情報を画面出力するプログラム
!	Type:UG
!	Input:UConf,UG
!	Output:UG
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:2017.11.10
!	Other:
!***********************************/
subroutine UTestVTK(UG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

        write(6,*) UG%Tri%Total
        do iCell=1, UG%GI%RealCells
            write(6,*) UG%Tri%Point(iCell,1),UG%Tri%Point(iCell,2),UG%Tri%Point(iCell,3)
            write(6,*) ""
        end do
        stop

return
end subroutine UTestVTK
