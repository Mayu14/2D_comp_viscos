!***********************************/
!	Name:エラーメッセージを吐くプログラム
!	Alias:ErrorMessage
!	Description:プログラムのいろいろな場所で不正な値を確認するとcallで呼ばれて飛んでくる，場所ごとに異なる引数を与えてデバッグを助けてくれるはず...だったのだが...
!	Type:
!	Input:iNumberOfError
!	Output:
!	Note:エラーメッセージ吐いてプログラムを止める．他の変数を参照できないため，エラー原因を探れない．
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.07
!	Other:
!***********************************/
subroutine ErrorMessage(iNumberOfError)
    implicit none
    integer, intent(in) :: iNumberOfError

    if(iNumberOfError == 1) then
        write(6,*) "density is 0 in Conserve2Primitive"
        stop

    else if(iNumberOfError == 2) then
        write(6,*) "Pressure is 0 in Conserve2Primitive"
        stop
    else if(iNumberOfError == 3) then
        write(6,*) "1D Flux become NaN"
        stop
    else if(iNumberOfError == 4) then
        write(6,*) "2D Flux become NaN"
        stop
    else if(iNumberOfError == 5) then
        write(6,*) "BoundaryCondition is not Apllied."
        stop
    else if(iNumberOfError == 6) then
        write(6,*) "ConservedQuantity is NaN"
        stop
    else if(iNumberOfError == 7) then
         write(6,*) "ConservedQuantity become below 0"
    end if

return
end subroutine ErrorMessage
