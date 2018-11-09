!***********************************/
!	Name:前処理のみで使う構造型の定義モジュール
!	Alias:StructVar_OnlyUsePreProcess
!	Description:凸包を作るときのみに用いる
!	Type:UnstructuredGrid
!	Input:UG%GM%Dimension,UG%xxx%Point,UG%GM%Width (xxxは使用する種類のセルデータ)
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.20
!	Update:-
!	Other:
!***********************************/
 module StructVar_OnlyUsePreProcess
implicit none
    type MakeConvexHull
        integer, allocatable :: OutlineCell(:) !0:InternalCell, 1:OutlineCell
    end type
 end module StructVar_OnlyUsePreProcess
