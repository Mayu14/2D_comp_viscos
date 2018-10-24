!***********************************/
!	Name:しばしば使いそうなループ変数の宣言省略用モジュール
!	Alias:LoopVar_Mod
!	Description:構造格子用，非構造格子用とりあえずもろもろ揃えてみた
!	Type:integer
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.06
!	Update:2017.11.10
!	Other:
!***********************************/
module LoopVar_Mod_mk2
implicit none
!ここで定義する変数はループ変数として以外の使い方はしない！
!ポインタは適宜「次元等による処理分け」「略記」の目的で使う
!その都合上，常にnullifyで初期化して，ポインターを指定し直して使う
    integer :: iLoop !普通のループ用変数

    integer, target :: iVariable !変数のループ用
    integer, pointer :: iVar

    integer, target :: iCenterX,iCenterY,iCenterZ !セル中心で定義される変数のループ用
    integer, pointer :: iCX,iCY,iCZ

    integer, target :: iFaceX,iFaceY,iFaceZ !セル界面で定義される変数のループ用
    integer, pointer :: iFX,iFY,iFZ

    integer, pointer :: iTX,iTY,iTZ !時間Time

    integer :: iXmax,iYmax,iZmax
    integer :: iSide
    integer :: iDim

    integer, target :: iPoint, iLocalPoint
    integer, pointer :: iPT,iLPT
    integer, target :: iEdge, iLocalEdge
    integer, pointer :: iED,iLED
    integer, target :: iCell
    integer, target :: iTriangleCell
    integer, pointer :: iTC

    integer :: iGrid

    integer, target :: iElement
    integer, pointer :: iELM

    integer :: iEndPoint1, iEndPoint2
    integer :: iVertex1, iVertex2, iVertex3
    integer :: iShareEdge !共有界面
    integer :: iAdjacentCell,iAdjacentEdge !隣接要素，隣接要素から見た共有界面の局所界面番号

    integer :: iFrontCell, iBackCell

end module LoopVar_Mod
