!***********************************/
!	Name:vtkの平面切り出しプログラム
!	Alias:Generate2DVTKFrom3D
!	Description:次元のvtkファイルから2次元の平面を切り出すプログラム
!	Type:d
!	Input:z平面を切り出したいvtkファイル
!	Output:z平面のみになったぺらっぺらのvtkファイル
!	Note:0番目の点のz座標を元に平面を切り出します．
!	Author:Akitaka Toyota
!	Date:2017.11.02
!	Update:
!	Other:四面体のみで構成された平板以外だとうまくいかないと思います．(z平面以外を切ったりするオプションを付ける予定は)ないです
!***********************************/
program Generate2DVTKFrom3D
    implicit none
    character(len=64) :: cReadVtkName,cWriteVtkName,cAnnotate
    integer(kind=8) :: iTotalPoints,iTotalCells,iTotalCellData
    integer(kind=8) :: iNewTotalPoints,iNewTotalCells,iNewTotalCellData
    integer(kind=8) :: iLoop,iTmp,iTmp2
    integer(kind=8) :: iTriNum,iTetraNum
    integer(kind=8), allocatable :: iOld2New(:) !点番号,古い点番号→新しい点番号 !0を返す=消去済み
    double precision, allocatable :: Coordinate(:,:) !点番号,xyz

    integer(8), allocatable :: iCells(:,:),iSwapCells(:,:),iNewCells(:,:) !要素番号，局所点番号で大域点番号を返す

    write(6,*) "Please Input the filename of 3D-Vtk what you want to convert 2D."
    read(5,*) cReadVtkName
    !cReadVtkName = "square2027.vtk"

    write(6,*) "Please Input the filename of 2D-Vtk what you want to generate."
    read(5,*) cWriteVtkName
    !cWriteVtkName = "toyotaMK2.vtk"


    open(unit = 1, file =trim(adjustl(cReadVtkName)), status = 'unknown')
        read(1,*) cAnnotate !header
        read(1,*) cAnnotate !Name
        read(1,*) cAnnotate !Ascii
        read(1,*) cAnnotate !DATASET
        read(1,*) cAnnotate, iTotalPoints !点の数

        allocate(iOld2New(iTotalPoints),Coordinate(iTotalPoints,3))
        iOld2New = -1 !iOld2Newに-1が残っている点を廃棄する
        Coordinate = 10.0d0**30 !変な値が入ってる時に判別するため

            read(1,*) Coordinate(1,1),Coordinate(1,2),Coordinate(1,3) !1点目の確保／z座標基準の確保
            iOld2New(1) = 1 !vtkの点番号0を1として扱う
            iTmp = 2
!座標情報読み込み
        do iLoop = 2, iTotalPoints
            read(1,*) Coordinate(iTmp,1),Coordinate(iTmp,2),Coordinate(iTmp,3) !2点目〜最終点までの座標を読み取る

            if(Coordinate(iTmp,3) == Coordinate(1,3)) then !1点目とz座標が同一であるときiTmpを1加算
                iOld2New(iLoop) = iTmp !iLoop番目の点の点番号を記録
                iTmp = iTmp+1 !削除しない点の総数の数え上げ
                !iNew2Old(iTmp) = iLoop
            end if !1点目とz座標が異なる場合はiTmpが増えない = 次のループで同じiTmpに値が上書きされる
        end do
        iNewTotalPoints = iTmp-1 !点の総数を記録

!セルデータ読み込み
        read(1,*) cAnnotate,iTotalCells,iTotalCellData !総セル数,各セルの頂点を構成する点の番号データの総数

        iTetraNum = -4*iTotalCells + iTotalCellData !四面体の数を求める
        iTriNum = 5*iTotalCells - iTotalCellData !三角形面の数を求める(連立方程式を直接法で解いてる)

        allocate(iCells(iTriNum,4),iSwapCells(iTriNum,4),iNewCells(iTriNum,4))
        iCells = -1 !セル番号は0始まりのため-1で初期化する
        iSwapCells = -35
        iNewCells = 0 !iCellsを+1してから格納するため合わせる

        do iLoop=1, iTetraNum
            read(1,*) cAnnotate,cAnnotate,cAnnotate,cAnnotate,cAnnotate
        end do
        do iLoop = iTetraNum+1, iTotalCells
            read(1,*) iCells(iLoop-iTetraNum,1),iCells(iLoop-iTetraNum,2),iCells(iLoop-iTetraNum,3),iCells(iLoop-iTetraNum,4) !三角形面を構成する(旧)点番号を記録
        end do
    close(1)

    !vtk用の点番号は0~Max-1までなので番号を1からに振り直す．
    iCells = iCells + 1

!セル番号の入れ替え
    iTmp = 0 !セル数
    iTmp2 = 0 !データ数
    do iLoop=1,iTriNum !全セルについて
        if(iOld2New(iCells(iLoop,2)) == -1) then
        else if(iOld2New(iCells(iLoop,3)) == -1) then
        else if(iOld2New(iCells(iLoop,4)) == -1) then
        else !面内すべての点が-1を返さないとき"
            iTmp = iTmp + 1
            iTmp2 = iTmp2 + 4
            iSwapCells(iTmp,1) = 4 !あとで-1するため
            iSwapCells(iTmp,2) = iOld2New(iCells(iLoop,2)) !セルを構成する点の番号を新たな番号へと変更する
            iSwapCells(iTmp,3) = iOld2New(iCells(iLoop,3))
            iSwapCells(iTmp,4) = iOld2New(iCells(iLoop,4))
        end if
    end do
!この時点でiTmpには消去対象ではない(残った)セルの数が入っている
    iNewTotalCells = iTmp
    iNewTotalCellData = iTmp2

    iSwapCells = iSwapCells - 1

!VTKを記録するまえに点番号を0始まりへ戻す
    iNewCells = iNewCells - 1

    open(unit = 1, file =trim(adjustl(cWriteVtkName)), status = 'unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,"('2D-Unstructured')")
        write(1,"('ASCII')")
        write(1,"('DATASET UNSTRUCTURED_GRID')")
        write(1,"('POINTS ',(1x,i14),' double')") iNewTotalPoints

        do iLoop=1, iNewTotalPoints
            write(1,"(3(1x,f22.17))") Coordinate(iLoop,1),Coordinate(iLoop,2),0.0d0
        end do
        write(1,*) ""

        write(1,"('CELLS ',2(1x,i14))") iNewTotalCells,iNewTotalCellData
        do iLoop=1, iNewTotalCells
                write(1,"(4(1x,i14))") iSwapCells(iLoop,1),iSwapCells(iLoop,2),iSwapCells(iLoop,3),iSwapCells(iLoop,4)
        end do
        write(1,*) ""

        write(1,"('CELL_TYPES ',(1x,i14))") iNewTotalCells
        do iLoop=1, iNewTotalCells
            write(1,"(1x,i1)") 5
        end do
    close(1)

stop
end program Generate2DVTKFrom3D

