!***********************************/
!	Name:非構造格子における共有辺を除去して辺番号を割り当てたり，格子の隣接関係を与えたりするプログラム（クソ）
!	Alias:UMakeEdgeNumber
!	Description:とりあえずやってることの概要は内部サブルーチンの上に注記で書く，絶対書き直す
!	Type:UnstructuredGrid
!	Input:UG%CD%Point,UG%GI%RealCells (xxxは使用する種類のセルデータ)
!	Output:UG%Line,UG%Tri
!	Note:クソクソアンドクソ(KKK)
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:2018.11.02
!	Other:現時点でハイブリッド格子・3次元格子への対応は未定(準備すらない)
!***********************************/
subroutine UMakeEdgeNumber(UG) !読み辛いことこの上ないがこれ以上見たくもないので修正はしない予定(不具合が出た場合はその限りではない)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG

    integer :: iOtherEdge,iExit
    integer, allocatable :: tmpC2C(:,:), tmpC2E(:,:) !Cell2Point(:,:) = UG%Tri%Point
    integer, allocatable :: tmpE2C(:,:,:), tmpE2P(:,:)
    integer, allocatable :: CopyE2P(:,:)
    integer, allocatable :: CheckBound(:)
    integer, allocatable :: ReplaceE(:)
    integer, allocatable :: StuffE(:)
    integer :: iDeleteNum,iTmpTotalEdge,iNewEdge,iBoundNum

    call Preprocessing
    call SortE2P
    call CheckOverlap
    call StuffAndOrder
    !call BoundaryTestOutput
    deallocate(tmpC2C,tmpC2E,tmpE2C,tmpE2P,CopyE2P,CheckBound,ReplaceE,StuffE)
return

contains
!辺の数を予測される最大の数で仮定して配列を取得(仮定：全セルが共有界面を持たない)
    subroutine Preprocessing
    implicit none
        iTmpTotalEdge=UG%GI%RealCells*3
        allocate(tmpC2C(UG%GI%RealCells,3),tmpC2E(UG%GI%RealCells,3),tmpE2C(iTmpTotalEdge,2,2),tmpE2P(iTmpTotalEdge,2))
        tmpC2C=0
        tmpE2C=0
        tmpE2P=0

        UG%Quad%Total=0 !今のところ三角形要素以外に対応する気ゼロ
        iTC => iTriangleCell
        iLED => iLocalEdge
        iED => iTriangleCell

        do iTriangleCell=1, UG%GI%RealCells !すべての三角形セルについて
            do iLocalEdge = 1, 3
                tmpC2E(iTC,iLED) = 3*(iTC-1) + iLED !借置き界面番号は3*(要素番号-1)+x (x=1~3) !要素番号iの全体辺番号を3(i-1)(i=1,3)で確保
                tmpE2C(tmpC2E(iTC,iLED),1,1) = iTC !とりあえず全部の界面は表として接していると仮定する
                tmpE2C(tmpC2E(iTC,iLED),1,2) = iLED !全体界面番号とセル番号から局所界面番号を返すための準備
            end do
            tmpE2P(tmpC2E(iTC,1),1) = UG%Tri%Point(iTC,2) !点番号を借置き界面番号と対応させる
            tmpE2P(tmpC2E(iTC,1),2) = UG%Tri%Point(iTC,3) !局所界面番号x=iの界面は局所点番号jとkの点から構成される(i,j,kは相異なる)
            tmpE2P(tmpC2E(iTC,2),1) = UG%Tri%Point(iTC,3)
            tmpE2P(tmpC2E(iTC,2),2) = UG%Tri%Point(iTC,1)
            tmpE2P(tmpC2E(iTC,3),1) = UG%Tri%Point(iTC,1)
            tmpE2P(tmpC2E(iTC,3),2) = UG%Tri%Point(iTC,2)
        end do
        !この時点でC2EとE2C,E2Pのデータ自体は完成しているからあとは間引くだけ
        return
    end subroutine Preprocessing

!___________________________________________\
!界面番号から点番号を与えるE2Pを複製し，保有する点番号を小さい順に並べ替える(あとで検索しやすくするため)
!ただし，E2Pにおける点の順番はvtkを再構成するために必要不可欠であるため，複製のみを編集し，使い捨てにする
    subroutine SortE2P
    implicit none
    !UG%Tri%Edgeの複製
        allocate(CopyE2P(iTmpTotalEdge,2))
        CopyE2P = tmpE2P

        iED => iEdge
        do iEdge=1, iTmpTotalEdge !点番号を小さい順に並び換え
            if(tmpE2P(iED,1) > tmpE2P(iED,2)) then
                CopyE2P(iED,1) = tmpE2P(iED,2)
                CopyE2P(iED,2) = tmpE2P(iED,1)
            end if
        end do
        return
    end subroutine SortE2P

!_________________________________\
!SortE2Pで並べ替えたデータを元に，共有界面を検索する
!共有界面が現れたときは下のif文内のルールで処理する
    subroutine CheckOverlap
    implicit none
        !異なる界面番号で重複定義されているセルを特定する
        allocate(ReplaceE(iTmpTotalEdge),CheckBound(iTmpTotalEdge))
        ReplaceE=0
        CheckBound = 0
        iDeleteNum=0
        iBoundNum = 0
        iED => iEdge

        do iEdge=1, iTmpTotalEdge !すべての借置界面番号について
            iExit = 0 !(終了条件を満たすとループ出る)
            if(ReplaceE(iEdge) /= 0) cycle
            if(iEdge == iTmpTotalEdge) then ! その界面が最後の界面であり，削除対象になっていない場合
                call MakeVirtualCell !その界面は外部境界である

            else    ! その界面が最後の界面でなければ
                do iOtherEdge=iEdge+1, iTmpTotalEdge         !他の界面と重複しないような組み合わせを取り，
                    if((ReplaceE(iOtherEdge) /= 0) .and. (iOtherEdge /= iTmpTotalEdge)) cycle       !選択した界面の除外がまだ決定していない&その界面が最後の界面じゃない場合は
                    if(CopyE2P(iED,1) == CopyE2P(iOtherEdge,1)) then    !それらの界面を構成するすべての点の点番号を比較して
                        if(CopyE2P(iED,2) == CopyE2P(iOtherEdge,2)) then !全部が一致していたならばそれらは同じ界面を指しているから
                            tmpC2C(tmpE2C(iED,1,1),TmpEdge2LocalEdge(iED)) = tmpE2C(iOtherEdge,1,1)  !隣接セル番号を結びつけるUG%Tri%Cellを作成し
                            tmpC2C(tmpE2C(iOtherEdge,1,1),TmpEdge2LocalEdge(iOtherEdge)) = tmpE2C(iED,1,1) !界面の両側のセル番号を記録
                            iDeleteNum = iDeleteNum+1 !重複界面のカウントを行う
                            if(tmpE2C(iED,1,1) < tmpE2C(iOtherEdge,1,1)) then !元のセル番号が小さな側を残し，
                                tmpE2C(iED,2,:) = tmpE2C(iOtherEdge,1,:) !元のセル番号が大きな側は面の裏側として記録する
                                ReplaceE(iED) = iED !置換用の番号として，残る方の界面は自身の番号を，
                                ReplaceE(iOtherEdge) = iED !置換される側の界面は置換後の番号をそれぞれ記録する
                            else
                                tmpE2C(iOtherEdge,2,:) = tmpE2C(iED,1,:)
                                ReplaceE(iED) = iOtherEdge
                                ReplaceE(iOtherEdge) = iOtherEdge
                            end if
                            iExit = 1 !このiEdgeについては終了条件を満たしたのでiExitを1にする
                        end if
                    end if

                    if(iOtherEdge == iTmpTotalEdge) then !最後のセルまで走査が到達してしまったとき
                        if(iExit == 0) then !まだ共有界面が見つかっていないならば
                            call MakeVirtualCell !その界面は外部境界であるから
                        end if
                    end if
                    if(iExit == 1) exit
                end do
            end if
        end do
!この時点ですべての面について走査が終わる
    UG%GI%Edges = iTmpTotalEdge-iDeleteNum !結局残ることになった界面は元の界面数-重複界面数
    if(UG%GM%Dimension == 2) UG%Line%Total = UG%GI%Edges
    UG%VC%Total = iBoundNum !外部境界のセル数
    UG%GI%AllCells = UG%GI%RealCells + UG%VC%Total
return
end subroutine CheckOverlap

!____________________________________________tmpE2C(i________
!CheckOverlapにて共有界面を存在持たないセルがあったときに呼び出される
!共有界面を持たない=外部境界に面しているセルであるため仮想格子を作成し，要素番号から相互に呼び出せるようにする
subroutine MakeVirtualCell
implicit none
    CheckBound(iED) = 1
    iBoundNum = iBoundNum+1  !界面は外部境界に面しているから，境界要素数をカウントし
    tmpE2C(iED,2,1) = UG%GI%RealCells + iBoundNum !界面の裏側として仮想格子の番号を与え
    tmpE2C(iED,1,2) = TmpEdge2LocalEdge(iED) !局所面番号を記録し，
    tmpE2C(iED,2,2) = 1 !仮想格子側の局所界面番号を1とする
    tmpC2C(tmpE2C(iED,1,1),TmpEdge2LocalEdge(iED)) = UG%GI%RealCells + iBoundNum !界面の属する実セルの隣接格子として仮想格子の番号を与え

    ReplaceE(iED) = iED !番号を置き換える際に，自身の番号で置き換えられるようにする
return
end subroutine MakeVirtualCell

!____________________________________________________
!重複しているとリストアップしたデータを元にStuff(詰めて)Order(並べる)
!先のCheckOverlapにて重複界面とみなされた仮置辺番号をReplace()に入れると変換後の辺番号が出てうる
!これを元に，変換後の番号=自分の番号である界面にのみ新しい番号を与えてC2CやE2Pなどのデータをコピーする
subroutine StuffAndOrder
implicit none
!番号を伝並べ替える
    allocate(UG%Tri%Cell(UG%GI%RealCells,3),UG%Tri%Edge(UG%GI%RealCells,3))
    allocate(UG%Line%Cell(UG%Line%Total,2,2),UG%Line%Point(UG%Line%Total,2)) !最終版としてジャストサイズの配列を用意する
    allocate(UG%VC%Cell(UG%GI%RealCells+1:UG%GI%AllCells,2),UG%VC%Edge(UG%GI%RealCells+1:UG%GI%AllCells))
    allocate(StuffE(iTmpTotalEdge))
    StuffE = 0

    iBoundNum = 0
    iNewEdge = 0
    do iEdge = 1, iTmpTotalEdge
        if(ReplaceE(iED) /= iED) cycle !置き換えられることが決まっているセルは無視する
        iNewEdge = iNewEdge + 1 !置換されずに残るセルがあったとき，1から詰めた新しい番号をiNewEdgeとして与えて
        StuffE(ReplaceE(iED)) = iNewEdge !詰めるための配列StuffEのうちReplaceEが存在する部分だけiNewEdgeを入れる(置き換えられることが確定している部分ではStuffE=0)
        !移し替え
        UG%Line%Cell(iNewEdge,:,:) = tmpE2C(iED,:,:)
        UG%Line%Point(iNewEdge,:) = tmpE2P(iED,:)

        if(CheckBound(iED) == 1) then
            iBoundNum = iBoundNum+1
            UG%Tri%Cell(UG%Line%Cell(iNewEdge,1,1),TmpEdge2LocalEdge(iED)) = UG%GI%RealCells + iBoundNum !実セルから仮想セルへ
            UG%VC%Edge(UG%GI%RealCells + iBoundNum) = iNewEdge !仮想セルから実セルへ
            UG%VC%Cell(UG%GI%RealCells + iBoundNum,:) = UG%Line%Cell(iNewEdge,1,:)
        end if
    end do

    do iCell=1, UG%GI%RealCells
        do iLocalEdge=1,3
            UG%Tri%Cell(iCell,iLED) = tmpC2C(iCell,iLED)
            UG%Tri%Edge(iCell,iLED) = StuffE(ReplaceE(tmpC2E(iCell,iLocalEdge)))
        end do
    end do
return
end subroutine StuffAndOrder

!仮置き辺番号における局所辺番号は3要素iに対して3(i-1)+1,3(i-1)+2,3(i-1)+3であったことから計算で求めることができる
!毎度毎度mod(iEdge+2,3)+1と書くと訳が分からなくなりそうだから関数を作ったが，こっちの名前の方が長くて分かりづらさが増した可能性
!記述ミスが起きないからいいけども...
function TmpEdge2LocalEdge(iEdge) result(iLocalEdge)
implicit none
    integer, intent(in) :: iEdge
    integer :: iLocalEdge

    iLocalEdge = mod(iEdge+2,3)+1
    return
end function TmpEdge2LocalEdge

!ためしに境界の格子のみを1,それ以外を0として出力してみるテスト用プログラム(ここに書くべきではないやろなぁ...)
subroutine BoundaryTestOutput
implicit none
double precision, allocatable :: iCounter4Bound(:)
   open(unit=1,file="GridTest2.vtk",status='unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,"('TestGrid')")
        write(1,"('ASCII')")
        write(1,"('DATASET UNSTRUCTURED_GRID')")
        write(1,"('POINTS ',(1x,i3),' double')") UG%GI%Points

        do iLoop=1, UG%GI%Points
            write(1,"(3(1x,f22.17))") UG%CD%Point(iLoop,1),UG%CD%Point(iLoop,2),UG%CD%Point(iLoop,3)
        end do
        write(1,*) ""

        write(1,"('CELLS ',2(1x,i4))") UG%GI%RealCells,UG%GI%RealCells*4
        do iLoop=1, UG%GI%RealCells
            write(1,"(4(1x,i4))") 3,UG%Tri%Point(iLoop,1)-1,UG%Tri%Point(iLoop,2)-1,UG%Tri%Point(iLoop,3)-1
        end do
        write(1,*) ""

        write(1,"('CELL_TYPES ',(1x,i4))") UG%GI%RealCells
        do iLoop=1, UG%GI%RealCells
            write(1,"(1x,i1)") 5
        end do
        write(1,*) ""

        write(1,"('POINT_DATA ',i7)") UG%GI%Points
        write(1,"('SCALARS point float')")
        write(1,"('LOOKUP_TABLE default')")
        allocate(iCounter4Bound(0:UG%GI%Points-1))
        iCounter4Bound = 0.0d0
        do iCell=1, UG%GI%RealCells
            do iLocalEdge=1,3
                if(UG%Tri%Cell(iCell,iLocalEdge) > UG%GI%RealCells) then
                    iCounter4Bound(UG%Line%Point(UG%Tri%Edge(iCell,iLocalEdge),1)) = 1.0d0
                    iCounter4Bound(UG%Line%Point(UG%Tri%Edge(iCell,iLocalEdge),2)) = 1.0d0
                end if
            end do
        end do
        do iPoint=0, UG%GI%Points-1
            write(1,"(1x,f22.17)") iCounter4Bound(iPoint)
        end do
    close(1)
    return
    end subroutine BoundaryTestOutput
end subroutine UMakeEdgeNumber
