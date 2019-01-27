!***********************************/
!	Name:衝突三角形の作る多角形における頂点並べ替えプログラム
!	Alias:SortVertexByAngle
!	Description:同時に共有面積の計算も実行しているため注意
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.07
!	Update:
!	Other:not parallelize
!***********************************/
    subroutine OSortVertexByAngle(CUV)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    type(CountUpVertex), intent(inout) :: CUV

    integer :: iPoint
    integer :: iYMinNumber
    integer :: iSort1, iSort2, iTmp, iExit
    integer, allocatable :: AccessOrder(:)
    double precision, allocatable :: AngleOfPointVector(:)


    allocate(AccessOrder(CUV%iTotalVertexNum+1))
    allocate(AngleOfPointVector(CUV%iTotalVertexNum))
    !minlocでy方向最小座標を持つ配列番号を特定
    !全成分から特定した配列番号の要素を引く
    iYMinNumber = minloc(CUV%AdditionVertex(1:6,2),1)
    do iTmp=1, CUV%iTotalVertexNum
        CUV%AdditionVertex(iTmp,:) = CUV%AdditionVertex(iTmp,:) - CUV%AdditionVertex(iYMinNumber,:)
    end do

    !xy平面における角度の値を求める
    do iPoint=1,CUV%iTotalVertexNum
        AngleOfPointVector(iPoint) = atan2(CUV%AdditionVertex(iPoint,2), CUV%AdditionVertex(iPoint,1))
    end do

    !角度が小さい順に点の順序の並び換えを行う...というより点の処理順序を決定する(内挿用点番号1→4番的な)
    !データ数が十分少ないため，挿入ソートを利用
    !内挿用点番号iとjの角度を比較し，
    !i>jならj-1として再びiとjの比較
    !j-1<iとなったところでiとjをスワップ
    !以降繰り返し
!挿入ソート
    do iTmp=1, CUV%iTotalVertexNum+1
        AccessOrder(iTmp)=iTmp
    end do

    do iSort1=2, CUV%iTotalVertexNum !2番目から最後の要素まで
        iTmp=iSort1 !現在の並べ替え対象をiTmp番とする
        if(AngleOfPointVector(iSort1-1) > AngleOfPointVector(iTmp)) then !iSort-1番とiSort番の値を比較し，iSort番の方が小さければ入れ替え発生
            iExit = 0 !離脱判定の初期化
            do iSort2 = iTmp, 2, -1 !iSort番目までの値を順々に右へずらし，A<iSortとなる初めての位置にデータを挿入する
                if(AngleOfPointVector(iSort2-1) > AngleOfPointVector(iTmp)) then
                    AccessOrder(iSort2) = AccessOrder(iSort2-1)
                else
                    AccessOrder(iSort2) = AccessOrder(iSort1)
                    iExit = 1
                end if
                if(iExit == 1) exit
            end do
        end if
    end do

    AccessOrder(CUV%iTotalVertexNum+1) = AccessOrder(1)

    CUV%SharedArea = 0.0d0
    do iPoint=1, CUV%iTotalVertexNum
        CUV%SharedArea = CUV%SharedArea &
            &   + 0.5d0*AbsCrossProduct(CUV%AdditionVertex(AccessOrder(iPoint),:),CUV%AdditionVertex(AccessOrder(iPoint+1),:))
    end do

    return
    end subroutine OSortVertexByAngle
