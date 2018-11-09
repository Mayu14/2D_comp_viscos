!***********************************/
!	Name:Target格子をグローバル格子に登録する
!	Alias:RegisterIntoGlobal
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:現状，面積の代わりにセル中心間距離の逆数を代入しているところに注意．さらに言えば，複数セル間に跨る要素はすべて中心が属する側に所属するものとして見なしている
!	Author:Akitaka Toyota
!	Date:2017.12.26
!	Update:
!	Other:
!***********************************/
    subroutine RegisterIntoGlobal(iGrid,GG,iGlobalCellNumber,Geom,CUV,iLoop,RIG)
    use StructVar_Mod
    use FrequentOperation
    use RegistVar_Mod4OMP
    implicit none

    integer, intent(in) :: iGrid
    type(GlobalGrid), intent(inout) :: GG
    integer, intent(inout) :: iGlobalCellNumber !なんと後で再利用する(内挿対象が見つからなくてG格子からの内蔵を行う際の登録部分で使う) !恐らく廃止
    type(Geometry), intent(inout) :: Geom
    type(CountUpVertex), intent(in) :: CUV

    integer, intent(inout) :: iLoop
    type(RegisterIntoGlobalWithOMP), intent(inout) :: RIG

    !type(DonorElement), pointer :: RIG%Temp
        !構造格子のデータの一部をGG%GASGにコピーして使う，というか前処理の段階でグローバル格子に対応するGG%GASGを作る

        !(非構造格子の座標r,Auxiliary格子の座標r0，構造格子負端からのセル数n，格子幅dxについて
        !r-r0 > n*dxなる関係があることから，nについて解くと以下のような式になる(※2次元用で書いたためzは常に固定値になる)

    !do iLoop = 0, 3 !セル中心+セル頂点に関するループ
    iLoop = 0
        !RIG%GlobalGridNum(1:2) = int((NSR%C%NextCoordsInG(iCell,1:2) - GG%GASG%Bound(1,1:2))/(GG%GASG%Width(1:2)))+1 !補助格子におけるセル中心位置
        call GetLocalNumber2D(CUV%TargetCoords(iLoop,1:2),GG%GASG%Width(1:2),GG%GASG%Bound(1,1:2),RIG%GlobalGridNum(1:2))
        RIG%GlobalGridNum(3) = 1

        !取扱いの利便性を上げるためグローバルセル番号(詳細は↓)に変換する．
        !(i,j,k) <=>!2D:(j-1)*(Lx)+i !3D:(k-1)*(Ly*Lx)+(j-1)*(jmax-1)+i

        iGlobalCellNumber = Local2Global(RIG%GlobalGridNum(1),RIG%GlobalGridNum(2),GG%GASG%CellNumber(1))

        !Global格子がLocal格子から参照された回数を記録 !このサブルーチンの外側で初期化しておく

        GG%GASG%ContentsAS2U(iGlobalCellNumber) = GG%GASG%ContentsAS2U(iGlobalCellNumber) + 1 !この場合Global 2 Localなので注意

        if(GG%GASG%ContentsAS2U(iGlobalCellNumber) == 1) then
            nullify(GG%GDE(iGlobalCellNumber)%DE) !リストの初期化
            Geom%Interpolated(RIG%GlobalGridNum(1),RIG%GlobalGridNum(2),1) = 1 !時間積分からの保護(globalを補間する格子はすべてn+1段階)
        end if
        !最後に構造格子に所属する非構造格子への対応を記録する

        !セル中心間距離
        !G格子のセル中心位置
        call GetCellCenterCoords(RIG%GlobalGridNum(1),RIG%GlobalGridNum(2),GG%GASG%Bound(1,1:2),GG%GASG%Width(1:2),RIG%CellCenterCoords(1:2))
        RIG%CellCenterCoords(3) = 0.0d0

        allocate(RIG%Temp) !新規ノードの確保
            RIG%Temp%Cell = CUV%iTargetNumber(0) !ノードへの記録
            RIG%Temp%Grid = iGrid
            RIG%Temp%SharedArea = 1.0d0/(sqrt(dot_product(&
                &   CUV%TargetCoords(iLoop,1:3)-RIG%CellCenterCoords(1:3),&
                &   CUV%TargetCoords(iLoop,1:3)-RIG%CellCenterCoords(1:3))))

            !if(iLoop == 0) then !セル中心のとき
            !    RIG%Temp%SharedArea = 0.4d0 * RIG%Temp%SharedArea !(セル中心が属する=2/5)
            !else !セル頂点のとき
            !    RIG%Temp%SharedArea = 0.2d0 * RIG%Temp%SharedArea !(セル中心が属する=1/5)
            !end if

            RIG%Temp%Next => GG%GDE(iGlobalCellNumber)%DE !前のノードへのリンク作成
            GG%GDE(iGlobalCellNumber)%DE => RIG%Temp !データの保存

    !end do

    return
    end subroutine RegisterIntoGlobal
