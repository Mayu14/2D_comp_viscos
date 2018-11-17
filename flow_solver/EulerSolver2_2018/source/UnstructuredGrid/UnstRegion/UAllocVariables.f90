!***********************************/
!	Name:非構造ソルバー用の配列確保ルーチン
!	Alias:UAllocVariables
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:MultipleOverset用に変更
!	Author:Akitaka Toyota
!	Date:2017.11.06?
!	Update:2018.02.03
!	Other:
!***********************************/
    subroutine UAllocVariables(UConf,UG,UCC,UCE)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(Configulation) :: UConf
    type(UnstructuredGrid) :: UG
    type(CellCenter) :: UCC
    type(CellEdge) :: UCE

        allocate(UG%Tri%Point(UG%Tri%Total,3)) !三角形要素の構成点番号は(要素数,3)
        allocate(UG%Tri%Edge(UG%Tri%Total,3))  !三角形要素の構成辺番号は(要素数,3)
        allocate(UG%Tri%Cell(UG%Tri%Total,3))  !三角形要素の隣接要素番号は(要素数,3)
        allocate(UG%Tri%Type(UG%Tri%Total)) !三角形要素の種類（与えるべき境界条件（内部境界条件を利用するため））

        allocate(UG%Tri%Belongs2Wall(UG%GI%RealCells)) !三角形要素に対する最近傍の物体表面辺番号    ! only-use Baldwin-Lomax model
        allocate(UG%Tri%Distance(UG%GI%RealCells)) !三角形要素に対する最近傍の物体表面辺番号からの界面垂直方向距離 ! only-use Baldwin-Lomax model

        allocate(UG%Line%Point(UG%Line%Total,2)) !辺要素の構成点番号は(要素数,2)
        allocate(UG%Line%Cell(UG%Line%Total,2,2)) !辺要素の隣接要素番号は(要素数,2,2) !最後の引数はセル番号と，相手セルから見た界面自身の局所界面番号の切替
        allocate(UG%VC%Cell(UG%GI%RealCells+1:UG%GI%AllCells,2)) !仮想要素の隣接要素番号は(要素数,2) !最後の引数は相手から見た共有辺の辺番号
        allocate(UG%VC%Edge(UG%GI%RealCells+1:UG%GI%AllCells))   !仮想要素の共有界面番号
        allocate(UG%CD%Point(UG%GI%Points,4)) !点の座標は(点の総数,3方向)
        allocate(UG%CD%Edge(UG%GI%Edges,3))   !界面の中心座標は(界面総数,3方向)
        allocate(UG%CD%Cell(UG%GI%AllCells,4)) !要素の中心座標は(要素数,3方向)
        allocate(UG%GM%Area(UG%GI%Edges))      !界面の面積は(界面の数)
        allocate(UG%GM%Bound(2,3))
        allocate(UG%GM%Volume(UG%GI%RealCells)) !要素の体積は(要素の数)
        allocate(UG%GM%Normal(UG%GI%Edges,3)) !界面の表向き法線ベクトルは(界面の数)
        allocate(UG%GM%Width(UG%GI%RealCells,3,3)) !セル中心から界面までの距離は(実格子の数,要素内の界面数,3方向)
        allocate(UG%GM%CellType(UG%GI%AllCells,1,1)) !セル番号を入れるとセルの条件境界が返る

        allocate(UG%GM%BC%InFlowVariable(5),UG%GM%BC%OutFlowVariable(5))
        allocate(UG%InscribedCircle(UG%GI%RealCells))
        if(UConf%UseMUSCL == 1) then!
            allocate(UG%GM%AverageWidth(UG%GI%RealCells)) !セル中心から界面までの平均距離は(実格子番号)だけ確保
        end if

        if(UConf%UseOverSet == 1) then
            allocate(UG%CH%PointNum(UG%CH%iTotal))
        end if

        if(UG%IO%iTotal /= 0) then
            allocate(UG%IO%PointNum(UG%IO%iTotal))
        end if

        iDim = UG%GM%Dimension
        allocate(UCC%ConservedQuantity(iDim+2,UG%GI%AllCells,1,1)) !保存変数は(変数の種類，要素番号,1,1)で確保
        allocate(UCC%PrimitiveVariable(iDim+2,UG%GI%AllCells,1,1)) !基礎変数は(変数の種類，要素番号,1,1)で確保
        allocate(UCE%NormalFluxDiff(iDim+2,1,1,1,UG%GI%Edges)) !法線方向流束は(変数の種類，1,1,1,大域面番号)で確保
        UCC%ConservedQuantity = 0.0d0
        UCC%PrimitiveVariable = 0.0d0
        UCE%NormalFluxDiff = 0.0d0

        if(UConf%UseLocalTimeStep == 1 .or. UConf%UseVariableTime == 1 .or. UConf%UseOverSet == 1) then
            allocate(UCC%TimeWidth(UG%GI%RealCells,1,1)) !セル中心における時間刻み幅は(実要素数,1,1)で確保
            !allocate(UCC%TmpTimeWidth(UG%GI%RealCells,1,1,3)) !セルの各境界における時間刻み幅は(実要素数,1,1,局所界面の数)で確保
            allocate(UCC%TmpTimeWidth(UG%GI%AllCells,1,1,1)) !セルの各境界における時間刻み幅は(実要素数,1,1,局所界面の数)で確保
            UCC%TmpTimeWidth = 0.0d0
        else
            allocate(UCC%TimeWidth(1,1,1))
        end if
            UCC%TimeWidth = 0.0d0

        if(UConf%UseRRK2 == 1) then
            allocate(UCC%RungeKuttaG1(iDim+2,UG%GI%AllCells,1,1)) !保存変数は(変数の種類，要素番号,1,1)で確保
            allocate(UCC%RungeKuttaG3(iDim+2,UG%GI%AllCells,1,1)) !保存変数は(変数の種類，要素番号,1,1)で確保
            allocate(UCC%PreviousQuantity(iDim+2,UG%GI%AllCells,1,1)) !保存変数は(変数の種類，要素番号,1,1)で確保
            UCC%RungeKuttaG1 = 0.0d0
            UCC%RungeKuttaG3 = 0.0d0
            UCC%PreviousQuantity = 0.0d0
        end if

        if(UConf%UseMUSCL == 1) then
            allocate(UCC%GradientOfVariable(iDim+2,UG%GI%RealCells,1,1,3)) !保存変数の勾配ベクトルは(変数の種類，実要素番号,1,1,3方向)で確保
            allocate(UCC%VariationOfVariable(2,iDim+2,UG%GI%RealCells,1,1)) !保存変数の勾配最大値・最小値は(1最大値,2最小値，変数の種類，実要素番号,1,1)で確保
            allocate(UCC%LimiterFunction(iDim+2,UG%GI%RealCells,1,1)) !変数・セル毎の勾配制限関数
            UCC%GradientOfVariable = 0.0d0
            UCC%VariationOfVariable = 0.0d0
            UCC%LimiterFunction = 0.0d0

            allocate(UCE%NormalGradient(iDim+2,1,1,2,UG%GI%Edges)) !∇q(rj-ri)=Δ_はセルが保有する局所界面の数だけ確保する
            allocate(UCE%RebuildQunatity(iDim+2,1,1,2,UG%GI%Edges)) !再構築物理量はセルが保有する局所界面の数だけ確保する
            allocate(UCE%TmpLimiterFunction(iDim+2,1,1,2,UG%GI%Edges)) !仮置勾配制限関数は界面の表裏ごとに定義 !(次元,1,1,表裏,面番号)
            UCE%NormalGradient = 0.0d0
            UCE%RebuildQunatity = 0.0d0
            UCE%TmpLimiterFunction = 0.0d0
        end if

        if(UConf%UseOverSet == 1) then
            allocate(UCC%InterpolatedQuantity(iDim+2,UG%GI%AllCells,1,1)) !n+1段階の自格子から内挿する際に全セルの処理が終わるまで上書きしないためのデータ退避場所
            allocate(UG%GM%Interpolated(UG%GI%AllCells,1,1)) !時間が進んでいる格子から内挿されたとき，この値を1にして時間積分の対象から外す
        end if

        if(UConf%TurbulenceModel /= 0) then
            ! 廃止予定
            allocate(UCC%AbsoluteVortisity(UG%GI%RealCells, 1, 1))
            allocate(UCC%TurbulenceViscosity(UG%GI%RealCells, 1, 1))
            allocate(UCC%StrainRateTensor(iDim, iDim, UG%GI%AllCells, 1, 1))    !uvw, xyz, icell, 1, 1
            allocate(UCC%LaminarViscosity(UG%GI%AllCells, 1, 1))

            allocate(UCE%AbsoluteVortisity(UG%GI%Edges, 1, 1))
            allocate(UCE%LaminarViscosity(UG%GI%Edges))
            allocate(UCE%TurbulenceViscosity(UG%GI%Edges, 1, 1))
            allocate(UCE%StrainRateTensor(iDim, iDim, UG%GI%Edges, 1, 1))
            allocate(UG%GM%BC%VW(UG%VC%Total - UG%GI%OutlineCells))
        end if
        return
    end subroutine UAllocVariables

