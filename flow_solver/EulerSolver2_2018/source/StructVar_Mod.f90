!***********************************/
!	Name:EulerEQの引数宣言を簡略化するため&スコープ制御のために作った構造型の宣言モジュール
!	Alias:StructVar_Mod
!	Description:
!	Type:いっぱい(現在14個収録)
!	Input:
!	Output:
!	Note:StructVarは構造型変数用(For Structured Variables)モジュールの略，構造格子用変数ではない
!	Author:Akitaka Toyota
!	Date:2017.10.06
!	Update:2017.12.19
!	Other:EulerEQ用に作ったけど，色々追加すればNSEQも対応できると思われる
!***********************************/
module StructVar_Mod
implicit none
!COMMON StructuredVariable
!__________________________________________________________
    type ResidualHistory
        integer :: iTime = 1
        integer :: iLastOutput = 1
        double precision, allocatable :: MaxResidual(:, :) !1~5:variable+6:Mix
        double precision, allocatable :: AveResidual(:, :) !1~5:variable+6:Mix
        double precision, allocatable :: TotalResidual(:, :)    !1~5:variable+6:Mix
    end type

    type CellCenter !セル中心(Center)にて定義される量
        double precision, allocatable :: ConservedQuantity(:,:,:,:) !保存変数(密度，運動量，全エネルギー)(非構造：変数番号，要素番号,1,1)(構造：変数番号，x,y,z)
        double precision, allocatable :: PreviousQuantity(:,:,:,:)  !時間積分処理に入る前の保存変数
        double precision, allocatable :: PastQuantity(:,:,:,:)      !NaN発生時刻から十分離れた時刻での保存変数
        double precision, allocatable :: RungeKuttaG1(:,:,:,:)       !ルンゲクッタ1段階目における保存変数の値
        double precision, allocatable :: RungeKuttaG3(:,:,:,:)       !ルンゲクッタ2段階目における保存変数の値

        double precision, allocatable :: PrimitiveVariable(:,:,:,:) !基本変数(密度，速度，全エネルギー) (構成は一緒)

        double precision, allocatable :: LaminarViscosity(:, :, :)
        double precision, allocatable :: EddyViscosity(:, :, :)
        double precision, allocatable :: VelocityNorm(:, :, :)  !(i,j,k)
        double precision, allocatable :: Temparature(:, :, :)    ! (i, j, k)
        double precision, allocatable :: StrainRateTensor(:, :, :, :, :)   !ひずみ速度テンソル(uvw, xyz, i, j, k)
        double precision, allocatable :: TemparatureGrad(:,:,:,:) !セル中心温度勾配(xyz,i,j,k)
        double precision, allocatable :: AbsoluteVortisity(:, :, :) ! (i, j, k)
        !double precision, allocatable :: MachNumber(:,:,:)  ! i,j,k

        double precision, allocatable :: GradientOfVariable(:,:,:,:,:) !各セルにおける基礎変数勾配(xyz方向)(変数番号，x座標,y座標,z座標,xyz方向),非構造：(変数番号，セル番号，1,1,xyz)
        double precision, allocatable :: VariationOfVariable(:,:,:,:,:) !各セルの隣接セルにおける基礎変数勾配の最大値(1で最大,2で最小)(1or2,変数,x,y,z)

        double precision, allocatable :: LimiterFunction(:,:,:,:) !セルにおける流束制限関数(xyz座標のみ構成は同じ)

        double precision, allocatable :: TmpTimeWidth(:,:,:,:) !構造xyz，面番号， 非構造：(要素番号,1,1,局所面番号)
        double precision, allocatable :: TimeWidth(:,:,:) !構造xyz，非構造:(要素番号,1,1)
        !for Overset
        double precision, allocatable :: InterpolatedQuantity(:,:,:,:) !補間したあとの保存変数(入れ方は一緒)
        integer :: iEndFlag !定常流計算を打ち切る判定
        integer :: ConvergeCondition = 1    ! 0:RMS, 1:Max
        type(ResidualHistory) :: RH

        double precision, allocatable :: debug(:, :)  ! debug
    end type

    type CellEdge !セル界面(Edge)にて定義される量
        double precision, allocatable :: NormalFluxDiff(:,:,:,:,:) !流束差の法線成分(FluxDifferenceInNormalDirection)(非構造：変数番号,1,1,1,面番号)(構造：変数番号，x,y,z,面番号(表裏))面番号はWest,East,South,North,Bottom,Topの順
        double precision, allocatable :: RebuildQunatity(:,:,:,:,:) !セル界面上での再構築物理量(構成は同じ)
        double precision, allocatable :: NormalGradient(:,:,:,:,:) !セル界面における法線方向の基礎変数勾配
        double precision, allocatable :: TmpLimiterFunction(:,:,:,:,:) !セル界面における流束制限関数(全界面の最小値をセル中心の流束制限関数として採用する)

        double precision, allocatable :: LaminarViscosity(:, :, :)
        double precision, allocatable :: EddyViscosity(:, :, :)
        double precision, allocatable :: StrainRateTensor(:, :, :, :, :)   !ひずみ速度テンソル(uvw, xyz, i, j, k)
        double precision, allocatable :: AbsoluteVortisity(:, :, :) ! (i, j, k)
        double precision, allocatable :: TemparatureGrad(:,:,:,:,:)   !(xyz,i,j,k,面番号)
    end type

    type ViscosityWall
        integer :: iGlobalEdge  ! 物体表面の局所境界辺番号→大域辺番号を返すやつ
        !integer :: iNumberOfMember ! その壁に所属する要素の総数
        integer :: iNumberOfMemberEdge ! その壁に所属する界面の総数
        !integer, allocatable :: iMemberCell(:)  ! その壁に所属する要素のセル番号を壁から近い順に格納(距離情報はセル側が保有)
        integer, allocatable :: iMemberEdge(:)  ! その壁に所属する界面の界面番号を壁から近い順に格納(距離情報は界面側が保有)
        double precision :: curvature   ! 物体表面の曲率
    end type

    type BoundaryCondition
        double precision, allocatable :: InFlowVariable(:) !流入する物理量のデータ(基礎変数形式)
        double precision, allocatable :: OutFlowVariable(:) !流入する物理量のデータ(基礎変数形式)
        integer :: iWallTotal = 0   ! 壁条件の界面総数
        type(ViscosityWall), allocatable :: VW(:)
    end type

    type Geometry !幾何的な情報
        double precision, allocatable :: Volume(:) !セル体積!セル番号,(構造格子の場合1)
        double precision, allocatable :: Area(:)  !セル界面面積!面番号(構造の場合xyz)非構造の場合(大域面番号)
        double precision, allocatable :: Width(:,:,:) !セル中心から界面までの距離!(構造の場合面番号(構造の場合xyz)，セル番号(構造の場合1)),(非構造(セル番号，局所面番号)
        integer :: Dimension !計算次元
    !only StructuredGrid
        double precision, allocatable :: Bound(:,:) !1下界，2上界，xyz方向
        integer, allocatable :: CellNumber(:) !構造のみ要素数!(構造の場合xyz方向)
        integer, allocatable :: ContentsG2L(:) !(GlobalCellNumber);Contents'(Local Cell's) Number include the GCell returns from GC-Number
        integer, allocatable :: WallType(:) !Boundary Condition !1:Io/Outflow 2:Wall 3:Reflect
    !only UnstructuredGrid
        double precision, allocatable :: Normal(:,:) !(界面番号,xyz)を入れると表側法線ベクトルの各成分が返る !非構造のみ
        double precision, allocatable :: AverageWidth(:) !格子i周りのセル間距離の平均(構造の場合1,非構造の場合セル数で確保する)
        integer, allocatable :: CellType(:,:,:) !構造の場合(i,j,k)，非構造の場合(セル番号,1,1) !セルの属する境界の種類を返す !0:実セル, 1:流出入境界, 2:壁, 3:反射境界，4:重合格子境界
        !これはセル中心量 or 幾何学量なのでは...？
        !今迄のプログラムに変更を加えずに...という要請を考えるとその辺に加えるべし
        integer, allocatable :: Interpolated(:,:,:) !xyz !unstructuredのとき(CellNum,1,1)
        type(BoundaryCondition) :: BC
    end type

    type Configulation !計算条件等の格納用
        integer :: SwitchProgram
        integer :: UseReadRegion = 1 !1でデータから格子情報読み出し
        integer :: UseResume = 0 !1でデータから初期条件読み出し
            character(len=256) :: ResumeFileName !初期条件を読み出す際のファイル名
            integer :: ResumeFileFormat = 2 !初期条件ファイルの形式(1:vtk,2:gnuplot用txt)
            integer :: ResumeNumber !計算を再開する際の計算ステップ
            double precision, allocatable :: ResumeInFlowVars(:) !流入する物理量のデータ(基礎変数形式)
            double precision, allocatable :: ResumeOutFlowVars(:) !流入する物理量のデータ(基礎変数形式)
        integer :: UseReadBoundaryCond = 0  !1でデータから境界条件読み出し
        integer :: UseVariableTime !1で各時間ごとの時間刻み幅が可変になる
        integer :: UseLocalTimeStep = 0 !1で局所時間刻みを使用
        integer :: UseSteadyCalc = 0 !1で定常計算
        integer :: UseFluxMethod = 1    ! 0でRoe法，1でSLAU2
        integer :: UseRRK2 = 0  !1で2段階有理ルンゲクッタ法による時間2次精度計算を実行
        integer :: UseMUSCL = 0   !1でRoeの流束差分離法による空間2次精度計算を実行
            integer :: KindLimiter = 1 !1:Venkatakrishnan, 2(未実装)
        integer :: UseOverSet = 0
        integer :: TurbulenceModel = 0  !0:invicid flow, 1:Baldwin-Lomax(RANS)
        integer :: UseSutherlandLaw = 1    !0:not use, 1:use
        integer :: UseJobParallel = 0   ! 0:not use, 1:use
        double precision :: dAttackAngle    ! 迎角(JobParallel用)
        character(len=256) :: cGridName  ! mayuファイル名(JobParallel用)
        character(len=256) :: cFileName  ! 出力vtk名(JobParallel用)
        character(len=256) :: cDirectory = ""   ! 出力先フルパス
        character(len=256) :: cLogName = "" ! 計算ログファイル名
        character(len=256) :: cCaseName = ""    ! 計算ケース名
        integer :: my_rank  ! for MPI
        integer :: my_thread    ! for OpenMP
        integer :: CalcEnv = 0  ! 0:lab, 1:supercomputer
        integer :: OutputStatus = 0 ! 0:Normal, 1:FinalStep, 2:ResumeStep
    end type

    type RoeAverage !Roe平均
        integer :: TypeOfQuantity !1:保存変数，2:基礎変数
        double precision, allocatable :: SideQuantity(:,:) !セル両側の物理量
        double precision, allocatable :: RoeAverage(:) !音速，速度(*次元)，エンタルピ
        double precision, allocatable :: GeneralFlux(:,:,:)
        integer :: Direction !1:x,2:y,3:z
    end type

    type RungeKuttaParameter !2段階ルンゲクッタ(2-order Rational Runge-Kutta Method)用のパラメータ
        double precision :: G11 !g1とg1の内積
        double precision :: G13 !g1とg3の内積
        double precision :: InvG33 !g3とg3の内積の逆数
    end type

    type FDSMatrix !Flux Difference Splitting Methodに関する量
        double precision, allocatable :: RightEigenMatrix(:,:)
        double precision, allocatable :: LeftEigenMatrix(:,:)
        double precision, allocatable :: AbsEigenValues(:,:)
    end type
!__________________________________________________________
!Only UnstructruedGrid
!Cellフィールドの使用方法
!界面のとき(界面要素番号，裏表)の形で隣接する要素番号を返す，Cell(iEdge,iDirection) = iCellsとなる
!要素のとき(要素番号，局所界面番号)の形で隣接要素番号を返す，即ちCell(iCell,iLocalEdge) = iAdjacentCell となる
!Edgeフィールドの使用方法
!要素のときのみ使用，(要素番号，局所界面番号)の形で全体界面番号を返す，即ちEdge(iCell,iLocalEdge) = iEdge
    type Line !辺データ格納用
        integer :: Total !辺要素の総数
        integer, allocatable :: Point(:,:) !(辺要素番号，局所点番号)で全体点番号を返すLine%Point(iEdge,iLocalPoint)=iPoint
        integer, allocatable :: Cell(:,:,:) !辺要素の隣接要素番号は(要素数,2,2) !最後の引数はセル番号と，相手セルから見た界面自身の局所界面番号の切替
        integer, allocatable :: Belongs2Wall(:)   !(大域界面番号)で最も近い壁境界の辺番号を返す
        double precision, allocatable :: Distance(:)    !(大域界面番号)で最も近い壁境界からの距離を返す
        double precision, allocatable :: Angle(:)   ! 辺中心から表裏のセルまでのベクトル同士がなす角
        integer, allocatable :: NextCell(:,:)    ! (辺要素番号，裏表)辺両側のセルをi, i+1としたとき，i-1, i+2に相当するセルのセル番号)
    end type

    type Triangle !三角形要素データ格納用
        integer :: Total !三角形要素の総数
        integer, allocatable :: Point(:,:) !(三角形要素番号，局所点番号)で全体点番号を返すTriangle%Point(iCell,iLocalPoint)=iPoint
        integer, allocatable :: Edge(:,:) !(三角形要素番号，局所辺番号)で全体面番号を返す
        integer, allocatable :: Cell(:,:) !!(面要素番号，裏表or局所界面番号)で隣接セル番号を返すTriangle%Cell(iEdge,iDirection) = iCells
        integer, allocatable :: Type(:) !(要素番号)でセルが属する境界の種類を返す
        integer, allocatable :: Belongs2Wall(:)   !(三角形要素番号)で最も近い壁境界の辺番号を返す
        double precision, allocatable :: Distance(:)    !(三角形要素番号)で最も近い壁境界からの距離を返す
    end type

    type Quadrilateral !四辺形要素データ格納用
        integer :: Total !四辺形要素の総数
        integer, allocatable :: Point(:,:) !(四辺形要素番号，局所点番号)で全体点番号を返すTriangle%Point(iCell,iLocalPoint)=iPoint
        integer, allocatable :: Edge(:,:) !(四辺形要素番号，局所辺番号)で全体面番号を返す
        integer, allocatable :: Cell(:,:) !!(面要素番号，裏表)で隣接セル番号を返すTriangle%Cell(iEdge,iDirection) = iCells
        integer, allocatable :: Type(:) !(要素番号)でセルが属する境界の種類を返す
    end type

!0122
    type InternalPointOfEdge !境界の界面についてのみ与える，もう1段pointer属性構造型を噛ませなければ配列宣言不可
        double precision, dimension(3) :: Coords
        double precision :: Distance
        type(InternalPointOfEdge), pointer :: Next
    end type

    !type AllocatableInternalPointOfEdge
    !    type(InternalPointOfEdge), pointer :: IPE
    !end type

!0123
    type CuttedEdge
        integer :: InternalPointNumber
        integer :: EdgeNumber, GridNumber !Input at Preprocess
        integer :: CrossPattern !1:Global 2:Local 3:Mixed !1の場合カットセル対象から除外(あるいはL格子同士でカットセル)
        integer :: GlobalNumber1, GlobalNumber2 !端点1,2の属するグローバルセルのセル番号
        double precision, allocatable :: Normal(:) !Unit Normal Vector
        double precision, allocatable :: EndPoint(:,:) !1,2, xyz !Input Everystep?
        double precision :: Length !sqrt(dot_product(EndPoint(2,1:3)-EndPoint(1,1:3),EndPoint(2,1:3)-EndPoint(1,1:3)))
        double precision, allocatable :: TmpInternalPoint(:) !xyz
        type(InternalPointOfEdge), pointer :: IPE
    end type
!0123

    type VirtualCell !仮想格子データ格納用
        integer :: Total !仮想格子の総数
        integer, allocatable :: Edge(:) !Edge(仮想格子番号)で仮想格子の共有面番号を返す
        integer, allocatable :: Cell(:,:) !Cell(仮想格子番号,1or2)で1:隣接実格子番号,2:隣接格子から見た共有界面の局所界面番号を返す
        integer, allocatable :: Type(:) !1:Outline, 2:MaterialSurface
        type(CuttedEdge), allocatable :: CutEdge(:)
    end type

    type Coordinate
        double precision, allocatable :: Cell(:,:) !(要素番号，1:x,2:y,3:z)
        double precision, allocatable :: Edge(:,:)  !(大域界面番号，1:x,2:y,3:z)
        double precision, allocatable :: Point(:,:) !(点番号，1:x,2:y,3:z)
    end type

    type GridInformation
        integer :: RealCells !すべての実格子数
        integer :: AllCells !実格子+仮想格子の総数
        integer :: Edges !すべての界面の数(形状に依らない)
        integer :: Points !すべての点の数免型
        integer :: OutlineCells !外周部の仮想セル数 !内部境界ではない仮想セルの総数
    end type

    type ConvexHull
        integer :: iTotal !外周を構成する点の総数
        integer, allocatable :: PointNum(:) !点番号(特に順序は気にしていない)
    end type

    type InternalObject
        integer :: iTotal !物体を構成する点の総数
        integer :: iMinNumber !=minval(PointNum)
        integer, allocatable :: PointNum(:) !整列された点番号
    end type

    type UnstructuredGrid !非構造格子関連情報はすべてここにまとめる？
    !CellType
        type(Line) :: Line
        type(Triangle) :: Tri
        type(Quadrilateral) :: Quad
        type(VirtualCell) :: VC
    !UnstructuredGridInformation
        type(Coordinate) :: CD
        type(GridInformation) :: GI
        type(Geometry) :: GM
        type(ConvexHull) :: CH
        type(InternalObject) :: IO !allocateされたりされなかったりしろ(適当)
!        type(BoundaryCondition) :: BC
        integer :: Number !格子の番号
        double precision, allocatable :: Overlap(:,:) !(1-2,Donor格子番号),第1引数(nとn+1段階の区別)，第2引数で0:接触なし,1:N+1下位格子と被りあり,2:内部に物体境界を持ちN+1上位格子と被りあり，3:N下位格子と被りありを返す
        double precision :: RoughRadius !大雑把に求めた半径
        double precision, allocatable :: InscribedCircle(:) !CFL条件および累積移動量の判定に用いるΔx(ここでは格子中の最長辺の長さをΔxとして採用する)
        integer :: InternalBoundary !格子が内部境界(物体)をもつとき1,ないとき0 !Body Number 1,2,3... Not have body = 0
        double precision :: InternalRadius !内部物体の半径 !円形領域のみで用いる
    end type
!__________________________________________________________
!for parallel calculate (only flow solver)
    type UpwindFluxWithOMP
        double precision, allocatable :: DeltaQuantity(:)
        double precision, allocatable :: RoeFlux(:)
    end type

    type GetGeneralFluxWithOMP
        double precision, allocatable :: Velocity(:,:)
        double precision, allocatable :: KinemeticEnergy(:),NormalVelocity(:),InverseDensity(:),TotalEnergy(:) !KinemeticEnergy,NormalP4O%GGF%Velocity,inverse density
        double precision, allocatable :: FluxJacobianMatrix(:,:,:)
    end type

    type GetRoeAverageWithOMP
        double precision, allocatable :: SqrtDensity(:),Velocity(:,:),Enthalpy(:)
        double precision :: InverseDensitySum
    end type

    type GetEigenMatrixWithOMP
        double precision :: CalcAssistB1,CalcAssistB2,InverseSSS !Sound Speed Squared
        double precision :: KinemeticEnergy, NormalVelocity
    end type

    type PrivateVar4OMP
        type(UpwindFluxWithOMP) :: UF
        type(GetGeneralFluxWithOMP) :: GGF
        type(GetRoeAverageWithOMP) :: GRA
        type(GetEigenMatrixWithOMP) :: GEM
        integer :: Cell1, Cell2
    end type

    type CalcTimeWidthWithOMP
        double precision :: LocalSoundSpeed !セルにおける局所的な音速
        double precision :: AbsoluteVelocity,AbsoluteVelocity2 !速度絶対値とその2乗
    end type
!__________________________________________________________
!以下重合格子法用
    type GravityCenterDisplacement !重心の変位量
        double precision, allocatable :: Translation(:) !x,y,z !並進移動量
        double precision, allocatable :: Rotation(:) !x,y,z !radian !回転移動量
    end type

    type TransformationMatrix !すべて4行4列で確保するが，1:3,1:3の指定をすれば回転行列として使用可能
        double precision, allocatable :: Current2Next(:,:)  !N段階(計算済み)の座標からN+1段階(未計算)の座標系への変換行列
        double precision, allocatable :: Top2Next(:,:)       !グローバル座標系からN+1段階(未計算)座標系への変換行列
        double precision, allocatable :: Top2Current(:,:)   !グローバル座標系からN段階(計算済み)座標系への変換行列
        double precision, allocatable :: Next2Current(:,:)  !N+1段階(未計算)の座標からN段階(計算済み)の座標系への変換行列
        double precision, allocatable :: Next2Top(:,:)       !N+1段階(未計算)の座標からグローバル座標系への変換行列
        double precision, allocatable :: Current2Top(:,:)   !N段階(計算済み)の座標からグローバル座標系への変換行列
    end type

    type RelativeCoordinate
        !double precision, allocatable :: Point(:,:) !(PointNumber,Direction)
        double precision, allocatable :: Cell(:,:) !(CellNumber,Direction)
        double precision, allocatable :: Point(:,:) !(PointNumber,Direction)
        double precision, allocatable :: GravityCenter(:,:) !(Time1or2,xyz) Time1=Initial,Time2=Now
    end type

    type Structured2Unstructured
        integer :: Cell !UnstGridNum
        type(Structured2Unstructured), pointer :: Next
    end type

    type AllocatableS2U
        type(Structured2Unstructured), pointer :: S2U
    end type

    type AuxiliaryStructuredGrid
        double precision, allocatable :: Bound(:,:) !1下界，2上界，xyz方向
        double precision, allocatable :: Width(:) !xyz
        integer, allocatable :: CellNumber(:) !xyz
        integer, allocatable :: RelatedU2AS(:,:) !
!        integer, allocatable :: RelatedAS2U(:,:) !(Global Number, ContentsNumber) UnstGridNum
        integer, allocatable :: ContentsAS2U(:) !ContentsNumber returns from PS-Number
        !integer :: Inclusion = 6!1セルあたり内包セル数
        integer :: TotalCell
        type(AllocatableS2U), allocatable :: AS2U(:) !(GlobalCellNum)%S2U%Cellで非構造格子の番号を返す，=>%S2U%Nextで次の要素
    end type

    type AnotherCoordinateRelation
        double precision, allocatable :: PresentCoordsInG(:,:)
        double precision, allocatable :: NextCoordsInG(:,:)
        !double precision, allocatable :: NextCoordsInN(:,:)
    end type

    type NextStepRelation
        type(AnotherCoordinateRelation) :: C !Cell
        type(AnotherCoordinateRelation) :: P !Point
    end type

    type OMPOverSet

    end type

    type DonorElement
        integer :: Grid
        integer :: Cell
        double precision :: SharedArea
        type(DonorElement), pointer :: Next
    end type

    type AllocatableDonorElement
        type(DonorElement), pointer :: DE
    end type

    type InfluenceSeparation
        integer :: InfluenceDepth
        integer, dimension(3) :: BoundaryAdjacent
    end type

    type MoveGrid
        type(GravityCenterDisplacement) :: GCD
        type(RelativeCoordinate) :: RC
        type(AuxiliaryStructuredGrid) :: ASG
        type(NextStepRelation) :: NSR
        type(TransformationMatrix) :: TM
        type(OMPOverSet) :: OOS
        type(AllocatableDonorElement), allocatable :: ADE(:) !(要素番号)で対応するセルのDE情報を返す
        type(InfluenceSeparation), allocatable :: IS(:)
        integer :: iStepFromLastMove
        double precision :: MaxInfluenceRange, LastInfluenceRange
    end type

    type CountUpVertex
        !Target
        integer, allocatable :: iTargetNumber(:), iDonorNumber(:) !(0:3) 0:Cell,1:3:Point
        double precision, allocatable :: TargetCoords(:,:), DonorCoords(:,:) !(0:3,3) 1st 0:Cell,1:3:Point,2nd xyz
        integer :: iTotalVertexNum !補間したい要素と，それに接触する別の要素とが作り出す多角形の頂点の総数(1,2,3,4,6のどれかが入る)
        double precision, allocatable :: AdditionVertex(:,:) !(局所点番号,1:3)で局所点番号ごとにxyz座標を返す
        double precision :: TargetVolume,DonorVolume !セル面積(体積)
        integer, allocatable :: AuxiliaryCellNumber(:) !検索で引っかかった要素を5番とし，1~9で周囲の補助格子のグローバル番号を返す
        double precision :: SharedArea !共有面積 !一時的な保存用(別の変数に渡すまで保持してその後破棄する)
        !for Point
        integer :: iLocalPoint, iDonorPoint, iInOrOut !三角形内の局所点番号，検査対象点の番号，検査対象点が三角形の中に含まれるかどうかの判定=1で中，=0で外
        !for Edge
        integer :: iEdgeP11, iEdgeP12, iEdgeP21, iEdgeP22 !交差判定を行う線分の端点(11,12が判定対象の線分両端，21,22が三角形を構成する線分の両端)
        integer :: iTotalDonor !Donor Cell Number
    end type

    type LocalInterpolate
        integer, dimension(2) :: iCenter !グローバルセルの指標
        double precision :: AreaSum !重みの総和(全セルの面積和または中心間距離の総和)
        double precision, allocatable :: Moment(:) !重み×物理量の総和
        double precision :: MinLimiterValue !内挿に用いるセルが持つ流束制限関数の最小値
        double precision :: Threshold !流束制限関数の閾値
        double precision, allocatable :: ArithmeticAverage(:) !内挿に用いるセルの物理量の算術平均
        integer :: SummationCellNumber !内挿に用いるセル数の総和
        double precision :: UpwindArea, DownwindArea !風上，風下それぞれの重みの総和
        double precision, allocatable :: UpwindMoment(:),DownwindMoment(:),TmpMoment(:) !風上，風下，一時的な重み物理量の総和
    end type

    type MotionWait
        type(GravityCenterDisplacement) :: Accumulate
        integer :: iMotionWait !次の格子移動までの必要ステップ数
        double precision, allocatable :: StepDistance(:)
        integer :: iEstimateWait ! = iStepFromLastMove
    end type

    type InternalPointOfCell !配列宣言しないのでこれでok
        double precision, dimension(3) :: Coords
        type(InternalPointOfCell), pointer :: Next
    end type

    type InternalEdgeOfCell
        integer :: EdgeNum, GridNum
        double precision :: CutRatio
        type(InternalEdgeOfCell), pointer :: Next
    end type

    type CutCell !グローバル格子の全要素に対して配列型で宣言，というかGG以下に入れる，こいつはallocatableで確保しとく
        integer :: UseCutCell
        integer :: UnsolvedFlag !セルマージ問題が未解決であるというフラグ
        integer :: AlivePointCheck !must Initialize every step
        integer, dimension(3) :: iAlivePoint !(1or2) 1,2:Start or End, return point number of G-CEll
        integer, dimension(2,3) :: AliveCoords
        integer :: NumberOfPoint !内在する点の総数
        integer :: NumberOfEdge !内在する界面の総数
        double precision :: CutVolume !切断後の面積
        integer :: LongestEdgeNumber
        type(InternalPointOfCell), pointer :: IPC !内点の座標に関する単連結リスト
        type(InternalEdgeOfCell), pointer :: IEC !内在界面に関する単連結リスト
        type(InternalEdgeOfCell), pointer :: IEC2
    end type
!0122

    type GlobalGrid
        type(AuxiliaryStructuredGrid) :: GASG !for Global Grid
        type(AllocatableDonorElement), allocatable :: GDE(:) !forGlobalGrid
        type(CutCell), allocatable :: CutC(:)
    end type

    type LocalGrids
        type(UnstructuredGrid), allocatable :: UG(:) !(格子番号)に対して対応する格子データを返す
        type(MoveGrid), allocatable :: MG(:) !(格子番号)に対し対応するデータを返す
        type(OMPOverSet), allocatable :: OOS(:) !上に同じ
        type(CellCenter), allocatable :: CC(:)
        type(CellCenter), allocatable :: CC4MB(:) !Cell Center for Moving Body (Influence Separation Method)
        type(CellEdge), allocatable :: CE(:)
        integer :: TotalGridsNumber !使用する格子数
        type(GlobalGrid) :: GG
        type(MotionWait), allocatable :: MW(:)
    end type

    type ComputeInfluenceAreaWithOMP
        double precision :: LocalSoundSpeed
        double precision, allocatable :: FlowVelcity(:)
        double precision, allocatable :: BodyVelocity(:)
        double precision, allocatable :: UnitDirectionVector(:)
        double precision :: Distance, ContourOfBody, NormalVelocity
    end type

    type CalculateMotionWaitWithOMP
        double precision :: LocalSoundSpeed !セルにおける局所的な音速
        double precision :: AbsoluteVelocity,AbsoluteVelocity2 !速度絶対値とその2乗
        double precision :: ShockPassingTime !衝撃波が通過するのに必要な時間
    end type

    type StopWatch
        double precision, allocatable :: LapTime(:,:)
    end type

    type ShockTubeTest
        double precision :: TimeStep
        double precision, allocatable :: Position(:)
        double precision, allocatable :: PrimitiveVariable(:,:) !1:Density, 2:XVelocity, 3:YVelocity 4:Pressure
        double precision, allocatable :: ErrorNorm(:)
    end type

    type TranslationTest
        double precision :: ErrorNorm
    end type

    type AeroCharacteristics
        double precision, allocatable :: coefficient(:, :)  !1:lift or drag(1or2), 2:timestep
        double precision, allocatable :: pressure_coefficient(:, :)    ! 1:Wall Number, 2:timestep
        double precision :: residual = 1000000
    end type

    type NumericalExperiment
        type(TranslationTest) :: TT
        type(ShockTubeTest) :: STT
        !type(AeroCharacteristics) :: AC
    end type

    type SurfaceEdge
        integer :: iGEdgeNum !大域辺番号
        integer :: iLEdgeNum !局所辺番号
        integer :: iyNormalSign !法線ベクトルy方向成分の符号
        double precision :: dxCoord ! 辺中心のx座標
    end type

    type DummySE
        type(SurfaceEdge), allocatable :: SE(:)
    end type

end module StructVar_Mod

