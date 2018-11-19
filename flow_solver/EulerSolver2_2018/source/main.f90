!***********************************/
!	Name:オイラー方程式を解くプログラム
!	Alias:EulerEQ1D
!	Description:構造格子用・非構造格子用を切り替え可能
!	Type:
!	Input:vtk or MAYUGrid
!	Output:vtk
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.06
!	Update:2017.11.10
!	Other:
!***********************************/
program EulerEQ1D
    use StructVar_Mod
    use ConstantVar_Mod
    use NumericalEx
    implicit none
    integer :: iStep,iStartStep,iStep4Plot
    double precision :: StartTime,EndTime,TotalTime
    double precision, allocatable :: WatchTime(:,:) !(1Boundary,2Flux,3Integral:1Start,2End,3Total)
    type(Configulation) :: Conf
    character :: cConfirm
    allocate(WatchTime(3,3))
    TotalTime = 0.0d0
    WatchTime = 0.0d0
    call ReadConfigulation(Conf)
    !iSwitch = 3
    !write(6,*) "Please input use CPU number."
    !read(5,*) CoreNumberOfCPU
    CoreNumberOfCPU = 4

    if(Conf%SwitchProgram == 0) then
        call StructEulerEQ(Conf)
    else if(Conf%SwitchProgram == 1) then
        call UnstructEulerEQ(Conf)
    !else if(iSwitch == 2) then
    !    call OverSetEuler
    else if(Conf%SwitchProgram == 3) then
        call MultipleOverSet(Conf)

    else if(Conf%SwitchProgram == 4) then
        call SteadUnstructEuler(Conf)

    else if(Conf%SwitchProgram == 5) then
        call JobParallelNS(Conf)
    end if
    stop

contains
    subroutine MultipleOverSet(MConf)
!$ use omp_lib
    implicit none
    type(LocalGrids) :: LG
    !type(GlobalGrid) :: GG
    type(Configulation), intent(in) :: MConf
    type(Geometry) :: Geom !for Global
    type(CellCenter) :: CC !for Global
    type(CellEdge) :: CE !for Global
    integer :: iGrid!, iTimeCombination, iCell
    type(StopWatch) :: SW
    integer :: iMoveNumber !初回移動時のみポインタ解放を防ぐ
    integer :: iMoved !その計算ステップでどれかつでも格子が移動したら1，そうでなければ0
    integer :: iBody=0 !Number of Body
    integer :: iStepFromLastMove
    integer :: iDonorGrid

    allocate(SW%LapTime(3,12))
    SW%LapTime = 0.0d0
!初回のみ回るプログラム
!計算用の設定を読み込んで

        allocate(LG%MG(GridNumber))
        allocate(LG%MW(GridNumber))
!        allocate(LG%OOS(GridNumber))
        allocate(LG%UG(GridNumber))
        allocate(LG%CC(GridNumber))
        allocate(LG%CE(GridNumber))

        LG%TotalGridsNumber = GridNumber
!すべての非構造格子の情報を読み込んで = この時点で格子の数が分かっている必要あり = 格子の数は計算用の設定に書くべき
    do iGrid = 1, GridNumber
        call UReadUnStrGrid(MConf,LG%CC(iGrid),LG%CE(iGrid),LG%UG(iGrid))
        LG%UG(iGrid)%Number = iGrid !データを入力する段階で格子の平均要素サイズが小さい順に並べ替えておくこと
        !LG%UG(iGrid)%DeltaXofCFLCondition = maxval(LG%UG(iGrid)%GM%Area) !かなり雑な仮定(すべての界面の長さがほぼ同一であり，衝撃波が各セルを通過する際は必ずその行程を要するという仮定)
        if(LG%UG(iGrid)%InternalBoundary /= 0) then
            iBody = iBody + 1
            LG%UG(iGrid)%InternalBoundary = iBody
        end if
        call URegistCellType(LG%UG(iGrid))
        call UReadInflowOutflowCondition(LG%UG(iGrid))
    end do

    allocate(LG%CC4MB(iBody))
    do iGrid = 1, GridNumber
        if(LG%UG(iGrid)%InternalBoundary /= 0) then
            call AllocCC4MB(MConf,LG%CC(iGrid),LG%CC4MB(LG%UG(iGrid)%InternalBoundary))
        end if
    end do
!グローバル格子の情報を読み込んで
    call SetRegion(MConf,Geom)
!構造格子用の前処理(配列確保〜初期条件代入)をして
    call Initialize(MConf,Geom,CC,CE)
!重合格子法用の前処理(配列の確保と，相対座標・補助格子の生成，初期位置の読み込み)をして
    call AllocVariables4MultipleOverset(LG,Geom)
    do iGrid = 1, GridNumber
        call OMakeRelativeCoordinate(LG%UG(iGrid),LG%MG(iGrid))
        call OMakeAuxiliaryStructuredGrid(LG%UG(iGrid),LG%MG(iGrid)%RC,LG%MG(iGrid)%ASG)
    end do
!非構造格子用の初期条件を代入して
    do iGrid = 1, GridNumber
        call MInitialize(MConf,LG%UG(iGrid),LG%MG(iGrid),LG%CC(iGrid),Geom)
        if(LG%UG(iGrid)%InternalBoundary /= 0) call MInitialize(MConf,LG%UG(iGrid),LG%MG(iGrid),LG%CC4MB(LG%UG(iGrid)%InternalBoundary),Geom)
    end do

!初回計算用のデータ修正
    do iGrid=1, GridNumber
        call DataFabrication4FirstStep(LG%UG(iGrid),LG%MG(iGrid),LG%CC(iGrid))
        if(LG%UG(iGrid)%InternalBoundary /= 0) call DataFabrication4FirstStep(LG%UG(iGrid),LG%MG(iGrid),LG%CC4MB(LG%UG(iGrid)%InternalBoundary))
        LG%MW(iGrid)%iEstimateWait = 10**4 !Example
    end do

!ループ演算部
        iMoveNumber = 0
        iStep4Plot = 0
        iStepFromLastMove = 0
    if(MConf%UseResume /= 1) then
        iStartStep = 1
    else
        write(6,*) "Please enter the Resume Number"
        read(5,*) iStartStep
        write(6,*) "Do you want to replace the iteration count",IterationNumber,"with the iteration count from",iStartStep,"[y/n]?"
        read(5,*) cConfirm
        if(cConfirm == "y") IterationNumber = IterationNumber + iStartStep
    end if
!初期条件の出力
        call Conserve2Primitive(Geom,CC)
        Geom%Interpolated = 0
        CC%LimiterFunction = 1.0d0
        call MOutputSTRUCT(Geom,iStartStep-1,CC)
        do iGrid=1, GridNumber
            call UConserve2Primitive(LG%UG(iGrid),LG%CC(iGrid))
            LG%CC(iGrid)%LimiterFunction = 1.0d0
            call MOutput(MConf,LG%UG(iGrid),LG%MG(iGrid)%NSR%P,LG%CC(iGrid),iStartStep-1,iGrid,LG%MG(iGrid))
            if(LG%UG(iGrid)%InternalBoundary /= 0) call MOutput4Object(LG%UG(iGrid),LG%MG(iGrid)%NSR%P,iStartStep-1,iGrid)
        end do

    do iStep = iStartStep, IterationNumber
        !write(6,*) iStep,"th"
!$ SW%LapTime(1,1) = omp_get_wtime()
    !重合格子法(全格子についての内部ループ)
        do iGrid = 1, GridNumber
            !LG%MG(iGrid)%iStepFromLastMove = LG%MG(iGrid)%iStepFromLastMove + 1 !it count another place
            if(iStep == 1) call ORefreshStoreDisplacement(LG%MW(iGrid)) !累積してきた移動量を初期化する
            !移動量の累積
            !LG%MW(iGrid)%iMotionWait = 1 !この式を入れると累積移動量による移動判定のみが有効になる
            call OStoreDisplacement(LG%MW(iGrid),iStep,iGrid) !iStep段階におけるグローバル変位量(並進，回転)を入力
            !call OMoveCheck(LG%MW(iGrid),LG%UG(iGrid)%DeltaXofCFLCondition) !累積移動量による格子移動条件を満たしているか判定を行う
            !call OMoveCheck(LG%MW(iGrid),LG%UG(iGrid)%GM%AverageWidth(1),sqrt(Geom%Area(1)**2+Geom%Area(2)**2)) !累積移動量による格子移動条件を満たしているか判定を行う
            !call OMoveCheck(LG%MW(iGrid),LG%UG(iGrid)%InscribedCircle(1),2.0d0*max(Geom%Area(1),Geom%Area(2),LG%UG(iGrid)%InscribedCircle(1)))
            call OMoveCheck(LG%MW(iGrid),LG%UG(iGrid)%InscribedCircle(1),2.0d0*LG%UG(iGrid)%InscribedCircle(1))
            !↑では回転移動量の計算用として1番セルの中心から界面までの距離の平均値を採用してみた(適当だしあとから変えたらinじゃねーの？)
        end do

        if(GridNumber /= 1) then
            do iGrid=1, GridNumber-1
                if(LG%UG(iGrid)%InternalBoundary /= 0) then
                    do iDonorGrid=iGrid+1,GridNumber
                        if(LG%UG(iDonorGrid)%InternalBoundary /= 0) then
                            call OMoveCheckOfInfluence(LG%MG(iGrid),LG%MG(iDonorGrid),LG%UG(iGrid),LG%MW(iGrid),LG%MW(iDonorGrid))
                        end if
                    end do
                end if
            end do
        end if

!$ SW%LapTime(2,1) = omp_get_wtime()
!$ SW%LapTime(3,1) = SW%LapTime(3,1) + SW%LapTime(2,1) - SW%LapTime(1,1)
!$ SW%LapTime(1,6) = omp_get_wtime()
        iMoved = 0
        do iGrid = 1, GridNumber
            !if(mod(iStep,2) == 0) LG%MW(iGrid)%iMotionWait = 0
            if(LG%MW(iGrid)%iMotionWait == 0 .and. iStep /= 1) then !動かす格子は条件を満たしているもののみ
                iMoved = 1
                write(6,*) iGrid,"th Grid Moving!",iStep4Plot,iStep
                !ローカル格子の変位量の入力
                call OGetDisplacement(LG%MG(iGrid),LG%MW(iGrid)) !累積した移動量をすべて突っ込む
                call ORefreshStoreDisplacement(LG%MW(iGrid)) !累積してきた移動量を初期化する
                !座標変換行列の生成
                call OMakeTransformationMatrix(LG%MG(iGrid),Geom)
                !次段階座標の計算
                call OComputeNextCoords(iStep,LG%UG(iGrid),LG%MG(iGrid),LG%CC(iGrid))

            else !移動していない格子について
                LG%MG(iGrid)%TM%Current2Top = LG%MG(iGrid)%TM%Next2Top
                LG%MG(iGrid)%TM%Top2Current = LG%MG(iGrid)%TM%Top2Next
            end if
        end do
        !ここで一旦ループ出る(上記作業のみなら処理順序の制約がないから並列化しやすい)
!$ SW%LapTime(2,6) = omp_get_wtime()
!$ SW%LapTime(3,6) = SW%LapTime(3,6) + SW%LapTime(2,6) - SW%LapTime(1,6)
!$ SW%LapTime(1,7) = omp_get_wtime()
        if(iStep == 1 .or. iMoved == 1) then !どれか1つの格子でも移動していたら，それ以外の格子で内挿に用いるセルが変化しない可能性は非常に低いためこれは移動があれば必ず回す
            !毎回最大数が変化するポインタの解放
            if(iMoveNumber /= 0) then !初回のみポインタ解放を行わない(確保していないため)
                do iGrid=1,GridNumber
                    call FreeMemory(LG%UG(iGrid)%GI%AllCells,LG%MG(iGrid))
                end do
                call MemoryReleaseOfGDE(LG%GG,Geom%CellNumber(1)*Geom%CellNumber(2))

            end if
            iMoveNumber = iMoveNumber + 1

            Geom%Interpolated = 0 !初期化 !各格子のInterpolatedの初期化はRegister4InterpolationMain内にて行う
            !すべての格子の現段階座標＆次段階座標が出そろっているので
            do iGrid = 1, GridNumber
                !内挿に用いる要素の特定および共有面積の保存

                call Register4InterpolationMain(LG,iGrid,Geom)
                !if(iMoved == 1) call MOutput(MConf,LG%UG(iGrid),LG%MG(iGrid)%NSR%P,LG%CC(iGrid),iStep-1,iGrid,LG%MG(iGrid)) !For Test

                !call Preregist4CutCell !カットセルで用いる新たな点座標等の格納
                !LG%MG(iGrid)%iStepFromLastMove = 1 !この処理はLocalInterpolate内部で行う
            end do
            !call MakeCutCellAndMergeCell !セルのカットとマージを実行する

        end if
!$ SW%LapTime(2,7) = omp_get_wtime()
!$ SW%LapTime(3,7) = SW%LapTime(3,7) + SW%LapTime(2,7) - SW%LapTime(1,7)
    !ここから処理順序の制約あり
    !全ローカル格子についてループ
    do iGrid=1, GridNumber
        !内挿計算
        !call InterpolationOfLocalGrid(iGrid,LG,CC,Geom,iStep)
        !call ShockConvergeLocalInterpolate(iGrid,LG,CC,Geom,iStep,MConf)
!$ SW%LapTime(1,8) = omp_get_wtime()
        call ControlLocalInterpolation(iGrid,LG,CC,Geom,iStep,MConf)
!$ SW%LapTime(2,8) = omp_get_wtime()
!$ SW%LapTime(3,8) = SW%LapTime(3,8) + SW%LapTime(2,8) - SW%LapTime(1,8)
!$ SW%LapTime(1,12) = omp_get_wtime()
        !運動量の座標変換(ローカル系での値に修正)
        call OTransferInflowMomentum(LG%UG(iGrid),LG%CC(iGrid),LG%MG(iGrid))
!$ SW%LapTime(2,12) = omp_get_wtime()
!$ SW%LapTime(3,12) = SW%LapTime(3,12) + SW%LapTime(2,12) - SW%LapTime(1,12)
        !ローカル格子でのEulerEQ
        call OUnstEuler(iStep,MConf,LG%UG(iGrid),LG%CC(iGrid),LG%CE(iGrid),SW)

        !運動量の座標変換(グローバル系での値に修正)
        call OTransferOutflowMomentum(LG%UG(iGrid),LG%CC(iGrid),LG%MG(iGrid))

        if(LG%UG(iGrid)%InternalBoundary /= 0) then
            call InfluenceSeparationMain(iStep,MConf,LG%MW(iGrid),LG%UG(iGrid),LG%MG(iGrid),LG%CC(iGrid),LG%CC4MB(LG%UG(iGrid)%InternalBoundary),LG%CE(iGrid),SW)

        end if

!$ SW%LapTime(1,10) = omp_get_wtime()
        !格子移動を行った = MotionWaitの値が0になっているとき，CFL条件からモーションウェイトを再計算する
        if(LG%MW(iGrid)%iMotionWait == 0) then
            call OCalculateMotionWait(LG%CC(iGrid),LG%UG(iGrid),FixedTimeStep,LG%MW(iGrid),Geom)
        end if
!$ SW%LapTime(2,10) = omp_get_wtime()
!$ SW%LapTime(3,10) = SW%LapTime(3,10) + SW%LapTime(2,10) - SW%LapTime(1,10)
        if(mod(iStep,OutputInterval) == 0) then
            iStep4Plot = iStep/OutputInterval
            call MOutput(MConf,LG%UG(iGrid),LG%MG(iGrid)%NSR%P,LG%CC(iGrid),iStep4Plot,iGrid,LG%MG(iGrid))
            if(LG%UG(iGrid)%InternalBoundary /= 0) call MOutput4Object(LG%UG(iGrid),LG%MG(iGrid)%NSR%P,iStep4Plot,iGrid)
        end if

        if(iMoved == 1) then
            call MOutput(MConf,LG%UG(iGrid),LG%MG(iGrid)%NSR%P,LG%CC(iGrid),iStep,iGrid,LG%MG(iGrid))
            if(LG%UG(iGrid)%InternalBoundary /= 0) call MOutput4Object(LG%UG(iGrid),LG%MG(iGrid)%NSR%P,iStep,iGrid)
            call TComputeResidualNorm(iStep,iGrid,LG%CC(iGrid),LG%UG(iGrid),LG%MG(iGrid),Geom)
        end if

        LG%MW(iGrid)%iMotionWait = LG%MW(iGrid)%iMotionWait - 1 !計算が1ステップ終了したのでウェイトを1減らす
        !残差ノルムの計算
        !call TComputeResidualNorm(iStep,iGrid,LG%CC(iGrid),LG%UG(iGrid),LG%MG(iGrid),Geom)

    end do
!$ SW%LapTime(1,11) = omp_get_wtime()
    !グローバル格子にて
        !内挿計算
        call InterpolationOfGlobalGrid(CC,Geom,LG)
!$ SW%LapTime(2,11) = omp_get_wtime()
!$ SW%LapTime(3,11) = SW%LapTime(3,11) + SW%LapTime(2,11) - SW%LapTime(1,11)
        !グローバル格子でのEulerEQ
        call OStructEuler(iStep,MConf,Geom,CC,CE,SW)

        if(mod(iStep,OutputInterval) == 0) then
            iStep4Plot = iStep / OutputInterval
            call MOutputSTRUCT(Geom,iStep4Plot,CC)
        end if

        if(iMoved == 1) then
            call MOutputSTRUCT(Geom,iStep,CC)
            call TComputeResidualNorm(iStep,0,CC,LG%UG(1),LG%MG(1),Geom)
        end if

        !if(iStep4Plot /= iMoveNumber) then
        !    iStep4Plot = iStep4Plot + 1
        !    call MOutputSTRUCT(Geom,iStep4Plot,CC)
        !    do iGrid=1, GridNumber
        !        call MOutput(MConf,LG%UG(iGrid),LG%MG(iGrid)%NSR%P,LG%CC(iGrid),iStep4Plot,iGrid)
        !    end do
        !end if

!        call TComputeResidualNorm(iStep,0,CC,LG%UG(1),LG%MG(1),Geom)
    !規定回数でループ終了
    end do
!$ SW%LapTime(1,1) = SW%LapTime(3,4)+SW%LapTime(3,5)
!$ SW%LapTime(2,1) = SW%LapTime(3,2)+SW%LapTime(3,3)
!$ SW%LapTime(3,1) = SW%LapTime(3,1)+SW%LapTime(3,6)+SW%LapTime(3,7)+SW%LapTime(3,8)+SW%LapTime(3,9)+SW%LapTime(3,10)+SW%LapTime(3,11)+SW%LapTime(3,12)
!$ SW%LapTime(1,4) = SW%LapTime(1,1) + SW%LapTime(2,1) + SW%LapTime(3,1)
!$ write(6,*) "Struct",SW%LapTime(1,1),SW%LapTime(1,1)/SW%LapTime(1,4),"%"
!$ write(6,*) "Unstruct",SW%LapTime(2,1),SW%LapTime(2,1)/SW%LapTime(1,4),"%"
!$ write(6,*) "Overset",SW%LapTime(3,1),SW%LapTime(3,1)/SW%LapTime(1,4),"%"
!$ write(6,*) "Total",SW%LapTime(1,4)

        do iGrid=1,GridNumber
            call FreeMemory(LG%UG(iGrid)%GI%AllCells,LG%MG(iGrid))
        end do
        call MemoryReleaseOfGDE(LG%GG,Geom%CellNumber(1)*Geom%CellNumber(2))

    return
    end subroutine


    subroutine UnstructEulerEQ(UConf)
    implicit none
    type(CellCenter) :: UCC
    type(CellEdge) :: UCE
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid) :: UG
    type(StopWatch) :: SW
    integer :: iLoop,iSplit !時間分割
    type(NumericalExperiment) :: NE
    integer :: ExactSolResolution
    allocate(SW%LapTime(3,7))

        call UReadUnStrGrid(UConf,UCC,UCE,UG)
        !call UTestAlloc(UCC,UCE,UG)
        call UInitialize(UConf,UG,UCC) !ReadInitial,MakeInternalBoundary
        iStep = 0
        call UReadInflowOutflowCondition(UG)
        call UOutput(UConf,UG,UCC,iStep)

        !max(UG%GM%Bound(2,1)-UG%GM%Bound(1,1),UG%GM%Bound(2,2)-UG%GM%Bound(1,2)/
        !ExactSolResolution = UG%GI%AllCells
        !call PreNumericalEx(ExactSolResolution,NE)

        iStartStep = 1
        UG%GM%Interpolated = 0
        write(6,*) "Start Calculating"
        write(6,*) CoreNumberOfCPU,"Parallel"
        UCC%iEndFlag = 0
        do iStep = iStartStep, IterationNumber

            if(UConf%UseLocalTimeStep == 0) then
                call UCheckCFL4FixedTime(UG,UCC,iSplit)
                FixedTimeStep = FixedTimeStep / dble(iSplit)
            else
                iSplit = 1
            end if

            do iLoop=1, iSplit
                call OUnstEuler(iStep,UConf,UG,UCC,UCE,SW)
            end do

            FixedTimeStep = DefaultTimeStep

            if(mod(iStep,OutputInterval) == 0) then
                if(UConf%UseRRK2 == 1) then
                    if(mod(iStep,2*OutputInterval) == 0) then
                        iStep4Plot = iStep / (2*OutputInterval)
                        call UOutput(UConf,UG,UCC,iStep4Plot)
                    end if
                else
                    iStep4Plot = iStep / OutputInterval
                    call UOutput(UConf,UG,UCC,iStep4Plot)
                end if
            end if
            if(UCC%iEndFlag == 2) exit

            !call NumericalExperimentation(iStep,UConf,NE,UCC,UG%GM,UG)

        end do

        call UOutput(UConf,UG,UCC,iStep)

    open(unit=1,file="ResultU/CalcInfoU.txt",status='unknown')
        write(1,"('この計算は2次元衝撃波管問題を非構造三角形格子にて')")
        write(1,"('総要素数 ',(1x,i6))") UG%GI%RealCells
        write(1,"('総界面数 ',(1x,i6))") UG%GI%Edges
        write(1,"('総節点数 ',(1x,i6))") UG%GI%Points
        write(1,"('時間刻み幅 ',(1x,f8.6))") FixedTimeStep
        write(1,"('空間 ',(1x,i2),'次精度')") UConf%UseMUSCL+1
        write(1,"('時間 ',(1x,i2),'次精度')") UConf%UseRRK2+1
        write(1,"('の下で無次元時間',(1x,i7),'まで解きました．')") iStep4Plot*OutputInterval
        write(1,*) ""
        !write(1,"('計算を続行したい場合はCalcConfigにて')")
        !write(1,"('UseResume=',(1x,i1))") 1
        !write(1,"('ResumeFileFormat=',(1x,i1))") 3
        !write(1,"('ResumeNumber=',(1x,i12))") iStep
        !write(1,*) "と設定のうえ，適宜反復回数を増やしてプログラムを再実行してください"
        !write(1,*) ""
        !write(1,*) "※その他の設定を変更した場合，計算を続行できるかは未検証です"
    close(1)


    return
    end subroutine UnstructEulerEQ

    subroutine StructEulerEQ(Conf)
    implicit none
    type(CellCenter) :: CC
    type(CellEdge) :: CE
    type(Configulation), intent(in) :: Conf
    type(Geometry) :: Geom
    type(NumericalExperiment) :: NE

    write(6,*) "Complete to read Configulation File"

    call SetRegion(Conf,Geom)
    write(6,*) "Complete to Calculate Region Data"

    !call TestRegionData(Geom)
    call PreNumericalEx(max(Geom%CellNumber(1),Geom%CellNumber(2)),NE)

    call Initialize(Conf,Geom,CC,CE)
    write(6,*) "Complete to input Initial Condition"

    if(Conf%UseResume /= 1) then
        iStartStep = 0
        if(Geom%Dimension == 1) call OutputGnuPlot(Geom,iStartStep,CC)
        if(Geom%Dimension >= 2) call OutputResult(Geom,iStartStep,CC)
        iStartStep = 1
        write(6,*) "Output Initial State into VTK"
    else if(Conf%UseResume == 1) then
        if(Conf%UseRRK2 == 0) iStartStep = (Conf%ResumeNumber-1)*OutputInterval+1
        if(Conf%UseRRK2 == 1) iStartStep = (Conf%ResumeNumber-1)*2*OutputInterval+1
    end if

    open(unit=101,file="EXP/DensityError.txt",status='unknown')
    open(unit=102,file="EXP/XVelocityError.txt",status='unknown')
    open(unit=103,file="EXP/YVelocityError.txt",status='unknown')
    open(unit=104,file="EXP/PressureError.txt",status='unknown')

    TotalTime = 0.0d0
    do iStep = iStartStep, IterationNumber
        call cpu_time(StartTime)
        !write(6,*) iStep,"th Calculate..."

        call SetBoundary(Geom,CC)
        !write(6,*) "Complete to input Boundary data"
        !call TestBoundary(Geom,CC)

        call CalcFlux(Conf,Geom,CC,CE)
        !write(6,*) "Complete to Calculate Flux"

        call TimeIntegral(Conf,Geom,CE,CC,iStep)
        call cpu_time(EndTime)
        TotalTime = TotalTime + (EndTime - StartTime)
        !write(6,*) "Complete to Integration"
        if(mod(iStep,OutputInterval) == 0) then
            if(Conf%UseRRK2 == 1) then
                if(mod(iStep,2*OutputInterval) == 0) then
                    iStep4Plot = iStep / (2*OutputInterval)
                    if(Geom%Dimension == 1) call OutputGnuPlot(Geom,iStep4Plot,CC)
                    if(Geom%Dimension >= 2) call OutputResult(Geom,iStartStep,CC)
                end if
            else
                iStep4Plot = iStep / OutputInterval
                if(Geom%Dimension == 1) call OutputGnuPlot(Geom,iStep4Plot,CC)
                if(Geom%Dimension >= 2) call OutputResult(Geom,iStep4Plot,CC)
    !        write(6,*) FixedTimeStep*float(iStep),"sec."
            end if
        end if

        call NumericalExperimentation(iStep,Conf,NE,CC,Geom)

    end do

    close(100)

    write(6,*) FixedTimeStep*float(iStep4Plot * OutputInterval),"sec."

    if(Geom%Dimension == 1) then
        if(iStep > 10) call MergeGnuPlot(Geom,iStep4Plot)
        call Output4Resume(Geom,iStep4Plot)
    end if

    open(unit=1,file="Result/CalcInfoS.txt",status='unknown')
        write(1,"('この計算は',(1x,i1),'次元衝撃波管問題を構造直交格子にて')") Geom%Dimension
        write(1,"('X方向格子数 ',(1x,i6))") Geom%CellNumber(1)
        if(Geom%CellNumber(2) /= 0) then
            write(1,"('Y方向格子数 ',(1x,i6))") Geom%CellNumber(2)
        end if
        write(1,"('時間刻み幅 ',(1x,f8.6))") FixedTimeStep
        write(1,"('空間 ',(1x,i2),'次精度')") Conf%UseMUSCL+1
        write(1,"('時間 ',(1x,i2),'次精度')") Conf%UseRRK2+1
        write(1,"('の下で無次元時間',(1x,i7),'まで解きました．')") iStep4Plot*OutputInterval
        write(1,"('計算にあたりcpuを使用した時間は ',(1x,f22.1),'秒でした')") TotalTime
        write(1,"('そのうち境界条件の代入に用いた時間は ',(1x,f22.1),'秒')") WatchTime(1,3)
        write(1,"('そのうち流束計算に用いた時間は ',(1x,f22.1),'秒')") WatchTime(2,3)
        write(1,"('そのうち時間積分に用いた時間は ',(1x,f22.1),'秒でした')") WatchTime(3,3)
        write(1,"('よって流束計算時間が全体に占める割合は ',(1x,f8.1),'%でした')") &
            & WatchTime(2,3)/TotalTime*100.0d0
        write(1,*) ""
        !write(1,"('計算を続行したい場合はCalcConfigにて')")
        !write(1,"('UseResume=',(1x,i1))") 1
        !write(1,"('ResumeFileFormat=',(1x,i1))") 3
        !write(1,"('ResumeNumber=',(1x,i12))") iStep
        !write(1,*) "と設定のうえ，適宜反復回数を増やしてプログラムを再実行してください"
        !write(1,*) ""
        !write(1,*) "※その他の設定を変更した場合，計算を続行できるかは未検証です"
    close(1)

    return
    end subroutine StructEulerEQ

    subroutine SteadUnstructEuler(UConf)
    implicit none
    type(CellCenter) :: UCC
    type(CellEdge) :: UCE
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid) :: UG
    type(StopWatch) :: SW
    integer :: iLoop,iSplit !時間分割
    !type(NumericalExperiment) :: NE
    !integer :: ExactSolResolution
    type(MoveGrid) :: MG
    type(MotionWait) :: MW
    type(AeroCharacteristics) :: UAC
    allocate(SW%LapTime(3,7))

        call UReadUnStrGrid(UConf,UCC,UCE,UG)
        !call UTestAlloc(UCC,UCE,UG)
        call UInitialize(UConf,UG,UCC) !ReadInitial,MakeInternalBoundary
        iStep = 0
        call UReadInflowOutflowCondition(UG)
        call UOutput(UConf,UG,UCC,iStep)

        !max(UG%GM%Bound(2,1)-UG%GM%Bound(1,1),UG%GM%Bound(2,2)-UG%GM%Bound(1,2)/

        iStartStep = 1
        UG%GM%Interpolated = 0
        write(6,*) "Start Calculating"
        write(6,*) CoreNumberOfCPU,"Parallel"
        UCC%iEndFlag = 0
        allocate(MW%StepDistance(1:3))
        allocate(UAC%coefficient(2, int(IterationNumber/OutputInterval)))
        do iStep = iStartStep, IterationNumber
            if(UConf%UseLocalTimeStep == 0) then
                call UCheckCFL4FixedTime(UG,UCC,iSplit)
                FixedTimeStep = FixedTimeStep / dble(iSplit)
            else
                iSplit = 1
            end if

            MW%StepDistance(1:3) = - UG%GM%BC%InFlowVariable(2:4)*FixedTimeStep
            call VelocityCorrection(1,UG,MG,UCC,MW,FixedTimeStep,99.0d0)

            do iLoop=1, iSplit
                call OUnstEuler(iStep,UConf,UG,UCC,UCE,SW)
            end do

            call VelocityCorrection(2,UG,MG,UCC,MW,FixedTimeStep,99.0d0)
            FixedTimeStep = DefaultTimeStep

            if(mod(iStep,OutputInterval) == 0) then
                if(UConf%UseRRK2 == 1) then
                    if(mod(iStep,2*OutputInterval) == 0) then
                        iStep4Plot = iStep / (2*OutputInterval)
                        call UOutput(UConf,UG,UCC,iStep4Plot)
                    end if
                else
                    iStep4Plot = iStep / OutputInterval
                    call UOutput(UConf,UG,UCC,iStep4Plot)
                    call UCalcAeroCharacteristics(UCC, UG, iStep4Plot, UAC)
                    write(6,*) UAC%coefficient(1, iStep4Plot), UAC%coefficient(2, iStep4Plot)
                end if
            end if

            if(UCC%iEndFlag == 2) exit

            !call NumericalExperimentation(iStep,UConf,NE,UCC,UG%GM,UG)

        end do

        call UOutput(UConf,UG,UCC,iStep)
        call UOutput_Characteristics(UConf, UAC)
    return
    end subroutine SteadUnstructEuler


end program EulerEQ1D
