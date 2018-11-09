!***********************************/
!	Name:内挿元となるドナー要素の検索・登録プログラム
!	Alias:Register4InterpolationMain
!	Description:ドナー要素の共有面積も計算する
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.19
!	Update:2017.12.19
!	Other:
!***********************************/
subroutine Register4InterpolationMain(LG,iGrid,Geom)
    use StructVar_Mod
!    use LoopVar_Mod
    use ConstantVar_Mod , only:CoreNumberOfCPU
    use RegistVar_Mod4OMP

    implicit none
    type(LocalGrids), intent(inout) :: LG
    type(Geometry), intent(inout) :: Geom
    integer, intent(in) :: iGrid
    type(CountUpVertex) :: CUV

    integer :: iDonorGrid, iTimeCombination
    integer :: iTargetCell
    integer :: iKindOfDonor, iFindDonor
    integer :: iGlobalCellNumber
    integer :: iLoop

    type(PrivateVarOfRegist) :: PVR
    allocate(PVR%RG2L%CenterDistanceVector(3),PVR%RG2L%iCenter(2),PVR%RG2L%iCenter2(2))
    allocate(PVR%RIG%CellCenterCoords(3),PVR%RIG%GlobalGridNum(3))
    allocate(PVR%RSDILG%MWDIAG%AuxiliaryGridNum(3),PVR%RSDILG%SSDIG%CenterDistanceVector(3))
!内挿用のDonor格子を見つけて共有面積を保存する
!    write(6,*) "Registoration of Interpolation-Cells"
        call AllocCUV !CUVの確保

        LG%UG(iGrid)%GM%Interpolated = 0 !時間積分をスキップするフラグの初期化
        LG%GG%GASG%ContentsAS2U = 0 !initialize

    !格子ごとのOverlap情報の作成
    do iDonorGrid = 1, LG%TotalGridsNumber
    !1.Lp格子と，n+1段階のLx格子(x<p)またはn段階のLy格子(p=<y)格子について接触判定を行う
    !物体を扱えないため：La<Li<Lbについて，LiとLaはN+1同士，LiとLbはN+1とNとで接触判定を行うため，1つとして同じ接触判定はない
    !n+1同士の接触判定は両側から見たものとして登録することで処理を減らすAtoB is BtoA

        if(iDonorGrid < iGrid) then
            iTimeCombination = 0 !N+1同士
        else
            iTimeCombination = 1 !N+1とN
        end if

        call OContactDetermination(LG%UG(iGrid),LG%UG(iDonorGrid),LG%MG(iGrid),LG%MG(iDonorGrid),iTimeCombination)

    end do

!$omp parallel num_threads(CoreNumberOfCPU),shared(LG,iGrid,Geom),firstprivate(iTargetCell,iKindOfDonor,iFindDonor,CUV,iGlobalCellNumber,PVR,iLoop,iTimeCombination)
!$omp do
!内挿用のDonor格子を見つけて共有面積を保存する
    do iTargetCell=1, LG%UG(iGrid)%GI%RealCells !Li格子の全実要素について

        call InputCUV(1,0,iTargetCell,LG%UG(iGrid),LG%MG(iGrid)%NSR,CUV,LG%MG(iGrid))
        !グローバル格子における位置判定+登録を行う
        call RegisterIntoGlobal(iGrid,LG%GG,iGlobalCellNumber,Geom,CUV,iLoop,PVR%RIG) !CutCell入れるなら消す？

        iFindDonor = 0
        do iKindOfDonor=1, 3
            if(iKindOfDonor <= 2 .or. iFindDonor == 0) then !上位格子+物体あり格子からの補間が終わった時点でまだ補間要素が見つかっていない時下位格子の検索を行う

                do iDonorGrid=1, LG%TotalGridsNumber !計算領域中の全格子について
                    if(iDonorGrid < iGrid) then
                        iTimeCombination = 0 !N+1同士
                    else
                        iTimeCombination = 1 !N+1とN
                    end if

                    if(LG%UG(iGrid)%Overlap(iTimeCombination+1,iDonorGrid) == iKindOfDonor) then !iKind=1で上位格子，iKind=2で物体ありの下位格子，iKind=3で下位格子
                        call RoughSearchDonorInLocalGrid(iTimeCombination,LG%MG(iDonorGrid),LG%UG(iDonorGrid),LG%MG(iGrid),CUV,iFindDonor,iLoop,PVR%RSDILG)
                        call InputCUV(1,0,iTargetCell,LG%UG(iGrid),LG%MG(iGrid)%NSR,CUV,LG%MG(iGrid)) !Donor格子が変わるたびにその座標系に合わせる必要があるため，毎回グローバル座標系の値に戻す
                    end if

                end do

            else !上位格子+物体あり下位格子の検索が終わった時点で補間用要素が見つかっているとき = 積分不要
                LG%UG(iGrid)%GM%Interpolated(iTargetCell,1,1) = 1
            end if
        end do

    !この時点で要素と被る可能性のあるすべてのローカル格子との接触判定が終わっている
    !そのため，この時点でiFindDonorが0であるならばグローバル格子からの補間が必要となる
        if(LG%UG(iGrid)%GM%CellType(iTargetCell,1,1) == 99) then
            call RegisterGlobal2Local(iGlobalCellNumber,LG%MG(iGrid),iFindDonor,1,CUV,Geom,9999,iLoop,PVR%RG2L)
            iFindDonor = 1

        end if


        if(iFindDonor == 0) then

            call RegisterGlobal2Local(iGlobalCellNumber,LG%MG(iGrid),iFindDonor,0,CUV,Geom,9999,iLoop,PVR%RG2L)

        end if

    end do
!$omp end do
!$omp end parallel
    deallocate(PVR%RG2L%TempDE,PVR%RIG%Temp,PVR%RSDILG%SSDIG%ASDCN%TempDE,PVR%RSDILG%TempS2U)

    do iTargetCell = LG%UG(iGrid)%GI%RealCells+1, LG%UG(iGrid)%GI%AllCells !すべての仮想セルについて
        !仮想セルのセル中心座標から下位格子のデータを拾ってくる
        !仮想セルの補間は当初想定していなかったために大きめの修正なしでCUVが使えない．そのため，暫定的にグローバルから拾うだけに留める
        CUV%TargetCoords(0,1:3) = LG%MG(iGrid)%NSR%C%NextCoordsInG(iTargetCell,1:3)
        CUV%DonorCoords(0,1:3) = LG%MG(iGrid)%NSR%C%NextCoordsInG(LG%UG(iGrid)%VC%Cell(iTargetCell,1),1:3)
        CUV%TargetCoords(1,1:3) = LG%UG(igrid)%GM%Normal(LG%UG(iGrid)%VC%Edge(iTargetCell),1:3) !Normal Vector On Boundary Surface

        call RegisterGlobal2Local(iGlobalCellNumber,LG%MG(iGrid),9999,2,CUV,Geom,iTargetCell,iLoop,PVR%RG2L)
    end do

    !if(LG%UG(iGrid)%InternalBoundary /= 0) then
    !    call ZeroFill
    !end if

    !VC Loop:To Search G-Cell what Overlap with OuterBoundary Of L-Cell
    !and nullify Outer Boundary G-Cell's Interpolation relationship
    call ExcludeImperfectOverlapped(LG%UG(iGrid),LG%MG(iGrid),LG%GG,Geom)
    call ReRegistBoundaryCell(LG%UG(iGrid),LG%MG(iGrid),Geom)

return
contains

    subroutine AllocCUV
    implicit none

        allocate(CUV%AdditionVertex(10,3))
        allocate(CUV%AuxiliaryCellNumber(9))
        allocate(CUV%DonorCoords(0:3,4))
        allocate(CUV%iDonorNumber(0:3))
        allocate(CUV%iTargetNumber(0:3))
        allocate(CUV%TargetCoords(0:3,4))
        CUV%TargetCoords(:,4) = 1.0d0

    return
    end subroutine AllocCUV

    subroutine ZeroFill !実際使おうとするとかなり面倒なので却下
    use FrequentOperation
    implicit none
    integer :: iCenterX, iCenterY, iGlobalNum, GetGlobalNumber2D
    type(DonorElement), pointer :: TempDE
    double precision, allocatable :: CenterCoords(:), DistanceVector(:)


    allocate(CenterCoords(3),DistanceVector(3))

        do iCenterY = 1, Geom%CellNumber(2)
            do iCenterX = 1, Geom%CellNumber(1)
                call GetCellCenterCoords(iCenterX,iCenterY,Geom%Bound(1,1:2),2.0d0*Geom%Width(1,1,1:2),CenterCoords(1:2))
                CenterCoords(3) = 0.0d0

                DistanceVector = CenterCoords(1:3) - LG%MG(iGrid)%RC%GravityCenter(2,1:3)

                if(LG%UG(iGrid)%InternalRadius > sqrt(dot_product(DistanceVector,DistanceVector))) then !物体の内部とかぶっているなら(クソ雑判定)
                    iGlobalNum = GetGlobalNumber2D(CenterCoords(1:2),2.0d0*Geom%Width(1,1,1:2),Geom%Bound(1,1:2),Geom%CellNumber(1))

                    Geom%Interpolated(iCenterX,iCenterY,1) = 1
                    nullify(LG%GG%GDE(iGlobalNum)%DE)
                    allocate(TempDE) !新規ノードの確保
                        TempDE%Cell = 0 !ノードへの記録
                        TempDE%Grid = -1 !物体の内側
                        TempDE%SharedArea = 1.0d0

                        TempDE%Next => LG%GG%GDE(iGlobalNum)%DE !前のノードへのリンク作成
                        LG%GG%GDE(iGlobalNum)%DE => TempDE !データの保存

                end if
            end do
        end do

        return
    end subroutine ZeroFill

end subroutine Register4InterpolationMain
