!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:AllocVariables4MultipleOverset
!	Description:
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.11.16
!	Update:2017.11.24
!	Other:
!***********************************/
subroutine AllocVariables4MultipleOverset(LG,Geom)
    use StructVar_Mod
    implicit none
    type(LocalGrids), intent(inout) :: LG
    type(Geometry), intent(inout) :: Geom
    integer :: iGrid

    write(6,*) "Alloc Variables for Overset"
!計算全体を通して1つだけ確保しておけばよいもの
    !call AllocCUV 必要メモリ的に大した大きさではないのでCountUpVertexは宣言と同時に確保，ルーチンにおける局所変数とします
!各格子ごとに必要なもの

    call AllocGG(LG%GG,Geom)
    do iGrid=1, LG%TotalGridsNumber
        allocate(LG%UG(iGrid)%Overlap(2,LG%TotalGridsNumber))
        call AllocMG(LG%MG(iGrid),LG%UG(iGrid))
        allocate(LG%MW(iGrid)%Accumulate%Translation(3),LG%MW(iGrid)%Accumulate%Rotation(3),LG%MW(iGrid)%StepDistance(0:3)) !0:angle
    end do



return
contains
    subroutine AllocMG(MG,UG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(MoveGrid), intent(inout) :: MG
    type(UnstructuredGrid), intent(in) :: UG

    !for AllocatableDonorElement
    allocate(MG%ADE(UG%GI%AllCells)) !総セル数

    !for AuxiliaryStructuredGrid
    allocate(MG%ASG%Width(3),MG%ASG%CellNumber(3)) !xyz
    allocate(MG%ASG%Bound(2,3)) !上界・下界，xyz

    !for GravityCenterDisplacement
    allocate(MG%GCD%Translation(3),MG%GCD%Rotation(3)) !重心位置における並進運動と回転運動

    !for NextstepRelation%AnotherCoordinates
    allocate(MG%NSR%C%PresentCoordsInG(UG%GI%AllCells,4)) !総セル数,xyzw
    allocate(MG%NSR%C%NextCoordsInG(UG%GI%AllCells,4)) !総セル数,xyzw
    allocate(MG%NSR%P%PresentCoordsInG(UG%GI%Points,4)) !総点数,xyzw
    allocate(MG%NSR%P%NextCoordsInG(UG%GI%Points,4)) !総点数,xyzw

    !for OpenMP
!    allocate(MG%OOS) !ここ未実装エリア

    !fpr RelativeCoordinates
    allocate(MG%RC%Cell(UG%GI%AllCells,4)) !総セル数, xyzw
    allocate(MG%RC%Point(UG%GI%Points,4)) !総点数, xyzw
    allocate(MG%RC%GravityCenter(2,4)) !初期・現在, xyzw

    !for TransformationMatrix
    allocate(MG%TM%Current2Next(4,4), MG%TM%Current2Top(4,4), MG%TM%Next2Current(4,4)) !座標変換行列，G2Ln,G2Ln+1,Ln2Ln+1の相互
    allocate(MG%TM%Next2Top(4,4),MG%TM%Top2Current(4,4), MG%TM%Top2Next(4,4))
        MG%TM%Current2Next = 0.0d0
        do iLoop=1,4
            MG%TM%Current2Next(iLoop,iLoop) = 1.0d0
        end do
        !単位行列で初期化
        MG%TM%Current2Top = MG%TM%Current2Next
        MG%TM%Next2Current= MG%TM%Current2Next
        MG%TM%Next2Top = MG%TM%Current2Next
        MG%TM%Top2Current = MG%TM%Current2Next
        MG%TM%Top2Next = MG%TM%Current2Next

    !for Influence Separation
    allocate(MG%IS(UG%GI%RealCells))
        MG%IS%InfluenceDepth = 0

    return
    end subroutine AllocMG

    subroutine AllocGG(GG, Geom)
    use StructVar_Mod
    implicit none
    type(GlobalGrid), intent(inout) :: GG
    type(Geometry), intent(in) :: Geom
    !Geomの情報が手に入ってからじゃないと回せないです

    !GG%GASG (AuxiliaryStructuredGrid)→グローバル格子から内挿元のローカル格子を辿るための変数

        allocate(GG%GASG%RelatedU2AS(1,1)) !使いません
        allocate(GG%GASG%AS2U(1)) !使いません ※内挿元の情報はDonorElementを用いた方が便利であるため

    allocate(GG%GASG%Width(3))
        GG%GASG%Width(1:2) = 2.0d0*Geom%Width(1,1,1:2)
        GG%GASG%Width(3) = 0.0d0

    allocate(GG%GASG%Bound(2,3)) !上界・下界,xyz
        GG%GASG%Bound(:,1:2) = Geom%Bound(:,1:2)
        GG%GASG%Bound(:,3) = 0.0d0

    allocate(GG%GASG%CellNumber(3)) !xyz
        GG%GASG%CellNumber(1:2) = Geom%CellNumber(1:2) !各方向のセル数
        GG%GASG%CellNumber(3) = 1
        GG%GASG%TotalCell = max(Geom%CellNumber(1),1)*max(Geom%CellNumber(2),1)*max(Geom%CellNumber(3),1) !セルの総数

    allocate(GG%GASG%ContentsAS2U(GG%GASG%TotalCell)) !総セル数
        GG%GASG%ContentsAS2U(:) = 0

    !GG%GDE (DonorElement)→実際のドナー情報(セル番号・格子番号・共有面積)をたどるためのデータ
    allocate(GG%GDE(GG%GASG%TotalCell))

    return
    end subroutine AllocGG



end subroutine AllocVariables4MultipleOverset
