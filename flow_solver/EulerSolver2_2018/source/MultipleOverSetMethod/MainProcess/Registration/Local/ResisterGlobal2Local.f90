!***********************************/
!	Name:構造格子からデータを補間するための登録
!	Alias:RegisterGlobal2Local
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.26
!	Update:2018.01.06
!	Other:add nullify pointer
!***********************************/
    subroutine RegisterGlobal2Local(iGlobalCellNumber,T_MG,iProcessType,BoundOption,CUV,Geom,iVirtualCell,iLoop,RG2L)
    use StructVar_Mod
    use ConstantVar_Mod
    use FrequentOperation
    use RegistVar_Mod4OMP
    implicit none

    integer, intent(in) :: iGlobalCellNumber
    type(MoveGrid), intent(inout) :: T_MG !Target's MoveGrid
    integer, intent(in) :: iProcessType !初めての内挿か，そうではないか
    integer, intent(in) :: BoundOption !0:一般のセル，1:境界のセル，2:仮想セル
    type(CountUpVertex), intent(in) :: CUV
    type(Geometry), intent(in) :: Geom
    integer, intent(in) :: iVirtualCell
    integer, intent(inout) :: iLoop
    type(RegisterGlobal2LocalWithOMP), intent(inout) :: RG2L

    !integer :: RG2L%iRegisterCellNumber, RG2L%iRegisterCellNumber2
    !integer, allocatable :: RG2L%iCenter(:), RG2L%iCenter2(:) !RG2L%iCenterX,RG2L%iCenterY,RG2L%iCenterZ
    !double precision, allocatable :: RG2L%CenterDistanceVector(:)
    !double precision :: RG2L%CenterDistance2,RG2L%CenterDistance

    !type(DonorElement), pointer :: RG2L%TempDE


    if(BoundOption == 0) then !境界セル以外のやつ
        allocate(RG2L%TempDE) !新規ノードの確保
        RG2L%TempDE%Cell = iGlobalCellNumber !DonorGrid's CellNumber !ノードへの記録
        RG2L%TempDE%Grid = 0 !iGlobalGrid !DonorGridNumber !グローバル格子の場合0を返す感じで
        RG2L%TempDE%SharedArea = 1.0d0
        nullify(RG2L%TempDE%Next)
        T_MG%ADE(CUV%iTargetNumber(0))%DE => RG2L%TempDE !データの保存

    else if(BoundOption == 1) then!境界セルのうちGenerate Cellに該当するやーつ
        !allocate(RG2L%iCenter(2))
        !allocate(RG2L%CenterDistanceVector(3))
!        do iLoop=0,3 !セル中心とセル頂点
            iLoop = 0
            call GetLocalNumber2D(CUV%TargetCoords(0,1:2),2.0d0*Geom%Width(1,1,1:2),Geom%Bound(1,1:2),RG2L%iCenter(1:2))
            RG2L%iRegisterCellNumber = Local2Global(RG2L%iCenter(1),RG2L%iCenter(2),Geom%CellNumber(1))

            call GetCellCenterCoords(RG2L%iCenter(1),RG2L%iCenter(2),Geom%Bound(1,1:2),2.0d0*Geom%Width(1,1,1:2),RG2L%CenterDistanceVector(1:2))
            RG2L%CenterDistanceVector(1:2) = CUV%TargetCoords(iLoop,1:2) - RG2L%CenterDistanceVector(1:2)
            RG2L%CenterDistanceVector(3) = 0.0d0

            !RG2L%CenterDistance2 = dot_product(RG2L%CenterDistanceVector,RG2L%CenterDistanceVector)
            RG2L%CenterDistance = sqrt(dot_product(RG2L%CenterDistanceVector,RG2L%CenterDistanceVector))

            allocate(RG2L%TempDE) !新規ノードの確保
            RG2L%TempDE%Cell = RG2L%iRegisterCellNumber !DonorGrid's CellNumber !ノードへの記録
            RG2L%TempDE%Grid = 0 !iGlobalGrid !DonorGridNumber !グローバル格子の場合0を返す感じで

            !if(iLoop == 0) then
                !RG2L%TempDE%SharedArea = 1.0d0/RG2L%CenterDistance2 *0.4d0
            !else
                !RG2L%TempDE%SharedArea = 1.0d0/RG2L%CenterDistance2 *0.2d0
            !end if

            RG2L%TempDE%SharedArea = 1.0d0/RG2L%CenterDistance !暫定

            if(iProcessType == 0) then !内挿対象がここまでで登録されていない(前のノードが存在しないとき)
                nullify(RG2L%TempDE%Next)
            else
                RG2L%TempDE%Next => T_MG%ADE(CUV%iTargetNumber(0))%DE !前のノードへのリンク作成
            end if
            T_MG%ADE(CUV%iTargetNumber(0))%DE => RG2L%TempDE !データの保存
        !end do

    else !for Virtual Cell
        !allocate(RG2L%iCenter(2),RG2L%iCenter2(2))
        !Virtual Cell Coords
            call GetLocalNumber2D(CUV%TargetCoords(0,1:2),2.0d0*Geom%Width(1,1,1:2),Geom%Bound(1,1:2),RG2L%iCenter(1:2))
            RG2L%iRegisterCellNumber = Local2Global(RG2L%iCenter(1),RG2L%iCenter(2),Geom%CellNumber(1))

            call GetLocalNumber2D(CUV%DonorCoords(0,1:2),2.0d0*Geom%Width(1,1,1:2),Geom%Bound(1,1:2),RG2L%iCenter2(1:2))
            RG2L%iRegisterCellNumber2 = Local2Global(RG2L%iCenter2(1),RG2L%iCenter2(2),Geom%CellNumber(1))

        !仮想セルと隣接実セルとが同じセルに存在していた場合
        if(RG2L%iRegisterCellNumber == RG2L%iRegisterCellNumber2) then
            !界面の単位法線ベクトルの値で面の方向が分かるから，と，(1,0)の
                if(atan2(CUV%TargetCoords(1,2),CUV%TargetCoords(1,1)) < -7.0d0/9.0d0*dPi .or. 7.0d0/9.0d0*dPi <= atan2(CUV%TargetCoords(1,2),CUV%TargetCoords(1,1))) then !RG2L%iCenter(1) - 1
                    RG2L%iCenter(1) = RG2L%iCenter(1) - 1
                else if(-7.0d0/9.0d0*dPi <= atan2(CUV%TargetCoords(1,2),CUV%TargetCoords(1,1)) < -5.0d0/9.0d0*dPi) then !RG2L%iCenter(1)-1, RG2L%iCenter(2)-1
                    RG2L%iCenter(1) = RG2L%iCenter(1)-1
                    RG2L%iCenter(2) = RG2L%iCenter(2)-1
                else if(-5.0d0/9.0d0*dPi <= atan2(CUV%TargetCoords(1,2),CUV%TargetCoords(1,1)) < -3.0d0/9.0d0*dPi) then !RG2L%iCenter(1), RG2L%iCenter(2)-1
                    RG2L%iCenter(2) = RG2L%iCenter(2)-1
                else if(-3.0d0/9.0d0*dPi <= atan2(CUV%TargetCoords(1,2),CUV%TargetCoords(1,1)) < -1.0d0/9.0d0*dPi) then !RG2L%iCenter(1)+1, RG2L%iCenter(2)-1
                    RG2L%iCenter(1) = RG2L%iCenter(1)+1
                    RG2L%iCenter(2) = RG2L%iCenter(2)-1
                else if(-1.0d0/9.0d0*dPi <= atan2(CUV%TargetCoords(1,2),CUV%TargetCoords(1,1)) < 1.0d0/9.0d0*dPi) then !RG2L%iCenter(1)+1, RG2L%iCenter(2)
                    RG2L%iCenter(1) = RG2L%iCenter(1)+1
                else if(1.0d0/9.0d0*dPi <= atan2(CUV%TargetCoords(1,2),CUV%TargetCoords(1,1)) < 3.0d0/9.0d0*dPi) then !RG2L%iCenter(1)+1, RG2L%iCenter(2)+1
                    RG2L%iCenter(1) = RG2L%iCenter(1)+1
                    RG2L%iCenter(2) = RG2L%iCenter(2)+1
                else if(3.0d0/9.0d0*dPi <= atan2(CUV%TargetCoords(1,2),CUV%TargetCoords(1,1)) < 5.0d0/9.0d0*dPi) then !RG2L%iCenter(1), RG2L%iCenter(2)+1
                    RG2L%iCenter(2) = RG2L%iCenter(2)+1
                else if(5.0d0/9.0d0*dPi <= atan2(CUV%TargetCoords(1,2),CUV%TargetCoords(1,1)) < 7.0d0/9.0d0*dPi) then !RG2L%iCenter(1)-1, RG2L%iCenter(2)+1
                    RG2L%iCenter(1) = RG2L%iCenter(1)-1
                    RG2L%iCenter(2) = RG2L%iCenter(2)+1
                end if
            RG2L%iRegisterCellNumber = Local2Global(RG2L%iCenter(1),RG2L%iCenter(2),Geom%CellNumber(1))

        end if

        allocate(RG2L%TempDE)
            RG2L%TempDE%Cell = RG2L%iRegisterCellNumber
            RG2L%TempDE%Grid = 0
            RG2L%TempDE%SharedArea = 1.0d0
            nullify(RG2L%TempDE%Next)
            T_MG%ADE(iVirtualCell)%DE => RG2L%TempDE !データの保存

    end if

    return
    end subroutine RegisterGlobal2Local
