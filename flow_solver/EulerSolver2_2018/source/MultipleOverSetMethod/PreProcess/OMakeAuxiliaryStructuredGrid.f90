!***********************************/
!	Name:ローカル格子に対応した補助構造格子の生成
!	Alias:OMakeAuxiliaryStructuredGrid
!	Description:ローカル格子移動計算を効率化するための補助格子を生成する
!	Type:
!	Input:UG,MG
!	Output:
!	Note:補助構造格子の各セル内にセル中心が存在するローカル格子を登録し，相互に検索できるようにする
!	Author:Akitaka Toyota
!	Date:2017.11.25
!	Update:2017.12.29
!	Other:適切なUG,MGさえ入れれば補助格子が返るため，複数格子用の変更は不要
!***********************************/
subroutine OMakeAuxiliaryStructuredGrid(UG,RC,ASG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(in) :: UG
    type(RelativeCoordinate), intent(in) :: RC
    type(AuxiliaryStructuredGrid), intent(inout) :: ASG

    double precision :: AverageVolume, TmpWidth
!    type(Structured2Unstructured), pointer :: StU
    write(6,*) "Make Auxiliary Grid for Overset"
    call MakeAuxiliaryGrid
    call CorrespondGrid
    !call CheckAuxiliaryGrid

    !do iCell=1, ASG%TotalCell

    !StU => ASG%AS2U(iCell)%S2U
    !if (associated(StU)) write(6,*) "Target",iCell
    !do while(associated(StU))
    !    write(6,*) StU%Cell
    !    StU => StU%Next
    !end do
    !end do
    !stop

return
contains

    subroutine MakeAuxiliaryGrid
    implicit none
    double precision :: SmallIncrement
    !for 2D-Grid

        SmallIncrement = 10.0d0**(-9)
    !この値は本当に適当で，0より大きい十分に小さな値であれば正直なんでもいい
    !この値の幅だけローカル格子からはみ出す無駄な領域が増えるため，計算中に消えない程度にできるだけ小さく取ればそれで十分

        ASG%Bound(1,1) = minval(RC%Cell(:,1),1) - SmallIncrement !まったく同じ座標にしてしまうとグローバル番号を検索する式が不正な値を取り，配列がオーバーフローしてしまう
        ASG%Bound(2,1) = maxval(RC%Cell(:,1),1) + SmallIncrement
        ASG%Bound(1,2) = minval(RC%Cell(:,2),1) - SmallIncrement
        ASG%Bound(2,2) = maxval(RC%Cell(:,2),1) + SmallIncrement
        ASG%Bound(1,3) = minval(RC%Cell(:,3),1) - SmallIncrement
        ASG%Bound(2,3) = maxval(RC%Cell(:,3),1) + SmallIncrement
        !width fix to integer
        AverageVolume = sum(UG%GM%Volume(:)) / dble(UG%GI%RealCells) !要素の平均体積を計算する
        TmpWidth = 2.0d0* sqrt(AverageVolume)

        ASG%CellNumber(1:2) = int((ASG%Bound(2,1:2) - ASG%Bound(1,1:2))/TmpWidth)
        ASG%CellNumber(3) = 1

        ASG%Width(:) = (ASG%Bound(2,:) - ASG%Bound(1,:))/dble(ASG%CellNumber(:))
        ASG%TotalCell = ASG%CellNumber(1) * ASG%CellNumber(2) * ASG%CellNumber(3)

        return
    end subroutine MakeAuxiliaryGrid


    subroutine CorrespondGrid
    use FrequentOperation
    implicit none
    !type(RelativeCoordinate), intent(in) :: RC
    !type(AuxiliaryStructuredGrid), intent(inout) :: ASG

    type(Structured2Unstructured), pointer :: Temp
    integer :: iGlobalCellNumber

    allocate(ASG%RelatedU2AS(UG%GI%AllCells,3)) !(非構造要素番号,xyz)で補助格子の座標軸に沿った番号を返す !この変数いらない説が濃厚 !Local変数として破棄して大丈夫なやーつっぽい
    allocate(ASG%ContentsAS2U(ASG%TotalCell)) !(補助格子のグローバル番号)で補助格子が保有する非構造要素数が返る
    allocate(ASG%AS2U(ASG%TotalCell)) !(補助格子のグローバル番号)%S2U%Cellで内包セル，%S2U%Nextで次の内包セル

        ASG%ContentsAS2U = 0

        do iCell=1, UG%GI%RealCells !すべての非構造実セルについて
            !Auxiliary格子の座標は配列の添字と対応しているため，
            !(非構造格子の座標r,Auxiliary格子の座標r0，構造格子負端からのセル数n，格子幅dxについて
            !r-r0 > n*dxなる関係があることから，nについて解くと以下のような式になる(※2次元用で書いたためzは常に固定値になる)

            call GetLocalNumber2D(RC%Cell(iCell,1:2),ASG%Width(1:2),ASG%Bound(1,1:2),ASG%RelatedU2AS(iCell,1:2))
            ASG%RelatedU2AS(iCell,3) = 1

            !取扱いの利便性を上げるためグローバルセル番号(詳細は↓)に変換する．
            !(i,j,k) <=>!2D:(j-1)*(Lx)+i !3D:(k-1)*(Ly*Lx)+(j-1)*(jmax-1)+i
            iGlobalCellNumber = Local2Global(ASG%RelatedU2AS(iCell,1),ASG%RelatedU2AS(iCell,2),ASG%CellNumber(1))

            !Auxiliary格子ijkが非構造格子から参照された回数を記録
            ASG%ContentsAS2U(iGlobalCellNumber) = ASG%ContentsAS2U(iGlobalCellNumber) + 1

            if(ASG%ContentsAS2U(iGlobalCellNumber) == 1) then
                nullify(ASG%AS2U(iGlobalCellNumber)%S2U) !リストの初期化
            end if
            !最後に構造格子に所属する非構造格子への対応を記録する
            allocate(Temp) !新規ノードの確保
            Temp%Cell = iCell !ノードへの記録
            Temp%Next => ASG%AS2U(iGlobalCellNumber)%S2U !前のノードへのリンク作成
            ASG%AS2U(iGlobalCellNumber)%S2U => Temp !データの保存

        end do

        return
    end subroutine CorrespondGrid

    subroutine CheckAuxiliaryGrid
    use LoopVar_Mod
    implicit none

    open(unit = 1, file ="Auxiliary1.vtk", status = 'unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,*) "AuxiliaryTest"
        write(1,"('ASCII')")
        write(1,"('DATASET UNSTRUCTURED_GRID')")
        write(1,"('POINTS ',(1x,i7),' double')") UG%GI%Points

        do iPoint=1, UG%GI%Points
            !write(1,"(3(1x,f22.17))") (UG%CD%Point(iPoint,iLoop),iLoop=1,3) !Adjust Point Number
            write(1,"(3(1x,f22.17))") UG%CD%Point(iPoint,1),UG%CD%Point(iPoint,2),UG%CD%Point(iPoint,3)+1 !for OverSet
        end do
        write(1,*) ""

        write(1,"('CELLS ',2(1x,i7))") UG%GI%RealCells,UG%GI%RealCells*4
        do iCell=1, UG%GI%RealCells
                write(1,"(4(1x,i7))") 3,(UG%Tri%Point(iCell,iLoop)-1,iLoop=1,3) !Adjusted Point Number
        end do
        write(1,*) ""

        write(1,"('CELL_TYPES ',(1x,i7))") UG%GI%RealCells
        do iCell=1, UG%GI%RealCells
            write(1,"(1x,i1)") 5
        end do

        write(1,"('CELL_DATA ',i7)") UG%GI%RealCells

        write(1,"('SCALARS AuxiliaryNumber float')")
        write(1,"('LOOKUP_TABLE default')")
        do iCell=1, UG%GI%RealCells
            write(1,"((2x,i6))") (ASG%RelatedU2AS(iCell,2)-1)*ASG%CellNumber(1) + ASG%RelatedU2AS(iCell,1)
            !write(6,*) (ASG%RelatedU2AS(iCell,2)-1)*ASG%CellNumber(1) + ASG%RelatedU2AS(iCell,1)
        end do
    close(1)

    open(unit=1,file="Auxiliary2.vtk",status='unknown')
        write(1,"('# vtk DataFile Version 3.0')")
        write(1,"('Auxiliary2')")
        write(1,"('ASCII')")
        write(1,"('DATASET RECTILINEAR_GRID')")
        write(1,"('DIMENSIONS ',3(1x,i5))") ASG%CellNumber(1)+1,ASG%CellNumber(2)+1,1

!座標情報
        write(1,"('X_COORDINATES ',i5,' float')") ASG%CellNumber(1)+1
            do iCenterX=1, ASG%CellNumber(1)+1
                write(1,"(f22.17)") ASG%Bound(1,1) + ASG%Width(1) * (dble(iCenterX)-1.0d0)
            end do

        write(1,"('Y_COORDINATES ',i5,' float')") ASG%CellNumber(2)+1
        do iCenterY=1, ASG%CellNumber(2)+1
            write(1,"(f22.17)") ASG%Bound(1,2) + ASG%Width(2) * (dble(iCenterY)-1.0d0)
        end do

        write(1,"('Z_COORDINATES ',i5,' float')") 1
        write(1,"(f22.17)") 0.0

        write(1,"('CELL_DATA ',i7)") (ASG%CellNumber(1)) * (ASG%CellNumber(2))

        write(1,"('SCALARS Contents float')")
        write(1,"('LOOKUP_TABLE default')")
            do iCenterY=1,ASG%CellNumber(2)
                do iCenterX=1,ASG%CellNumber(1)
                    write(1,"((2x,i6))") ASG%ContentsAS2U((iCenterY-1)*ASG%CellNumber(1) + iCenterX)
                end do
            end do

        write(1,"('SCALARS GlobalNumber float')")
        write(1,"('LOOKUP_TABLE default')")
            do iCenterY=1,ASG%CellNumber(2)
                do iCenterX=1,ASG%CellNumber(1)
                    write(1,"((2x,i6))") (iCenterY-1)*ASG%CellNumber(1) + iCenterX
                end do
            end do

        close(1)

    return
    end subroutine CheckAuxiliaryGrid


end subroutine OMakeAuxiliaryStructuredGrid
