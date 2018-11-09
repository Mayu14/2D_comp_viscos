!***********************************/
!	Name:セル座標のローカル座標系への変換
!	Alias:OMakeRelativeCoordinate
!	Description:ローカル格子の重心を基準とした座標系でセル座標を書き改める.格子点座標は出力時のみ必要で，それは初期位置に累計座標変換行列を作用させることで得られるから割愛
!	Type:
!	Input:
!	Output:
!	Note:第1ステップにおける格子初期位置の設定および，円近似接触判定に用いる重心円半径の計算もここで行っている．
!	Author:Akitaka Toyota
!	Date:2017.11.16
!	Update:2018.01.15
!	Other:initial position
!***********************************/
subroutine OMakeRelativeCoordinate(UG,MG)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none
    type(UnstructuredGrid), intent(inout) :: UG
    type(MoveGrid), intent(inout) :: MG

    write(6,*) "Change Coordinate Global into Gravity Center"
    call SetGravityCenter

    do iCell = 1, UG%GI%AllCells
        MG%RC%Cell(iCell,1:3) = UG%CD%Cell(iCell,1:3) - MG%RC%GravityCenter(1,1:3)
        MG%RC%Cell(iCell,4) = 1.0d0
    end do

    do iPoint = 1, UG%GI%Points
        MG%RC%Point(iPoint,1:3) = UG%CD%Point(iPoint,1:3) - MG%RC%GravityCenter(1,1:3)
        MG%RC%Point(iPoint,4) = 1.0d0
    end do

    call InputPresentCoords4FirstStep

    call RoughRadiusOfLocalGrid

return
contains
    subroutine SetGravityCenter
    implicit none
!read region data
!this is sample GravityCenter
!作成した格子の座標(相対座標に変換する前)の重心位置であることに注意
    !write(6,*) ubound(MG%RC%GravityCenter)
    MG%RC%GravityCenter(:,1) = 0.0d0
    MG%RC%GravityCenter(:,2) = 0.0d0
    MG%RC%GravityCenter(:,3) = 0.0d0
    MG%RC%GravityCenter(:,4) = 1.0d0

    return
    end subroutine SetGravityCenter


    subroutine InputPresentCoords4FirstStep
    implicit none
    double precision, allocatable :: Transration(:),Rotation(:)
    type(Geometry) :: TmpGM
    type(CellCenter) :: TmpCC
    allocate(Transration(3),Rotation(3))
    !Transration(1) = 5.0d0
    !Transration(2) = 0.0d0
    !格子の初期位置はここで変更する
    Transration = 0.0d0
    Rotation = 0.0d0
    write(6,*) "please input initial position(x,y,r)"
    read(5,*) Transration(1),Transration(2),Rotation(3)

    MG%GCD%Translation = Transration
    MG%GCD%Rotation = Rotation
    TmpGM%Dimension = 2
    call OMakeTransformationMatrix(MG,TmpGM) !Geom ha muimi
    call OComputeNextCoords(99,UG,MG,TmpCC)

    MG%NSR%C%PresentCoordsInG = MG%NSR%C%NextCoordsInG
    MG%NSR%P%PresentCoordsInG = MG%NSR%P%NextCoordsInG

    return
    end subroutine

    subroutine RoughRadiusOfLocalGrid
    implicit none

        UG%RoughRadius = sqrt(maxval(sum(MG%RC%Point(:,1:3)**2,2))) !重心からの座標成分の内積が最大のものの平方根 = 重心からもっとも遠い点を含む円の半径

    return
    end subroutine RoughRadiusOfLocalGrid

end subroutine OMakeRelativeCoordinate
