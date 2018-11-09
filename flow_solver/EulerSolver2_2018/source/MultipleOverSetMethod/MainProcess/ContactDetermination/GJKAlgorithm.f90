!***********************************/
!	Name:GJKアルゴリズムによるL格子同士の詳細な接触判定プログラム
!	Alias:GJKAlgorithm
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.19
!	Update:2017.12.19
!	Other:
!***********************************/
subroutine GJKAlgorithm(iKind,MG1,MG2,UG1,UG2,iToF)
    use StructVar_Mod
!    use LoopVar_Mod
    implicit none
    integer, intent(in) :: iKind !0のときN+1段階同士の接触判定，1のときN+1とN段階での接触判定
    type(MoveGrid), intent(in) :: MG1, MG2
    type(UnstructuredGrid), intent(in) :: UG1, UG2
    integer, intent(out) :: iToF
!1.Lp格子と，n+1段階のLx格子(x<p)またはn段階のLy格子(p=<y)格子について接触判定を行う
!L1格子は2~N，L2格子は3~Nのように接触判定を行う

    integer :: iPoint
    integer :: iMax1,iMin1,iMax2,iMin2,iMax3,iMin3
    integer :: iDeleted !1回目の探索で接触判定ができなかったときに，どのサポート写像を消去するかの判断基準
    double precision, allocatable :: Array4SerachPoint1(:), Array4SerachPoint2(:)
    double precision, allocatable :: SupportMap(:,:) !サポート写像(1~3番，xyz成分) 座標値を返す
    double precision, allocatable :: SupportVector3(:) !xyz成分
    double precision, allocatable :: NormalVector23(:) !xyz成分
    integer :: iTestCount = 0
    allocate(SupportMap(3,3))
    allocate(SupportVector3(3), NormalVector23(3))
    allocate(Array4SerachPoint1(UG1%CH%iTotal),Array4SerachPoint2(UG2%CH%iTotal))

    iDeleted = 0 !初期化
    !次段階と次段階の格子の接触判定を行うケース？
    !1.適当にベクトルp0を設定する((1,0,0)とかでいい，むしろ1,0なら内積いらないし処理が楽)
    !2.La格子，Lb格子の輪郭点列すべてについて，p0との内積を取り，次のように置く
    do iPoint=1, UG1%CH%iTotal
        !x成分の値が最大となる輪郭要素を探しだし，その点の座標を取得する
        Array4SerachPoint1(iPoint) = MG1%NSR%P%NextCoordsInG(UG1%CH%PointNum(iPoint),1)
    end do

    if(iKind == 0) then !N+1 to N+1

        do iPoint=1,UG2%CH%iTotal
            Array4SerachPoint2(iPoint) = MG2%NSR%P%NextCoordsInG(UG2%CH%PointNum(iPoint),1)
        end do

    else !N+1 to N

        do iPoint=1,UG2%CH%iTotal
            Array4SerachPoint2(iPoint) = MG2%NSR%P%PresentCoordsInG(UG2%CH%PointNum(iPoint),1)
        end do
    end if

        iMax1 = maxloc(Array4SerachPoint1,1) !※最大値が存在する配列上の番号を整数値として受け取っている
        iMin1 = minloc(Array4SerachPoint2,1)
        iMin2 = minloc(Array4SerachPoint1,1)
        iMax2 = maxloc(Array4SerachPoint2,1)

        iToF = -1
    do

        iTestCount = iTestCount + 1
        !(1,0,0)ベクトルとの内積を取ったものとしてサポート写像を計算
        if(iKind == 0) then !N+1 to N+1
            SupportMap(1,:) = &
            &   MG1%NSR%P%NextCoordsInG(UG1%CH%PointNum(iMax1),:) &
            & - MG2%NSR%P%NextCoordsInG(UG2%CH%PointNum(iMin1),:)

            SupportMap(2,:) = &
            &   MG1%NSR%P%NextCoordsInG(UG1%CH%PointNum(iMin2),:) &
            & - MG2%NSR%P%NextCoordsInG(UG2%CH%PointNum(iMax2),:)
        else !N+1 to N
            SupportMap(1,:) = &
            &   MG1%NSR%P%NextCoordsInG(UG1%CH%PointNum(iMax1),:) &
            & - MG2%NSR%P%PresentCoordsInG(UG2%CH%PointNum(iMin1),:)

            SupportMap(2,:) = &
            &   MG1%NSR%P%NextCoordsInG(UG1%CH%PointNum(iMin2),:) &
            & - MG2%NSR%P%PresentCoordsInG(UG2%CH%PointNum(iMax2),:)
        end if

    !3.P0とP1との内積が正であるとき，凸包の定義よりミンコフスキー差は原点を含まない = 衝突していない→exit
        if(dot_product(SupportMap(1,:),SupportMap(2,:)) > 0.0d0) iToF = 0
        if(iToF == 0) exit
        !if(UG1%Number == 2 .and. UG2%Number == 1) then

    !4.線分P0P1に垂直でかつ，原点へ向かう方向の(P1との内積が負の)ベクトルをp2とする
        !P0P1との内積0条件より法線ベクトル自体は以下の通り定まる
        SupportVector3(1) =  SupportMap(2,2) - SupportMap(1,2)
        SupportVector3(2) = -SupportMap(2,1) + SupportMap(1,1)
        SupportVector3(3) =  0.0d0

        if(dot_product(SupportMap(2,:),SupportVector3) > 0.0d0) SupportVector3 = -SupportVector3

    !5.p2方向のサポート写像を取る
        do iPoint=1, UG1%CH%iTotal
            Array4SerachPoint1(iPoint) = dot_product(SupportVector3,MG1%NSR%P%NextCoordsInG(UG1%CH%PointNum(iPoint),:))
        end do

        if(iKind == 0) then !N+1 to N+1
            do iPoint=1, UG2%CH%iTotal
                Array4SerachPoint2(iPoint) = dot_product(SupportVector3,MG2%NSR%P%NextCoordsInG(UG2%CH%PointNum(iPoint),:))
            end do
        else !N+1 to N
            do iPoint=1, UG2%CH%iTotal
                Array4SerachPoint2(iPoint) = dot_product(SupportVector3,MG2%NSR%P%PresentCoordsInG(UG2%CH%PointNum(iPoint),:))
            end do
        end if

        iMax3 = maxloc(Array4SerachPoint1,1)
        iMin3 = minloc(Array4SerachPoint2,1)

        if(iKind == 0) then !N+1 to N+1
            SupportMap(3,:) = &
            &   MG1%NSR%P%NextCoordsInG(UG1%CH%PointNum(iMax3),:) &
            & - MG2%NSR%P%NextCoordsInG(UG2%CH%PointNum(iMin3),:)
        else !N+1 to N
            SupportMap(3,:) = &
            &   MG1%NSR%P%NextCoordsInG(UG1%CH%PointNum(iMax3),:) &
            & - MG2%NSR%P%PresentCoordsInG(UG2%CH%PointNum(iMin3),:)
        end if

    !p2とP2との内積が負であるとき，凸包の定義よりミンコフスキー差は原点を含まない = 衝突していない
        if(dot_product(SupportMap(3,:),SupportVector3) < 0) iToF = 0
        if(iToF == 0) exit

    !6.p0,p1,p2のそれぞれについてベクトルpiと，線分PjPkの原点方向法線との内積が負であるものを1つ選んで消去する
        NormalVector23(1) =  SupportMap(3,2) - SupportMap(2,2)
        NormalVector23(2) = -SupportMap(3,1) + SupportMap(2,1)
        NormalVector23(3) =  0.0d0
        if(dot_product(SupportMap(3,:),NormalVector23) > 0.0d0) NormalVector23 = -NormalVector23 !原点方向に向ける

        if(dot_product(SupportMap(1,:),NormalVector23) < 0) then !サポート写像1を消去
            SupportMap(1,:) = SupportMap(3,:)
            iDeleted = 1
        end if

        if(iDeleted == 1) cycle

        NormalVector23(1) =  SupportMap(1,2) - SupportMap(3,2)
        NormalVector23(2) = -SupportMap(1,1) + SupportMap(3,1)
        NormalVector23(3) =  0.0d0
        if(dot_product(SupportMap(1,:),NormalVector23) > 0.0d0) NormalVector23 = -NormalVector23 !原点方向に向ける

        if(dot_product(SupportMap(2,:),NormalVector23) < 0) then !サポート写像2を消去
            SupportMap(2,:) = SupportMap(3,:)
            iDeleted = 1
        end if
        if(iDeleted == 1) cycle

        NormalVector23(1) =  SupportMap(2,2) - SupportMap(1,2)
        NormalVector23(2) = -SupportMap(2,1) + SupportMap(1,1)
        NormalVector23(3) =  0.0d0
        if(dot_product(SupportMap(2,:),NormalVector23) > 0.0d0) NormalVector23 = -NormalVector23 !原点方向に向ける

        if(dot_product(SupportMap(3,:),NormalVector23) < 0) then !サポート写像3を消去
            write(6,*) "ERROR! SupportMapping3 is DELETED."
        end if

        if(iDeleted == 0) iToF = 1
        if(iToF == 1) exit

    end do

return
contains
    subroutine OutputCrossImage
    implicit none
    character(len=16) :: cGJKTest
    integer :: iLoop
       write(cGJKTest, '("GJKTest/Test.txt")')

    open(unit=1,file=cGJKTest,status='unknown')
        do iPoint=1, UG1%CH%iTotal
            write(1,*) (MG1%NSR%P%NextCoordsInG(iPoint,iLoop),iLoop=1,3)
        end do

        do iPoint=1, UG2%CH%iTotal
            write(1,*) (MG2%NSR%P%PresentCoordsInG(iPoint,iLoop),iLoop=1,3)
        end do

    close(1)
    return
    end subroutine OutputCrossImage

end subroutine GJKAlgorithm
