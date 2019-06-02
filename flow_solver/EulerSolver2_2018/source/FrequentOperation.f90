!***********************************/
!	Name:頻繁に使うわりに組込関数が存在しないうえ毎回書くのが面倒な計算を1行で実行するためのサブルーチン・関数集
!	Alias:FrequentOperation
!	Description:
!	Type:double precision, allocatable(3)
!	Input:double precision, allocatable(3)
!	Output:double precision, allocatable(3) or double precision
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.11.10
!	Update:-
!	Other:
!***********************************/
module FrequentOperation
implicit none

public CrossProduct
public AbsCrossProduct
public CheckPointInternalTriangle
public CheckEdgeIntersect
public GetEdgeIntersectionCoords
public GetEdgeIntersectionCoordsWithTriangle
public AbsVector
public CheckBackOrFront
public GetLengthBetweenEdge
public SaveLog

contains
    subroutine CrossProduct(A,B,C) !外積計算
    implicit none
        double precision, intent(in) :: A(:), B(:) !入力：3次元列ベクトル2組
        double precision, intent(out) :: C(:)       !出力：3次元列ベクトル
            C(1) = A(2)*B(3) - A(3)*B(2)
            C(2) = A(3)*B(1) - A(1)*B(3)
            C(3) = A(1)*B(2) - A(2)*B(1)
    return
    end subroutine CrossProduct


    function AbsCrossProduct(A,B) result(R) !外積の絶対値
    implicit none
        double precision, intent(in) :: A(:), B(:) !入力：3次元列ベクトル2組
        double precision :: R                       !出力：スカラー倍精度実数
        double precision, dimension(3) :: C

        call CrossProduct(A,B,C)
        R = sqrt(C(1)**2 + C(2)**2 + C(3)**2)
    return
    end function AbsCrossProduct

    subroutine GetCellCenterCoords(LocalNumX,LocalNumY,Bound,Width,Point)
    implicit none
    integer, intent(in) :: LocalNumX, LocalNumY !局所点番号
    double precision, intent(in) :: Bound(:) !座標下端
    double precision, intent(in) :: Width(:) !界面幅Δx,Δy
    double precision, intent(out) :: Point(:) !セル中心端点座標

        Point(1) = Bound(1) + Width(1)*dble(LocalNumX-0.5d0)
        Point(2) = Bound(2) + Width(2)*dble(LocalNumY-0.5d0)

    return
    end subroutine GetCellCenterCoords

    subroutine GetEdgeEndPoint(LocalNumX,LocalNumY,Bound,Width,Point)
    implicit none
    integer, intent(in) :: LocalNumX, LocalNumY !局所点番号
    double precision, intent(in) :: Bound(:) !座標下端
    double precision, intent(in) :: Width(:) !界面幅Δx,Δy
    double precision, intent(out) :: Point(:) !界面端点座標

        Point(1) = Bound(1) + Width(1)*dble(LocalNumX-1)
        Point(2) = Bound(2) + Width(2)*dble(LocalNumY-1)

    return
    end subroutine GetEdgeEndPoint


    subroutine GetLocalNumber2D(Coords,Width,Bound,LocalNum)
    use StructVar_Mod
    implicit none

    double precision, intent(in) :: Coords(:), Width(:), Bound(:)
    integer, intent(out) :: LocalNum(:)

        LocalNum(1:2) = int((Coords(1:2) - Bound(1:2))/(Width(1:2)))+1 !補助格子におけるセル中心位置

    return
    end subroutine GetLocalNumber2D

    subroutine SwapInteger(A,B)
    implicit none
    integer, intent(inout) :: A,B
    integer :: Tmp
        Tmp = A
        A = B
        B = Tmp
        return
    end subroutine SwapInteger


    subroutine SwapDoublePrecision(A,B)
    implicit none
    double precision, intent(inout) :: A,B
    double precision :: Tmp
        Tmp = A
        A = B
        B = Tmp
        return
    end subroutine SwapDoublePrecision

    subroutine SwapDbleArray(A,B)
    implicit none
    double precision, intent(inout) :: A(:),B(:)
    double precision, allocatable :: Tmp(:)
        allocate(Tmp(ubound(A,1)))
        Tmp = A
        A = B
        B = Tmp
        return
    end subroutine SwapDbleArray


    subroutine GetPointCoordsFromCellCoords(Corner,CellCoords,HalfWidth,PointCoords)
    implicit none
        integer, intent(in) :: Corner
        double precision, intent(in) :: CellCoords(:), HalfWidth(:)
        double precision, intent(out) :: PointCoords(:)

        if(Corner == 1) then !Left Top
            PointCoords(1) = CellCoords(1) - HalfWidth(1)
            PointCoords(2) = CellCoords(2) + HalfWidth(2)

        else if(Corner == 2) then !Right Top
            PointCoords(1) = CellCoords(1) + HalfWidth(1)
            PointCoords(2) = CellCoords(2) + HalfWidth(2)

        else if(Corner == 3) then !Right Bottom
            PointCoords(1) = CellCoords(1) + HalfWidth(1)
            PointCoords(2) = CellCoords(2) - HalfWidth(2)

        else if(Corner == 4) then !Left Bottom
            PointCoords(1) = CellCoords(1) - HalfWidth(1)
            PointCoords(2) = CellCoords(2) - HalfWidth(2)

        end if

        return
    end subroutine GetPointCoordsFromCellCoords


    subroutine Global2Local(GlobalNumber,ixCellNumber,LocalNum) !グローバル番号，X方向の分割数，局所番号を返す
    implicit none

        integer, intent(in) :: GlobalNumber,ixCellNumber
        integer, intent(out) :: LocalNum(:)

        LocalNum(1) = GlobalNumber - int(GlobalNumber/ixCellNumber)*ixCellNumber
        if(LocalNum(1) /= 0) then
            LocalNum(2) = int(GlobalNumber/ixCellNumber) + 1
        else
            LocalNum(1) = ixCellNumber
            LocalNum(2) = int(GlobalNumber/ixCellNumber)
        end if

    return
    end subroutine Global2Local

    function Local2Global(iLocalX,iLocalY,ixCellNumber) result(iGlobalNumber)
    implicit none
        integer, intent(in) :: iLocalX, iLocalY, ixCellNumber
        integer :: iGlobalNumber

            iGlobalNumber = ixCellNumber*(iLocalY-1)+iLocalX

    return
    end function


    function CheckPointInternalTriangle(Point,Vertex1, Vertex2, Vertex3) result(TorF)
    implicit none
        double precision, intent(in) :: Point(:), Vertex1(:), Vertex2(:), Vertex3(:)
        integer :: TorF
        double precision, dimension(3,3) :: VertexOnRelatedCoords
        double precision :: AB, BC, CA, AA, BB

        VertexOnRelatedCoords(1,:) = Vertex1(:) - Point(:)
        VertexOnRelatedCoords(2,:) = Vertex2(:) - Point(:)
        VertexOnRelatedCoords(3,:) = Vertex3(:) - Point(:)


        AB = dot_product(VertexOnRelatedCoords(1,:),VertexOnRelatedCoords(2,:))
        BC = dot_product(VertexOnRelatedCoords(2,:),VertexOnRelatedCoords(3,:))
        CA = dot_product(VertexOnRelatedCoords(3,:),VertexOnRelatedCoords(1,:))
        AA = dot_product(VertexOnRelatedCoords(1,:),VertexOnRelatedCoords(1,:))

        if(AB*BC-AA*BC < 0.0d0) then
            TorF = 0
        else
            BB = dot_product(VertexOnRelatedCoords(2,:),VertexOnRelatedCoords(2,:))
            if(AB*BC-CA*BB < 0.0d0) then
                TorF = 0
            else
                TorF = 1
            end if
        end if

        return
    end function CheckPointInternalTriangle


    subroutine GetEdgeIntersectionCoords(Vertex11,Vertex12, Vertex21, Vertex22, IntersectCoords)
    implicit none
        double precision, intent(in) :: Vertex11(:), Vertex12(:), Vertex21(:), Vertex22(:)
        double precision, intent(out) :: IntersectCoords(:)

        double precision :: Ksi, Eta, Delta
        double precision :: Lambda, Mu !Coefficient of Line's Gradient

            !Check Intersect
            Ksi = (Vertex22(2)-Vertex21(2))*(Vertex22(1)-Vertex11(1)) - (Vertex22(1)-Vertex21(1))*(Vertex22(2)-Vertex11(2))
            Eta = (Vertex12(1)-Vertex11(1))*(Vertex22(2)-Vertex11(2)) - (Vertex12(2)-Vertex11(2))*(Vertex22(1)-Vertex11(1))
            Delta = (Vertex12(1)-Vertex11(1))*(Vertex22(2)-Vertex21(2)) - (Vertex12(2)-Vertex11(2))*(Vertex22(1)-Vertex21(1))
            if(Delta == 0.0d0) then
                write(6,*) Ksi,Eta,Delta
                stop
            end if
            Lambda = Ksi / Delta
            Mu = Eta / Delta

            IntersectCoords(1:2) = Vertex11(1:2) + Lambda*(Vertex12(1:2)-Vertex11(1:2))
            IntersectCoords(3) = 0.0d0

        return
    end subroutine GetEdgeIntersectionCoords


    function CheckEdgeIntersect(Vertex11,Vertex12, Vertex21, Vertex22) result(TorF)
    implicit none
        double precision, intent(in) :: Vertex11(:), Vertex12(:), Vertex21(:), Vertex22(:)
        integer :: TorF

        double precision :: Ksi, Eta

        !Check Intersect
        TorF = 0 !初期化
        Ksi = (Vertex11(1) - Vertex12(1))*(Vertex21(2)-Vertex11(2)) + (Vertex11(2)-Vertex12(2))*(Vertex11(1)-Vertex21(1))
        Eta = (Vertex11(1) - Vertex12(1))*(Vertex22(2)-Vertex11(2)) + (Vertex11(2)-Vertex12(2))*(Vertex11(1)-Vertex22(1))
        if(Ksi*Eta < 0.0d0) TorF = 1

        return
    end function CheckEdgeIntersect


    subroutine GetEdgeIntersectionCoordsWithTriangle(EdgeVertex1,EdgeVertex2,TriVertex1,TriVertex2,TriVertex3,IntersectCoords)
    implicit none
        double precision, intent(in) :: EdgeVertex1(:), EdgeVertex2(:)
        double precision, intent(in) :: TriVertex1(:), TriVertex2(:), TriVertex3(:)
        double precision, intent(out) :: IntersectCoords(:)
        integer :: iEndFlag

        iEndFlag = 0
        iEndFlag = CheckEdgeIntersect(EdgeVertex1,EdgeVertex2,TriVertex1,TriVertex2)
        if(iEndFlag == 1) then
            call GetEdgeIntersectionCoords(EdgeVertex1,EdgeVertex2,TriVertex1,TriVertex2,IntersectCoords)
        end if

        if(iEndFlag == 0) then
            iEndFlag = CheckEdgeIntersect(EdgeVertex1,EdgeVertex2,TriVertex2,TriVertex3)
            if(iEndFlag == 1) then
                call GetEdgeIntersectionCoords(EdgeVertex1,EdgeVertex2,TriVertex2,TriVertex3,IntersectCoords)
            end if
        end if

        if(iEndFlag == 0) then
            iEndFlag = CheckEdgeIntersect(EdgeVertex1,EdgeVertex2,TriVertex1,TriVertex3)
            if(iEndFlag == 1) then
                call GetEdgeIntersectionCoords(EdgeVertex1,EdgeVertex2,TriVertex1,TriVertex3,IntersectCoords)
            else  !Not intersect
                write(6,*) "Error Occured"
                stop
            end if
        end if

        return
    end subroutine GetEdgeIntersectionCoordsWithTriangle

    function AbsVector(vector) result(length)
        double precision, intent(in) :: vector(:)
        double precision :: length
        length = sqrt(dot_product(vector, vector))
        return
    end function AbsVector

    subroutine CheckBackOrFront(iMainCell, iAdjCell, iBackOrFront12, iBackOrFrontm11)
        implicit none
        integer, intent(in) :: iMainCell, iAdjCell
        integer, intent(out) :: iBackOrFront12
        double precision, intent(out) :: iBackOrFrontm11

        if(iMainCell > iAdjCell) then ! Main側の数字が大きい = Main側が表
            iBackOrFront12 = 1  ! GM%Width用
            iBackOrFrontm11 = -1.0d0    ! GM%Normalの向き調整用
        else    ! Main側の数字が小さい = Main側が裏
            iBackOrFront12 = 2
            iBackOrFrontm11 = 1.0d0
        end if
        return
    end subroutine CheckBackOrFront

    subroutine GetLengthBetweenEdge(UG, iEdge, iFrontCell, iBackCell, FrontLength, BackLength)
        use StructVar_Mod
        implicit none
        type(UnstructuredGrid), intent(in) :: UG
        integer, intent(in) :: iEdge
        integer, intent(out) :: iFrontCell, iBackCell
        double precision, intent(out) :: FrontLength, BackLength
            iFrontCell = UG%Line%Cell(iEdge,1,1)
            iBackCell =  UG%Line%Cell(iEdge,2,1)    !iBackCell >  iFrontCell

            FrontLength = AbsVector(UG%GM%Width(iFrontCell, UG%Line%Cell(iEdge,1,2), :))
            if(iBackCell > UG%GI%RealCells) then
                BackLength = FrontLength
            else
                BackLength = AbsVector(UG%GM%Width(iBackCell, UG%Line%Cell(iEdge,2,2), :))
            end if
        return
    end subroutine GetLengthBetweenEdge

    subroutine SaveLog(logMessage, UConf)
        use StructVar_Mod
        implicit none
        character(len=256), intent(in) :: logMessage
        type(Configulation), intent(in) :: UConf

        open(9999, file=UConf%cLogName, status="unknown", position="append")
            write(9999,*) trim(adjustl(logMessage))
        close(9999)

        return
    end subroutine SaveLog

    function get_2d_determinant(A) result(det)
        implicit none
        double precision, intent(in) :: A(:,:)
        double precision :: det
        det = A(1,1) * A(2,2) - A(1,2) * A(2,1)
        return
    end function get_2d_determinant

    subroutine set_2d_inverse_matrix(A, singular)
        implicit none
        double precision, intent(inout) :: A(:,:)
        double precision, dimension(2,2) :: tmpA
        double precision :: det, machine_eps = 0.0000000000001
        logical :: singular

        det = get_2d_determinant(A)
        !write(6,*) det
        if(abs(det) > machine_eps) then
            singular = .false.
            tmpA = A
            A(1,1) = tmpA(2,2) / det
            A(1,2) = - tmpA(1,2) / det
            A(2,1) = - tmpA(2,1) / det
            A(2,2) = tmpA(1,1) / det
            !write(6,*) tmpA
            !write(6,*) A
        else
            singular = .true.
        end if

        return
    end subroutine set_2d_inverse_matrix

end module FrequentOperation

