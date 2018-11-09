!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:CountUpVertexFromPoint
!	Description:セルiの補間に用いるセルjの頂点Pjがセルiの中に含まれるか判定して数え上げる
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.12.05
!	Update:
!	Other:not parallelize
!***********************************/
    subroutine OCountUpVertexFromEdge(CUV)
    use StructVar_Mod
    implicit none

    type(CountUpVertex), intent(inout) :: CUV

    integer :: iMyEdge, iDonorEdge
    integer :: iTorF

    do iMyEdge=1,3
        do iDonorEdge=1,3

            call CheckEdgeIntersectTriangle(CUV%TargetCoords(mod(iMyEdge,3)+1,1:3),CUV%TargetCoords(mod(iMyEdge+1,3)+1,1:3), &
                &   CUV%DonorCoords(mod(iDonorEdge,3)+1,1:3), CUV%DonorCoords(mod(iDonorEdge+1,3)+1,1:3),   &
                &   iTorF,CUV%AdditionVertex(CUV%iTotalVertexNum+1,1:3))

            if(iTorF == 1) then
                CUV%iTotalVertexNum = CUV%iTotalVertexNum + 1
            end if
        end do
    end do

    if(CUV%iTotalVertexNum < 6) then    !使ってない部分はめっちゃでかい値で埋めとく
        CUV%AdditionVertex(CUV%iTotalVertexNum+1:6,1:3) = 10.0d0**14 !minloc判定用
    end if

    return
contains

    subroutine CheckEdgeIntersectTriangle(Vertex11,Vertex12, Vertex21, Vertex22, TorF,IntersectCoords)
    implicit none
        double precision, intent(in) :: Vertex11(:), Vertex12(:), Vertex21(:), Vertex22(:)
        integer, intent(out) :: TorF
        double precision, intent(out) :: IntersectCoords(:)

        double precision :: Ksi, Eta, Delta
        double precision :: Lambda, Mu !Coefficient of Line's Gradient

    !Check Intersect
        TorF = 0 !初期化
        IntersectCoords = 0.0d0
        Ksi = (Vertex11(1) - Vertex12(1))*(Vertex21(2)-Vertex11(2)) + (Vertex11(2)-Vertex12(2))*(Vertex11(1)-Vertex21(1))
        Eta = (Vertex11(1) - Vertex12(1))*(Vertex22(2)-Vertex11(2)) + (Vertex11(2)-Vertex12(2))*(Vertex11(1)-Vertex22(1))
        if(Ksi*Eta < 0.0d0) TorF = 1

        if(TorF == 1) then !上の2行で交差判定は終わっているため，TorF=0であればここは通らずに出る
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

        end if

        return
    end subroutine CheckEdgeIntersectTriangle

    end subroutine OCountUpVertexFromEdge
