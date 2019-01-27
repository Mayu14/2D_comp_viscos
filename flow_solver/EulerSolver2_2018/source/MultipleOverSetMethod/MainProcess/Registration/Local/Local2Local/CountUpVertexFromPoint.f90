!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:CountUpVertexFromPoint
!	Description:セルiの補間に用いるセルjの頂点Pjがセルiの中に含まれるか判定して数え上げる
!	Type:Geom,CellCenter
!	Input:Geom,CC
!	Output:CC%PrimitiveVariable
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.05
!	Update:
!	Other:not parallelize
!***********************************/
    subroutine OCountUpVertexFromPoint(CUV)
    use StructVar_Mod
    implicit none

    type(CountUpVertex), intent(inout) :: CUV
    integer :: iLocalPoint

    CUV%iTotalVertexNum = 0


    do iLocalPoint=1, 3

        CUV%iInOrOut = CheckPointInternalTriangle(CUV%DonorCoords(iLocalPoint,1:3), &
        &   CUV%TargetCoords(1,1:3),CUV%TargetCoords(2,1:3),CUV%TargetCoords(3,1:3))

        if(CUV%iInOrOut == 1) then
            CUV%iTotalVertexNum = CUV%iTotalVertexNum + 1
            CUV%AdditionVertex(CUV%iTotalVertexNum,:) = CUV%DonorCoords(iLocalPoint,1:3)
        end if
    end do

    return
contains

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

end subroutine OCountUpVertexFromPoint
