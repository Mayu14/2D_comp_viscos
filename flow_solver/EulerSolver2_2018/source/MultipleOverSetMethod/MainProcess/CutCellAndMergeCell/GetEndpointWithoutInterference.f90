!***********************************/
!	Name:実際にG格子と接触している界面を取り出すプログラム
!	Alias:GetEndpointWithoutInterference
!	Description:端点を登録する
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.23
!	Update:
!	Other:not parallelize
!***********************************/
    subroutine GetEndpointWithoutInterference(LG,iGrid, iVirtualCell, iAdjacentCell)
    use StructVar_Mod
    use FrequentOperation
    implicit none

    type(LocalGrids), intent(inout) :: LG
    integer, intent(in) :: iGrid, iVirtualCell, iAdjacentCell

    type(DonorElement), pointer :: DE
    integer :: iSide, iLoop
    integer, allocatable :: Check(:), iLocalPoint(:)

    double precision, allocatable :: NewEndPoint(:,:,:) !iSide,xyz
    double precision, allocatable :: IndexOfEndPoint(:,:) !iSide,iGrid

    allocate(Check(2),iLocalPoint(3))
    allocate(NewEndPoint(2,3,LG%TotalGridsNumber))

        Check = 0
        NewEndPoint(1,:,:) = 10.0d0**6 !Very Small Number
        NewEndPoint(2,:,:) = 10.0d0**6 !Very Small Number

        DE => LG%MG(iGrid)%ADE(iAdjacentCell)%DE !rootへ移動

        do while(associated(DE)) !次のレコードがなくなるまで続行
            !座標値を取り出して，点の内外判定プログラムへ送る
            if(DE%Grid /= 0) then !G格子の要素は除外する

                iLocalPoint(:) = LG%UG(iGrid)%Tri%Point(DE%Cell,:) !内挿元セルの3頂点の頂点番号

                do iSide =1, 2 !接触判定したい界面の両側の端点について
                    !3角形の中か外かの判定用関数に放り込む 0で外，1で内
                    Check(iSide) = CheckPointInternalTriangle(LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(iSide,1:3),&
                        &   LG%MG(DE%Grid)%NSR%P%NextCoordsInG(iLocalPoint(1),1:3),&
                        &   LG%MG(DE%Grid)%NSR%P%NextCoordsInG(iLocalPoint(2),1:3),&
                        &   LG%MG(DE%Grid)%NSR%P%NextCoordsInG(iLocalPoint(3),1:3))



                    if(Check(iSide) == 1) then
                        call GetEdgeIntersectionCoordsWithTriangle(&
                            &   LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(1,1:3),LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(2,1:3), &
                            &   LG%MG(DE%Grid)%NSR%P%NextCoordsInG(iLocalPoint(1),1:3),&
                            &   LG%MG(DE%Grid)%NSR%P%NextCoordsInG(iLocalPoint(2),1:3),&
                            &   LG%MG(DE%Grid)%NSR%P%NextCoordsInG(iLocalPoint(3),1:3),&
                            &   NewEndPoint(iSide,1:3,DE%Grid))
                    !2つ以上の格子がぐちゃぐちゃに絡み合ってるときにNewEndPointの値が狂う可能性あり(2つの格子のことしか考えていないため) 01/23
                    !複数のNewEndPointを計算して，端点間距離が最小になる組み合わせを選べば正しい点が得られる
                    end if
                    Check = 0
                end do

            end if

        end do

        !もともとのEndPoint2点が作る距離ベクトル(向きは1→2)と，NewEndPointとの内積を取り，
        !iSide=1は正の最大値を，iSide=2は負の最大値を取るものを採用する
        IndexOfEndPoint = 0.0d0
        do iLoop = 1, LG%TotalGridsNumber
            do iSide = 1, 2
                if(NewEndPoint(iSide,1,iLoop) < 10.0d0**6) then
                    IndexOfEndPoint(iSide,iLoop) = dot_product(NewEndPoint(iSide,1:3,iLoop),&
                        &   LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(2,1:3)-LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(1,1:3))

                end if
            end do
        end do

        LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(1,1:3) = NewEndPoint(1,1:3,maxloc(IndexOfEndPoint(1,:),1))
        LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(2,1:3) = NewEndPoint(2,1:3,minloc(IndexOfEndPoint(1,:),1))

        LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%Length = &
        & sqrt(dot_product(LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(2,1:3)-LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(1,1:3),&
        &                    LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(2,1:3)-LG%UG(iGrid)%VC%CutEdge(iVirtualCell)%EndPoint(1,1:3)))


    return
    end subroutine GetEndpointWithoutInterference
