!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:MatchingWithDonorInAuxiliaryGrid
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
subroutine MatchingWithDonorInAuxiliaryGrid(D_MG,CUV,MWDIAG)
    use StructVar_Mod
    use FrequentOperation
    use RegistVar_Mod4OMP
    implicit none
    type(MoveGrid), intent(in) :: D_MG !Donor
    type(CountUpVertex), intent(inout) :: CUV

    type(MatchingWithDonorInAuxiliaryGridWithOMP), intent(inout) :: MWDIAG

    !integer, allocatable :: MWDIAG%AuxiliaryGridNum(:)
    !allocate(MWDIAG%AuxiliaryGridNum(3))

    !Auxiliary格子の座標は配列の添字と対応しているため，
    !(非構造格子の座標r,Auxiliary格子の座標r0，構造格子負端からのセル数n，格子幅dxについて
    !r-r0 > n*dxなる関係があることから，nについて解くと以下のような式になる(※2次元用で書いたためzは常に固定値になる)
    !write(6,*) CUV%TargetCoords(0,1),D_MG%ASG%Bound(1,1),D_MG%ASG%Width(1)
        call GetLocalNumber2D(CUV%TargetCoords(0,1:2),D_MG%ASG%Width(1:2),D_MG%ASG%Bound(1,1:2),MWDIAG%AuxiliaryGridNum(1:2))
        MWDIAG%AuxiliaryGridNum(3) = 1

        CUV%AuxiliaryCellNumber(5) = MakeASGNumber(D_MG%ASG,MWDIAG%AuxiliaryGridNum(1),MWDIAG%AuxiliaryGridNum(2),0,0) !(x,y)

        !write(6,*) CUV%TargetCoords(0,1),CUV%TargetCoords(0,2)
        !write(6,*) D_MG%ASG%Bound(1,1), D_MG%ASG%Bound(1,2)
        !write(6,*) D_MG%ASG%Width(1),D_MG%ASG%Width(2)
        !write(6,*) MWDIAG%AuxiliaryGridNum(1),MWDIAG%AuxiliaryGridNum(2)
        !write(6,*) CUV%AuxiliaryCellNumber(5)
        !stop

        if(CUV%AuxiliaryCellNumber(5) == 0) then !要素が存在する補助格子位置に該当する要素が存在しない(補助格子と格子の形状が大きく異なる時に発生しうる)
            CUV%AuxiliaryCellNumber(:) = 0 !すべての補助格子を無効とする
        else !補助格子中心に要素が存在するなら周辺格子の番号を求めておく

            CUV%AuxiliaryCellNumber(1) = MakeASGNumber(D_MG%ASG,MWDIAG%AuxiliaryGridNum(1),MWDIAG%AuxiliaryGridNum(2),-1,-1) !(x-1,y-1)
            CUV%AuxiliaryCellNumber(2) = MakeASGNumber(D_MG%ASG,MWDIAG%AuxiliaryGridNum(1),MWDIAG%AuxiliaryGridNum(2),0,-1)     !(x,y-1)
            CUV%AuxiliaryCellNumber(3) = MakeASGNumber(D_MG%ASG,MWDIAG%AuxiliaryGridNum(1),MWDIAG%AuxiliaryGridNum(2),1,-1) !(x+1,y-1)

            CUV%AuxiliaryCellNumber(4) = MakeASGNumber(D_MG%ASG,MWDIAG%AuxiliaryGridNum(1),MWDIAG%AuxiliaryGridNum(2),-1,0) !(x-1,y)

            CUV%AuxiliaryCellNumber(6) = MakeASGNumber(D_MG%ASG,MWDIAG%AuxiliaryGridNum(1),MWDIAG%AuxiliaryGridNum(2),1,0) !(x+1,y)

            CUV%AuxiliaryCellNumber(7) = MakeASGNumber(D_MG%ASG,MWDIAG%AuxiliaryGridNum(1),MWDIAG%AuxiliaryGridNum(2),-1,1)   !(x-1,y+1)
            CUV%AuxiliaryCellNumber(8) = MakeASGNumber(D_MG%ASG,MWDIAG%AuxiliaryGridNum(1),MWDIAG%AuxiliaryGridNum(2),0,1)       !(x,y+1)
            CUV%AuxiliaryCellNumber(9) = MakeASGNumber(D_MG%ASG,MWDIAG%AuxiliaryGridNum(1),MWDIAG%AuxiliaryGridNum(2),1,1)   !(x+1,y+1)

        end if

return
end subroutine MatchingWithDonorInAuxiliaryGrid
