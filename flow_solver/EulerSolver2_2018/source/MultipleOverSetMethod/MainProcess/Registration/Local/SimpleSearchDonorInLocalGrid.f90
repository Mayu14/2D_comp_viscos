!***********************************/
!	Name:内挿targetセルをドナー候補の座標系に変換する
!	Alias:SimpleSearchDonorInLocalGrid
!	Description:iTime=0でドナー格子のN+1座標系へ写す，iTime=1でN座標系へ
!	Type:
!	Input:
!	Output:
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.12.05
!	Update:2017.01.06
!	Other:addition of CUV%iTotalDonor
!***********************************/
    subroutine SimpleSearchDonorInLocalGrid(D_NSR,D_UG,iTime,iUCell,CUV,T_MG,iFindDonor,D_MG,SSDILG)
    use StructVar_Mod
    use RegistVar_Mod4OMP
    implicit none
    type(NextStepRelation), intent(in) :: D_NSR
    type(UnstructuredGrid), intent(in) :: D_UG
    integer, intent(in) :: iTime, iUCell !=TempS2U%Cell

    type(CountUpVertex), intent(inout) :: CUV
    type(MoveGrid), intent(inout) :: T_MG
    integer, intent(inout) :: iFindDonor
    type(SimpleSearchDonorInLocalGridWithOMP), intent(inout) :: SSDILG
    integer :: BoundaryBoxMethod
    doubleprecision :: m=1.0d0
    !double precision, allocatable :: SSDILG%CenterDistanceVector(:)
    !double precision :: SSDILG%CenterDistance2

    type(MoveGrid), intent(in) :: D_MG !!

    !allocate(SSDILG%CenterDistanceVector(3))

    !第1：距離の２乗和と面積和
    !接触判定を行う相手セルをCUVに登録
        call InputCUV(2,iTime,iUCell,D_UG,D_NSR,CUV,D_MG)

        if(BoundaryBoxMethod(CUV) == 1) then !補助格子に登録された要素がtarget要素と接触してそうな距離にある(=極端に離れていない)とき1
            CUV%iTotalDonor = CUV%iTotalDonor + 1

            SSDILG%CenterDistanceVector(1:3) = CUV%DonorCoords(0,1:3) - CUV%TargetCoords(0,1:3)
            SSDILG%CenterDistance2 = dot_product(SSDILG%CenterDistanceVector,SSDILG%CenterDistanceVector)

            if(SSDILG%CenterDistance2 /= 0.0d0) then
                CUV%SharedArea = (1.0d0/SSDILG%CenterDistance2)**m !ユークリッド距離の２乗の逆数

            else
                CUV%SharedArea = 10.0d0**15 !前の位置と現在の位置が一緒 = 内挿値を自分のセルの影響のみにしてやる(意図的に桁落ちを起こす)
            end if

            call AddSharedDonorCellNumber(D_UG,T_MG,CUV,SSDILG%ASDCN)
            iFindDonor = 1
                !end if

        end if

    return
    end subroutine SimpleSearchDonorInLocalGrid
