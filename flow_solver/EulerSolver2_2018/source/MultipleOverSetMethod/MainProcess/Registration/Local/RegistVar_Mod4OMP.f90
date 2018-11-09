!***********************************/
!	Name:
!	Alias:RegistVar_Mod4OMP
!	Description:
!	Type:いっぱい(現在14個収録)
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.02.09
!	Update:
!	Other:
!***********************************/
module RegistVar_Mod4OMP
use StructVar_Mod
implicit none

    type RegisterIntoGlobalWithOMP
        !iLoop
        integer, allocatable :: GlobalGridNum(:) !3
        double precision, allocatable :: CellCenterCoords(:) !3
        type(DonorElement), pointer :: Temp
    end type

        type MatchingWithDonorInAuxiliaryGridWithOMP
            integer, allocatable :: AuxiliaryGridNum(:) !3
        end type

            type AddSharedDonorCellNumberWithOMP
                type(DonorElement), pointer :: TempDE
            end type

        type SimpleSearchDonorInLocalGridWithOMP
            double precision, allocatable :: CenterDistanceVector(:) !3
            double precision :: CenterDistance2
            type(AddSharedDonorCellNumberWithOMP) :: ASDCN
        end type

    type RoughSearchDonorInLocalGridWithOMP
        !integer :: iLoop !旧iAuxiliary → iLoopに変更
        type(Structured2Unstructured), pointer :: TempS2U
        type(MatchingWithDonorInAuxiliaryGridWithOMP) :: MWDIAG
        type(SimpleSearchDonorInLocalGridWithOMP) :: SSDIG
    end type

    type RegisterGlobal2LocalWithOMP
        !iLoop
        integer :: iRegisterCellNumber
        integer :: iRegisterCellNumber2
        integer, allocatable :: iCenter(:), iCenter2(:) !iCenterX,iCenterY,iCenterZ
        double precision, allocatable :: CenterDistanceVector(:) !3
        double precision :: CenterDistance2,CenterDistance
        type(DonorElement), pointer :: TempDE
    end type

    type PrivateVarOfRegist
        type(RegisterIntoGlobalWithOMP) :: RIG
        type(RoughSearchDonorInLocalGridWithOMP) :: RSDILG
        type(RegisterGlobal2LocalWithOMP) :: RG2L
    end type

contains
    function MakeASGNumber(ASG,xNum, yNum, xCorrect, yCorrect) result(iGlobalNum)
    use StructVar_Mod
    use FrequentOperation
    implicit none
        integer, intent(in) :: xNum, yNum   !x,y方向位置
        integer, intent(in) :: xCorrect, yCorrect !x,y方向修正量
        type(AuxiliaryStructuredGrid), intent(in) :: ASG
        integer :: iGlobalNum


        iGlobalNum = Local2Global(xNum+xCorrect,yNum+yCorrect,ASG%CellNumber(1))

        !上で得られた値が不正な値を取るとき以下のどれかに引っかかって出力が0に変更される
        if(xNum > ASG%CellNumber(1)) then
            iGlobalNum = 0
        else if(xNum < 1) then !x方向位置が補助格子のx負境界より小さい
            iGlobalNum = 0
        else if(yNum > ASG%CellNumber(2)) then !y方向位置が補助格子のy正境界より大きい
            iGlobalNum = 0
        else if(yNum < 1) then !y方向位置が補助格子のy負境界より大きい
            iGlobalNum = 0

        else if(iGlobalNum < 1) then !Global number become minus
            iGlobalNum = 0
        else if(iGlobalNum > ASG%TotalCell) then !Global cell Become overflow
            iGlobalNum = 0
        !ここまでの条件はそもそも補助格子が存在しないケース
        else if(ASG%ContentsAS2U(iGlobalNum) == 0) then !補助格子自体は存在するが，その補助格子内に要素が存在しないケース
            iGlobalNum = 0
        end if

    return
    end function MakeASGNumber

end module RegistVar_Mod4OMP

