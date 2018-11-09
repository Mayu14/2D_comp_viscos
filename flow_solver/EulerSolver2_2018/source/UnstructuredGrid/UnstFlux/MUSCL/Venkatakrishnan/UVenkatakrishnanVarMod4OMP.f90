!***********************************/
!	Name:
!	Alias:UVenkatakrishnanVar_Mod4OMP
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.02.09
!	Update:
!	Other:
!***********************************/
module UVenkatakrishnanVar_Mod4OMP
use StructVar_Mod
implicit none

    type UFindVariationWithOMP
        integer :: iAdjacentCell
        double precision, allocatable :: Variation(:,:) !5vairable,2(max or min)
    end type

    type UGetNormalGradientWithOMP
        integer :: iFrontCell
        integer :: iBackCell
    end type

    type UVenkatakrishnanWithOMP
        double precision :: VE2 !Venkatakrishnan's epsilon squared
        integer, pointer :: iMaxMin
    end type

    type PrivateVarOfUVenkatakrishnan
        type(UFindVariationWithOMP) :: UFV
        type(UGetNormalGradientWithOMP) :: UGNG
        type(UVenkatakrishnanWithOMP) :: UV
    end type

end module UVenkatakrishnanVar_Mod4OMP
