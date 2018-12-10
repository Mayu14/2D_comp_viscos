!***********************************/
!	Name:界面にてひずみ速度テンソルを計算するためのプログラム
!	Alias:UGetStrainRateTensor_edge
!	Description:
!	Type:CellEdge
!	Input:Configulation,Geometory,CellCenter,CellEdge
!	Output:CE%NormalFluxDiff
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:1次元と2次元,1次精度と2次精度のみに対応
!***********************************/
subroutine UGetStrainRateTensor_edge(UConf, UG,UCC,UCE)
    use StructVar_Mod
    use LoopVar_Mod
    use FrequentOperation
    use ConstantVar_Mod, ci => ImaginaryNumber
    implicit none
    type(Configulation), intent(in) :: UConf
    type(UnstructuredGrid), intent(in) :: UG
    type(CellCenter), intent(inout) :: UCC
    type(CellEdge), intent(inout) :: UCE
    integer :: iFrontLocalEdge, iBackLocalEdge
    double precision, allocatable :: velocity_dif(:)
    double precision :: length
    allocate(velocity_dif(2))

    if(UConf%UseMUSCL == 1) then
        do iEdge = 1, UG%GI%Edges
            call UCentralDifferencePrepareAroundFace(UG, iEdge, iFrontCell, iFrontLocalEdge, iBackCell, iBackLocalEdge, length)

            velocity_dif(1:2) = UCE%RebuildQunatity(2:3, 1, 1, 2, iEdge) - UCE%RebuildQunatity(2:3, 1, 1, 1, iEdge) ! 表 - 裏

            UCE%StrainRateTensor(:, 1, iEdge, 1, 1) = velocity_dif / length * UG%GM%Normal(iEdge, 1)
            UCE%StrainRateTensor(:, 2, iEdge, 1, 1) = velocity_dif / length * (UG%GM%Normal(iEdge, 2))
            UCE%AbsoluteVortisity(iEdge, 1, 1) = abs(UCE%StrainRateTensor(1, 2, iEdge, 1, 1) &
                                             & - UCE%StrainRateTensor(2, 1, iEdge, 1, 1))
        end do

    else
        do iEdge = 1, UG%GI%Edges
            call UCentralDifferencePrepareAroundFace(UG, iEdge, iFrontCell, iFrontLocalEdge, iBackCell, iBackLocalEdge, length)

            velocity_dif(1:2) = UCC%PrimitiveVariable(2:3, iFrontCell, 1, 1) - UCC%PrimitiveVariable(2:3, iBackCell, 1, 1) ! 表 - 裏

            UCE%StrainRateTensor(:, 1, iEdge, 1, 1) = velocity_dif / length * UG%GM%Normal(iEdge, 1)
            UCE%StrainRateTensor(:, 2, iEdge, 1, 1) = velocity_dif / length * (UG%GM%Normal(iEdge, 2))
            UCE%AbsoluteVortisity(iEdge, 1, 1) = abs(UCE%StrainRateTensor(1, 2, iEdge, 1, 1) &
                                             & - UCE%StrainRateTensor(2, 1, iEdge, 1, 1))
        end do

    end if


    return
end subroutine UGetStrainRateTensor_edge
