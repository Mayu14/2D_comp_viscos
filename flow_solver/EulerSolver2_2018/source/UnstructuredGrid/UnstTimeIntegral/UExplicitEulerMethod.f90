!***********************************/
!	Name:オイラー陽解法による時間積分プログラム
!	Alias:ExplicitEulerMethod
!	Description:時間1次精度(基本これで十分)
!	Type:Geom,CellEdge,CellCenter
!	Input:Geom,CC,CE
!	Output:CC%ConservedQuantity
!	Note:1~3次元すべてに対応済み
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.12.26
!	Other:
!***********************************/
subroutine UExplicitEulerMethod(UG,UCE,UCC)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod
implicit none
    type(UnstructuredGrid),intent(in) :: UG
    type(CellEdge), intent(in) :: UCE
    type(CellCenter), intent(inout) :: UCC
    double precision,allocatable :: Flux(:,:)
    double precision, allocatable :: DeltaQuantity(:)

    double precision :: Residual
    Residual = 0.0d0
    allocate(Flux(UG%GM%Dimension+2,3))

    if(UCC%iEndFlag == 0) then !unsteady

        do iCell = 1, UG%GI%RealCells !すべての実セルについて

            if(UG%GM%Interpolated(iCell,1,1) == 0) then !n+1段階の格子からの内挿を行っていない(=積分済みでない)とき0

                do iLocalEdge = 1, 3 !法線ベクトルの方向に関する場合分け
                    iEdge = UG%Tri%Edge(iCell,iLocalEdge)
                    if(UG%Tri%Cell(iCell, iLocalEdge) > iCell) then
                        Flux(:,iLocalEdge) = - UCE%NormalFluxDiff(:,1,1,1,iEdge) * UG%GM%Area(iEdge)
                    else
                        Flux(:,iLocalEdge) = + UCE%NormalFluxDiff(:,1,1,1,iEdge) * UG%GM%Area(iEdge)
                    end if
                end do
                !EEM
                UCC%ConservedQuantity(:,iCell,1,1) = UCC%ConservedQuantity(:,iCell,1,1) - 1.0d0/UG%GM%Volume(iCell)&
                        &   * UCC%TimeWidth(iTC,1,1)*sum(Flux,2)
            end if

        end do

    else !steady
        allocate(DeltaQuantity(5))
        Residual = 0.0d0

        do iCell = 1, UG%GI%RealCells
            do iLocalEdge = 1, 3 !法線ベクトルの方向に関する場合分け
                iEdge = UG%Tri%Edge(iCell,iLocalEdge)
                if(UG%Tri%Cell(iCell, iLocalEdge) > iCell) then
                    Flux(:,iLocalEdge) = - UCE%NormalFluxDiff(:,1,1,1,iEdge) * UG%GM%Area(iEdge)
                else
                    Flux(:,iLocalEdge) = + UCE%NormalFluxDiff(:,1,1,1,iEdge) * UG%GM%Area(iEdge)
                end if
            end do

            DeltaQuantity =  - 1.0d0/UG%GM%Volume(iCell) * UCC%TimeWidth(iTC,1,1)*sum(Flux,2)
            !EEM
            UCC%ConservedQuantity(:,iCell,1,1) = UCC%ConservedQuantity(:,iCell,1,1) + DeltaQuantity

            !RMS
            UCC%RH%AveResidual(1:5, UCC%RH%iTime) = UCC%RH%AveResidual(1:5, UCC%RH%iTime) + abs(DeltaQuantity)
            UCC%RH%AveResidual(6, UCC%RH%iTime) = UCC%RH%AveResidual(6, UCC%RH%iTime) + 0.25d0*dot_product(DeltaQuantity,DeltaQuantity)

            !Max Residual
            UCC%RH%MaxResidual(1, UCC%RH%iTime) = max(UCC%RH%AveResidual(1, UCC%RH%iTime), abs(DeltaQuantity(1)))
            UCC%RH%MaxResidual(2, UCC%RH%iTime) = max(UCC%RH%AveResidual(2, UCC%RH%iTime), abs(DeltaQuantity(2)))
            UCC%RH%MaxResidual(3, UCC%RH%iTime) = max(UCC%RH%AveResidual(3, UCC%RH%iTime), abs(DeltaQuantity(3)))
            UCC%RH%MaxResidual(4, UCC%RH%iTime) = max(UCC%RH%AveResidual(4, UCC%RH%iTime), abs(DeltaQuantity(4)))
            UCC%RH%MaxResidual(5, UCC%RH%iTime) = max(UCC%RH%AveResidual(5, UCC%RH%iTime), abs(DeltaQuantity(5)))

        end do

        UCC%RH%TotalResidual(:, UCC%RH%iTime) = UCC%RH%AveResidual(:, UCC%RH%iTime) ! 平均化する前の値
        UCC%RH%AveResidual(:, UCC%RH%iTime) = UCC%RH%AveResidual(:, UCC%RH%iTime) / dble(UG%GI%RealCells)
        UCC%RH%MaxResidual(6, UCC%RH%iTime) = maxval(UCC%RH%MaxResidual(1:5, UCC%RH%iTime))

        if(UCC%ConvergeCondition == 0) then !RMS
            Residual = UCC%RH%AveResidual(6, UCC%RH%iTime)
        else if(UCC%ConvergeCondition == 1) then    !MAX
            Residual = UCC%RH%MaxResidual(6, UCC%RH%iTime)
        else if(UCC%ConvergeCondition == 2) then    ! grid total
            Residual = UCC%RH%TotalResidual(6, UCC%RH%iTime)
        end if

        UCC%RH%iTime = UCC%RH%iTime + 1
        !write(6,*) Converge_tolerance, sqrt(Residual/dble(UG%GI%RealCells))
        if(UCC%ConvergeCondition == 0) then
            Residual = sqrt(Residual/dble(UG%GI%RealCells))
            write(6,*) "RMS Residual : ", Residual

        else if(UCC%ConvergeCondition == 1) then
            Residual = sqrt(Residual)
            if(DetailedReport > 3) write(6,*) "MAX Residual : ", Residual

        else if(UCC%ConvergeCondition == 2) then
            if(DetailedReport > 3) write(6,*) "Total Residual : ", Residual

        end if

        !if(Residual <= Converge_tolerance) then
            !if(DetailedReport > 3) write(6,*) "converge... ", Converge_tolerance, Residual
            !UCC%iEndFlag = 2
        !end if

    end if

return
end subroutine UExplicitEulerMethod
