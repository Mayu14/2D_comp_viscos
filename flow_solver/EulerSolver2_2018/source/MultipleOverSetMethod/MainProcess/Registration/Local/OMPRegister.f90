!***********************************/
!	Name:内挿元となるドナー要素の検索・登録プログラム
!	Alias:Register4InterpolationMain
!	Description:ドナー要素の共有面積も計算する
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.12.19
!	Update:2017.12.19
!	Other:
!***********************************/
subroutine Register4InterpolationMain(LG,iGrid,Geom)
    use StructVar_Mod
!    use LoopVar_Mod
    use RegistVar_Mod4OMP
!$ use omp_lib
    implicit none
    type(LocalGrids), intent(inout) :: LG
    type(Geometry), intent(inout) :: Geom
    integer, intent(in) :: iGrid

    type(CountUpVertex) :: CUV
    integer :: iDonorGrid, iTimeCombination
    integer :: iTargetCell
    integer :: iKindOfDonor, iFindDonor
    integer :: iGlobalCellNumber
    integer :: iLoop
    type(PrivateVarOfRegist) :: PVR
    allocate(PVR%RG2L%CenterDistanceVector(3),PVR%RG2L%iCenter(3),PVR%RG2L%iCenter2(3))
    allocate(PVR%RIG%CellCenterCoords(3),PVR%RIG%GlobalGridNum(3))
    allocate(PVR%RSDILG%MWDIAG%AuxiliaryGridNum(3),PVR%RSDILG%SSDIG%CenterDistanceVector(3))
!内挿用のDonor格子を見つけて共有面積を保存する

!$omp parallel num_threads(CoreNumberOfCPU),shared(LG,iGrid,Geom),firstprivate(iTargetCell,iKindOfDonor,iFindDonor,CUV,iGlobalCellNumber,PVR,iLoop)
    do iTargetCell=1, LG%UG(iGrid)%GI%RealCells !Li格子の全実要素について

        call InputCUV(1,0,iTargetCell,LG%UG(iGrid),LG%MG(iGrid)%NSR,CUV,LG%MG(iGrid))
        !グローバル格子における位置判定+登録を行う
        call RegisterIntoGlobal(iGrid,LG%GG,iGlobalCellNumber,Geom,CUV,iLoop,PVR%RIG) !CutCell入れるなら消す？

        iFindDonor = 0
        do iKindOfDonor=1, 3
            if(iKindOfDonor <= 2 .or. iFindDonor == 0) then !上位格子+物体あり格子からの補間が終わった時点でまだ補間要素が見つかっていない時下位格子の検索を行う

                do iDonorGrid=1, LG%TotalGridsNumber !計算領域中の全格子について

                    if(LG%UG(iGrid)%Overlap(iTimeCombination+1,iDonorGrid) == iKindOfDonor) then !iKind=1で上位格子，iKind=2で物体ありの下位格子，iKind=3で下位格子
                        call RoughSearchDonorInLocalGrid(iTimeCombination,LG%MG(iDonorGrid),LG%UG(iDonorGrid),LG%MG(iGrid),CUV,iFindDonor)
                        call InputCUV(1,0,iTargetCell,LG%UG(iGrid),LG%MG(iGrid)%NSR,CUV,LG%MG(iGrid)) !Donor格子が変わるたびにその座標系に合わせる必要があるため，毎回グローバル座標系の値に戻す
                    end if

                end do

            else !上位格子+物体あり下位格子の検索が終わった時点で補間用要素が見つかっているとき = 積分不要
                LG%UG(iGrid)%GM%Interpolated(iTargetCell,1,1) = 1
            end if
        end do

    !この時点で要素と被る可能性のあるすべてのローカル格子との接触判定が終わっている
    !そのため，この時点でiFindDonorが0であるならばグローバル格子からの補間が必要となる
        if(LG%UG(iGrid)%GM%CellType(iTargetCell,1,1) == 99) then
            call RegisterGlobal2Local(iGlobalCellNumber,LG%MG(iGrid),iFindDonor,1,CUV,Geom,9999)
            iFindDonor = 1
        end if


        if(iFindDonor == 0) then
            call RegisterGlobal2Local(iGlobalCellNumber,LG%MG(iGrid),iFindDonor,0,CUV,Geom,9999)
        end if

    end do
return
end subroutine Register4InterpolationMain
