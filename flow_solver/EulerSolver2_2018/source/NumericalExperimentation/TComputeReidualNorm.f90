!***********************************/
!	Name:厳密値からの残差ノルムを計算する
!	Alias:
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:
!	Date:
!	Update:
!	Other:
!***********************************/
subroutine TComputeResidualNorm(iStep,iTargetGrid,CC,UG,MG,Geom)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none

    integer, intent(in) :: iStep,iTargetGrid
    type(CellCenter),intent(in) :: CC
    type(UnstructuredGrid), intent(in) :: UG
    type(MoveGrid), intent(in) :: MG
    type(Geometry), intent(in) :: Geom

    integer :: iInitialCondition !MShockTubeICのICと統一すること
    double precision :: midValue !中央座標の位置
    double precision :: ResidualNorm, ExactValue, TotalVolume

    character(len=64) :: cDirectory,cFileName
    character(len=32) :: cGrid

    iInitialCondition=2 !ここで調節

    if(iStep == 1) then
        write(cGrid,*) iTargetGrid
        cDirectory = "EXP/"!//trim(adjustl(cGrid))//"/" !UConf%SaveDirectiry
        cFileName = trim(adjustl(cDirectory))//"ResidualNormHistory"//trim(adjustl(cGrid))//".txt"
        open(unit = iTargetGrid+100, file =trim(adjustl(cFileName)), status = 'unknown')
    end if

    ResidualNorm = 0.0d0
    TotalVolume = 0.0d0

    if(iTargetGrid == 0) then !構造格子の場合
        iXmax = Geom%CellNumber(1)
        iYmax = max(1,Geom%CellNumber(2)) !セル総数と1の大きい方がiYmax
        iZmax = max(1,Geom%CellNumber(3))

        do iCenterZ = 1, iZmax
            do iCenterY = 1, iYmax
                do iCenterX = 1, iXmax
                    if(iInitialCondition == 1) then
                        if(iCenterX <= iXmax/2) then
                            ExactValue = 1.0d0
                        else
                            ExactValue = 0.1d0
                        end if
                    else !if(iInitialCondition == 2) then
                        if(iCenterY <= iYmax/2) then
                            ExactValue = 1.0d0
                        else
                            ExactValue = 0.1d0
                        end if
                    end if
                    ResidualNorm = ResidualNorm + sqrt((CC%ConservedQuantity(1,iCenterX,iCenterY,iCenterZ) - ExactValue)**2)
                end do
            end do
            if(iDim <= 2) exit
        end do

        ResidualNorm = ResidualNorm / (iXmax*iYmax*iZmax)

    else !非構造格子の場合

        if(iInitialCondition == 2) then !H:0<=Y<0.5 L:0.5<Y<=1.0
            midValue = 0.5d0*(Geom%Bound(1,2)+Geom%Bound(2,2))
        else !H:0<=X<0.5 L:0.5<X<=1.0
            midValue = 0.5d0*(Geom%Bound(1,1)+Geom%Bound(2,1))
        end if

        do iCell=1, UG%GI%RealCells
            if(MG%RC%Cell(iCell,iInitialCondition) < midValue) then !H side
                ExactValue = 1.0d0
            else
                ExactValue = 0.1d0
            end if
                ResidualNorm = ResidualNorm + sqrt((CC%ConservedQuantity(1,iCell,1,1) - ExactValue)**2)*UG%GM%Volume(iCell)
                TotalVolume = TotalVolume + UG%GM%Volume(iCell)
        end do
        ResidualNorm = ResidualNorm / TotalVolume

    end if

    write(iTargetGrid+100,"(3(1x,f22.17))") ResidualNorm,maxval(CC%ConservedQuantity(1,:,1,1)),minval(CC%ConservedQuantity(1,:,1,1))

return
end subroutine TComputeResidualNorm
