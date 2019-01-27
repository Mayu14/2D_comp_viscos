!***********************************/
!	Name:CFL条件に基づく局所最大時間刻みを求めるプログラム
!	Alias:CalcTimeWidth
!	Description:局所時間刻み幅を使うのか統一時間刻みを使うのかは外部ルーチンで制御
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.02.01
!	Update:
!	Other:
!***********************************/
subroutine CalcTimeWidth(Geom,CC)
use StructVar_Mod
use ConstantVar_Mod, Gamma => SpecificOfHeatRatio, CFL => CourantFriedrichsLewyCondition
implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CalcTimeWidthWithOMP) :: CTW
    integer :: iCenterX, iCenterY
    double precision :: AverageCellWidthHalf

    AverageCellWidthHalf = min(Geom%Width(1,1,1),Geom%Width(1,1,2)) !sqrtと迷ったが，歪んだ直交格子を考えればminが正しいことが分かる


!$omp parallel num_threads(CoreNumberOfCPU),shared(Geom,CC,AverageCellWidthHalf),private(iCenterX,iCenterY,CTW)
!$omp do
    do iCenterY=1, Geom%CellNumber(2)
        do iCenterX=1, Geom%CellNumber(1)
            if(Geom%Interpolated(iCenterX,iCenterY,1) == 0) then
                CTW%AbsoluteVelocity2 = dot_product(CC%ConservedQuantity(2:3,iCenterX,iCenterY,1),CC%ConservedQuantity(2:3,iCenterX,iCenterY,1)) &
                                        &   / (CC%ConservedQuantity(1,iCenterX,iCenterY,1)**2) !速度L2ノルム

                if((Gamma*Gmin1*(CC%ConservedQuantity(4,iCenterX,iCenterY,1) / CC%ConservedQuantity(1,iCenterX,iCenterY,1) &
                                    & - 0.5d0 * CTW%AbsoluteVelocity2)) < 0.0d0) then
                        write(6,*) iCenterX,iCenterY
                        write(6,*) CC%ConservedQuantity(1,iCenterX,iCenterY,1), CC%ConservedQuantity(4,iCenterX,iCenterY,1)
                        write(6,*) CC%ConservedQuantity(2,iCenterX,iCenterY,1), CC%ConservedQuantity(3,iCenterX,iCenterY,1)
                end if

                CTW%LocalSoundSpeed = sqrt(Gamma*Gmin1*(CC%ConservedQuantity(4,iCenterX,iCenterY,1) / CC%ConservedQuantity(1,iCenterX,iCenterY,1) &
                                    & - 0.5d0 * CTW%AbsoluteVelocity2)) !局所音速

                CTW%AbsoluteVelocity = sqrt(CTW%AbsoluteVelocity2) !速度絶対値

                CC%TmpTimeWidth(iCenterX,iCenterY,1,1) = CFL*2.0d0*AverageCellWidthHalf / (CTW%LocalSoundSpeed + CTW%AbsoluteVelocity) !衝撃波が通過する速度

            else
                CC%TmpTimeWidth(iCenterX,iCenterY,1,1) = FixedTimeStep
            end if
        end do
    end do
!$omp end do
!$omp end parallel

    do iCenterY = 1, Geom%CellNumber(2)
        CC%TmpTimeWidth(0,iCenterY,1,1) = CC%TmpTimeWidth(1,iCenterY,1,1)
        CC%TmpTimeWidth(Geom%CellNumber(1)+1,iCenterY,1,1) = CC%TmpTimeWidth(Geom%CellNumber(1),iCenterY,1,1)
    end do

    do iCenterX = 0, Geom%CellNumber(1)+1
        CC%TmpTimeWidth(iCenterX,0,1,1) = CC%TmpTimeWidth(iCenterY,1,1,1)
        CC%TmpTimeWidth(iCenterX,Geom%CellNumber(2)+1,1,1) = CC%TmpTimeWidth(iCenterY,Geom%CellNumber(2),1,1)
    end do

!$omp parallel num_threads(CoreNumberOfCPU),shared(Geom,CC),private(iCenterX,iCenterY)
!$omp do
    do iCenterY=1, Geom%CellNumber(2)
        do iCenterX=1, Geom%CellNumber(1)

            CC%TimeWidth(iCenterX,iCenterY,1) = min(CC%TmpTimeWidth(iCenterX,iCenterY,1,1),&
                    &   CC%TmpTimeWidth(iCenterX-1,iCenterY,1,1), CC%TmpTimeWidth(iCenterX+1,iCenterY,1,1), &
                    &   CC%TmpTimeWidth(iCenterX,iCenterY-1,1,1), CC%TmpTimeWidth(iCenterX,iCenterY+1,1,1))

        end do
    end do
!$omp end do
!$omp end parallel

return
end subroutine CalcTimeWidth
