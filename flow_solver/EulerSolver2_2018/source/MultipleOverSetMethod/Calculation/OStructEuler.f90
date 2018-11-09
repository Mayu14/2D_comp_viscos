!***********************************/
!	Name:非構造格子から構造格子へのデータ移動
!	Alias:OOverSetU2S
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
subroutine OStructEuler(iStep,OConf,Geom,CC,CE,SW)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod
!$ use omp_lib
    implicit none
    integer, intent(in) :: iStep
    type(Configulation), intent(in) :: OConf
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CellEdge), intent(inout) :: CE
    type(StopWatch), intent(inout) :: SW
    integer :: iSplit !時間分割

        call CheckCFL4FixedTime(Geom,CC,iSplit)
        FixedTimeStep = FixedTimeStep / dble(iSplit)
        !write(6,*) "Split number of StructEEM",iSplit
        do iLoop=1, iSplit

            call SetBoundary(Geom,CC)

!$ SW%LapTime(1,4) = omp_get_wtime()
            call CalcFlux(OConf,Geom,CC,CE)
!$ SW%LapTime(2,4) = omp_get_wtime()
!$ SW%LapTime(3,4) = SW%LapTime(3,4) + SW%LapTime(2,4) - SW%LapTime(1,4)
!$ SW%LapTime(1,5) = omp_get_wtime()
            call TimeIntegral(OConf,Geom,CE,CC,iStep)
!$ SW%LapTime(2,5) = omp_get_wtime()
!$ SW%LapTime(3,5) = SW%LapTime(3,5) + SW%LapTime(2,5) - SW%LapTime(1,5)
        end do

        FixedTimeStep = DefaultTimeStep
return
end subroutine OStructEuler
