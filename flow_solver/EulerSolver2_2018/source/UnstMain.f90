program UEulerEQ1D
    use StructVar_Mod
    use ConstantVar_Mod
    implicit none

    type(CellCenter) :: UCC
    type(CellEdge) :: UCE
    type(Configulation) :: UConf
    type(Geometry) :: UGeom
    integer :: iStep,iStartStep,iStep4Plot
    double precision :: StartTime,EndTime,TotalTime

!    call UReadConfigulation(Conf)

    call USetRegion(UConf,UGeom)

    !call UTestRegionData(Geom)

!    call UInitialize(Conf,Geom,CC,CE)

!    if(Conf%UseResume == 0) then
!        iStartStep = 0
!        if(Geom%Dimension == 1) call OutputGnuPlot(Geom,iStartStep,CC)
!        if(Geom%Dimension >= 2) call OutputResult(Conf,Geom,iStartStep,CC)
!        iStartStep = 1
!    else if(Conf%UseResume == 1) then
!        if(Conf%UseRRK2 == 0) iStartStep = (Conf%ResumeNumber-1)*OutputInterval+1
!        if(Conf%UseRRK2 == 1) iStartStep = (Conf%ResumeNumber-1)*2*OutputInterval+1
!    end if

!    TotalTime = 0.0d0
!    do iStep = iStartStep, IterationNumber
!        call cpu_time(StartTime)
        !write(6,*) iStep,"th Calculate..."

!        call SetBoundary(Geom,CC)
        !write(6,*) "Complete to input Boundary data"
        !call TestBoundary(Geom,CC)

!        call CalcFlux(Conf,Geom,CC,CE)
        !write(6,*) "Complete to Calculate Flux"

 !       call TimeIntegral(Conf,Geom,CE,CC,iStep)
!        call cpu_time(EndTime)
!        TotalTime = TotalTime + (EndTime - StartTime)
        !write(6,*) "Complete to Integration"
        !if(mod(iStep,OutputInterval) == 0) then
         !   if(Conf%UseRRK2 == 1) then
          !      if(mod(iStep,2*OutputInterval) == 0) then
           !         iStep4Plot = iStep / (2*OutputInterval)
            !        if(Geom%Dimension == 1) call OutputGnuPlot(Geom,iStep4Plot,CC)
             !       if(Geom%Dimension >= 2) call OutputResult(Conf,Geom,iStep4Plot,CC)
              !  end if
            !else
             !   iStep4Plot = iStep / OutputInterval
              !  if(Geom%Dimension == 1) call OutputGnuPlot(Geom,iStep4Plot,CC)
               ! if(Geom%Dimension >= 2) call OutputResult(Conf,Geom,iStep4Plot,CC)
    !        write(6,*) FixedTimeStep*float(iStep),"sec."
            !end if
        !end if

   !end do
    !write(6,*) FixedTimeStep*float(iStep4Plot * OutputInterval),"sec."

    !if(Geom%Dimension == 1) then
     !   if(iStep > 10) call MergeGnuPlot(Geom,iStep4Plot)
      !  call Output4Resume(Geom,iStep4Plot)
    !end if

    !open(unit=1,file="CalcInfo.txt",status='unknown')
     !   write(1,"('この計算は',(1x,i1),'次元衝撃波管問題を')") Geom%Dimension
      !  write(1,"('格子数 ',(1x,i6))") Geom%CellNumber(1)
       ! write(1,"('時間刻み幅 ',(1x,f8.6))") FixedTimeStep
        !write(1,"('空間 ',(1x,i2),'次精度')") Conf%UseMUSCL+1
!        write(1,"('時間 ',(1x,i2),'次精度')") Conf%UseRRK2+1
 !       write(1,"('の下で無次元時間',(1x,i7),'まで解きました．')") iStep4Plot*OutputInterval
  !      write(1,"('計算にあたりcpuを使用した時間は ',(1x,f22.1),'秒でした')") TotalTime
   !     write(1,*) ""
    !    write(1,"('計算を続行したい場合はCalcConfigにて')")
     !   write(1,"('UseResume=',(1x,i1))") 1
      !  write(1,"('ResumeFileFormat=',(1x,i1))") 3
       ! write(1,"('ResumeNumber=',(1x,i12))") iStep
        !write(1,*) "と設定のうえ，適宜反復回数を増やしてプログラムを再実行してください"
!        write(1,*) ""
 !       write(1,*) "※その他の設定を変更した場合，計算を続行できるかは未検証です"
  !  close(1)

stop
end program
