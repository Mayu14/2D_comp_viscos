!***********************************/
!	Name:Venkatakrishnanリミッターを求めるプログラム
!	Alias:Venkatakrishnan
!	Description:既に求めてあるΔmax,Δmin,Δ_を元に各格子における流束制限関数を求める
!	Type:CellCenter
!	Input:Geometory,CC%VariationOfVariable,CC%GradientOfVariable,CE%NormalGradient
!	Output:CC%LimiterFunction
!	Note:法線ベクトルについての考慮は式自体に内包してあるのでいらない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2018.01.11
!	Other:
!***********************************/
subroutine Venkatakrishnan(Geom,CC,CE)
    use StructVar_Mod
    use LoopVar_Mod
    use UVenkatakrishnanVar_Mod4OMP
    use ConstantVar_Mod, VPK => VenkatakrishnanParameterK
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(inout) :: CC
    type(CellEdge), intent(inout) :: CE

    type(UVenkatakrishnanWithOMP) :: UV

        UV%VE2 = (VPK * 2.0d0*min(Geom%Width(1,1,1) , Geom%Width(1,1,2)))**3
        nullify(iCX,iCY,iCZ,iFX,iFY,iFZ,iVar)
        iVar => iVariable
        iFX => iFaceX
        iCX => iFaceX
        iFY => iFaceY
        iCY => iFaceY
        iFZ => iFaceZ
        iCZ => iFaceZ
    if(Geom%Dimension == 1) then
!端のセルではセル周囲変動量(Variation)を定義できないため計算不可であるため1.0(制限関数の最大値)を与える
!なんなら端だけ処理を分けてもいいけれど，セル中心の制限関数を計算するときのmin関数で影響が消えるはずなので問題ない...？ !問題出てたみたいです(01/09)
        CE%TmpLimiterFunction(:,0,1,1,1) = 1.0d0

        do iFaceX = 1, Geom%CellNumber(1)
            do iSide=1,2
            do iVariable=1,3
            if(CE%NormalGradient(iVar,iFX,1,1,iSide) > 0.0d0) then
                UV%iMaxMin => iOne
            else if(CE%NormalGradient(iVar,iFX,1,1,iSide) <= 0.0d0) then
                UV%iMaxMin => iTwo
            end if

            CE%TmpLimiterFunction(iVar,iFX,1,1,iSide) = &
                &   (CC%VariationOfVariable(UV%iMaxMin,iVar,iCX,1,1)**2 + UV%VE2 &
                    &   + 2.0d0 * CE%NormalGradient(iVar,iFX,1,1,iSide)*CC%VariationOfVariable(UV%iMaxMin,iVar,iCX,1,1)) &
                & / (CC%VariationOfVariable(UV%iMaxMin,iVar,iCX,1,1)**2  + 2.0d0*CE%NormalGradient(iVar,iFX,1,1,iSide)**2&
                    &   + CC%VariationOfVariable(UV%iMaxMin,iVar,iCX,1,1)*CE%NormalGradient(iVar,iFX,1,1,iSide) + UV%VE2)
            end do
            end do
        end do

        iCX => iCenterX
        iFX => iCenterX

        do iCenterX = 1, Geom%CellNumber(1)-1
            CC%LimiterFunction(:,iCX,1,1) = min(CE%TmpLimiterFunction(:,iFX-1,1,1,2),CE%TmpLimiterFunction(:,iFX,1,1,1),&
                &   CE%TmpLimiterFunction(:,iFX,1,1,2),CE%TmpLimiterFunction(:,iFX+1,1,1,1))
        end do
            CC%LimiterFunction(:,iCX,1,1) = min(CE%TmpLimiterFunction(:,iFX-1,1,1,1), &
                &   CE%TmpLimiterFunction(:,iFX-1,1,1,2),CE%TmpLimiterFunction(:,iFX,1,1,1))

    else if(Geom%Dimension == 2) then
        CE%TmpLimiterFunction = 1.0d0
!!$omp parallel num_threads(CoreNumberOfCPU),shared(Geom,CE,CC),firstprivate(iFaceX,iFaceY,iCenterX,iCenterY,iVariable,iLoop,UV)
!!$omp do
        do iFaceY = 1, Geom%CellNumber(2) !セル中心についてのループ
            do iFaceX = 1, Geom%CellNumber(1)
                do iVariable=1, Geom%Dimension+2 !変数についてのループ
                    do iLoop=1, 4 !界面についてのループ

                        if(CE%NormalGradient(iVariable,iFaceX,iFaceY,1,iLoop) == 0.0d0) then
                            CE%TmpLimiterFunction(iVariable,iFaceX,iFaceY,1,iLoop) = 1.0d0

                        else ! NormalGradient /= 0.0d0

                            if(CE%NormalGradient(iVariable,iFaceX,iFaceY,1,iLoop) > 0.0d0) then
                                UV%iMaxMin => iOne !Δmax
                            else !if(CE%NormalGradient(iVar,iFX,iFY,1,iLoop) < 0.0d0) then
                                UV%iMaxMin => iTwo !Δmin
                            end if

                            CE%TmpLimiterFunction(iVariable,iFaceX,iFaceY,1,iLoop) = &
                                &   (CC%VariationOfVariable(UV%iMaxMin,iVariable,iFaceX,iFaceY,1)**2 + UV%VE2 &
                                    &   + 2.0d0 * CE%NormalGradient(iVariable,iFaceX,iFaceY,1,iLoop)*CC%VariationOfVariable(UV%iMaxMin,iVariable,iFaceX,iFaceY,1)) &
                                & / (CC%VariationOfVariable(UV%iMaxMin,iVariable,iFaceX,iFaceY,1)**2  + 2.0d0*CE%NormalGradient(iVariable,iFaceX,iFaceY,1,iLoop)**2&
                                    &   + CC%VariationOfVariable(UV%iMaxMin,iVariable,iFaceX,iFaceY,1)*CE%NormalGradient(iVariable,iFaceX,iFaceY,1,iLoop) + UV%VE2)

                        end if
                    end do
                end do
            end do
        end do
!!$omp end do
!!$omp end parallel
!!$omp barrier

!Make Imaginary Cell's Limiter

!!$omp parallel num_threads(CoreNumberOfCPU),shared(Geom,CE,CC),private(iFaceX,iFaceY,iCenterX,iCenterY)
!!$omp do
        do iCenterY =1 ,Geom%CellNumber(2)
            do iCenterX = 1, Geom%CellNumber(1)
               CC%LimiterFunction(:,iCenterX,iCenterY,1) = &
                    &   min(CE%TmpLimiterFunction(:,iCenterX,iCenterY,1,1), CE%TmpLimiterFunction(:,iCenterX-1,iCenterY,1,2) &
                    &       ,CE%TmpLimiterFunction(:,iCenterX,iCenterY,1,2), CE%TmpLimiterFunction(:,iCenterX+1,iCenterY,1,1) &
                    &       ,CE%TmpLimiterFunction(:,iCenterX,iCenterY,1,3), CE%TmpLimiterFunction(:,iCenterX,iCenterY-1,1,4) &
                    &       ,CE%TmpLimiterFunction(:,iCenterX,iCenterY,1,4), CE%TmpLimiterFunction(:,iCenterX,iCenterY+1,1,3))
            end do
        end do

!!$omp end do
!!$omp end parallel
        nullify(iCX,iCY,iFX,iFY)

    else if(Geom%Dimension == 3) then
        CE%TmpLimiterFunction(:,0,:,:,1) = 1.0d0
        CE%TmpLimiterFunction(:,:,0,:,2) = 1.0d0
        CE%TmpLimiterFunction(:,:,:,0,3) = 1.0d0

        do iFaceZ =1, Geom%CellNumber(3)
            do iFaceY = 1, Geom%CellNumber(2)
                do iFaceX = 1, Geom%CellNumber(1)
                    do iVariable = 1, Geom%Dimension+2
                        do iLoop = 1, Geom%Dimension

            if(CE%NormalGradient(iVar,iFX,iFY,iFZ,iLoop) > 0.0d0) then
                UV%iMaxMin => iOne
            else if(CE%NormalGradient(iVar,iFX,iFY,iFZ,iLoop) <= 0.0d0) then
                UV%iMaxMin => iTwo
            end if

            CE%TmpLimiterFunction(iVar,iFX,iFY,iFZ,iLoop) = &
                &   (CC%VariationOfVariable(UV%iMaxMin,iVar,iCX,iCY,iCZ)**2 + UV%VE2 &
                    &   + 2.0d0 * CE%NormalGradient(iVar,iFX,iFY,iFZ,iLoop)*CC%VariationOfVariable(UV%iMaxMin,iVar,iCX,iCY,iCZ)) &
                & / (CC%VariationOfVariable(UV%iMaxMin,iVar,iCX,iCY,iCZ)**2  + 2.0d0*CE%NormalGradient(iVar,iFX,iFY,iFZ,iLoop)**2&
                    &   + CC%VariationOfVariable(UV%iMaxMin,iVar,iCX,iCY,iCZ)*CE%NormalGradient(iVar,iFX,iFY,iFZ,iLoop) + UV%VE2)

                        end do
                    end do
                end do
            end do
        end do

        iCX => iCenterX
        iCY => iCenterY
        iFX => iCenterX
        iFY => iCenterY
        iCZ => iCenterZ
        iFZ => iCenterZ
        do iCenterZ = 1, Geom%CellNumber(3)
            do iCenterY = 1 ,Geom%CellNumber(2)
                do iCenterX = 1, Geom%CellNumber(1)
                CC%LimiterFunction(:,iCX,iCY,iCZ) = &
                    &   min(CE%TmpLimiterFunction(:,iFX-1,iFY,iFZ,1),CE%TmpLimiterFunction(:,iFX,iFY,iFZ,1),&
                    &       CE%TmpLimiterFunction(:,iFX,iFY-1,iFZ,2),CE%TmpLimiterFunction(:,iFX,iFY,iFZ,2),&
                    &       CE%TmpLimiterFunction(:,iFX,iFY,iFZ-1,3),CE%TmpLimiterFunction(:,iFX,iFY,iFZ,3))
                    end do
            end do
        end do

    end if
    return
end subroutine Venkatakrishnan
