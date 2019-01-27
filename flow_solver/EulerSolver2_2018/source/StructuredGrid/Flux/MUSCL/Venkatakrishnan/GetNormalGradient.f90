!***********************************/
!	Name:VenkatakrishnanリミッターのΔ_を求めるプログラム
!	Alias:GetNormalGradient
!	Description:セル内の基礎変数勾配にセル中心，界面中心間の位置ベクトルとの要素積を取る
!	Type:CellEdge
!	Input:Geometory,CC%GradientOfVariable,
!	Output:CE%NormalGradient
!	Note:法線ベクトルについての考慮は式自体に内包してあるのでいらない
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2018.01.11
!	Other:
!***********************************/
subroutine GetNormalGradient(Geom,CC,CE)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, only:CoreNumberOfCPU
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(in) :: CC
    type(CellEdge), intent(inout) :: CE
    nullify(iCX,iCY,iCZ,iFX,iFY,iFZ,iVar)
    iFX => iFaceX
    iFY => iFaceY
    iFZ => iFaceZ
    iCX => iFaceX

    if(Geom%Dimension == 1) then
        CE%NormalGradient(:,0,1,1,2) = CC%GradientOfVariable(:,1,1,1,1)*Geom%Width(1,1,1)*(-1.0d0)
        CE%NormalGradient(:,Geom%CellNumber(1),1,1,1) = CC%GradientOfVariable(:,Geom%CellNumber(1),1,1,1)*Geom%Width(1,1,1)
        do iFaceX = 1,Geom%CellNumber(1)-1
            CE%NormalGradient(:,iCX,1,1,1) = CC%GradientOfVariable(:,iCX,1,1,1)*Geom%Width(1,1,1)*(-1.0d0) !要素west界面
            CE%NormalGradient(:,iCX,1,1,2) = CC%GradientOfVariable(:,iCX,1,1,1)*Geom%Width(1,1,1) !要素east界面
        end do

!以下書き直し？NormalGradientの計算方法について
    else if(Geom%Dimension == 2) then

            CE%NormalGradient = 0.0d0 !仮想格子のzerofillも兼ねて
!!$omp parallel num_threads(CoreNumberOfCPU),shared(Geom,CE,CC),firstprivate(iFaceX,iFaceY)
!!$omp do
        do iFaceY = 1, Geom%CellNumber(2) !01/11
            do iFaceX = 1, Geom%CellNumber(1)
                CE%NormalGradient(:,iFaceX,iFaceY,1,1) = CC%GradientOfVariable(:,iFaceX,iFaceY,1,1)*Geom%Width(1,1,1)*(-1.0d0) !West界面
                CE%NormalGradient(:,iFaceX,iFaceY,1,2) = CC%GradientOfVariable(:,iFaceX,iFaceY,1,1)*Geom%Width(1,1,1) !East界面
                CE%NormalGradient(:,iFaceX,iFaceY,1,3) = CC%GradientOfVariable(:,iFaceX,iFaceY,1,2)*Geom%Width(1,1,2)*(-1.0d0) !South界面
                CE%NormalGradient(:,iFaceX,iFaceY,1,4) = CC%GradientOfVariable(:,iFaceX,iFaceY,1,2)*Geom%Width(1,1,2) !North界面
            end do
        end do
!!$omp end do
!!$omp end parallel

    else if(Geom%Dimension == 3) then
        do iFaceZ = 1, Geom%CellNumber(3)
            do iFaceY = 1, Geom%CellNumber(2)
        CE%NormalGradient(:,0,iFY,iFZ,2) = CC%GradientOfVariable(:,1,iFY,iFZ,1)*Geom%Width(1,1,1)*(-1.0d0)
        CE%NormalGradient(:,Geom%CellNumber(1),iFY,iFZ,1) = CC%GradientOfVariable(:,Geom%CellNumber(1),iFY,iFZ,1)*Geom%Width(1,1,1)
            end do
        end do

        do iFaceZ = 1, Geom%CellNumber(3)
            do iFaceX = 1, Geom%CellNumber(1)
        CE%NormalGradient(:,iFX,0,iFZ,4) = CC%GradientOfVariable(:,iFX,1,iFZ,2)*Geom%Width(1,1,2)*(-1.0d0)
        CE%NormalGradient(:,iFX,Geom%CellNumber(2),1,3) = CC%GradientOfVariable(:,iFX,Geom%CellNumber(2),iFZ,1)*Geom%Width(1,1,2)
            end do
        end do

        do iFaceX = 1, Geom%CellNumber(1)
            do iFaceY = 1, Geom%CellNumber(2)
        CE%NormalGradient(:,iFX,iFY,0,6) = CC%GradientOfVariable(:,iFX,iFY,1,3)*Geom%Width(1,1,3)*(-1.0d0)
        CE%NormalGradient(:,iFX,iFY,Geom%CellNumber(3),5) = CC%GradientOfVariable(:,iFX,iFY,Geom%CellNumber(3),3)*Geom%Width(1,1,3)
            end do
        end do

        do iFaceZ = 1, Geom%CellNumber(3)
            do iFaceY = 1, Geom%CellNumber(2)
                do iFaceX = 1, Geom%CellNumber(1)
                    CE%NormalGradient(:,iFX,iFY,iFZ,1) = CC%GradientOfVariable(:,iFX,iFY,iFZ,1)*Geom%Width(1,1,1)*(-1.0d0)
                    CE%NormalGradient(:,iFX,iFY,iFZ,2) = CC%GradientOfVariable(:,iFX,iFY,iFZ,1)*Geom%Width(1,1,1)
                    CE%NormalGradient(:,iFX,iFY,iFZ,3) = CC%GradientOfVariable(:,iFX,iFY,iFZ,2)*Geom%Width(1,1,2)*(-1.0d0)
                    CE%NormalGradient(:,iFX,iFY,iFZ,4) = CC%GradientOfVariable(:,iFX,iFY,iFZ,2)*Geom%Width(1,1,2)
                    CE%NormalGradient(:,iFX,iFY,iFZ,5) = CC%GradientOfVariable(:,iFX,iFY,iFZ,3)*Geom%Width(1,1,3)*(-1.0d0)
                    CE%NormalGradient(:,iFX,iFY,iFZ,6) = CC%GradientOfVariable(:,iFX,iFY,iFZ,3)*Geom%Width(1,1,3)
                end do
            end do
        end do
    end if
    return
end subroutine GetNormalGradient
