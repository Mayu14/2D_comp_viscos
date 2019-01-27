!***********************************/
!	Name:流束制限関数を求いて，界面にTVDを満足する基礎変数を再構築するためのプログラム
!	Alias:GetQuantityOnSurface
!	Description:
!	Type:CellEdge
!	Input:Configulation,Geometory,CellCenter,CellEdge
!	Output:CE%RebuildQuantity
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2018.02.09
!	Other:
!***********************************/
subroutine GetQuantityOnSurface(Geom,CC,CE)
    use StructVar_Mod
    use LoopVar_Mod
    use ConstantVar_Mod, only:CoreNumberOfCPU
    implicit none
    type(Geometry), intent(in) :: Geom
    type(CellCenter), intent(in) :: CC
    type(CellEdge), intent(inout) :: CE

    nullify(iVar,iCX,iFX,iCY,iFY,iCZ,iFZ)
    iCX => iFaceX
    iFX => iFaceX
    iCY => iFaceY
    iFY => iFaceY
    iCZ => iFaceZ
    iFZ => iFaceZ

    if(Geom%Dimension == 1) then
        CE%RebuildQunatity(:,0,1,1,2) = CC%PrimitiveVariable(:,0,1,1) &
            & + CC%LimiterFunction(:,1,1,1)*CE%NormalGradient(:,0,1,1,2) !0番の境界における勾配のみ逆方向に計算したため，その補正

        CE%RebuildQunatity(:,Geom%CellNumber(1)+1,1,1,1) = CC%PrimitiveVariable(:,Geom%CellNumber(1)+1,1,1) &
            & - CC%LimiterFunction(:,Geom%CellNumber(1),1,1)*CE%NormalGradient(:,Geom%CellNumber(1),1,1,1) !境界条件入れるべきか？

        do iFaceX = 1, Geom%CellNumber(1)
            CE%RebuildQunatity(:,iFX,1,1,1) = CC%PrimitiveVariable(:,iCX,1,1) &
                & - CC%LimiterFunction(:,iCX,1,1)*CE%NormalGradient(:,iCX,1,1,1)

            CE%RebuildQunatity(:,iFX,1,1,2) = CC%PrimitiveVariable(:,iCX,1,1) &
                & + CC%LimiterFunction(:,iCX,1,1)*CE%NormalGradient(:,iCX,1,1,2)
        end do


    else if(Geom%Dimension == 2) then
        do iFaceY = 1, Geom%CellNumber(2)
            CE%RebuildQunatity(:,0,iFY,1,2) = CC%PrimitiveVariable(:,0,iCY,1)
            CE%RebuildQunatity(:,Geom%CellNumber(1)+1,iFY,1,1) = CC%PrimitiveVariable(:,Geom%CellNumber(1)+1,iCY,1)
        end do

        do iFaceX = 1, Geom%CellNumber(1)
            CE%RebuildQunatity(:,iFX,0,1,4) = CC%PrimitiveVariable(:,iFX,0,1)
            CE%RebuildQunatity(:,iFX,Geom%CellNumber(2)+1,1,3) = CC%PrimitiveVariable(:,iCX,Geom%CellNumber(2)+1,1)
        end do

!$omp parallel num_threads(CoreNumberOfCPU),shared(Geom,CE,CC),private(iFaceX,iFaceY)
!$omp do
        do iFaceY = 1, Geom%CellNumber(2)
            do iFaceX = 1, Geom%CellNumber(1)
                CE%RebuildQunatity(:,iFaceX,iFaceY,1,1) = CC%PrimitiveVariable(:,iFaceX,iFaceY,1) &
                    & - CC%LimiterFunction(:,iFaceX,iFaceY,1)*CE%NormalGradient(:,iFaceX,iFaceY,1,1)

                CE%RebuildQunatity(:,iFaceX,iFaceY,1,2) = CC%PrimitiveVariable(:,iFaceX,iFaceY,1) &
                    & + CC%LimiterFunction(:,iFaceX,iFaceY,1)*CE%NormalGradient(:,iFaceX,iFaceY,1,2)

                CE%RebuildQunatity(:,iFaceX,iFaceY,1,3) = CC%PrimitiveVariable(:,iFaceX,iFaceY,1) &
                    & - CC%LimiterFunction(:,iFaceX,iFaceY,1)*CE%NormalGradient(:,iFaceX,iFaceY,1,3)

                CE%RebuildQunatity(:,iFaceX,iFaceY,1,4) = CC%PrimitiveVariable(:,iFaceX,iFaceY,1) &
                    & + CC%LimiterFunction(:,iFaceX,iFaceY,1)*CE%NormalGradient(:,iFaceX,iFaceY,1,4)

            end do
        end do
!$omp end do
!$omp end parallel

    else if(Geom%Dimension == 3) then
        do iFaceZ = 1, Geom%CellNumber(3)
            do iFaceY = 1, Geom%CellNumber(2)
                CE%RebuildQunatity(:,0,iFY,iFZ,2) = CC%PrimitiveVariable(:,0,iCY,iCZ) &
                    & + CC%LimiterFunction(:,1,iCY,iCZ)*CE%NormalGradient(:,0,iFY,iFZ,2) !0番の境界における勾配のみ逆方向に計算したため，その補正

                CE%RebuildQunatity(:,Geom%CellNumber(1)+1,iFY,iFZ,1) = CC%PrimitiveVariable(:,Geom%CellNumber(1)+1,iCY,iCZ) &
                & - CC%LimiterFunction(:,Geom%CellNumber(1),iCY,iCZ)*CE%NormalGradient(:,Geom%CellNumber(1),iFY,iFZ,1) !境界条件入れるべきか？
            end do
        end do

        do iFaceZ = 1, Geom%CellNumber(3)
            do iFaceX = 1, Geom%CellNumber(1)
                CE%RebuildQunatity(:,iFX,0,iFZ,4) = CC%PrimitiveVariable(:,iCX,0,iCZ) &
                    & + CC%LimiterFunction(:,iCX,1,iCZ)*CE%NormalGradient(:,iFX,0,iFZ,4)

                CE%RebuildQunatity(:,iFX,Geom%CellNumber(2)+1,iFZ,3) = CC%PrimitiveVariable(:,iCX,Geom%CellNumber(2)+1,iCZ) &
                & - CC%LimiterFunction(:,iCX,Geom%CellNumber(2),iCZ)*CE%NormalGradient(:,iFX,Geom%CellNumber(2),iFZ,3) !境界条件入れるべきか？
            end do
        end do

        do iFaceY = 1, Geom%CellNumber(2)
            do iFaceX = 1, Geom%CellNumber(1)
                CE%RebuildQunatity(:,iFX,iFY,0,6) = CC%PrimitiveVariable(:,iCX,iCY,0) &
                    & + CC%LimiterFunction(:,iCX,iCY,1)*CE%NormalGradient(:,iFX,iFY,1,6)

                CE%RebuildQunatity(:,iFX,iFY,Geom%CellNumber(3)+1,5) = CC%PrimitiveVariable(:,iCX,iCY,Geom%CellNumber(3)+1) &
                & - CC%LimiterFunction(:,iCX,iCY,Geom%CellNumber(3))*CE%NormalGradient(:,iFX,iFY,Geom%CellNumber(3),5) !境界条件入れるべきか？
            end do
        end do


        do iFaceZ = 1, Geom%CellNumber(3)
            do iFaceY = 1, Geom%CellNumber(2)
                do iFaceX = 1, Geom%CellNumber(1)
                CE%RebuildQunatity(:,iFX,iFY,iFZ,1) = CC%PrimitiveVariable(:,iCX,iCY,iCZ) &
                    & - CC%LimiterFunction(:,iCX,iCY,iCZ)*CE%NormalGradient(:,iFX,iFY,iFZ,1)

                CE%RebuildQunatity(:,iFX,iFY,iFZ,2) = CC%PrimitiveVariable(:,iCX,iCY,iCZ) &
                    & + CC%LimiterFunction(:,iCX,iCY,iCZ)*CE%NormalGradient(:,iFX,iFY,iFZ,2)

                CE%RebuildQunatity(:,iFX,iFY,iFZ,3) = CC%PrimitiveVariable(:,iCX,iCY,iCZ) &
                    & - CC%LimiterFunction(:,iCX,iCY,iCZ)*CE%NormalGradient(:,iFX,iFY,iFZ,3)

                CE%RebuildQunatity(:,iFX,iFY,iFZ,4) = CC%PrimitiveVariable(:,iCX,iCY,iCZ) &
                    & + CC%LimiterFunction(:,iCX,iCY,iCZ)*CE%NormalGradient(:,iFX,iFY,iFZ,4)

                CE%RebuildQunatity(:,iFX,iFY,iFZ,5) = CC%PrimitiveVariable(:,iCX,iCY,iCZ) &
                    & - CC%LimiterFunction(:,iCX,iCY,iCZ)*CE%NormalGradient(:,iFX,iFY,iFZ,5)

                CE%RebuildQunatity(:,iFX,iFY,iFZ,6) = CC%PrimitiveVariable(:,iCX,iCY,iCZ) &
                    & + CC%LimiterFunction(:,iCX,iCY,iCZ)*CE%NormalGradient(:,iFX,iFY,iFZ,6)
                end do
            end do
        end do
    end if
return
end subroutine GetQuantityOnSurface
