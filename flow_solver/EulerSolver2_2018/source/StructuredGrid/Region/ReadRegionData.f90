!***********************************/
!	Name:StrGridから計算領域の設定を読み込むプログラム
!	Alias:ReadRegionData
!	Description:計算次元と計算領域の分割数，領域の上限下限を読み取る
!	Type:Geom
!	Input:StrGrid(外部入力)
!	Output:Geom%Dimension,Geom%Bound,Geom%CellNumber
!	Note:
!	Author:Akitaka Toyota
!	Date:2017.10.26
!	Update:2017.11.09
!	Other:
!***********************************/
subroutine ReadRegionData(Geom)
    use StructVar_Mod
    use LoopVar_Mod
    implicit none

    type(Geometry),intent(inout) :: Geom
    character :: cAnnotate
    double precision :: Swap

    open(unit=1, file='StrGrid',status='unknown')
    read(1,*) cAnnotate,Geom%Dimension

    allocate(Geom%CellNumber(3))
    allocate(Geom%Bound(2,Geom%Dimension))
    allocate(Geom%WallType(2*Geom%Dimension))
    allocate(Geom%BC%InFlowVariable(Geom%Dimension+2),Geom%BC%OutFlowVariable(Geom%Dimension+2))
!ここでいろいろallocate

    Geom%CellNumber = 0.0d0
    Geom%Bound=0.0d0

    do iLoop = 1, Geom%Dimension
        read(1,*) cAnnotate,Geom%CellNumber(iLoop)
    end do
        do iLoop = Geom%Dimension+1, 3
            read(1,*) cAnnotate,cAnnotate
        end do

    do iLoop = 1, Geom%Dimension
        read(1,*) cAnnotate,Geom%Bound(1,iLoop),Geom%Bound(2,iLoop)
    end do
        do iLoop = Geom%Dimension+1,3
            read(1,*) cAnnotate,cAnnotate,cAnnotate
        end do

    do iLoop = 1,  Geom%Dimension*2
        read(1,*) cAnnotate, Geom%WallType(iLoop)
    end do
        do iLoop = Geom%Dimension*2+1, 6
            read(1,*) cAnnotate,cAnnotate
        end do

    read(1,*) cAnnotate !Inflow Condition or OutFlow Condition
        iLoop=0
        do iVariable = 1, 5
            if(Geom%Dimension+1 < iVariable .and. iVariable < 5) then
                read(1,*) cAnnotate, cAnnotate
            else
                iLoop = iLoop + 1
                read(1,*) cAnnotate, Geom%BC%InFlowVariable(iLoop)
            end if
        end do

    read(1,*) cAnnotate !Inflow Condition or OutFlow Condition
        iLoop=0
        do iVariable = 1, 5
            if(Geom%Dimension+1 < iVariable .and. iVariable < 5) then
                read(1,*) cAnnotate, cAnnotate
            else
                iLoop = iLoop + 1
                read(1,*) cAnnotate, Geom%BC%OutFlowVariable(iLoop)
            end if
        end do

    close(1)

    call InOutFlowVariableCorrection

    do iLoop = 1, Geom%Dimension
        if(Geom%Bound(1,iLoop) >= Geom%Bound(2,iLoop)) then
            Swap = Geom%Bound(1,iLoop)
            Geom%Bound(1,iLoop) = Geom%Bound(2,iLoop)
            Geom%Bound(2,iLoop) = Swap
        end if
    end do

return
contains
    subroutine InOutFlowVariableCorrection
    use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
    implicit none
        Geom%BC%InFlowVariable(Geom%Dimension+2) = Geom%BC%InFlowVariable(1)/Gamma !無次元化
        Geom%BC%InFlowVariable(Geom%Dimension+2) = &
        &   Geom%BC%InFlowVariable(Geom%Dimension+2)/Gmin1 + 0.5d0 * Geom%BC%InFlowVariable(1) &
        &   * dot_product(Geom%BC%InFlowVariable(2:Geom%Dimension+1),Geom%BC%InFlowVariable(2:Geom%Dimension+1))

        Geom%BC%InFlowVariable(2:Geom%Dimension+1) = Geom%BC%InFlowVariable(2:Geom%Dimension+1)*Geom%BC%InFlowVariable(1) !速度→運動量

        Geom%BC%OutFlowVariable(Geom%Dimension+2) = Geom%BC%OutFlowVariable(1)/Gamma !圧力の無次元化
        Geom%BC%OutFlowVariable(Geom%Dimension+2) = &
        &   Geom%BC%OutFlowVariable(Geom%Dimension+2)/Gmin1 + 0.5d0 * Geom%BC%OutFlowVariable(1) &
        &   * dot_product(Geom%BC%OutFlowVariable(2:Geom%Dimension+1),Geom%BC%OutFlowVariable(2:Geom%Dimension+1)) !圧力→エネルギー

        Geom%BC%OutFlowVariable(2:Geom%Dimension+1) = Geom%BC%OutFlowVariable(2:Geom%Dimension+1)*Geom%BC%OutFlowVariable(1) !速度→運動量

    return
    end subroutine

end subroutine ReadRegionData
