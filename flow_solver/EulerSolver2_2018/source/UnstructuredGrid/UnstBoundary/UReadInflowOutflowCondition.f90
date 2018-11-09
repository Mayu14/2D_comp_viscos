!***********************************/
!	Name:流入条件・流出条件の読み込み
!	Alias:UReadInflowOutflowCondition
!	Description:
!	Type:
!	Input:
!	Output:
!	Note:
!	Author:Akitaka Toyota
!	Date:2018.01.20
!	Update:
!	Other:
!***********************************/
subroutine UReadInflowOutflowCondition(UG)
use StructVar_Mod
use LoopVar_Mod
use ConstantVar_Mod, Gamma => SpecificOfHeatRatio
implicit none
    type(UnstructuredGrid), intent(inout) :: UG
    character(len=64) :: cFileName
    character :: cAnnotate

    !write(6,*) "Please input the name of In/Outflow condition's datafile"
    !read(5,*) cFileName
    cFileName = "InOutFlowCondition"

    open(unit=1,file=trim(adjustl(cFileName)),status='unknown')
        read(1,*) cAnnotate
        do iLoop=1,5
            read(1,*) cAnnotate, UG%GM%BC%InFlowVariable(iLoop)
        end do
            UG%GM%BC%InFlowVariable(5) = UG%GM%BC%InFlowVariable(1)/Gamma !無次元化
            !UG%GM%BC%InFlowVariable(5) = &
            !&   UG%GM%BC%InFlowVariable(5)/Gmin1 + 0.5d0 * UG%GM%BC%InFlowVariable(1) &
            !&   * dot_product(UG%GM%BC%InFlowVariable(2:4),UG%GM%BC%InFlowVariable(2:4))

            !UG%GM%BC%InFlowVariable(2:4) = UG%GM%BC%InFlowVariable(2:4)*UG%GM%BC%InFlowVariable(1) !速度→運動量

        read(1,*) cAnnotate
        do iLoop=1,5
            read(1,*) cAnnotate, UG%GM%BC%OutFlowVariable(iLoop)
        end do
            UG%GM%BC%OutFlowVariable(5) = UG%GM%BC%OutFlowVariable(1)/Gamma !圧力の無次元化
            !UG%GM%BC%OutFlowVariable(5) = &
            !&   UG%GM%BC%OutFlowVariable(5)/Gmin1 + 0.5d0 * UG%GM%BC%OutFlowVariable(1) &
            !&   * dot_product(UG%GM%BC%OutFlowVariable(2:4),UG%GM%BC%OutFlowVariable(2:4)) !圧力→エネルギー

            !UG%GM%BC%OutFlowVariable(2:4) = UG%GM%BC%OutFlowVariable(2:4)*UG%GM%BC%OutFlowVariable(1) !速度→運動量
    close(1)
!write(6,*) (UG%GM%BC%InFlowVariable(iLoop),iLoop=1,5)
!write(6,*) (UG%GM%BC%OutFlowVariable(iLoop),iLoop=1,5)
!stop

return
end subroutine UReadInflowOutflowCondition
